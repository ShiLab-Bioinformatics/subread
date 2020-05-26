/***************************************************************

   The Subread software package is free software package: 
   you can redistribute it and/or modify it under the terms
   of the GNU General Public License as published by the 
   Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   Subread is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty
   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   
   See the GNU General Public License for more details.

   Authors: Drs Yang Liao and Wei Shi

  ***************************************************************/
  
  
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>


#include "subread.h"
#include "sublog.h"
#include "core.h"
#include "input-files.h"
#include "sorted-hashtable.h"

#include "core-indel.h"
#include "core-junction.h"

static struct option long_options[] =
{
	{"memory-optimisation",  required_argument, 0, 0},
	{"paired-end", no_argument,0,'p'},
	{0, 0, 0, 0}
};

int load_global_context_forindel(global_context_t * context)
{
	char tmp_fname [MAX_FILE_NAME_LENGTH];

	warning_file_limit();
	context -> input_reads.avg_read_length = 200;//guess_reads_density_format(context->config.first_read_file , 0, NULL, NULL, NULL);
	if(context -> input_reads.avg_read_length<0 )context -> input_reads.avg_read_length = 250;
	if(context -> input_reads.avg_read_length<0 ||geinput_open_sam(context->config.first_read_file, &context->input_reads.first_read_file, context -> input_reads.is_paired_end_reads))
	{
		SUBREADprintf("Unable to open '%s' as a read file. Please check if it exists, you have the permission to read it, and it is in the FASTQ/FASTA/PLAIN format.\n", context->config.first_read_file);
		return -1;
	}


	struct stat ginp1_stat;
	stat(context->config.first_read_file , &ginp1_stat);
	context->input_reads.first_read_file_size = ginp1_stat.st_size;

	sprintf(tmp_fname, "%s.reads", context->config.index_prefix);
	if(!does_file_exist(tmp_fname))
	{
		SUBREADprintf("Unable top open index '%s'. Please make sure that the correct prefix is specified and you have the permission to read these files. For example, if there are files '/opt/my_index.reads', '/opt/my_index.files' and etc, the index prefix should be specified as '/opt/my_index' without any suffix. \n", context->config.index_prefix);
		return -1;
	}


	if(context->config.space_type == GENE_SPACE_COLOR)
		sprintf(tmp_fname, "%s.00.c.tab", context->config.index_prefix);
	else
		sprintf(tmp_fname, "%s.00.b.tab", context->config.index_prefix);
	if(!does_file_exist(tmp_fname))
	{
		SUBREADprintf("Your reads are in the %s space but the index was not built in the same space. Unable to precess the reads.\n", context->config.space_type == GENE_SPACE_COLOR?"color":"base");
		return -1;
	}

	context->index_block_number = 0; 
	while(1)
	{
		sprintf(tmp_fname, "%s.%02d.%c.tab", context->config.index_prefix, context->index_block_number, context->config.space_type == GENE_SPACE_COLOR?'c':'b');
		if(!does_file_exist(tmp_fname))break;
		context->index_block_number ++;
	}

	context->current_index_block_number = 0;
	load_offsets(&context->chromosome_table, context->config.index_prefix);

	
	context->all_processed_reads = 0;
	context->all_mapped_reads = 0;

	memset( context->all_value_indexes , 0 , 100 * sizeof(gene_value_index_t));

	return 0;
}

int load_indexes_forindel(global_context_t * global_context)
{
	int block_no;
	for(block_no = 0; block_no< global_context->index_block_number; block_no++)
	{
		char tmp_fname[MAX_FILE_NAME_LENGTH];
		sprintf(tmp_fname, "%s.%02d.%c.array", global_context->config.index_prefix, block_no,  global_context->config.space_type == GENE_SPACE_COLOR?'c':'b');
		if(gvindex_load(&global_context -> all_value_indexes[block_no], tmp_fname)) return -1;
	}
	return 0;
}

int transfer_SAM_to_assembly(global_context_t * global_context)
{
	char linebuf[max(2*MAX_READ_LENGTH+300,3000)];
	int allreads=0,mapped=0;
	char mate_chro[MAX_CHROMOSOME_NAME_LEN+1];
	unsigned int mate_pos;
	int flags_mate;
	char cigar_mate[EXON_MAX_CIGAR_LEN+1];
	int is_warning_shown = 0;

	HashTable * local_reassembly_pileup_files;

	print_in_box(80,0,0,"Transforming the sam file into: %s", global_context->config.temp_file_prefix);

	local_reassembly_pileup_files = HashTableCreate(100);
	HashTableSetDeallocationFunctions(local_reassembly_pileup_files, NULL, NULL);
	HashTableSetKeyComparisonFunction(local_reassembly_pileup_files, my_strcmp);
	HashTableSetHashFunction(local_reassembly_pileup_files ,HashTableStringHashFunction);

	while(1)
	{
		unsigned int read_anchor_position =0xffffffff;
		char read_name[MAX_READ_NAME_LEN+1];
		int flags;
		char chro_name[MAX_CHROMOSOME_NAME_LEN+1];
		unsigned int pos, pair_dist;
		char cigar[EXON_MAX_CIGAR_LEN+1];
		char read_text[MAX_READ_LENGTH+1];
		char qual_text[MAX_READ_LENGTH+1];
		int read_len, is_repeated, mapping_quality, is_anchor_certain=1;
		int pos_delta = 0;

		if(feof((FILE *)global_context -> input_reads.first_read_file.input_fp))break;

		read_text[0]=0;
		qual_text[0]=0;
		geinput_readline(&global_context -> input_reads.first_read_file, linebuf,0);
		int res = parse_SAM_line(linebuf, read_name,& flags, chro_name, & pos, cigar, & mapping_quality, & pair_dist, read_text , qual_text, & read_len, & is_repeated);
		int cigar_cursor = 0;
		int firstM = 1,xx=0;

		if(res == 0&& (flags &1) && (!global_context -> input_reads.is_paired_end_reads))
		{
			if(!is_warning_shown)
			{
				SUBREADprintf("Warning: the input seems to be paired end. You may like to specify '-p' option.\n");
				is_warning_shown=1;	
			}
		}

		if(res==0)
		{
			for(; cigar[cigar_cursor]; cigar_cursor++)
			{
				char nch = cigar[cigar_cursor]; 
				if(nch>='0'&&nch<='9')
				{
					xx=xx*10+(nch-'0');
				}else
				{
					if(nch=='M') firstM=0;
					else if(nch=='S')
					{
						if(firstM)
							pos_delta = xx;
					}
					xx=0;
				}
			}

			pos -= pos_delta;
		}

		if(res == 1) {chro_name[0]='*'; chro_name[1]=0;}
		//printf("RES=%d\tMAPPED=%s,%u\n", res, chro_name, pos);

		if(res == 0)	// mapped
		{
			read_anchor_position = linear_gene_position(&global_context->chromosome_table , chro_name, pos) - 1;
			//printf("MAPPED_LINEAR=%u\n", read_anchor_position);

			mapped++;
		}
		else if(res == 1 && global_context -> input_reads.is_paired_end_reads)	// unmapped
		{
			is_anchor_certain=0;
			if(mate_chro[0])
			{
				if(mate_chro[0]!='*')
				{
					unsigned int mate_linear_pos = linear_gene_position(&global_context->chromosome_table , mate_chro, mate_pos) - 1;
					read_anchor_position = mate_linear_pos + global_context -> config.expected_pair_distance * ((((flags_mate & SAM_FLAG_FIRST_READ_IN_PAIR) && (flags_mate & SAM_FLAG_REVERSE_STRAND_MATCHED))|| ((flags_mate & SAM_FLAG_SECOND_READ_IN_PAIR)&& (flags_mate & SAM_FLAG_REVERSE_STRAND_MATCHED))) ?-1:1);
				}
				//printf("RECOVERED 1: %u - %s ; LEN=%d ; YOU_ARE_FIRST=%d\n%s\n%s\n", read_anchor_position,  read_name, read_len, flags_mate & SAM_FLAG_FIRST_READ_IN_PAIR, read_text, qual_text);
			}
			else
			{
				char read_text_null[MAX_READ_LENGTH+1];
				char qual_text_null[MAX_READ_LENGTH+1];
				int X_read_len=0;

				geinput_readline_back(&global_context -> input_reads.first_read_file, linebuf);
				res = parse_SAM_line(linebuf, read_name,& flags_mate, mate_chro, & mate_pos, cigar_mate, & mapping_quality, & pair_dist, read_text_null , qual_text_null, & X_read_len, & is_repeated);
				if(res==0)
				{
					unsigned int mate_linear_pos = linear_gene_position(&global_context->chromosome_table , mate_chro, mate_pos) - 1;
					read_anchor_position = mate_linear_pos + global_context -> config.expected_pair_distance * ((((flags_mate & SAM_FLAG_FIRST_READ_IN_PAIR) && (flags_mate & SAM_FLAG_REVERSE_STRAND_MATCHED))|| ((flags_mate & SAM_FLAG_SECOND_READ_IN_PAIR)&& (flags_mate & SAM_FLAG_REVERSE_STRAND_MATCHED))) ?-1:1);
				}
				//printf("RECOVERED 2: %u - %s ; LEN=%d ; YOU_ARE_FIRST=%d\n%s\n%s\n", read_anchor_position,  read_name, read_len, flags_mate & SAM_FLAG_FIRST_READ_IN_PAIR, read_text, qual_text);

				//read_anchor_position += (global_context -> config.expected_pair_distance/3) *((1.0*(read_anchor_position % 3751)  - 3751/2) / 3751.);
			}
			if(read_anchor_position<0xffff0000)
			{
				if((flags_mate & SAM_FLAG_REVERSE_STRAND_MATCHED) == (flags & SAM_FLAG_REVERSE_STRAND_MATCHED))
				{
					reverse_read(read_text, read_len, global_context->config.space_type);
					reverse_quality(qual_text, read_len);
				}
			}
		}

		//printf("R=%s; RL=%d; anchor=%u\n", read_name, read_len, read_anchor_position);

		if(read_anchor_position<0xffff0000 && read_len >=2)
		{
			//if(!is_anchor_certain)printf("Uncertain Anchor at %u = %s\n", read_anchor_position, read_name);
			//printf("INSERT LOCAL: %u\n", read_anchor_position);
			write_local_reassembly(global_context, local_reassembly_pileup_files, read_anchor_position, read_name , read_text , qual_text , read_len, is_anchor_certain);
		}
		if(allreads %2==1)
			mate_chro[0] = 0;
		else
		{
			strcpy(mate_chro, chro_name);
			mate_pos = pos ;
			flags_mate = flags;
			
		}
		allreads++;

	}
	print_in_box(80,0,0,"Processed %d reads; %d mapped.\n", allreads, mapped);

	destroy_pileup_table(local_reassembly_pileup_files);	

	return 0;
}
int print_configuration_forindel(global_context_t * global_context)
{
	SUBREADputs("");
	print_subread_logo();
	SUBREADputs("");


	print_in_box(80,1,1,"subindel setting");
	print_in_box(80,0,1,"");
	print_in_box(80,0,0,"                Input file : %s", global_context -> config.first_read_file );
	print_in_box(80,0,0,"               Output file : %s", global_context -> config.output_prefix);
	print_in_box(80,0,0,"                Index name : %s", global_context -> config.index_prefix);
	print_in_box(80,0,1,"");
	print_in_box(80,0,0,"      Maximum indel length : %d", global_context->config.max_indel_length);
	print_in_box(80,0,0,"           Paired-end mode : %s", global_context -> input_reads.is_paired_end_reads?"yes":"no");

	if(global_context -> input_reads.is_paired_end_reads)
		print_in_box(80,0,0,"  Expected Paired distance : %d", global_context->config.expected_pair_distance);

	print_in_box(80,0,1,"");
	print_in_box(80,2,1,"");
	SUBREADputs("");
	print_in_box(80,1,1,"Running");
	print_in_box(80,0,1,"");


	return 0;
}

int print_summary(global_context_t * global_context)
{
	double timepass = miltime() - global_context->start_time;
 
	print_in_box(80,0,1,"");
	print_in_box(80,2,1,"");

	SUBREADputs("");
	print_in_box(80,1,1,"Finished");
	print_in_box(80,0,1,"");
	print_in_box(80,0,0,"            De novo indels : %u", global_context -> all_indels);
	print_in_box(80,0,0,"                 Time cost : %.1f minutes.", timepass/60);
	print_in_box(80,0,1,"");
	print_in_box(80,2,1,"");

	return 0;
}

void print_usage_subindel()
{
	printf("\nSubIndel: detect short and long indels\n");
	printf("\nUsage:\n  subindel -i <SAM file> -g <subread index> -o <output VCF> {-d <expected fragment distance>} {-I <max indel length>} {--paired-end}\n\n");
	printf("Example:\n  subindel -i my_paired_end_reads.SAM -g my_index -o my_result -d 300 -I 200 --paired-end \n\n");
}


int main(int argc, char ** argv)
{

	int c;
	int option_index = 0;	
	int ret=0;
	global_context_t * global_context;

	struct timeval xtime; 
	gettimeofday(&xtime,NULL);
	myrand_srand(time(NULL)^xtime.tv_usec);



	global_context = (global_context_t*)malloc(sizeof(global_context_t));
	init_global_context(global_context);
	global_context->config.entry_program_name = CORE_PROGRAM_SUBINDEL;

	optind = 1;
	opterr = 1;
	optopt = 63;


	while ((c = getopt_long (argc, argv, "pi:g:o:I:d:?v", long_options, &option_index)) != -1)
	{
		switch(c)
		{
			case 'i':
				strncpy(global_context->config.first_read_file, optarg, MAX_FILE_NAME_LENGTH-1);
			break;
			case 'p':
				global_context -> input_reads.is_paired_end_reads=1;
			break;
			case 'g':
				strncpy(global_context->config.index_prefix, optarg, MAX_FILE_NAME_LENGTH-1);
			break;
			case 'o':
				strncpy(global_context->config.output_prefix, optarg, MAX_FILE_NAME_LENGTH-1);
			break;
			case 'd':
				global_context->config.expected_pair_distance = atoi(optarg);
			break;
			case 'v':
				core_version_number("subindel");
				return 0;
			break;
			case 'I':
				global_context->config.max_indel_length = atoi(optarg);
				if(global_context->config.max_indel_length>=16)
				{
					global_context->config.reassembly_subread_length = 16;
					global_context->config.reassembly_window_multiplex = 3;
					global_context->config.reassembly_start_read_number = 3;
					global_context->config.reassembly_tolerable_voting = 0;
					global_context->config.reassembly_window_alleles = 3;
					global_context->config.reassembly_key_length = 28;
				}

			break;
			default:
			case 'h':
				print_usage_subindel();
				return -1;
			break;
		}
	}

	if(argc<3)
	{
		print_usage_subindel();
		return -1;
	}
	ret = print_configuration_forindel(global_context);

	warning_file_type(global_context->config.first_read_file,global_context->config.is_BAM_input?FILE_TYPE_BAM:FILE_TYPE_SAM);

	ret = ret || load_global_context_forindel(global_context);
	ret = ret || load_indexes_forindel(global_context);
	ret = ret || transfer_SAM_to_assembly(global_context);
	ret = ret || init_indel_tables(global_context);
	ret = ret || finalise_long_insertions_by_hashtable(global_context);
	ret = ret || write_indel_final_results(global_context);
	ret = ret || destroy_indel_module(global_context);
	ret = ret || print_summary(global_context);

	free(global_context);
	return ret;
}
