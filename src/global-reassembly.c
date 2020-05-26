/***************************************************************

   The Subread and Rsubread software packages are free
   software packages:
 
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
 


#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


#ifndef MAKE_STANDALONE
  #include <R.h>
#endif

#include <zlib.h>
#include <math.h>
#include <pthread.h>
#include <getopt.h>
#include "subread.h"
#include "sambam-file.h"
#include "gene-algorithms.h"
#include "input-files.h"
#include "long-hashtable.h"
#include "HelperFunctions.h"


static struct option GRA_long_options[] =
{
	{"in", required_argument, 0, 'i'},
	{"BAMinput", no_argument, 0, 'b'},
	{"out", required_argument, 0, 'o'},
	{"threads", required_argument, 0, 'T'},
	{"reportLength", required_argument, 0, 'L'},
	{"requiredVotes", required_argument, 0, 'V'},
	{"extensionVotes", required_argument, 0, 'v'},
	{"reverseUnmapped", no_argument, 0, 'R'},
	{"phred64", no_argument, 0, '6'},
	{"20bpMM", required_argument, 0, '2'},
	{"hugeMemory", no_argument, 0, 'H'},
	{"tmpDir",required_argument, 0, 't'},
	{"trimQuality", required_argument, 0, 'Q'},
	{"ignorePairedNames", no_argument, 0, 'C'},
	{0, 0, 0, 0}
};

#define GRA_FRAG_USED_BEFORE 1
#define GRA_FRAG_USED_NOW 2

typedef struct{

	short length_1;
	short replicate_count_1;
	//short sequencing_quality_1;
	short possible_extensions_to_read_1;

	short length_2;
	short replicate_count_2;
	//short sequencing_quality_2;
	short possible_extensions_to_read_2;
	short overlapping;

	unsigned int flags;

	unsigned long long list_file_pos;
	unsigned long long replicate_representitive_1;
	unsigned long long replicate_representitive_2;

} GRA_frag_properties_t;


#define GRA_DEBUG_ASSEMBLING 4

typedef struct {
	char input_file_name [MAX_FILE_NAME_LENGTH];
	char input_file_name2 [MAX_FILE_NAME_LENGTH];
	char output_file_name [MAX_FILE_NAME_LENGTH];
	int is_BAM_input;
	int total_threads;
	int shortest_read_allowed;
	int shortest_startread_allowed;
	int subread_extract_step;
	int read_trim_base_quality;
	int min_overlapping_bases_between_2ends;
	int min_extension_votes;
	int min_overlap_votes;
	int subread_length;
	int maximum_mismatch_in20bp;
	int reverse_unmapped_reads;
	int phred_offset;
	int check_read_names;
	unsigned int maximum_contig_length;
	unsigned int minimum_contig_length;

	unsigned long long debug_flags;

	// ====== Global Runtime ======
	int no_more_assembly;
	unsigned long long all_fragments;
	unsigned long long all_assemblies;
	unsigned long long marked_reported_fragments;
	unsigned long long start_read_pointer;
	unsigned long long all_possible_extensions;
	double start_time;
	int fd_fragproperties;
	FILE * fp_fragproperties;
	FILE * fp_fraglist;
	FILE * fp_output;
	GRA_frag_properties_t * frag_properties_array;
	lnhash_t assembling_index;

	// ====== Global System ======
	char tmp_file_prefix[MAX_FILE_NAME_LENGTH];
	unsigned int random_seeds[8];
	pid_t system_pid;

	// ====== Current Contig ======
	GRA_frag_properties_t * start_fragment_props;
	unsigned long long current_frag_no;
	int current_frag_reversed_to_start;
	int use_huge_memory;

	char * contig_str;
	unsigned int contig_len;
	unsigned int contig_buffer_size;

} GRA_global_context_t;

void GRA_print_usage()
{
	SUBREADputs("\nUsage:");
	SUBREADputs("\n  ./globalReassembly -i input.{SAM|BAM|FQ} -I input2.FQ -o output.txt\\\n                     {--BAMinput} {--reverseUnmapped} {--phred64} {--reportLength <int>}\n                    {--requiredVotes <int>} {--extensionVotes <int>} {--trimQuality <int>\n                    {--hugeMemory}}");
}

void GRA_init_extension_part(GRA_global_context_t * global_context)
{
	if(global_context->contig_str)
	{
		free(global_context->contig_str);
		global_context->contig_str = NULL;
	}
	global_context->contig_len = 0;
	global_context->contig_buffer_size = 0;

}

int GRA_add_new_extension_part(GRA_global_context_t * global_context, char * bases, int len, int add_after_tail)
{
	if(global_context -> maximum_contig_length <= global_context->contig_len+len)
		return -1;

	if(global_context->contig_str==NULL)
	{
		global_context -> contig_buffer_size = max(MAX_READ_LENGTH * 1.5, len+1);
		global_context -> contig_str = malloc(global_context -> contig_buffer_size);
		global_context -> contig_str[0]=0;
	}
	else if(global_context -> contig_buffer_size<global_context->contig_len+len+1)
	{
		global_context -> contig_buffer_size *= 2;
		global_context -> contig_str = realloc(global_context -> contig_str, global_context -> contig_buffer_size);
	}

	if(add_after_tail)
	{
		memcpy(global_context -> contig_str + global_context->contig_len, bases, len);
		*(global_context -> contig_str + global_context->contig_len+len)=0;
	}
	else
	{
		int xk1;
		for(xk1=global_context->contig_len; xk1>=0; xk1--) // include the last '0' AFTER the string
			*(global_context -> contig_str + xk1 + len) = *(global_context -> contig_str +xk1);
		memcpy(global_context -> contig_str, bases, len);
	}
	
	global_context -> contig_len += len;
	return 0;
}

unsigned int GR_subread2int(GRA_global_context_t * global_context, const char * key)
{
	unsigned int i, ret = 0;

	for (i=0; i< global_context -> subread_length ; i++)
	{
			ret = ret << 2;
			ret |= color2int (key[i]);
	}

	return ret;
}

int GRA_check_config (GRA_global_context_t * global_context)
{
	if(global_context -> input_file_name[0]==0 || global_context -> output_file_name[0]==0)
	{
		SUBREADputs("ERROR: You have to specify the input file and the output file.");
		return 1;
	}

	SUBREADputs("");
	SUBREADputs("");
	SUBREADputs(" ============== Global Reassembly of Reads ==============");
	SUBREADputs("");
	SUBREADprintf("\n\tInput file1\t:\t%s\n", global_context -> input_file_name);
	if(global_context -> input_file_name2[0])
		SUBREADprintf("\tInput file2\t:\t%s\n", global_context -> input_file_name2);
	SUBREADprintf("\tOutput file\t:\t%s\n", global_context -> output_file_name);
	if(global_context -> input_file_name2[0])
		SUBREADputs("\tFile Format\t:\tFASTQ");
	else{
		if(global_context -> is_BAM_input)
			SUBREADputs("\tFile Format\t:\tBAM");
		else
			SUBREADputs("\tFile Format\t:\tSAM");
	}

	SUBREADputs("");

	if(strcmp(global_context -> output_file_name, "//STDERR") == 0)
	{
		global_context -> fp_output = stderr;
	}
	if(strcmp(global_context -> output_file_name, "//STDOUT") == 0)
	{
		global_context -> fp_output = stdout;
	}
	else
	{
		global_context -> fp_output = fopen(global_context -> output_file_name, "w");
	}

	return 0;
}


void GRA_init_context(GRA_global_context_t * global_context)
{
	char mac_rand[13];
	FILE * system_random_fp; 
	memset(global_context, 0 , sizeof(GRA_global_context_t));

	global_context -> system_pid = getpid();
	system_random_fp = fopen("/dev/urandom","rb");
	assert(system_random_fp);
	fread(global_context -> random_seeds, 8 * sizeof(int), 1, system_random_fp);
	fclose(system_random_fp);

	mac_or_rand_str(mac_rand);
	sprintf(global_context -> tmp_file_prefix, "GRAtmp-%06d-%s", global_context -> system_pid, mac_rand);

	// best for simulation = 40
	global_context -> shortest_read_allowed = 40; 
	global_context -> shortest_startread_allowed = 60; 
	global_context -> min_overlapping_bases_between_2ends = 24;

	global_context -> total_threads = 1;

	//#warning "================ REMOVE ' + 1' FROM THE NEXT LINE ===================="
	global_context -> maximum_mismatch_in20bp = 0 + 1;
	global_context -> phred_offset = 33;
	global_context -> read_trim_base_quality = 30 + global_context -> phred_offset; 

	//#warning "================ REMOVE ' - 3' FROM THE NEXT LINE ===================="
	global_context -> subread_length = 16 - 3;
	global_context -> subread_extract_step = 5;


	//#warning "================ REMOVE ' - 2' FROM THE NEXT LINE ===================="
	global_context -> min_overlap_votes = 4 - 2;


	//#warning "================ REMOVE ' - 1' FROM THE NEXT LINE ===================="
	global_context -> min_extension_votes = 3 - 1;
	global_context -> contig_str = NULL;
	global_context -> start_time = miltime();
	global_context -> maximum_contig_length = 50000000;


	global_context -> minimum_contig_length = 251;
	global_context -> reverse_unmapped_reads = 0;
	global_context -> check_read_names = 1;



	global_context -> debug_flags = 0*GRA_DEBUG_ASSEMBLING;
}

int GRA_scan_best_overlap(GRA_global_context_t * global_context, char * piece, char * read, int rlen)
{
	int xk1;
	int best_loc = -1;
	for(xk1=0; xk1 < rlen - global_context -> min_overlapping_bases_between_2ends; xk1++)	
	{
		if(memcmp(piece, read + xk1, global_context -> min_overlapping_bases_between_2ends)==0)
		{
			if(best_loc<0)
				best_loc=xk1;
			else return -1;
		}
	}
	return best_loc;
}



// r1 is a zero-terminating string
// r2 is not a zero-terminating string.
// r1 and r2 are on opposite strands.
int calc_2r_overlapping(GRA_global_context_t * global_context, char * r1, char * r2, int r2len)
{

	int r1len = strlen(r1);
	if(r2len < global_context -> shortest_read_allowed || r1len < global_context -> shortest_read_allowed) return 0;
	char r2_rev[100];

	int r2_offset = r2len>100?r2len-100:0;
	memcpy(r2_rev, r2 + r2_offset, r2len - r2_offset);

	r2len = r2len - r2_offset;
	reverse_read(r2_rev, r2len, GENE_SPACE_BASE); 

	int scan_in_r2 = GRA_scan_best_overlap(global_context, r1 + r1len - global_context -> min_overlapping_bases_between_2ends, r2_rev, r2len);
	int scan_in_r1 = GRA_scan_best_overlap(global_context, r2_rev, r1, r1len);

	if(scan_in_r2 >0 && scan_in_r1 >0)
	{
		int overlapping_bases = scan_in_r2 +  global_context -> min_overlapping_bases_between_2ends;
		if(overlapping_bases == r1len - scan_in_r1)
		{
			return overlapping_bases;
		}
	}
	return 0;

}

#define SKIP_LINE {char tmpchr; while (1){tmpchr = fgetc(current_fp); if(tmpchr=='\n') break;}}

int GRA_trim_reads(GRA_global_context_t * global_context)
{
	int line_buff_len = MAX_READ_LENGTH * 2 + MAX_READ_NAME_LEN + MAX_CHROMOSOME_NAME_LEN * 3 + 500; 
	char tmp_fname[MAX_FILE_NAME_LENGTH],
		* read1_tmp = malloc(line_buff_len),
		* line_buff = malloc(line_buff_len),
		* in_sequence = malloc(MAX_READ_LENGTH),
		* in_quality_string = malloc(MAX_READ_LENGTH),
		* in_readname = malloc(MAX_READ_NAME_LEN);
	FILE * fp_fraglist, * fp_fragproperties;

	sprintf(tmp_fname, "%s.fragproperties", global_context -> tmp_file_prefix);
	fp_fragproperties = fopen(tmp_fname,"wb");

	sprintf(tmp_fname, "%s.fraglist", global_context -> tmp_file_prefix);
	fp_fraglist = fopen(tmp_fname,"wb");

	if(!(fp_fragproperties && fp_fraglist)){
		SUBREADprintf("Unable to create temporary files in the currect directory!\n");
		return 1;
	}


	SamBam_FILE * reader = NULL;
	FILE * fq1fp = NULL;
	FILE * fq2fp = NULL;

	if(global_context -> input_file_name2[0])
	{
		fq1fp = fopen(global_context -> input_file_name,"r");
		fq2fp = fopen(global_context -> input_file_name2,"r");
		if((!fq1fp)||(!fq2fp))
		{
			SUBREADprintf("Unable to open input files '%s' and '%s' as FASTQ files!\n", 
				global_context -> input_file_name,
				global_context -> input_file_name2
				);
			return 1;
		}

	}
	else
	{
		reader = SamBam_fopen(global_context -> input_file_name,
			global_context->is_BAM_input?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM);

		if(!reader)
		{
			SUBREADprintf("Unable to open input file '%s' as a %s file!\n", 
				global_context -> input_file_name,
				global_context->is_BAM_input?"BAM":"SAM"
				);
			return 1;
		}
	}
	unsigned long long short_frag_number = 0;
	unsigned long long frag_number = 0, overlapping_frags = 0;
	unsigned long readname_hash = 0;
	GRA_frag_properties_t frag_props;
	memset(&frag_props, 0, sizeof(frag_props));

	int  in_flags_1= -1, is_second_read = 0;
	while(1)
	{
		char * readret;
		char in_cigar[100];
		int xk1, in_flags = 0, in_mapq = 0, in_isrepeated = 0, in_rl = 0;
		unsigned int in_pos=0, in_tlen=0;
		char in_chro[MAX_CHROMOSOME_NAME_LEN];

		if(reader)
		{
			readret = SamBam_fgets(reader, line_buff, line_buff_len, 1); 
			if(!readret) break;
			if(line_buff[0]=='@') continue;

			parse_SAM_line(line_buff, in_readname, &in_flags, in_chro, &in_pos, in_cigar, &in_mapq, &in_tlen, in_sequence , in_quality_string, &in_rl, &in_isrepeated);
		}
		else
		{
			FILE * current_fp = is_second_read?fq2fp:fq1fp;
			readret = fgets(in_readname , MAX_READ_NAME_LEN-1, current_fp);
			if(readret)
				in_readname[strlen(in_readname)-1]=0;
			else	break;
			fgets(in_sequence , MAX_READ_LENGTH-1, current_fp);
			SKIP_LINE;
			fgets(in_quality_string , MAX_READ_LENGTH-1, current_fp);

			in_rl = strlen(in_sequence) - 1;
			in_sequence[in_rl]=0;
			in_quality_string[in_rl]=0;
			if(global_context -> phred_offset == 64)
			{
				for(xk1=0;xk1<in_rl;xk1++)in_quality_string[xk1]-=31;	 // phred 64 -> 33
			}

			in_flags = is_second_read?141:77;
			is_second_read = !is_second_read;
		}

		for(xk1=0; in_readname[xk1]; xk1++)
			if(in_readname[xk1]=='/'){in_readname[xk1]=0; break;
			}
		

		if(in_flags & SAM_FLAG_SECOND_READ_IN_PAIR)
		{
			if(global_context -> check_read_names)
				assert(readname_hash == fc_chro_hash(in_readname));
		}
		else{
			unsigned long new_readname_hash = fc_chro_hash(in_readname); 
			readname_hash = new_readname_hash;
		}


		int best_qual_window_value = 0, current_qual_sum = 0;
		short best_qual_window_to_head = -1;

		// xk1 is the last base being added into the window.
		for(xk1 = 0 ; xk1 < in_rl; xk1++)
		{
			current_qual_sum += in_quality_string[xk1];
			if(xk1>=10)
			{
				if(current_qual_sum > best_qual_window_value)
				{
					best_qual_window_value = current_qual_sum;
					best_qual_window_to_head = xk1 - 10;
				}
				current_qual_sum -= in_quality_string[xk1 - 10];
			}
		}

		char * new_read_start = in_sequence;
		char * new_qual_start = in_quality_string;
		for(xk1 = best_qual_window_to_head; xk1>=0; xk1--)
		{
			if(in_quality_string[xk1] < global_context -> read_trim_base_quality)
			{
				new_read_start = in_sequence + xk1;
				new_qual_start = in_quality_string + xk1;
				break;
			}
		}

		for(xk1 = best_qual_window_to_head; xk1< in_rl; xk1++)
		{
			if(in_quality_string[xk1] < global_context -> read_trim_base_quality)
			{
				in_sequence[xk1]=0;
				in_quality_string[xk1]=0;
				break;
			}
		}


		int new_rl = strlen(new_read_start);
		short min_quality_value = 9999;
		for(xk1 = 0; xk1<new_rl; xk1++)
		{
			min_quality_value = min(min_quality_value, new_qual_start[xk1]);
		}
		
		if((in_flags & SAM_FLAG_REVERSE_STRAND_MATCHED) || (global_context -> reverse_unmapped_reads && (in_flags & SAM_FLAG_UNMAPPED) && (in_flags & SAM_FLAG_MATE_UNMATCHED) && (in_flags & SAM_FLAG_SECOND_READ_IN_PAIR)))
		{
			reverse_read(new_read_start, new_rl, GENE_SPACE_BASE); 
			reverse_quality(new_qual_start, new_rl);
		}

		unsigned long long pos_fraglist = ftello(fp_fraglist);
		fwrite(new_read_start, new_rl, 1, fp_fraglist);
		fputc('\n', fp_fraglist);
		fwrite(new_qual_start, new_rl, 1, fp_fraglist);
		fputc('\n', fp_fraglist);

		if(in_flags & SAM_FLAG_SECOND_READ_IN_PAIR)
		{
			//frag_props . sequencing_quality_2 = min_quality_value;
			frag_props . length_2 = new_rl;
			frag_props . replicate_count_2 = 0;
			frag_props . flags = 0;
			frag_props . overlapping = calc_2r_overlapping(global_context, read1_tmp, new_read_start, new_rl);
			if(global_context -> min_overlapping_bases_between_2ends <= frag_props . overlapping )
				overlapping_frags ++;
			fwrite(&frag_props, sizeof(frag_props), 1, fp_fragproperties);
			if(frag_props . length_1 < global_context -> shortest_read_allowed || frag_props . length_2 < global_context -> shortest_read_allowed)short_frag_number++;

			memset(&frag_props, 0, sizeof(frag_props));
		}
		else
		{
			//frag_props . sequencing_quality_1 = min_quality_value;
			frag_props . length_1 = new_rl;
			frag_props . replicate_count_1 = 0;
			frag_props . list_file_pos = pos_fraglist;
			strcpy(read1_tmp, new_read_start);
			read1_tmp[new_rl] = 0;
			in_flags_1 = in_flags;
		}

		// the frag list is simply the concatination if SEQ_1, QUAL_1, SEQ_2, QUAL_2, ... 
	
		if(in_flags & SAM_FLAG_SECOND_READ_IN_PAIR) 
			frag_number ++;
	}
	global_context -> all_fragments = frag_number;

	SUBREADprintf("Finished read trimming. %llu fragments were processed; %llu fragments were too short.\nThere are %llu fragments satisfying the overlapping requirement.\n\n", 
		global_context -> all_fragments, short_frag_number, overlapping_frags);

	fclose(fp_fraglist);
	fclose(fp_fragproperties);

	free(in_readname);
	free(in_sequence);
	free(in_quality_string);
	free(line_buff);
	if(reader)
		SamBam_fclose(reader);
	if(fq1fp)fclose(fq1fp);
	if(fq2fp)fclose(fq2fp);

	assert(!(frag_props . length_1 || frag_props . length_2));

	return 0;
}

unsigned long long get_frag_number_part(unsigned long long q, int * is_second_read, int * is_reversed, int * mapped_pos)
{
	unsigned long long pos_part = q & 0x3fff;
	unsigned long long info_part = q >> 14;
	if(pos_part > 0x2000)
	{
		pos_part -= 0x4000;
		info_part ++;
	}
	if(mapped_pos)
	{
		(*is_second_read) = (info_part & 0x2)?1:0;
		(*is_reversed) = (info_part & 0x1)?1:0;
		(*mapped_pos) = pos_part;
	}
	return info_part >> 2;
}

void get_frag_sequence(GRA_global_context_t* global_context, GRA_frag_properties_t * prop, int is_read2, char * seq_buffer)
{
	unsigned long long fpos = prop -> list_file_pos;
	if(is_read2) fpos += prop -> length_1 *2 + 2;
	int seq_len = is_read2 ? prop -> length_2:prop -> length_1;
	unsigned long long cpos = ftello(global_context -> fp_fraglist);

	//SUBREADprintf("FPOS=%llu, CPOS=%llu\n", fpos, cpos);
	fseeko(global_context -> fp_fraglist , fpos, SEEK_SET);
	fread(seq_buffer, seq_len,1 , global_context -> fp_fraglist);
	seq_buffer[seq_len]=0;
	fseeko(global_context -> fp_fraglist , cpos, SEEK_SET);

	assert(!strchr(seq_buffer,'|'));

	//if(is_read2)SUBREADprintf("QFSR=%s\n", seq_buffer);
	//else SUBREADprintf("QFSA=%s\n", seq_buffer);
}

void search_read_extension_number(GRA_global_context_t* global_context, lnhash_vote_record_t* vote_rec, unsigned long long frag_number, GRA_frag_properties_t * query_prop , int is_second_read)
{
	if(get_frag_number_part(vote_rec -> head_position, NULL, NULL, NULL) == frag_number) return;

	int mapp_target = (vote_rec -> head_position & 0x3fff);

	if(vote_rec -> coverage_start < global_context -> subread_extract_step && mapp_target >0 && mapp_target < 1300)
	{
		int next_is_second=0, next_is_reversed=0, next_pos=0;
		GRA_frag_properties_t * frag_prop = global_context -> frag_properties_array + get_frag_number_part(vote_rec -> head_position, & next_is_second, &next_is_reversed, &next_pos); 
		int target_rlen = next_is_second?frag_prop -> length_2:frag_prop -> length_1;
		int rlen = is_second_read?query_prop -> length_2:query_prop -> length_1;


		int should_have_votes = target_rlen - vote_rec -> coverage_start ;
		should_have_votes = min(rlen, should_have_votes);
		should_have_votes = (should_have_votes-15 - global_context -> subread_extract_step) / global_context -> subread_extract_step + 1;
		if(vote_rec -> num_votes >= should_have_votes && vote_rec -> num_votes >= global_context -> min_extension_votes + 2)
		{
			if(is_second_read)
				query_prop -> possible_extensions_to_read_2 ++;
			else
				query_prop -> possible_extensions_to_read_1 ++;

			global_context->all_possible_extensions++;
		}
	}
}

int search_replicate_representitive(GRA_global_context_t* global_context, lnhash_vote_record_t * vote_rec, unsigned long long frag_number, GRA_frag_properties_t * query_prop, int query_is_read2, unsigned long long * representitiva_id, int * target_is_read2)
{
	if(get_frag_number_part(vote_rec -> head_position, NULL, NULL, NULL) == frag_number) return 0;

	int rlen = query_is_read2?query_prop -> length_2:query_prop -> length_1;
	int mapp_target = (vote_rec -> head_position & 0x3fff);
	if(mapp_target<1300 && vote_rec -> coverage_start < global_context -> subread_extract_step && vote_rec -> coverage_end > rlen - global_context -> subread_extract_step - 1 && (vote_rec -> coverage_end - 15 - vote_rec -> coverage_start - global_context -> subread_extract_step) / global_context -> subread_extract_step + 1 <= vote_rec -> num_votes)
	{
		int next_is_second=0, next_is_reversed=0, next_pos=0;
		GRA_frag_properties_t * frag_prop = global_context -> frag_properties_array + get_frag_number_part(vote_rec -> head_position, & next_is_second, &next_is_reversed, &next_pos); 
		int target_rlen = next_is_second?frag_prop -> length_2:frag_prop -> length_1;

		if(target_rlen >= rlen + mapp_target &&  vote_rec -> coverage_end - vote_rec -> coverage_start + mapp_target <= target_rlen)
		{
			/*
			char read_seq_query[MAX_READ_LENGTH], read_seq_target[MAX_READ_LENGTH];
			int target_is_reversed = (vote_rec -> head_position & 0x4000)?1:0, xk1, mm=0;

			get_frag_sequence(global_context, frag_prop,  next_is_second, read_seq_target);
			get_frag_sequence(global_context, query_prop, query_is_read2 , read_seq_query );

			if(target_is_reversed) reverse_read(read_seq_target , target_rlen, GENE_SPACE_BASE);

			for(xk1=0; read_seq_query[xk1]; xk1++)
				mm+= (read_seq_target[xk1 + mapp_target] == 'N' || read_seq_target[xk1 + mapp_target] != read_seq_query[xk1]);
			*/
			int mm = 0;
			if(!mm)
			{
				(*representitiva_id) = get_frag_number_part(vote_rec -> head_position, NULL, NULL, NULL);
				(*target_is_read2) = next_is_second;

				if(rlen == target_rlen)
					return (*representitiva_id) < frag_number;
				return 1;
			}
		}
	}
	return 0;
}

int GRA_build_read_index(GRA_global_context_t* global_context)
{
	unsigned long long report_step = max(1, global_context -> all_fragments/20);
	FILE * fp_fragproperties, * fp_fraglist;
	int is_seed_frags;
	GRA_frag_properties_t frag_prop;

	lnhash_create(&global_context -> assembling_index, 666661);

	char tmp_fname[MAX_FILE_NAME_LENGTH],
		* in_sequence_1 = malloc(MAX_READ_LENGTH),
		* in_quality_string_1 = malloc(MAX_READ_LENGTH),
		* in_sequence_2 = malloc(MAX_READ_LENGTH),
		* in_quality_string_2 = malloc(MAX_READ_LENGTH);


	lnhash_t first_index;
	lnhash_create(&first_index, 666661);

	sprintf(tmp_fname, "%s.fragproperties", global_context -> tmp_file_prefix);
	fp_fragproperties = fopen(tmp_fname,"rb");
	global_context -> fp_fragproperties = fp_fragproperties;

	sprintf(tmp_fname, "%s.fraglist", global_context -> tmp_file_prefix);
	fp_fraglist = fopen(tmp_fname,"rb");
	global_context -> fp_fraglist = fp_fraglist;

	for(is_seed_frags = 1; is_seed_frags>=0; is_seed_frags--)
	{
		unsigned long long frag_number = 0;
		fseek(fp_fragproperties, 0, SEEK_SET);

		while(1)
		{
			int nrecord = fread(&frag_prop, sizeof(frag_prop), 1, fp_fragproperties);
			if(nrecord<1) break;


			int is_good_seed = 0;
			if(frag_prop.length_1 >=  global_context -> shortest_startread_allowed && frag_prop.length_2 >=  global_context -> shortest_startread_allowed )
				is_good_seed = 1;
			if(frag_prop.replicate_representitive_1 !=0 && frag_prop.replicate_representitive_2 != 0)
				is_good_seed = 0;
			if(frag_prop.overlapping < global_context -> min_overlapping_bases_between_2ends)
				is_good_seed = 0;
	
			if(is_good_seed == is_seed_frags && frag_prop.length_1 >= global_context -> shortest_read_allowed && frag_prop.length_2 >= global_context -> shortest_read_allowed)
			{

				fseeko(fp_fraglist, frag_prop.list_file_pos, SEEK_SET);
				fread(in_sequence_1, frag_prop.length_1 , 1, fp_fraglist);
				fgetc(fp_fraglist);
				fread(in_quality_string_1, frag_prop.length_1 , 1, fp_fraglist);
				fgetc(fp_fraglist);
				fread(in_sequence_2, frag_prop.length_2 , 1, fp_fraglist);
				fgetc(fp_fraglist);
				fread(in_quality_string_2, frag_prop.length_2 , 1, fp_fraglist);
				fgetc(fp_fraglist);

				int is_second_read, xxi, current_len;
				char * current_seq;
				for(xxi = 0; xxi<2; xxi++)
				{
					is_second_read = xxi;

					current_seq = is_second_read?in_sequence_2:in_sequence_1;
					current_len = is_second_read?frag_prop.length_2:frag_prop.length_1;
					unsigned long long read_pos_code = (frag_number<<16) | (is_second_read << 15);
					int extract_cursor;
					int current_repeated = frag_prop.replicate_representitive_2 && frag_prop.replicate_representitive_1;

					if(!current_repeated)
					{
						for(extract_cursor=0; extract_cursor< current_len-15; extract_cursor += global_context -> subread_extract_step)
						{

							char * subread_txt = current_seq + extract_cursor;
							unsigned int subread_int = GR_subread2int(global_context, subread_txt); 

							lnhash_insert(&global_context -> assembling_index, subread_int, read_pos_code | extract_cursor);
						}			
					}
				}

			}

			if(frag_number % report_step == 0)
				SUBREADprintf("  Building reassembly index, %.1f%% finished in %.1f minutes.\n", frag_number*100./global_context -> all_fragments, (miltime() - global_context -> start_time)/60);

			frag_number++;
		}
	}
	fclose(fp_fragproperties);
	free(in_sequence_1);
	free(in_sequence_2);
	free(in_quality_string_1);
	free(in_quality_string_2);


	lnhash_resort(&global_context -> assembling_index);
	SUBREADputs("\nThe reassembly index is built.\n");

	return 0;
}

int GRA_mark_read_extension(GRA_global_context_t* global_context)
{
	unsigned long long report_step = max(1, global_context -> all_fragments/20);
	FILE  * fp_fraglist;
	GRA_frag_properties_t * frag_prop;
	char tmp_fname[MAX_FILE_NAME_LENGTH],
		* in_sequence_1 = malloc(MAX_READ_LENGTH),
		* in_quality_string_1 = malloc(MAX_READ_LENGTH),
		* in_sequence_2 = malloc(MAX_READ_LENGTH),
		* in_quality_string_2 = malloc(MAX_READ_LENGTH);


	lnhash_t first_index;
	lnhash_create(&first_index, 666661);

	sprintf(tmp_fname, "%s.fraglist", global_context -> tmp_file_prefix);
	fp_fraglist = fopen(tmp_fname,"rb");
	global_context -> fp_fraglist = fp_fraglist;

	unsigned long long frag_number = 0;

	while(1)
	{
		if(frag_number >=  global_context -> all_fragments) break;
		frag_prop = global_context -> frag_properties_array + frag_number;

		if(frag_prop -> length_1 >= global_context -> shortest_read_allowed && frag_prop -> length_2 >= global_context -> shortest_read_allowed)
		{

			fseeko(fp_fraglist, frag_prop -> list_file_pos, SEEK_SET);
			fread(in_sequence_1, frag_prop -> length_1 , 1, fp_fraglist);
			fgetc(fp_fraglist);
			fread(in_quality_string_1, frag_prop -> length_1 , 1, fp_fraglist);
			fgetc(fp_fraglist);
			fread(in_sequence_2, frag_prop -> length_2 , 1, fp_fraglist);
			fgetc(fp_fraglist);
			fread(in_quality_string_2, frag_prop -> length_2 , 1, fp_fraglist);
			fgetc(fp_fraglist);

			int is_second_read, xxi, current_len;
			char * current_seq;
			for(xxi = 0; xxi<2; xxi++)
			{
				is_second_read = xxi;

				current_seq = is_second_read?in_sequence_2:in_sequence_1;
				current_len = is_second_read?frag_prop -> length_2:frag_prop -> length_1;
				unsigned long long read_pos_code = (frag_number<<16) | (is_second_read << 15);
				int extract_cursor;

				for(extract_cursor=0; extract_cursor< current_len-15; extract_cursor += global_context -> subread_extract_step)
				{
					char * subread_txt = current_seq + extract_cursor;
					unsigned int subread_int = GR_subread2int(global_context, subread_txt); 

					lnhash_insert(&first_index, subread_int, read_pos_code | extract_cursor);
				}
			}
		}
		if(frag_number % report_step == 0)
			SUBREADprintf("  Building first index, %.1f%% finished in %.1f minutes.\n", frag_number*100./global_context -> all_fragments, (miltime() - global_context -> start_time)/60);
		frag_number ++;

	}

	SUBREADprintf("\nThe index has been built.\n\n");

	fseeko(fp_fraglist,0, SEEK_SET);
	frag_number = 0;

	lnhash_resort(&first_index);

	lnhash_vote_t * votes = malloc(sizeof(lnhash_vote_t));
	unsigned long long all_rep_found = 0;


	while(1)
	{
		if(frag_number >=  global_context -> all_fragments) break;
		frag_prop = global_context -> frag_properties_array + frag_number;

		if(frag_prop -> length_1 >= global_context -> shortest_read_allowed && frag_prop -> length_2 >= global_context -> shortest_read_allowed)
		{
			fseeko(fp_fraglist, frag_prop -> list_file_pos, SEEK_SET);
			fread(in_sequence_1, frag_prop -> length_1 , 1, fp_fraglist);
			fgetc(fp_fraglist);
			fread(in_quality_string_1, frag_prop -> length_1 , 1, fp_fraglist);
			fgetc(fp_fraglist);
			fread(in_sequence_2, frag_prop -> length_2 , 1, fp_fraglist);
			fgetc(fp_fraglist);
			fread(in_quality_string_2, frag_prop -> length_2 , 1, fp_fraglist);
			fgetc(fp_fraglist);

			int is_second_read, xxi, current_len;
			char * current_seq;
			for(xxi = 0; xxi<2; xxi++)
			{
				lnhash_init_votes(votes);
				is_second_read = xxi;

				current_seq = is_second_read?in_sequence_2:in_sequence_1;
				current_len = is_second_read?frag_prop -> length_2:frag_prop -> length_1;
				int extract_cursor;

				for(extract_cursor=0; extract_cursor< current_len-15; extract_cursor ++)
				{
					char * subread_txt = current_seq + extract_cursor;
					unsigned int subread_int = GR_subread2int(global_context, subread_txt); 

					lnhash_query(&first_index, votes, subread_int, extract_cursor);
				}

				int vote_X, vote_Y;

				for(vote_Y=0; vote_Y<LNHASH_VOTE_ARRAY_HEIGHT; vote_Y++)
					for(vote_X=0; vote_X<votes -> vote_record_items[vote_Y]; vote_X++)
					{
						//unsigned long long representitiva_id = 0;
						//int target_is_second_read = 0;
						lnhash_vote_record_t * vote_1_rec = &(votes -> vote_records[vote_Y][vote_X]);
						if(vote_1_rec -> num_votes < global_context -> min_extension_votes)continue;
						//int is_rep_found = search_replicate_representitive(global_context, vote_1_rec, frag_number, frag_prop , is_second_read , &representitiva_id, &target_is_second_read);
						/*if(is_rep_found)
						{
							//if(is_second_read) 	frag_prop -> replicate_representitive_2 = representitiva_id + 1;
							//else			frag_prop -> replicate_representitive_1 = representitiva_id + 1;

							GRA_frag_properties_t * target_prop = global_context -> frag_properties_array + representitiva_id;	

							if(target_is_second_read)	target_prop -> replicate_count_2 = min(65500,  target_prop -> replicate_count_2+1);
							else				target_prop -> replicate_count_1 = min(65500,  target_prop -> replicate_count_1+1);
							all_rep_found +=1;
						}*/

						search_read_extension_number(global_context, vote_1_rec, frag_number, frag_prop , is_second_read);
					}
			}
		}

		//if(frag_prop.replicate_representitive_2||frag_prop.replicate_representitive_1)
		//	SUBREADprintf("REP=%llu  ,  %llu\n", frag_prop.replicate_representitive_1, frag_prop.replicate_representitive_2);

		if(frag_number % report_step == 0)
		{
			SUBREADprintf("  Searching for extensions, %.1f%% finished in %.1f minutes. %llu extensions were found.\n", frag_number*100./global_context -> all_fragments, (miltime() - global_context -> start_time)/60, global_context -> all_possible_extensions);
		}
		frag_number ++;
	}
	free(votes);

	SUBREADprintf("\n%llu replicates were marked.\n\n", all_rep_found);

	fclose(fp_fraglist);
	free(in_sequence_1);
	free(in_sequence_2);
	free(in_quality_string_1);
	free(in_quality_string_2);

	lnhash_destroy(&first_index);


	return 0;

} 

void GRA_add_output_contig(GRA_global_context_t* global_context)
{
	int xk1, is_rev;

	lnhash_vote_t * votes = malloc(sizeof(lnhash_vote_t));

	//SUBREADprintf("===========\n");

	lnhash_init_votes(votes);
	for(is_rev = 0; is_rev < 2; is_rev++)
	{

		
		for(xk1 = 0; xk1 < global_context -> contig_len - 15; xk1+=1)
		{
			char *subread_txt = global_context -> contig_str + xk1;
			unsigned int subread_int = GR_subread2int(global_context, subread_txt); 
			lnhash_query(&global_context -> assembling_index, votes, subread_int, xk1);
		}
		if(0==is_rev) reverse_read(global_context -> contig_str, global_context -> contig_len, GENE_SPACE_BASE);

	}

	int examined_vote_item = 0;

	lnhash_vote_record_t * sorted_votes = NULL;
	int used_votes_items = sorted_voting_table_offset(votes, &sorted_votes,5);

	long long int last_frag_no = -1;
	int curr_frag_first_rec_index = -1;
	
	for(examined_vote_item = 0; ; examined_vote_item++)
	{
		lnhash_vote_record_t * rec = sorted_votes + examined_vote_item;
		
		unsigned long long this_frag_no = 0;
		if(examined_vote_item < used_votes_items)
			this_frag_no = get_frag_number_part(rec -> head_position, NULL, NULL, NULL);

		if(this_frag_no != last_frag_no || examined_vote_item >= used_votes_items)
		{

			if(last_frag_no >=0)
			{
				int x1, x2;
				GRA_frag_properties_t * target_prop = global_context -> frag_properties_array + last_frag_no;
				int is_marked = 0;
				for(x1 = curr_frag_first_rec_index; x1 < examined_vote_item; x1++)
				{
					int map_to_read2=0, map_to_reversed=0, map_to_pos=0;
					get_frag_number_part(sorted_votes[x1].head_position, &map_to_read2 , &map_to_reversed, &map_to_pos);

					int target_1_votes = (target_prop -> length_1 - 15 - global_context -> subread_extract_step) / global_context -> subread_extract_step; 
					for(x2 = x1 + 1; x2 < examined_vote_item; x2++)
					{

						int  map_to_read22=0, map_to_reversed2=0, map_to_pos2=0;
						get_frag_number_part(sorted_votes[x2].head_position, &map_to_read22 , &map_to_reversed2, &map_to_pos2);
						int target_2_votes = (target_prop -> length_2 - 15 - global_context -> subread_extract_step) / global_context -> subread_extract_step; 

						if(map_to_read22 != map_to_read2)
						{
							int total_curr_votes = sorted_votes[x2].num_votes + sorted_votes[x1].num_votes;
							if(total_curr_votes >= target_1_votes + target_2_votes + 1)
							{
								if(0==(target_prop -> flags & GRA_FRAG_USED_BEFORE))
								{
									global_context -> marked_reported_fragments ++;
									target_prop -> flags |= GRA_FRAG_USED_BEFORE; 
									is_marked = 1;
								}
							}
						}
						if(is_marked)break;
					}
					if(is_marked)break;
				}

			}
			curr_frag_first_rec_index = examined_vote_item;
			last_frag_no  = this_frag_no;
		}

		if(examined_vote_item >= used_votes_items)break;
	}
	free(sorted_votes);
	/*
	for(vote_Y = 0; vote_Y < LNHASH_VOTE_ARRAY_HEIGHT ; vote_Y++) 
		for(vote_X = 0; vote_X < votes->vote_record_items[vote_Y] ; vote_X++) 
		{
			int vote2_X, vote2_Y, map_to_read2=0, map_to_reversed=0, map_to_pos=0;
			lnhash_vote_record_t * voterec_1 = &(votes -> vote_records[vote_Y][vote_X]);


		//	SUBREADprintf("[%s] Votes=%d; HPos=%llX; CovStart=%d; OffsetHpos=%llX\n", is_rev?"REV":"STR", voterec_1->num_votes , voterec_1 -> head_position  , voterec_1 -> coverage_start , voterec_1 -> head_position + voterec_1 -> coverage_start);

			if(voterec_1->num_votes < 5) continue; // ||  voterec_1 -> coverage_start <=  global_context -> subread_extract_step) continue;
			unsigned long long next_frag_no = get_frag_number_part(voterec_1 -> head_position + voterec_1 -> coverage_start, &map_to_read2 , &map_to_reversed, &map_to_pos);

			for(vote2_Y = 0; vote2_Y < LNHASH_VOTE_ARRAY_HEIGHT ; vote2_Y++) 
				for(vote2_X = 0; vote2_X < votes->vote_record_items[vote2_Y] ; vote2_X++) 
				{
					int  map_to_read22=0, map_to_reversed2=0, map_to_pos2=0;
					lnhash_vote_record_t * voterec_2 = &(votes -> vote_records[vote2_Y][vote2_X]);	
					if(voterec_2 == voterec_1 || voterec_2->num_votes < 5) continue; //  ||  voterec_2 -> coverage_start <=  global_context -> subread_extract_step) continue;
					unsigned long long next_frag_no2 = get_frag_number_part(voterec_2 -> head_position + voterec_2 -> coverage_start, &map_to_read22 , &map_to_reversed2, &map_to_pos2);
					if(next_frag_no2 == next_frag_no && map_to_read2 != map_to_read22)
					{
						GRA_frag_properties_t * target_prop = global_context -> frag_properties_array + next_frag_no2;
						int target_1_votes = (target_prop -> length_1 - 15 - global_context -> subread_extract_step) / global_context -> subread_extract_step; 
						int target_2_votes = (target_prop -> length_2 - 15 - global_context -> subread_extract_step) / global_context -> subread_extract_step; 
						int contig_1_votes = map_to_read2?voterec_2->num_votes:voterec_1 -> num_votes;
						int contig_2_votes = map_to_read2?voterec_1->num_votes:voterec_2 -> num_votes;

						//SUBREADprintf("FINE_EXAMINE: FRAG %llX\n ", next_frag_no);
						if(contig_1_votes + contig_2_votes >= target_1_votes + target_2_votes + 1)
						{
						//	SUBREADprintf("Marked read [%s]: %d>=%d (len %d)  ;  %d>=%d (len %d)\n", is_rev?"REV":"STR", contig_1_votes, target_1_votes - 1, target_prop -> length_1 , contig_2_votes, target_2_votes - 1, target_prop -> length_2);
							if(0==(target_prop -> flags & GRA_FRAG_USED_BEFORE))
							{
								global_context -> marked_reported_fragments ++;
								target_prop -> flags |= GRA_FRAG_USED_BEFORE; 
							}
						}
					}
				}

			examined_vote_item++;
			SUBREADprintf("Examined %d votes\n", examined_vote_item);
		}
	*/
	free(votes);
}


int GRA_write_current_assembly(GRA_global_context_t* global_context)
{
	if(global_context -> minimum_contig_length <= global_context -> contig_len)
	{
		fprintf(global_context -> fp_output, ">SEQ%llu_LEN%d\n%s\n", global_context -> all_assemblies++ , global_context -> contig_len, global_context -> contig_str);
		fflush(global_context -> fp_output);
		GRA_add_output_contig(global_context);
	}
	free(global_context -> contig_str);
	global_context -> contig_len = 0;
	global_context -> contig_buffer_size = 0;
	global_context -> contig_str = NULL;
	return 0;
}



typedef struct
{
	unsigned long long frag_number;

	int map_1_pos;
	int map_2_pos;

	int read_1_mapped_targ_2;
	int read_2_mapped_targ_2;

	int read_1_votes;
	int read_2_votes;

	int read_1_coverage_start;
	int read_1_coverage_end;

	int read_2_coverage_start;
	int read_2_coverage_end;
} GRA_votes_pairs;


int GRA_reassemble_reads_LR(GRA_global_context_t* global_context, int to_right)
{
	char * seq_buf_1 = malloc(MAX_READ_LENGTH);
	char * seq_buf_2 = malloc(MAX_READ_LENGTH);
	char * target_seq_buf = malloc(MAX_READ_LENGTH);
	char * best_bases_added = malloc(MAX_READ_LENGTH);
	lnhash_vote_t * votes1 = malloc(sizeof(lnhash_vote_t));
	lnhash_vote_t * votes2 = malloc(sizeof(lnhash_vote_t));
	GRA_votes_pairs * examine_vote_pairs = malloc(sizeof(GRA_votes_pairs) * LNHASH_VOTE_ARRAY_WIDTH * LNHASH_VOTE_ARRAY_HEIGHT);
	
	while(1)
	{
		GRA_frag_properties_t * current_prop = global_context -> frag_properties_array + global_context -> current_frag_no;
		int is_current_read2_extension = (global_context -> current_frag_reversed_to_start ^ to_right);

		int is_second_read, xxi, current_len;
		for(xxi = 0; xxi<2; xxi++)
		{
			is_second_read = xxi;
			char * seq_buf = is_second_read?seq_buf_2:seq_buf_1;
			lnhash_vote_t * votes = is_second_read?votes2:votes1;
			lnhash_init_votes(votes);

			get_frag_sequence(global_context, current_prop, is_second_read, seq_buf);

			current_len = is_second_read?current_prop -> length_2:current_prop -> length_1;
			int extract_cursor;

			for(extract_cursor=0; extract_cursor< current_len-15; extract_cursor ++)
			{
				char * subread_txt = seq_buf + extract_cursor;
				unsigned int subread_int = GR_subread2int(global_context, subread_txt); 

				lnhash_query(&global_context -> assembling_index, votes, subread_int, extract_cursor);
			}
		}


		// examine v1 and v2.
		//
		int examine_vote_pair_count = 0;
		memset(examine_vote_pairs, 0, sizeof(GRA_votes_pairs) * LNHASH_VOTE_ARRAY_WIDTH * LNHASH_VOTE_ARRAY_HEIGHT);

		lnhash_vote_record_t * vote1_sorted, * vote2_sorted;

		int used_votes_1 = sorted_voting_table(votes1, &vote1_sorted, global_context -> min_overlap_votes);
		int used_votes_2 = sorted_voting_table(votes2, &vote2_sorted, global_context -> min_overlap_votes);

		int xxj = 0, test_pair;
		int vote_X, vote_Y;
		long long int last_frag_no = -1;

		int xxi_max_votes = -1;
		int xxi_max_index = -1;
		int xxj_max_votes = -1;
		int xxj_max_index = -1;

		if(0)
		for(xxi = 0; xxi < 2; xxi++){
			
			lnhash_vote_record_t * cur_vote_sorted  =xxi?vote2_sorted:vote1_sorted;
			int cur_used_votes = xxi?used_votes_2:used_votes_1;
			for(xxj = 0; xxj < cur_used_votes; xxj++)
			{
				long long int current_frag_no = get_frag_number_part(cur_vote_sorted[xxj].head_position, NULL, NULL, NULL);
				SUBREADprintf("C=%llu;;  %d : FRAG=%llu RAW=%llu V=%d\n", global_context -> current_frag_no, xxi, current_frag_no, cur_vote_sorted[xxj].head_position, cur_vote_sorted[xxj].num_votes);
			}

			lnhash_vote_t * votes = is_second_read?votes2:votes1;
			for(vote_Y = 0; vote_Y < LNHASH_VOTE_ARRAY_HEIGHT ; vote_Y++) 
				for(vote_X = 0; vote_X < votes->vote_record_items[vote_Y] ; vote_X++) 
				{
					lnhash_vote_record_t * rec = &(votes->vote_records[vote_Y][vote_X]);
					long long int current_frag_no = get_frag_number_part(rec -> head_position, NULL, NULL, NULL);
					SUBREADprintf("R=%llu;;  %d : FRAG=%llu RAW=%llu V=%d\n", global_context -> current_frag_no, xxi, current_frag_no, rec -> head_position, rec->num_votes);
				}
		}

		xxj = 0;

		for(xxi = 0; xxi < used_votes_1; xxi++)
		{
			long long int current_frag_no = get_frag_number_part(vote1_sorted[xxi].head_position, NULL, NULL, NULL);
			if(current_frag_no != last_frag_no || xxi == used_votes_1 - 1)
			{
				if(last_frag_no>=0){
					while(xxj < used_votes_2){
						long long int current_grag_xxj = get_frag_number_part(vote2_sorted[xxj].head_position, NULL, NULL, NULL);
						if(current_grag_xxj > last_frag_no) break;
						if(current_grag_xxj == last_frag_no) 
						{
							//SUBREADprintf("CHIT! C = %llu\n", current_grag_xxj);
							if(xxj_max_votes < vote2_sorted[xxj].num_votes)
							{
								xxj_max_votes = vote2_sorted[xxj].num_votes;
								xxj_max_index = xxj;
							}
						}
						xxj ++;
					}

					
					if(xxi_max_votes>0 && xxj_max_votes>0 && last_frag_no != global_context -> current_frag_no){
						
						int map_to_read2 = 0, map_to_reversed = 0, map_to_pos = 0;

						examine_vote_pairs[examine_vote_pair_count].frag_number = last_frag_no;

						get_frag_number_part(vote1_sorted[xxi_max_index].head_position,&map_to_read2 , &map_to_reversed, &map_to_pos);

						examine_vote_pairs[examine_vote_pair_count].map_1_pos		  = map_to_pos;
						examine_vote_pairs[examine_vote_pair_count].read_1_votes 	  = xxi_max_votes;
						examine_vote_pairs[examine_vote_pair_count].read_1_mapped_targ_2  = map_to_read2;
						examine_vote_pairs[examine_vote_pair_count].read_1_coverage_start = vote1_sorted[xxi_max_index].coverage_start;
						examine_vote_pairs[examine_vote_pair_count].read_1_coverage_end	  = vote1_sorted[xxi_max_index].coverage_end;
		

						get_frag_number_part(vote2_sorted[xxj_max_index].head_position,&map_to_read2 , &map_to_reversed, &map_to_pos);

						examine_vote_pairs[examine_vote_pair_count].map_2_pos		  = map_to_pos;
						examine_vote_pairs[examine_vote_pair_count].read_2_votes 	  = xxj_max_votes;
						examine_vote_pairs[examine_vote_pair_count].read_2_mapped_targ_2  = map_to_read2;
						examine_vote_pairs[examine_vote_pair_count].read_2_coverage_start = vote2_sorted[xxj_max_index].coverage_start;
						examine_vote_pairs[examine_vote_pair_count].read_2_coverage_end	  = vote2_sorted[xxj_max_index].coverage_end;
						if(map_to_read2 != examine_vote_pairs[examine_vote_pair_count].read_1_mapped_targ_2)
							examine_vote_pair_count ++;
						else	memset(examine_vote_pairs + examine_vote_pair_count, 0, sizeof(GRA_votes_pairs));
					}
				}

				last_frag_no = current_frag_no;
				xxi_max_votes = -1;
				xxi_max_index = -1;
				xxj_max_votes = -1;
				xxj_max_index = -1;
			}

			if(xxi_max_votes < vote1_sorted[xxi].num_votes)
			{
				xxi_max_votes = vote1_sorted[xxi].num_votes;
				xxi_max_index = xxi;
			}


		}


		/*
		for(is_second_read = 0; is_second_read<2; is_second_read++)
		{
			lnhash_vote_t * votes = is_second_read?votes2:votes1;
			for(vote_Y = 0; vote_Y < LNHASH_VOTE_ARRAY_HEIGHT ; vote_Y++) 
				for(vote_X = 0; vote_X < votes->vote_record_items[vote_Y] ; vote_X++) 
				{
					lnhash_vote_record_t * voterec = &(votes -> vote_records[vote_Y][vote_X]);
					int map_to_read2 = 0, map_to_reversed = 0, map_to_pos = 0;
					unsigned long long next_frag_no = get_frag_number_part(voterec -> head_position,&map_to_read2 , &map_to_reversed, &map_to_pos);

					//SUBREADprintf("VV[%llu]: next_frag_no=%llX; frag=%llu; V=%d; is_second_read=%d; targ_second=%d; targ_rev=%d\n", global_context -> current_frag_no, voterec -> head_position, next_frag_no , voterec -> num_votes , is_second_read , map_to_read2, map_to_reversed);

					if(global_context -> current_frag_no == next_frag_no)continue;
					if(voterec -> num_votes < global_context -> min_overlap_votes) continue;
					if(map_to_reversed) continue;

					if(is_second_read)
					{
						for(test_pair =0; test_pair < examine_vote_pair_count; test_pair++)
						{
							if(examine_vote_pairs [test_pair].frag_number == next_frag_no && 
								voterec -> num_votes > examine_vote_pairs[test_pair].read_2_votes)
							{
								examine_vote_pairs[test_pair].map_2_pos = map_to_pos;
								examine_vote_pairs[test_pair].read_2_votes = voterec -> num_votes;
								examine_vote_pairs[test_pair].read_2_mapped_targ_2 = map_to_read2;
								examine_vote_pairs[test_pair].read_2_coverage_start = voterec -> coverage_start;
								examine_vote_pairs[test_pair].read_2_coverage_end = voterec -> coverage_end;
							}
						}
					}
					else
					{
						int is_found = 0;

						for(test_pair =0; test_pair < examine_vote_pair_count; test_pair++)
						{
							if(examine_vote_pairs [test_pair].frag_number == next_frag_no) 
							{
								if(voterec -> num_votes > examine_vote_pairs[test_pair].read_1_votes)
								{
									examine_vote_pairs[test_pair].map_1_pos = map_to_pos;
									examine_vote_pairs[test_pair].read_1_votes = voterec -> num_votes;
									examine_vote_pairs[test_pair].read_1_mapped_targ_2 = map_to_read2;
									examine_vote_pairs[test_pair].read_1_coverage_start = voterec -> coverage_start;
									examine_vote_pairs[test_pair].read_1_coverage_end = voterec -> coverage_end;
								}
								is_found = 1;
							}
						}

						if(!is_found)
						{
							examine_vote_pairs[examine_vote_pair_count].map_1_pos = map_to_pos;
							examine_vote_pairs[examine_vote_pair_count].read_1_votes = voterec -> num_votes;
							examine_vote_pairs[examine_vote_pair_count].read_1_mapped_targ_2 = map_to_read2;
							examine_vote_pairs[examine_vote_pair_count].frag_number = next_frag_no;
							examine_vote_pairs[examine_vote_pair_count].read_1_coverage_start = voterec -> coverage_start;
							examine_vote_pairs[examine_vote_pair_count].read_1_coverage_end = voterec -> coverage_end;
							examine_vote_pair_count++;
						}
					}
				}
		}
		*/

		if(global_context -> debug_flags & GRA_DEBUG_ASSEMBLING)
			SUBREADprintf("HITs=%d\n", examine_vote_pair_count);
		int max_vote_score = 0, best_bases_len = 0, xk1;
		int next_current_frag_reversed_to_start=0;
		unsigned long long next_current_frag_no=0;

		for(test_pair =0; test_pair < examine_vote_pair_count; test_pair++)
		{
			//SUBREADprintf("Examing %d / %d pairs.\n", test_pair, examine_vote_pair_count);
			GRA_votes_pairs * tpair = &(examine_vote_pairs[test_pair]);

			if(tpair -> read_1_votes<global_context -> min_overlap_votes || tpair -> read_2_votes<global_context -> min_overlap_votes){
				//SUBREADprintf("SKIP 01\n");
				continue;
			}
			if(tpair -> read_1_mapped_targ_2 + tpair -> read_2_mapped_targ_2 != 1){
				//SUBREADprintf("SKIP 02\n");
				continue;
			}

			int is_target_read2_used = is_current_read2_extension?tpair -> read_2_mapped_targ_2:tpair -> read_1_mapped_targ_2;
			int current_rlen = is_current_read2_extension ? current_prop->length_2: current_prop->length_1; 
			int current_mate_rlen = is_current_read2_extension ? current_prop->length_1: current_prop->length_2; 
			int current_map_pos = is_current_read2_extension ? tpair->map_2_pos:tpair->map_1_pos;

			if(current_map_pos <= 0){
				//SUBREADprintf("SKIP 03\n");
				continue;	// the progress should be long enough in each step.
			}


			GRA_frag_properties_t * target_prop = global_context -> frag_properties_array + tpair -> frag_number;
			if((target_prop -> flags & GRA_FRAG_USED_NOW) || (target_prop -> flags & GRA_FRAG_USED_BEFORE) )
			{
				//SUBREADprintf("SKIP 0A\n");
				continue;
			}

			int target_rlen = (is_target_read2_used) ? target_prop -> length_2:target_prop -> length_1;
			int target_mate_rlen = (is_target_read2_used) ? target_prop -> length_1:target_prop -> length_2;
			int target_mate_coverage_start = is_current_read2_extension ? tpair -> read_1_coverage_start: tpair -> read_2_coverage_start;
			int target_mate_coverage_end   = is_current_read2_extension ? tpair -> read_1_coverage_end  : tpair -> read_2_coverage_end  ;
			int target_mate_vpos = is_current_read2_extension ? tpair -> map_1_pos: tpair -> map_2_pos ;


			if(target_mate_coverage_start > global_context -> subread_extract_step + max(0, -target_mate_vpos))
			{
				//SUBREADprintf("SKIP 04\n");
				continue;
			}
			if(target_mate_coverage_end < min(target_mate_rlen - target_mate_vpos, current_mate_rlen) - global_context -> subread_extract_step - 16)
			{
				//SUBREADprintf("SKIP 05 : %d < min(%d, %d) - %d\n", target_mate_coverage_end, target_mate_rlen - target_mate_vpos, current_mate_rlen , global_context -> subread_extract_step - 16);
				continue;
			}

			int mate_mm = 0;
			if(1)
			{
				char target_mate_seq[1201];
				char * current_mate_seq = is_current_read2_extension?seq_buf_1:seq_buf_2;

				get_frag_sequence(global_context, target_prop, !is_target_read2_used, target_mate_seq);

				int current_mate_match_start = max(-target_mate_vpos,0);
				int target_mate_match_start = max(target_mate_vpos,0);
				int mates_overlapping_len = min(current_mate_rlen + min(0, target_mate_vpos), target_mate_rlen - max(0, target_mate_vpos));

				for(xk1 = 0; xk1 < mates_overlapping_len; xk1++)
					mate_mm += (*(current_mate_seq + current_mate_match_start + xk1) != *(target_mate_seq + xk1 + target_mate_match_start));
				if(mate_mm > 1)
				{
					/*
					char spaces[200];
					for(xk1=0;xk1<200;xk1++) spaces[xk1]=' ';
					spaces[abs(target_mate_vpos)] = 0;

					SUBREADprintf("%s%s\n", target_mate_vpos <0?"":spaces, current_mate_seq);
					SUBREADprintf("%s%s\n", target_mate_vpos >0?"":spaces, target_mate_seq);
					

					SUBREADprintf("SKIP 0B :%d in %d \n", mate_mm, mates_overlapping_len);
					*/
					continue;
				}
			}

			int next_reversed_map = to_right ^ is_target_read2_used;
			int current_read_votes = is_current_read2_extension?tpair -> read_2_votes:tpair -> read_1_votes;
			//int current_mate_votes = is_current_read2_extension?tpair -> read_1_votes:tpair -> read_2_votes;
			char * current_seq_buf = is_current_read2_extension?seq_buf_2:seq_buf_1;
			int mm=0, target_votes=0, new_bases_len;
			char * new_bases_added;
			
			if( current_read_votes <  global_context -> min_extension_votes)
			{
				//SUBREADprintf("SKIP 06\n");
				continue;
			}

			int max_overlap_end_on_current = min(current_rlen  , target_rlen - current_map_pos);
			get_frag_sequence(global_context , target_prop , is_target_read2_used , target_seq_buf);

			for(xk1 = 0; xk1 < max_overlap_end_on_current; xk1++)
				mm += (*(current_seq_buf + xk1) != *(target_seq_buf  + xk1 + current_map_pos));

			int max_mm = 0;
			if(max_overlap_end_on_current > 20) max_mm = global_context -> maximum_mismatch_in20bp;

			if(mm>max_mm)
			{
				//SUBREADprintf("SKIP 07 : %d > %d\n", mm, max_mm);
				continue;
			}

			target_votes = current_read_votes;

			if(global_context -> debug_flags & GRA_DEBUG_ASSEMBLING)
			{
				char spaces[MAX_FILE_NAME_LENGTH];
				SUBREADprintf(" === %s {%d} V=%d, MM=%d === NEXT_REV=%d \n", to_right?">>>":"<<<", is_current_read2_extension?2:1, current_read_votes, mm, next_reversed_map);

				for(xk1=0;xk1< current_map_pos; xk1++)
					spaces[xk1]=' ';
				spaces[xk1]=0;

				SUBREADprintf("Q=%s%s\n", spaces,current_seq_buf);
				SUBREADprintf("T=%s\n",target_seq_buf );
				SUBREADputs("");
			}

			new_bases_added = target_seq_buf;
			new_bases_len   = current_map_pos;


			int target_replicate_count = is_target_read2_used?target_prop -> possible_extensions_to_read_2 : target_prop -> possible_extensions_to_read_1;
			//int target_mate_replicate_count = is_target_read2_used?target_prop -> possible_extensions_to_read_1 : target_prop -> possible_extensions_to_read_2;
			int target_replicate_score = 0;

			//SUBREADprintf("REPLICATE_COUNT=%d\n", target_replicate_count);

			if(target_replicate_count > 20)
				target_replicate_score = 3000000;
			else if(target_replicate_count > 10)
				target_replicate_score = 2000000;
			else if(target_replicate_count > 5)
				target_replicate_score = 1000000;
			else 
				target_replicate_score = target_replicate_count*100000;

			//target_replicate_score = 0;

			int target_votes_score = current_read_votes * 0 * 100;

			int extension_score = new_bases_len * 30;

			int mm_score = 1000000 - mm?1000000:0;

			int max_support_score = 0, is_reversed, nVotes = 0;

			if(0 && new_bases_len <800)
			{
				char new_tailK_bases[1000];
				int old_tail_len = min(200, global_context -> contig_len);
				strncpy(new_tailK_bases, global_context -> contig_str + global_context -> contig_len - old_tail_len, 999);
				new_tailK_bases[old_tail_len] = 0;
				strncat(new_tailK_bases, new_bases_added, 999);
				new_tailK_bases[old_tail_len + new_bases_len] = 0;


				int cursor;
				for(is_reversed = 0; is_reversed<2; is_reversed++)
				{
					lnhash_init_votes(votes1);
					for(cursor = 0; cursor < old_tail_len + new_bases_len - 15; cursor ++)
					{
						char * subread = new_tailK_bases + cursor;
						unsigned int subread_key = GR_subread2int(global_context, subread); 

						lnhash_query(&global_context -> assembling_index, votes1, subread_key, cursor);
					}

					int vote_X2, vote_Y2;
				
					for(vote_Y = 0; vote_Y < LNHASH_VOTE_ARRAY_HEIGHT ; vote_Y++) 
						for(vote_X = 0; vote_X < votes1->vote_record_items[vote_Y] ; vote_X++) 
						{
							int rec1_read2, rec1_reversed, rec1_pos;
							lnhash_vote_record_t * voterec1 = &(votes1 -> vote_records[vote_Y][vote_X]);
							unsigned long long next_frag_no1 = get_frag_number_part(voterec1 -> head_position, &rec1_read2, &rec1_reversed, &rec1_pos);
							for(vote_Y2 = 0; vote_Y2 < LNHASH_VOTE_ARRAY_HEIGHT ; vote_Y2++) 
								for(vote_X2 = 0; vote_X2 < votes1->vote_record_items[vote_Y2] ; vote_X2++) 
								{
									int rec2_read2, rec2_reversed, rec2_pos;
									lnhash_vote_record_t * voterec2 = &(votes1 -> vote_records[vote_Y2][vote_X2]);
									unsigned long long next_frag_no2 = get_frag_number_part(voterec2 -> head_position, &rec2_read2, &rec2_reversed, &rec2_pos);

									if(next_frag_no2==next_frag_no1)
									{
										int min_start = min(voterec1 -> coverage_start, voterec2 -> coverage_start);
										int max_end = max(voterec1 -> coverage_end, voterec2 -> coverage_end);

										if(is_reversed)
										{
											int tmp = min_start;
											min_start = old_tail_len + new_bases_len - max_end;
											max_end = old_tail_len + new_bases_len - tmp;
										}

										assert(min_start < max_end);

										//SUBREADprintf("COV = %d ~ %d ; SPLIT AT %d\n", min_start, max_end , old_tail_len);
										if(min_start < old_tail_len - 16 && max_end >= min(old_tail_len + 16,  old_tail_len + new_bases_len - 1)) max_support_score = min(max_support_score +  500, 0x70000000);
									}
								}
							nVotes ++;
						}

					if(0==is_reversed) reverse_read(new_tailK_bases, old_tail_len + new_bases_len, GENE_SPACE_BASE);
				}
			}

			//SUBREADprintf("nSup=%d\tAllVotes=%d\tEXT=%d bases\n", max_support_score / 5000000, nVotes, new_bases_len);

			int target_vote_score = target_replicate_score + target_votes_score + mm_score + extension_score + max_support_score;

			if( max_vote_score < target_vote_score)
			{
				max_vote_score = target_vote_score;
				next_current_frag_no = tpair -> frag_number;
				next_current_frag_reversed_to_start = next_reversed_map;

				memcpy(best_bases_added , new_bases_added , new_bases_len);
				best_bases_len = new_bases_len;
			}

		}

		//SUBREADprintf("MAX_V = %d\n", max_vote_score);
		if(max_vote_score < 1)break;

		if(to_right)//is_current_read2_extension)
			reverse_read(best_bases_added, best_bases_len, GENE_SPACE_BASE);
		//SUBREADprintf("%s EXTEND , USING %s , ADD SEQ %d bases , THIS IS %s , NEXT IS %s\n", to_right?"RIGHT":"LEFT", is_current_read2_extension?"READ_2":"READ_1", best_bases_len,  global_context -> current_frag_reversed_to_start ? "REV_STRAND":"SAME_STRAND", next_current_frag_reversed_to_start?"REV_STRAND":"SAME_STRAND");
		int add_status = GRA_add_new_extension_part(global_context, best_bases_added, best_bases_len, to_right);
		if(add_status<0)break;
		if(add_status==0) 
		{

			global_context -> current_frag_no = next_current_frag_no;
			global_context -> current_frag_reversed_to_start = next_current_frag_reversed_to_start;

			if(global_context -> debug_flags & GRA_DEBUG_ASSEMBLING)
				SUBREADprintf("  === %s === REV=%d    FRAG=%llu === \n", to_right?">>>":"<<<" , global_context -> current_frag_reversed_to_start ,  global_context -> current_frag_no );
		}
		GRA_frag_properties_t * target_prop = global_context -> frag_properties_array + next_current_frag_no;
		target_prop -> flags |= GRA_FRAG_USED_NOW; 
		free(vote1_sorted);
		free(vote2_sorted);
	}
	free(best_bases_added);
	free(examine_vote_pairs);
	free(seq_buf_1);
	free(seq_buf_2);
	free(target_seq_buf);
	free(votes2);
	free(votes1);
	return 0;
}

int GRA_reassemble_reads(GRA_global_context_t* global_context)
{
	GRA_init_extension_part(global_context);

	char * seq_buf_1 = malloc(MAX_READ_LENGTH);
	get_frag_sequence(global_context, global_context -> start_fragment_props, 0 , seq_buf_1);
	GRA_add_new_extension_part(global_context , seq_buf_1 , global_context -> start_fragment_props -> length_1 - global_context -> start_fragment_props -> overlapping, 0);

	global_context -> current_frag_reversed_to_start = 0;
	global_context -> current_frag_no = global_context -> start_fragment_props - global_context -> frag_properties_array;
	int ret =  GRA_reassemble_reads_LR(global_context,0);

	get_frag_sequence(global_context, global_context -> start_fragment_props, 1 , seq_buf_1);
	reverse_read(seq_buf_1 , global_context -> start_fragment_props -> length_2, GENE_SPACE_BASE);
	GRA_add_new_extension_part(global_context , seq_buf_1 , global_context -> start_fragment_props -> length_2, 1);
	free(seq_buf_1);

	global_context -> current_frag_reversed_to_start = 0;
	global_context -> current_frag_no = global_context -> start_fragment_props - global_context -> frag_properties_array;
	return ret || GRA_reassemble_reads_LR(global_context,1);
}

void GRA_destroy_context(GRA_global_context_t* global_context)
{
	
	char tmp_fname[100];
	lnhash_destroy(&global_context -> assembling_index);
	if(global_context -> fd_fragproperties)
		close(global_context -> fd_fragproperties);
	else free(global_context -> frag_properties_array);

	if(stderr != global_context -> fp_output && stdout != global_context -> fp_output)
		fclose(global_context -> fp_output);


	sprintf(tmp_fname, "%s.fraglist", global_context -> tmp_file_prefix);
	unlink(tmp_fname);
	sprintf(tmp_fname, "%s.fragproperties", global_context -> tmp_file_prefix);
	unlink(tmp_fname);

	SUBREADprintf("\nRead assembly finished in %.1f minutes.\n\n", (miltime() - global_context -> start_time)/60);
}


int GRA_memory_map(GRA_global_context_t* global_context)
{
	char tmp_fname[MAX_FILE_NAME_LENGTH];

	sprintf(tmp_fname, "%s.fragproperties", global_context -> tmp_file_prefix);

	if(global_context -> use_huge_memory)
	{
		global_context -> fd_fragproperties = 0;

		global_context -> frag_properties_array = malloc(global_context -> all_fragments * sizeof(GRA_frag_properties_t));
		FILE * fp = fopen(tmp_fname,"r");
		fread(global_context -> frag_properties_array, sizeof(GRA_frag_properties_t), global_context -> all_fragments , fp);
		fclose(fp);
		SUBREADputs("\nThe property file is load into memory.\n");
	}
	else
	{
		global_context -> fd_fragproperties = open(tmp_fname,O_RDWR);
		global_context -> frag_properties_array = mmap64(NULL, global_context -> all_fragments * sizeof(GRA_frag_properties_t),
			PROT_WRITE|PROT_READ,  MAP_PRIVATE, global_context -> fd_fragproperties, 0);
		SUBREADputs("\nThe property file is mapped to memory.\n");
	}

	assert(global_context -> frag_properties_array);


	return 0;
}


int GRA_find_start_read(GRA_global_context_t* global_context)
{
	unsigned long long frag_number, possible_starts = 0, possible_used=0;
	unsigned long long report_step = max(1, global_context -> all_fragments*1./1000);

	for(frag_number = global_context -> start_read_pointer; frag_number < global_context -> all_fragments ; frag_number++)
	{
		if(frag_number % report_step == 0)
			SUBREADprintf("  Assembling contigs, %.1f%% (%llu) finished in %.1f minutes. %llu fragments were disabled.\n", frag_number*100./global_context -> all_fragments, frag_number, (miltime() - global_context -> start_time)/60, global_context -> marked_reported_fragments);

		GRA_frag_properties_t * frag_prop = global_context -> frag_properties_array + frag_number;

		if(frag_prop -> replicate_representitive_1 !=0 && frag_prop -> replicate_representitive_2 != 0)
			continue;
		if(frag_prop -> overlapping < global_context -> min_overlapping_bases_between_2ends) continue;

		possible_starts ++;

		if((frag_prop->flags & GRA_FRAG_USED_BEFORE) || (frag_prop->flags & GRA_FRAG_USED_NOW)){
			possible_used ++;
			continue;
		}

		if(frag_prop->length_1 >=  global_context -> shortest_startread_allowed && frag_prop->length_2 >=  global_context -> shortest_startread_allowed )
		{
			frag_prop -> flags |= GRA_FRAG_USED_NOW;
			global_context -> start_read_pointer = frag_number + 1;
			global_context -> start_fragment_props = frag_prop;
			return 0;
		} 
	}

	global_context -> no_more_assembly = 1;

	return 0;
}

int main(int argc, char ** argv)
{

	int c;
	int option_index = 0;
	GRA_global_context_t * global_context = malloc(sizeof(GRA_global_context_t));

	GRA_init_context(global_context);
	optind=0;
	opterr=1;
	optopt=63;
	int tmp_trim_qual = -1;

	if(argc<2)
	{
		GRA_print_usage();
		return -1;
	}

	char * tmptmp = malloc(MAX_FILE_NAME_LENGTH);
	while ((c = getopt_long (argc, argv, "2:t:Q:V:L:i:I:o:bHCR6", GRA_long_options, &option_index)) != -1)
		switch(c)
		{
			case '2':
				global_context -> maximum_mismatch_in20bp = atoi(optarg);
			break;
			case 'H':
				global_context -> use_huge_memory = 1;
			break;
			case 'R':
				global_context -> reverse_unmapped_reads = 1;
			break;
			case 'i':
				strcpy(global_context -> input_file_name, optarg);
			break;
			case 'I':
				strcpy(global_context -> input_file_name2, optarg);
			break;
			case 'b':
				global_context -> is_BAM_input = 1;
			break;
			case 'o':
				strcpy(global_context -> output_file_name, optarg);
			break;
			case 'L':
				global_context -> minimum_contig_length = atoi(optarg);
			break;
			case '6':
				global_context -> phred_offset = 64;
			break;
			case 'V':
				global_context -> min_overlap_votes = atoi(optarg);
			break;
			case 'v':
				global_context -> min_extension_votes = atoi(optarg);
			break;
			case 't':
				strcpy(tmptmp, optarg);
				strcat(tmptmp, "/");
				strcat(tmptmp, global_context -> tmp_file_prefix);
				strcpy(global_context -> tmp_file_prefix, tmptmp);
			break;
			case 'T':
				global_context -> total_threads = atoi(optarg);
			break;
			case 'Q':
				tmp_trim_qual = atoi(optarg);
			break;
			case 'C':
				global_context -> check_read_names = 0;
			break;

			default:
				GRA_print_usage();
				return -1;
		}
	
	free(tmptmp);

	int ret = 0;

	if(tmp_trim_qual >= 0)
		global_context -> read_trim_base_quality = tmp_trim_qual + global_context -> phred_offset;

	ret = ret || GRA_check_config(global_context);
	ret = ret || GRA_trim_reads(global_context);
	ret = ret || GRA_memory_map(global_context);
	//ret = ret || GRA_mark_read_extension(global_context);
	ret = ret || GRA_build_read_index(global_context);

	while(!ret){
		ret = ret || GRA_find_start_read(global_context);
		if(global_context -> no_more_assembly)
			break;

		ret = ret || GRA_reassemble_reads(global_context);
		ret = ret || GRA_write_current_assembly(global_context);
	}

	GRA_destroy_context(global_context);
	free(global_context);

	return ret;

}
