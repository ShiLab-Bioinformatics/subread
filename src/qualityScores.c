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
#include <string.h>
#include <zlib.h>
#include <getopt.h>
#include <ctype.h>
#include "subread.h"
#include "core.h"
#include "sambam-file.h"
#include "input-files.h"

//static int ERROR_PROB_INT_TABLE[] = { 1000000 , 794328 , 630957 , 501187 , 398107 , 316227 , 251188 , 199526 , 158489 , 125892 , 100000 , 79432 , 63095 , 50118 , 39810 , 31622 , 25118 , 19952 , 15848 , 12589 , 10000 , 7943 , 6309 , 5011 , 3981 , 3162 , 2511 , 1995 , 1584 , 1258 , 1000 , 794 , 630 , 501 , 398 , 316 , 251 , 199 , 158 , 125 , 100 , 79 , 63 , 50 , 39 , 31 , 25 , 19 , 15 , 12 , 10 , 7 , 6 , 5 , 3 , 3 , 2 , 1 , 1 , 1 , 1 , 0 , 0 , 0 ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


typedef struct {
	char * input_name;
	char * output_name;
	int input_file_type;
	int phred_offset;
	long long needed_reads;
	int sam_end;
	int max_read_length;
	char * IO_line_buff;


	FILE * result_fp;
	union{
		SamBam_FILE * sambam_reader;
		gzFile gzfq_reader;
		FILE * fq_reader; 
	};

	unsigned long long scored_reads;
	unsigned long long line_counter;
	unsigned long long scored_bases;
	unsigned long long total_phred_sum;
	unsigned long long total_errorprob_sum;
	unsigned long long * base_phred_sum;
	unsigned long long * base_errorprob_sum;

	int quality_offset_warning_shown;
} qualscore_context;

int finalise_qs_context(qualscore_context * qs_context, int prev_ret)
{
	fclose(qs_context -> result_fp);
	free(qs_context-> IO_line_buff);

	SUBREADprintf("\n");

	if(prev_ret)
	{
		SUBREADprintf("No results were generated.\n");
	}
	else
	{
		if(qs_context -> input_file_type == FILE_TYPE_FASTQ)
			fclose(qs_context -> fq_reader);
		else if(qs_context -> input_file_type == FILE_TYPE_SAM || qs_context -> input_file_type == FILE_TYPE_BAM)
			SamBam_fclose(qs_context -> sambam_reader);
		else if(qs_context -> input_file_type == FILE_TYPE_GZIP_FASTQ)
			gzclose(qs_context -> gzfq_reader);
		else return 1;

		//double avg_phred = qs_context -> total_phred_sum*1./ qs_context->scored_bases;
		//double avg_prob = qs_context -> total_errorprob_sum*1./1000000/qs_context->scored_bases;

		SUBREADprintf("Completed successfully. Quality scores for %llu reads (equally spaced in the file) are returned.\n", qs_context -> scored_reads);

		if(qs_context-> quality_offset_warning_shown)
			SUBREADprintf("However, the Phred score offset (%d) seemed to be wrong. The quality scores can be meaningless.\n", qs_context -> phred_offset);
	//	else
	//		SUBREADprintf("Average Phred score is %.1f. Average wrong sequence probability is %.2f%%.\n", avg_phred, avg_prob*100);
	}

	SUBREADprintf("\n");
	return 0;
}

int add_read_scores(qualscore_context * qs_context, char * qual_text, unsigned int total_cols) {
	int xk1, nch, out_cursor = 0, aval = 1;

	qs_context -> IO_line_buff[0]=0;
	for(xk1=0; xk1 < total_cols; xk1++){
		nch = -1;
		if(aval) nch = qual_text[xk1];	
		if(nch == 0){
			aval = 0;
			nch=-1;
		}

		if(nch>0){
			nch -= qs_context -> phred_offset;
			if(nch <=0 || nch >= 65)
			{
				if(!qs_context-> quality_offset_warning_shown)
				{
					SUBREADprintf("Warning: the Phred score offset (%d) seems wrong : it ended up with a phred score of %d.\n",  qs_context -> phred_offset, nch);
					qs_context-> quality_offset_warning_shown=1;
				}
			}
			out_cursor+=sprintf(qs_context->IO_line_buff+out_cursor,"%d,", nch);
		} else {
			out_cursor+=sprintf(qs_context->IO_line_buff+out_cursor,"NA,");
		}
	}
	if(out_cursor>0)
		qs_context->IO_line_buff[out_cursor-1]='\n';//remove the last comma and change it to a linebreak.
	fwrite(qs_context->IO_line_buff, 1, out_cursor, qs_context -> result_fp);
	qs_context -> scored_reads ++;

	return 0;
}

int qs_next_qual(qualscore_context * qs_context, char * qual_buff)
{
	int ret = 0, need_reverse=0;
	if(qs_context -> input_file_type == FILE_TYPE_FASTQ || qs_context -> input_file_type == FILE_TYPE_GZIP_FASTQ)
	{
		char * retstr;

		if(qs_context -> input_file_type == FILE_TYPE_FASTQ)
			retstr = fgets_noempty(qual_buff, 2999, qs_context -> fq_reader);
		else
			retstr = gzgets_noempty(qs_context -> fq_reader, qual_buff, 2999);

		if(!retstr) ret = 1;
		else
		{
			if(retstr[0]!='@')
			{
				SUBREADprintf("ERROR: The input fastq file has a wrong format.\n");
				return -1;
			}
			qs_context -> line_counter ++;


			if(qs_context -> input_file_type == FILE_TYPE_FASTQ)
				retstr = fgets_noempty(qual_buff, 2999, qs_context -> fq_reader);
			else
				retstr = gzgets_noempty(qs_context -> fq_reader, qual_buff, 2999);

			if(qs_context -> input_file_type == FILE_TYPE_FASTQ)
				retstr = fgets_noempty(qual_buff, 2999, qs_context -> fq_reader);
			else
				retstr = gzgets_noempty(qs_context -> fq_reader, qual_buff, 2999);

			if(!retstr) ret = 1;
			else 
			{
				if(retstr[0]!='+')
				{
					SUBREADprintf("ERROR: The input fastq file has a wrong format.\n");
					return -1;
				}
				qs_context -> line_counter += 2;
			}

			if(qs_context -> input_file_type == FILE_TYPE_FASTQ)
				retstr = fgets_noempty(qual_buff, 2999, qs_context -> fq_reader);
			else
				retstr = gzgets_noempty(qs_context -> fq_reader, qual_buff, 2999);
			
			if(!retstr) ret = 1;
		}
	}
	else// sam or bam
	{
		while(1)
		{
			char * retstr = SamBam_fgets(qs_context -> sambam_reader , qs_context -> IO_line_buff, 5999, qual_buff?1:0);
			if(!retstr) ret = 1;
			else
			{
				if(qs_context -> IO_line_buff[0]!='@')
				{
					char * tmpptr;
					int flags = 0;
					char * nstr = strtok_r(qs_context -> IO_line_buff, "\t", &tmpptr);	// name
					if(!nstr) continue;

					nstr = strtok_r(NULL, "\t", &tmpptr);	// name
					if(nstr) flags = atoi(nstr);
					else continue; 

					if(0&& qs_context -> sam_end>0 && !(flags &SAM_FLAG_PAIRED_TASK ))
					{
						SUBREADprintf("ERROR: The SAM/BAM input file contains single-end reads, please do not specify paired-end specific arguments (e.g., first-end or second-end).\n");
						ret = 1;
						break;
					}

					nstr = strtok_r(NULL, "\t", &tmpptr);	// chro 
					if(!nstr) continue; 

					nstr = strtok_r(NULL, "\t", &tmpptr);	// pos 
					if(!nstr) continue; 

					nstr = strtok_r(NULL, "\t", &tmpptr);	// mapq 
					if(!nstr) continue; 

					nstr = strtok_r(NULL, "\t", &tmpptr);	// cigar 
					if(!nstr) continue; 

					nstr = strtok_r(NULL, "\t", &tmpptr);	// matechr 
					if(!nstr) continue; 

					nstr = strtok_r(NULL, "\t", &tmpptr);	// matepos
					if(!nstr) continue; 

					nstr = strtok_r(NULL, "\t", &tmpptr);	// tlen 
					if(!nstr) continue; 

					nstr = strtok_r(NULL, "\t", &tmpptr);	// read_text 
					if(!nstr) continue; 

					nstr = strtok_r(NULL, "\t", &tmpptr);	// qual_text 
					if(!nstr) continue; 

					if((qs_context -> sam_end == 1 && (!(flags & SAM_FLAG_SECOND_READ_IN_PAIR))) || (qs_context -> sam_end == 2 && (flags & SAM_FLAG_SECOND_READ_IN_PAIR))  || (qs_context -> sam_end == 0))
					{
						if((flags & 256)==0 )
						{
							if(flags & SAM_FLAG_REVERSE_STRAND_MATCHED) need_reverse=1;
							strcpy(qual_buff, nstr);
							break;
						}
					}
				}
			}
			if(ret) break;
		}
	}

	int qblen = strlen(qual_buff);
	if(qual_buff[qblen-1]=='\n'){
		qual_buff[qblen-1]=0;
		qblen--;
	}
	if(need_reverse)
		reverse_quality(qual_buff, qblen);

	return ret;
}

int rewind_qs_file(qualscore_context * qs_context)
{
	if(qs_context -> input_file_type == FILE_TYPE_FASTQ)
		fclose(qs_context -> fq_reader);
	else if(qs_context -> input_file_type == FILE_TYPE_SAM || qs_context -> input_file_type == FILE_TYPE_BAM)
		SamBam_fclose(qs_context -> sambam_reader);
	else if(qs_context -> input_file_type == FILE_TYPE_GZIP_FASTQ)
		gzclose(qs_context -> gzfq_reader);

	if(qs_context -> input_file_type == FILE_TYPE_FASTQ)
		qs_context -> fq_reader = f_subr_open(qs_context->input_name, "r");
	else if(qs_context -> input_file_type == FILE_TYPE_SAM || qs_context -> input_file_type == FILE_TYPE_BAM)
		qs_context -> sambam_reader = SamBam_fopen(qs_context->input_name, (qs_context -> input_file_type == FILE_TYPE_BAM)?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM); 
	else if(qs_context -> input_file_type == FILE_TYPE_GZIP_FASTQ)
		qs_context -> gzfq_reader = gzopen(qs_context->input_name, "rb");


	return 0;
}

int qs_str_input_type(char * strtype)
{
	if(strcmp(strtype,"SAM")==0)
		return FILE_TYPE_SAM;
	else if(strcmp(strtype, "BAM")==0)
		return FILE_TYPE_BAM;
	else if(strcmp(strtype, "FASTQ")==0)
		return FILE_TYPE_FASTQ;
	else if(strcmp(strtype, "GZFASTQ")==0)
		return FILE_TYPE_GZIP_FASTQ;

	return -1;
}


int start_qs_context(qualscore_context * qs_context)
{
	int ret =0;
	SUBREADprintf("\nqualityScores %s\n\n", SUBREAD_VERSION);
	
	if(qs_context -> input_file_type == FILE_TYPE_FASTQ)
		qs_context -> fq_reader = f_subr_open(qs_context->input_name, "r");
	else if(qs_context -> input_file_type == FILE_TYPE_SAM || qs_context -> input_file_type == FILE_TYPE_BAM)
		qs_context -> sambam_reader = SamBam_fopen(qs_context->input_name, (qs_context -> input_file_type == FILE_TYPE_BAM)?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM); 
	else if(qs_context -> input_file_type == FILE_TYPE_GZIP_FASTQ)
		qs_context -> gzfq_reader = gzopen(qs_context->input_name, "rb");
	else ret = 1;

	qs_context -> result_fp = f_subr_open(qs_context->output_name, "w");
	if(!qs_context -> result_fp)
	{
		ret = 1;
		SUBREADprintf("ERROR: cannot open output file '%s'.\n", qs_context->output_name);
	}
	else 
	{
		if(ret || !qs_context -> fq_reader)
		{
			SUBREADprintf("ERROR: cannot open input file '%s'.\n", qs_context->input_name);
			SUBREADprintf("       The file does not exist, or it is not in the %s format.\n", qs_context -> input_file_type == FILE_TYPE_FASTQ?"FASTQ":(qs_context -> input_file_type == FILE_TYPE_SAM?"SAM":(qs_context -> input_file_type == FILE_TYPE_BAM?"BAM":"gzipped FASTQ")));
			ret=1;
		}
		else
		{
			qs_context -> IO_line_buff = malloc(6000);
			if(!qs_context -> IO_line_buff){
				SUBREADprintf("ERROR: cannot allocate 3KB memory??\n");
				ret=1;
			}
		}
	}


	if(qs_context->phred_offset!=33 && (qs_context -> input_file_type == FILE_TYPE_SAM || qs_context -> input_file_type == FILE_TYPE_BAM))
		SUBREADprintf("WARNING: your input file is in SAM or BAM format. In most cases, the Phred score offset in SAM or BAM is 33.\n");

	if(!ret)
		SUBREADprintf("Scan the input file...\n");

	return ret;
}

int retrieve_scores(char ** input, int *offset_pt, long long *size, int * sam_whichend, char **input_file_type, char ** output_sc){
  
	qualscore_context qs_context;
	memset(&qs_context, 0, sizeof(qualscore_context));

	// input_name is a FastQ file, a GZ FastQ file, a SAM file or a BAM file.
	qs_context.input_name = input[0];
	// output_sc_name is the final result file used in donwstream analysis.
	// the output_sc file has a format like
	// 40,40,40,40,40,....,30,30,10,3,3,3,3	   { elements = read length } 
	// 40,35,35,30,30,....,30,30,10,3,3,3,3    { elements = read length } 
	// .... {total lines = needed_reads +- 1}
	qs_context.output_name = output_sc[0];

	// input_file_type can be "BAM", "SAM", "FASTQ" or "GZFASTQ"
	char *input_file_type_str = *input_file_type;

	qs_context.input_file_type = qs_str_input_type(input_file_type_str);

	// if the input file is SAM or BAM, sam_end=0 : both ends; sam_end=1 : first end; sam_end=2 : second end.
	qs_context.sam_end = *sam_whichend;

	// offset is 33 or 64
	qs_context.phred_offset = *offset_pt;

	// total number of sampled reads
	qs_context.needed_reads = *size;

	double start_time =  miltime();

	int ret = start_qs_context(&qs_context);

	if(!ret)
	{

		unsigned long long read_number = 0;
		unsigned int all_read_max_len = 0;

		char *linebuff = malloc(3000);

		while (1)
		{
			int reti = qs_next_qual(&qs_context, linebuff);
			if(reti<0) ret = 1; 
			if(reti) break;
			read_number++;
			int rlen = strlen(linebuff);
			if(linebuff[rlen-1]=='\n') rlen--;
			all_read_max_len = max(all_read_max_len, rlen);
			if(read_number % 1000000 == 0)
				SUBREADprintf("  %llu reads have been scanned in %.1f seconds.\n", read_number, miltime() -  start_time);
		}

		if(!ret)
		{
			if(read_number<1)
				SUBREADprintf("Warning: the input file is empty%s.\n", (qs_context.sam_end ==2)?", or it does not contain any second-end read":"");
			else
			{

				rewind_qs_file(&qs_context);

				double sampling_rate = (read_number-0.1)/qs_context.needed_reads;
				 
				double next_sample_number = 0;
				sampling_rate = max(1.0, sampling_rate);

				SUBREADprintf("Totally %llu reads were scanned; the sampling interval is %d.\nNow extract read quality information...\n", read_number, (int)(sampling_rate));

				read_number = 0;

				while(1)
				{
					int reti = qs_next_qual(&qs_context, linebuff);
					if(reti) break;
					if(read_number >= next_sample_number - 0.0000000001)
					{
						add_read_scores(&qs_context, linebuff, all_read_max_len);
						next_sample_number += sampling_rate;
					}
					read_number++;

					if(read_number % 1000000 == 0)
						SUBREADprintf("  %llu reads have been scanned in %.1f seconds.\n", read_number, miltime() -  start_time);
				}
			}
		}

		free(linebuff);
	}

	finalise_qs_context(&qs_context, ret);
	return ret;
}


static struct option qs_long_options[] =
{
	{"BAMinput",  no_argument, 0, '9'},
	{"SAMinput",  no_argument, 0, '8'},
	{"gzFASTQinput",  no_argument, 0, '7'},
	{"FASTQinput",  no_argument, 0, '6'},
	{"first-end",  no_argument, 0, '1'},
	{"second-end",  no_argument, 0, '2'},
	{"both-ends",  no_argument, 0, '0'},
	{"counted-reads",  required_argument, 0, 'n'},
	{"phred-offset",  required_argument, 0, 'P'},
	{0, 0, 0, 0}
};

void qualscore_usage()
{
	SUBREADprintf("\nqualityScore Version %s\n\n", SUBREAD_VERSION);
	SUBREADputs("  Retrieve Phred score for read bases\n");
	SUBREADputs("Usage:");
	SUBREADputs("");
	SUBREADputs("  ./qualityScores [options] -i <input_file> -o <output_file>");
	SUBREADputs("");
	SUBREADputs("Required arguments:");
	SUBREADputs("");
	SUBREADputs("  -i <string>  Name of input file including read data. The default format is");
	SUBREADputs("               Fastq.");
	SUBREADputs("");
	SUBREADputs("  -o <string>  Name of output file that is a text file including Phred scores");
	SUBREADputs("               for each read base.");
	SUBREADputs("");
	SUBREADputs("Optional arguments:");
	SUBREADputs("");
	SUBREADputs("  --gzFASTQinput Input file is in gzipped Fastq format.");
	SUBREADputs("");
	SUBREADputs("  --BAMinput     Input file is in BAM format.");
	SUBREADputs("");
	SUBREADputs("  --SAMinput     Input file is in SAM format.");
	SUBREADputs("");
	SUBREADputs("  --first-end    Use only first reads in paired-end data. Only applicable for");
	SUBREADputs("                 paired-end BAM/SAM input.");
	SUBREADputs("");
	SUBREADputs("  --second-end   Use only second reads in paired-end data. Only applicable for");
	SUBREADputs("                 paired-end BAM/SAM input.");
	SUBREADputs("");
	SUBREADputs("  --counted-reads <int> Total number of reads to be extracted from the input");
	SUBREADputs("                 file. 10,000 by default.");
	SUBREADputs("");
	SUBREADputs("  --phred-offset <33|64> refer to subread aligner.");
	SUBREADputs("");

}

#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int main_qualityScores(int argc, char ** argv)
#endif
{

	int ret;
	int c;
	long long needed_reads = 10000ll;
	int option_index = 0 , offset_pt = 33, sam_end = 0;
	char in_name[MAX_FILE_NAME_LENGTH];
	char out_name[MAX_FILE_NAME_LENGTH];
	char * input_format = "FASTQ";
	char * in_nameptr = in_name;
	char * out_nameptr = out_name;

	in_name[0]=out_name[0]=0;

	sam_end = 1;	// default: first-end only
	optind = 0;
	opterr = 1;
	optopt = 63;


	while((c = getopt_long (argc, argv, "n:i:o:P:12987", qs_long_options, &option_index)) != -1)
	{
		switch(c){
			case 'n':
				needed_reads = atoll(optarg);
				break;
			case 'P':
				offset_pt = optarg[0]=='6'?64:33;
				break;
			case '9':
				input_format="BAM";
				break;
			case '8':
				input_format="SAM";
				break;
			case '7':
				input_format="GZFASTQ";
				break;
			case '6':
				input_format="FASTQ";
				break;
			case '0':
				sam_end = 0;
				break;
			case '1':
				sam_end = 1;
				break;
			case '2':
				sam_end = 2;
				break;
			case 'i':
				strcpy(in_name, optarg);
				break;
			case 'o':
				strcpy(out_name, optarg);
				break;
			case 0:	// long names
				break;
			default:
				qualscore_usage();
				return 0;
		}
	}

	if(in_name[0]==0 || out_name[0]==0)
	{
		qualscore_usage();
		return 0;
	}

	ret = retrieve_scores(&in_nameptr, &offset_pt, &needed_reads, &sam_end, &input_format, &out_nameptr);

	return ret;
}

