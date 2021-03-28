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

/***************************************************************

   The ASCII Art used in this file was generated using FIGlet and
   the big font, contributed by Glenn Chappell to FIGlet.
  
  ***************************************************************/
  
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <sys/time.h>
#include <getopt.h>
#include <sys/types.h>
#ifndef __MINGW32__
#include <sys/resource.h>
#endif
#include <unistd.h>
#include <sys/stat.h>
#include <locale.h>
#include <ctype.h>



#include "subread.h"
#include "sublog.h"
#include "core.h"
#include "input-files.h"
#include "input-blc.h"
#include "sorted-hashtable.h"
#include "HelperFunctions.h"

#include "core-bigtable.h"
#include "core-indel.h"
#include "core-junction.h"

static struct option long_options[] =
{
	{"memory-optimisation",  required_argument, 0, 0},
	{0, 0, 0, 0}
};

int (*progress_report_callback)(int, int, int);

void core_version_number(char * program)
{
	SUBREADprintf("\n%s v%s\n\n" , program, SUBREAD_VERSION);
}

int is_valid_float(char * optarg, char * optname){
	int xk1=0;
	while(1){
		int nch = optarg[xk1++];
		if(!nch){
			if(xk1 == 1){
				SUBREADprintf("Value for argumant %s-%s is missing.\n", optname[1]?"-":"", optname);
				return 0;
			}
			break;
		}
		if(( nch!='-' || xk1 > 1 ) && nch != '.' && !isdigit(nch)){
			SUBREADprintf("Value for argumant %s-%s is not a valid number: '%s'\n", optname[1]?"-":"", optname, optarg);
			return 0;
		}
	}
	return 1;
}

int is_valid_digit(char * optarg, char * optname){
	int xk1=0;
	while(1){
		int nch = optarg[xk1++];
		if(!nch){
			if(xk1 == 1){
				SUBREADprintf("Value for argumant %s-%s is missing.\n", optname[1]?"-":"", optname);
				return 0;
			}
			break;
		}
		if(( nch!='-' || xk1 > 1 ) && !isdigit(nch)){
			SUBREADprintf("Value for argumant %s-%s is not a valid number: '%s'\n", optname[1]?"-":"", optname, optarg);
			return 0;
		}
	}
	return 1;
}

int is_valid_digit_range(char * optarg, char * optname, int min, int max_inc){
	if(!is_valid_digit( optarg, optname )) return 0;

	int tv = atoi(optarg);
	if(tv < min || tv > max_inc){
		SUBREADprintf("Value for argumant %s-%s is out of range: %d to %d\n", optname[1]?"-":"", optname, min, max_inc);
		return 0;
	}
	return 1;
}

void warning_file_limit()
{
	#ifndef __MINGW32__
	struct rlimit limit_st;
	getrlimit(RLIMIT_NOFILE, & limit_st);

	{
		if(min(limit_st.rlim_cur , limit_st.rlim_max) < 400)
		{
			print_in_box(80,0,0,"WARNING This operation needs to open many files at same time,");
			print_in_box(80,0,0,"	but your OS only allows to open %d files.", min(limit_st.rlim_cur , limit_st.rlim_max));
			print_in_box(80,0,0,"	You can use command 'ulimit -n 500' to raise this limit");
			print_in_box(80,0,0,"	to 500, or the program may become very slow.");
			print_in_box(80,0,0,"");
		}
	}
	#endif
}

int exec_cmd(char * cmd, char * outstr, int out_limit){
	FILE* pipe = popen(cmd, "r");
	if (!pipe) return -1;

	outstr[0]=0;
	char * linebuf = malloc(3000);
	int out_ptr = 0;

	while (!feof(pipe)) {
		if (fgets(linebuf, 128, pipe) != NULL){
			if(out_ptr + strlen(linebuf) < out_limit){
				strcat(outstr+out_ptr, linebuf);
				out_ptr += strlen(linebuf);
			}
		}
	}
	pclose(pipe);
	free(linebuf);
	return out_ptr;
}

void print_in_box(int line_width, int is_boundary, int options, char * pattern,...)
{
	int put_color_for_colon, is_center, is_wrapped;
	va_list args;
	va_start(args , pattern);
	char * content, *out_line_buff;

	put_color_for_colon = (options & PRINT_BOX_NOCOLOR_FOR_COLON)?0:1;
	is_center = (options & PRINT_BOX_CENTER)?1:0;
	is_wrapped = (options & PRINT_BOX_WRAPPED)?1:0;
	
	content= malloc(1200);
	int content_len = vsprintf(content, pattern, args);
	out_line_buff= malloc(1200);
	out_line_buff[0]=0;;

	if(is_wrapped){
		int seg_i = 0;
		for(seg_i=0; seg_i < content_len; seg_i += line_width-7){
			strcpy(out_line_buff, content + seg_i);
			out_line_buff[line_width-7] = 0;
			
			print_in_box(line_width, is_boundary, options & (~PRINT_BOX_WRAPPED), out_line_buff);
		}
	}else{
		int is_R_code,x1,content_len = strlen(content), state, txt_len, is_cut = 0, real_lenwidth;

		is_R_code = 1;
		#ifdef MAKE_STANDALONE
			is_R_code = 0;
		#endif

		if(content_len>0&&content[content_len-1]=='\r'){
			content_len--;
			content[content_len] = 0;
		}

		if(content_len>0&&content[content_len-1]=='\n'){
			content_len--;
			content[content_len] = 0;
		}

		state = 0;
		txt_len = 0;
		real_lenwidth = line_width;
		for(x1 = 0; content [x1]; x1++)
		{
			char nch = content [x1];
			if(nch == CHAR_ESC)
				state = 1;
			if(state){
				real_lenwidth --;
			}else{
				txt_len++;
				
				if(txt_len == 80 - 6)
				{
					is_cut = 1;
				} 
			}

			if(nch == 'm' && state)
				state = 0;
		}

		if(is_cut)
		{
			state = 0;
			txt_len = 0;
			for(x1 = 0; content [x1]; x1++)
			{
				char nch = content [x1];
				if(nch == CHAR_ESC)
					state = 1;
				if(!state){
					txt_len++;
					if(txt_len == 80 - 9)
					{
						strcpy(content+x1, "\x1b[0m ...");
						content_len = line_width - 4;
						content_len = 80 - 4;
						line_width = 80;
						break;
					} 
				}
				if(nch == 'm' && state)
					state = 0;
			}
		}

		if(content_len==0 && is_boundary)
		{
			strcat(out_line_buff,is_boundary==1?"//":"\\\\");
			for(x1=0;x1<line_width-4;x1++)
				strcat(out_line_buff,"=");
			strcat(out_line_buff,is_boundary==1?"\\\\":"//");
			sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "%s", out_line_buff);

			free(content);
			free(out_line_buff);
			return;
		}
		else if(is_boundary)
		{
			int left_stars = (line_width - content_len)/2 - 1;
			int right_stars = line_width - content_len - 2 - left_stars;
			strcat(out_line_buff,is_boundary==1?"//":"\\\\");
			for(x1=0;x1<left_stars-2;x1++) strcat(out_line_buff,"=");
			sprintf(out_line_buff+strlen(out_line_buff),"%c[36m", CHAR_ESC);
			sprintf(out_line_buff+strlen(out_line_buff)," %s ", content);
			sprintf(out_line_buff+strlen(out_line_buff),"%c[0m", CHAR_ESC);
			for(x1=0;x1<right_stars-2;x1++) strcat(out_line_buff,"=");
			strcat(out_line_buff,is_boundary==1?"\\\\":"//");
			sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "%s", out_line_buff);

			free(content);
			free(out_line_buff);
			return;
		}

		int right_spaces, left_spaces;	
		if(is_center)
			left_spaces = (line_width - content_len)/2-2;
		else
			left_spaces = 1;

		right_spaces = line_width - 4 - content_len- left_spaces; 

		char spaces[81];
		memset(spaces , ' ', 80);
		spaces[0]='|';
		spaces[1]='|';
		spaces[80]=0;

		//sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,"||");
		
		//for(x1=0;x1<left_spaces;x1++) sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO," ");

		spaces[left_spaces+2] = 0;
		strcat(out_line_buff,spaces);

		if(is_R_code)
		{
			strcat(out_line_buff,content);
		}
		else
		{
			int col1w=-1;
			for(x1=0; content[x1]; x1++)
			{
				if(content[x1]==':')
				{
					col1w=x1;
					break;
				}
			}
			if(col1w>0 && col1w < content_len-1 && put_color_for_colon)
			{
				content[col1w+1]=0;
				strcat(out_line_buff,content);
				strcat(out_line_buff," ");
				sprintf(out_line_buff+strlen(out_line_buff),"%c[36m", CHAR_ESC);
				strcat(out_line_buff,content+col1w+2);
				sprintf(out_line_buff+strlen(out_line_buff),"%c[0m", CHAR_ESC);
			}
			else
				strcat(out_line_buff,content);
		}
	//	for(x1=0;x1<right_spaces - 1;x1++) sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO," ");
		
		memset(spaces , ' ', 80);
		spaces[79]='|';
		spaces[78]='|';
		
		right_spaces = max(1,right_spaces);
		sprintf(out_line_buff+strlen(out_line_buff)," %c[0m%s", CHAR_ESC , spaces + (78 - right_spaces + 1));
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, out_line_buff);
	}
	free(out_line_buff);
	free(content);
}




int show_summary(global_context_t * global_context)
{

#ifndef MAKE_STANDALONE
	// SAVE summary into files. This is only for Rsubread currently.
	char sumname[MAX_FILE_NAME_LENGTH+30];
	sprintf(sumname, "%s.summary", global_context->config.output_prefix);
	FILE * sumfp = fopen(sumname,"w");
	#ifdef __MINGW32__
	fprintf(sumfp, "Total_%s\t%I64u\n", global_context->input_reads.is_paired_end_reads?"fragments":"reads" , global_context -> all_processed_reads);
	#else
	fprintf(sumfp, "Total_%s\t%llu\n", global_context->input_reads.is_paired_end_reads?"fragments":"reads" , global_context -> all_processed_reads);
	#endif
	fprintf(sumfp, "Mapped_%s\t%u\n", global_context->input_reads.is_paired_end_reads?"fragments":"reads" , global_context -> all_mapped_reads);
	fprintf(sumfp, "Uniquely_mapped_%s\t%u\n", global_context->input_reads.is_paired_end_reads?"fragments":"reads" , global_context -> all_uniquely_mapped_reads);
	fprintf(sumfp, "Multi_mapping_%s\t%u\n", global_context->input_reads.is_paired_end_reads?"fragments":"reads" , global_context -> all_multimapping_reads);
	fprintf(sumfp, "Unmapped_%s\t%u\n", global_context->input_reads.is_paired_end_reads?"fragments":"reads" , global_context -> all_unmapped_reads);

	if(global_context->input_reads.is_paired_end_reads){
		#ifdef __MINGW32__
		fprintf(sumfp, "Properly_paired_fragments\t%I64u\n",global_context -> all_correct_PE_reads);
		#else
		fprintf(sumfp, "Properly_paired_fragments\t%llu\n",global_context -> all_correct_PE_reads);
		#endif
		fprintf(sumfp, "Singleton_fragments\t%u\n", global_context -> not_properly_pairs_only_one_end_mapped);
		fprintf(sumfp, "More_than_one_chr_fragments\t%u\n", global_context -> not_properly_pairs_different_chro);
		fprintf(sumfp, "Unexpected_strandness_fragments\t%u\n", global_context -> not_properly_different_strands);
		fprintf(sumfp, "Unexpected_template_length\t%u\n", global_context -> not_properly_pairs_TLEN_wrong);
		fprintf(sumfp, "Inversed_mapping\t%u\n", global_context -> not_properly_pairs_wrong_arrangement);
	}
	
	if(global_context->config.entry_program_name == CORE_PROGRAM_SUBJUNC && ( global_context -> config.prefer_donor_receptor_junctions || !(global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection)))
		fprintf(sumfp, "Junctions\t%u\n", global_context -> all_junctions);
	fprintf(sumfp, "Indels\t%u\n", global_context -> all_indels);


	fclose(sumfp);
#endif

	if(progress_report_callback)
	{
	        long long int all_reads_K = global_context -> all_processed_reads / 1000;
	        float mapped_reads_percentage = global_context -> all_mapped_reads * 1./global_context -> all_processed_reads;
	        if(global_context->input_reads.is_paired_end_reads) mapped_reads_percentage/=2;
	        progress_report_callback(10, 900000, (int) (miltime()-global_context->start_time));
	        progress_report_callback(10, 900010, (int) all_reads_K);
	        progress_report_callback(10, 900011, (int) (10000.*mapped_reads_percentage));
	}

	print_in_box(80,0,1,"  ");
	print_in_box(89,0,1,"  %c[36mCompleted successfully.%c[0m", CHAR_ESC, CHAR_ESC);
	print_in_box(80,0,1,"  ");
	print_in_box(80,2,1,"  ");
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "");
	print_in_box(80, 1,1,"  Summary");
	print_in_box(80, 0,1,"  ");

	#ifdef __MINGW32__
	if(global_context->input_reads.is_paired_end_reads)
		print_in_box(80, 0,0,"            Total fragments : %I64d" , global_context -> all_processed_reads);
	else
		print_in_box(80, 0,0,"                Total reads : %I64d" , global_context -> all_processed_reads);

	print_in_box(81, 0,0,"                     Mapped : %u (%.1f%%%%)", global_context -> all_mapped_reads,  global_context -> all_mapped_reads*100.0 / global_context -> all_processed_reads);
	print_in_box(80, 0,0,"            Uniquely mapped : %u",  global_context -> all_uniquely_mapped_reads);
	print_in_box(80, 0,0,"              Multi-mapping : %u",  global_context -> all_multimapping_reads);
	print_in_box(80, 0,1,"      ");
	print_in_box(80, 0,0,"                   Unmapped : %u",  global_context -> all_unmapped_reads);
	if(global_context->input_reads.is_paired_end_reads){
		print_in_box(80, 0,1,"      ");
		print_in_box(80, 0,0,"            Properly paired : %I64d", global_context -> all_correct_PE_reads);
		print_in_box(80, 0,0,"        Not properly paired : %I64d", global_context -> all_mapped_reads -  global_context -> all_correct_PE_reads);
		print_in_box(80, 0,0,"                  Singleton : %u", global_context -> not_properly_pairs_only_one_end_mapped);
		print_in_box(80, 0,0,"                   Chimeric : %u", global_context -> not_properly_pairs_different_chro);
		print_in_box(80, 0,0,"      Unexpected strandness : %u", global_context -> not_properly_different_strands);
	 	print_in_box(80, 0,0," Unexpected fragment length : %u", global_context -> not_properly_pairs_TLEN_wrong);
	 	print_in_box(80, 0,0,"      Unexpected read order : %u", global_context -> not_properly_pairs_wrong_arrangement);
	}

	print_in_box(80, 0,1,"      ");

	if(global_context->config.output_prefix[0])
	{
	        if(global_context->config.entry_program_name == CORE_PROGRAM_SUBJUNC && ( global_context -> config.prefer_donor_receptor_junctions || !(global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection)))
	                print_in_box(80, 0,0,"                  Junctions : %u", global_context -> all_junctions);
	        if((global_context-> config.do_fusion_detection || global_context-> config.do_long_del_detection))
	                print_in_box(80, 0,0,"                    Fusions : %u", global_context -> all_fusions);
	        print_in_box(80, 0,0,"                     Indels : %u", global_context -> all_indels);
	}
	#else
	if(global_context->input_reads.is_paired_end_reads)
		print_in_box(80, 0,0,"            Total fragments : %'llu" , global_context -> all_processed_reads);
	else
		print_in_box(80, 0,0,"                Total reads : %'llu" , global_context -> all_processed_reads);

	print_in_box(81, 0,0,"                     Mapped : %'u (%.1f%%%%)", global_context -> all_mapped_reads,  global_context -> all_mapped_reads*100.0 / global_context -> all_processed_reads);
	print_in_box(80, 0,0,"            Uniquely mapped : %'u",  global_context -> all_uniquely_mapped_reads);
	print_in_box(80, 0,0,"              Multi-mapping : %'u",  global_context -> all_multimapping_reads);
	print_in_box(80, 0,1,"      ");
	print_in_box(80, 0,0,"                   Unmapped : %'u",  global_context -> all_unmapped_reads);
	if(global_context->input_reads.is_paired_end_reads){
		print_in_box(80, 0,1,"      ");
		print_in_box(80, 0,0,"            Properly paired : %'llu", global_context -> all_correct_PE_reads);
		print_in_box(80, 0,0,"        Not properly paired : %'llu", global_context -> all_mapped_reads -  global_context -> all_correct_PE_reads);
		print_in_box(80, 0,0,"                  Singleton : %'u", global_context -> not_properly_pairs_only_one_end_mapped);
		print_in_box(80, 0,0,"                   Chimeric : %'u", global_context -> not_properly_pairs_different_chro);
		print_in_box(80, 0,0,"      Unexpected strandness : %'u", global_context -> not_properly_different_strands);
	 	print_in_box(80, 0,0," Unexpected fragment length : %'u", global_context -> not_properly_pairs_TLEN_wrong);
	 	print_in_box(80, 0,0,"      Unexpected read order : %'u", global_context -> not_properly_pairs_wrong_arrangement);
	}

	print_in_box(80, 0,1,"      ");

	if(global_context->config.output_prefix[0])
	{
	        if(global_context->config.entry_program_name == CORE_PROGRAM_SUBJUNC && ( global_context -> config.prefer_donor_receptor_junctions || !(global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection)))
	                print_in_box(80, 0,0,"                  Junctions : %'u", global_context -> all_junctions);
	        if((global_context-> config.do_fusion_detection || global_context-> config.do_long_del_detection))
	                print_in_box(80, 0,0,"                    Fusions : %'u", global_context -> all_fusions);
	        print_in_box(80, 0,0,"                     Indels : %'u", global_context -> all_indels);
	}

	#endif

	if(global_context -> is_phred_warning)
	{
		print_in_box(80, 0,1,"      ");
		print_in_box(80,0,0,"                    WARNING : Phred offset (%d) incorrect?", global_context->config.phred_score_format == FASTQ_PHRED33?33:64);
	}
	print_in_box(80, 0,1,"      ");
	print_in_box(80, 0,0,"               Running time : %.1f minutes", (miltime()-global_context->start_time)*1./60);

	if( global_context->input_reads.is_paired_end_reads && global_context -> config.reported_multi_best_reads<2 && global_context -> expected_TLEN_read_numbers < READPAIRS_FOR_CALC_EXPT_TLEN){
	    print_in_box(80, 0,1,"      ");
		print_in_box(80,0,0,"  NOTE : No enough read-pairs to derive expected fragment length.");
	}
/*
	print_in_box(80, 0,0,"    Running time 0 : %.2f minutes", global_context->timecost_load_index/60);
	print_in_box(80, 0,0,"    Running time 1 : %.2f minutes", global_context->timecost_voting/60);
	print_in_box(80, 0,0,"    Running time 2 : %.2f minutes", global_context->timecost_before_realign/60);
	print_in_box(80, 0,0,"    Running time 3 : %.2f minutes", global_context->timecost_for_realign/60);
*/
	print_in_box(80, 0,1,"");
	print_in_box(80, 2,1,"");
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "");

	return 0;
}



void show_progress(global_context_t * global_context, thread_context_t * thread_context, unsigned int current_read_no, int task)
{

	// Read_chunk_start is the file_offset of the very first read in the entire file.
	// current_circle_start_position_file1 is the file_offset of the first read in this 5-million read chunk (or whatever the chunk size is)


	if(global_context->config.scRNA_input_mode){
		if(task == 10){
			char minchr[10];
			float min_value = (miltime() - global_context -> start_time)*1./60;
			if(min_value < 9.91)
				sprintf(minchr, "%.01f", min_value);
			else sprintf(minchr, "% 3d", (int)min_value);
			print_in_box(80,0,0,"   processed % 3d million input reads in %s minutes", current_read_no/1000000, minchr);
		}
		//SUBREADprintf("TASK=%d, RNO=%u\n", task, current_read_no);
		return;
	}

	if(thread_context&&thread_context->thread_id)
	{
		SUBREADputs("show_progress can only be called by thread#0\n");
		return;
	}

	gene_input_t * ginp1 = &global_context->input_reads.first_read_file;

	unsigned long long ginp1_file_pos = geinput_file_offset(ginp1);

	if(task == STEP_VOTING)
	{
		unsigned long long real_read_number = global_context -> all_processed_reads + current_read_no;
		if(real_read_number>1000)
			global_context -> input_reads . avg_read_length = (ginp1_file_pos - ginp1 -> read_chunk_start) * 1./real_read_number ;
	}

	unsigned long long total_file_size = global_context -> input_reads.first_read_file_size;
	unsigned long long guessed_all_reads = total_file_size / global_context -> input_reads . avg_read_length;

	unsigned long long current_block_start_file_offset = global_context -> current_circle_start_abs_offset_file1;

	unsigned long long guessed_this_chunk_reads = (total_file_size - current_block_start_file_offset) / global_context -> input_reads . avg_read_length ;
	if(guessed_this_chunk_reads > global_context ->config.reads_per_chunk) guessed_this_chunk_reads = global_context ->config.reads_per_chunk;

	unsigned long long guessed_all_reads_before_this_chunk = current_block_start_file_offset / global_context -> input_reads . avg_read_length ;
	unsigned long long reads_finished_in_this_chunk = (ginp1_file_pos - current_block_start_file_offset) / global_context -> input_reads . avg_read_length;

	int is_thred_step_running = global_context->config.is_third_iteration_running ? 1:0;

	if(task != STEP_VOTING)
		reads_finished_in_this_chunk = (ginp1_file_pos - current_block_start_file_offset) / global_context -> input_reads . avg_read_length;

	unsigned long long finished_steps = guessed_all_reads_before_this_chunk * (global_context -> index_block_number * 2 + 1 + is_thred_step_running);
	

	// add steps for voting
	if(task == STEP_VOTING)
		finished_steps += guessed_this_chunk_reads * global_context -> current_index_block_number * 2;
	else if(task == STEP_ITERATION_TWO)
		finished_steps += guessed_this_chunk_reads * global_context -> index_block_number * 2;
	else if(task > STEP_ITERATION_TWO)
		finished_steps += guessed_this_chunk_reads *(global_context -> index_block_number * 2+1);

	if(task == STEP_VOTING)
		finished_steps += reads_finished_in_this_chunk*2;
	else	finished_steps += reads_finished_in_this_chunk;

	unsigned long long guessed_all_steps = guessed_all_reads * (global_context -> index_block_number *2 + 1 + is_thred_step_running);

	float finished_rate = finished_steps*1./guessed_all_steps;
	float reads_per_second = 0;

	if(task == STEP_VOTING)
		reads_per_second = finished_rate *1.*guessed_all_reads / (miltime() - global_context -> align_start_time);
	else
		reads_per_second = finished_rate *1.*guessed_all_reads / (miltime() - global_context -> start_time);

	if(current_read_no>1000 && !progress_report_callback)
	{
		char minchr[10];
		float min_value = (miltime() - global_context -> start_time)*1./60;
		if(min_value < 9.91)
			sprintf(minchr, "%.01f", min_value);
		else sprintf(minchr, "% 3d", (int)min_value);

	//	print_in_box(81,0,0, "%4d%%%% completed, %s mins elapsed, total=%dk %s, rate=%2.1fk/s\r", (int)(finished_rate*100), minchr,(int)(guessed_all_reads*1./1000), global_context -> input_reads.is_paired_end_reads?"frags":"reads", reads_per_second/1000);
		print_in_box(81,0,0, "%4d%%%% completed, %s mins elapsed, rate=%2.1fk %s per second\r", (int)(finished_rate*100), minchr, reads_per_second/1000, global_context -> input_reads.is_paired_end_reads?"fragments":"reads");
	}

	if(progress_report_callback)
	{
		progress_report_callback(10 ,task,(int)(10000*finished_rate));
		progress_report_callback(20 ,task,(int)(guessed_all_reads/1000));
	}
}


/*
int Xmain(int argc , char ** argv);

int main(int argc, char ** argv)
{
	int xk1=0;
	for(xk1=0; xk1<4; xk1++)
		Xmain(argc ,  argv);
	return 0;
}
*/


int parse_opts_core(int argc , char ** argv, global_context_t * global_context)
{
	int c;
	int option_index = 0;	

	optind = 1;
	opterr = 1;
	optopt = 63;

	while ((c = getopt_long (argc, argv, "ExsS:L:AHd:D:n:m:p:P:R:r:i:l:o:T:Q:I:t:B:b:Q:FcuUfM?", long_options, &option_index)) != -1)
	{
		switch(c)
		{
			case 'Q':
				global_context->config.multi_best_reads = atoi(optarg); 
				break;
			case 'H':
				global_context->config.use_hamming_distance_break_ties = 1;
				break;
			case 's':
				global_context->config.downscale_mapping_quality = 1;
				break;
			case 'M':
				global_context->config.do_big_margin_filtering_for_reads = 1;
				global_context->config.report_multi_mapping_reads = 0;
				break;

			case 'A':
				global_context->config.report_sam_file = 0;
				break;
			case 'E':
				global_context->config.max_mismatch_exonic_reads = 200;
				global_context->config.max_mismatch_junction_reads = 200;

				break;
			case 'f':
				global_context->config.max_mismatch_exonic_reads = 200;
				global_context->config.max_mismatch_junction_reads = 200;
				global_context->config.do_fusion_detection = 1;
				global_context->config.minimum_subread_for_first_read = 1;
				global_context->config.minimum_subread_for_second_read = 1;
				global_context->config.total_subreads = 28;
				global_context->config.report_no_unpaired_reads = 0;
				global_context->config.limited_tree_scan = 0;
				global_context->config.use_hamming_distance_in_exon = 1;
				break;
			case 'x':
				global_context->config.max_mismatch_exonic_reads = 10;
				global_context->config.max_mismatch_junction_reads = 1;
				global_context->config.ambiguous_mapping_tolerance = 39;
				global_context->config.extending_search_indels = 0;

				global_context->config.do_breakpoint_detection = 1;
				global_context->config.total_subreads = 14;
				global_context->config.minimum_subread_for_first_read = 3;
				global_context->config.minimum_subread_for_second_read = 1;
				global_context->config.high_quality_base_threshold = 990000;
				global_context->config.do_big_margin_filtering_for_junctions = 1;
				global_context->config.report_no_unpaired_reads = 0;
				global_context->config.limited_tree_scan = 1;
				global_context->config.use_hamming_distance_in_exon = 0;
				break;
			case 'S':
				global_context->config.is_first_read_reversed = optarg[0]=='r'?1:0;
				global_context->config.is_second_read_reversed = optarg[0]=='f'?0:1;
				break;
			case 'U':
				global_context->config.report_no_unpaired_reads = 1;
				break;
			case 'u':
				global_context->config.report_multi_mapping_reads = 0;
				break;
			case 'b':
				global_context->config.is_methylation_reads = 1;
				break;
			case 'D':
				global_context->config.maximum_pair_distance = atoi(optarg);
				break;
			case 'd':
				global_context->config.minimum_pair_distance = atoi(optarg);
				break;
			case 'n':
				global_context->config.total_subreads = atoi(optarg);
				break;
			case 'm':
				global_context->config.minimum_subread_for_first_read = atoi(optarg);
				break;
			case 'T':
				global_context->config.all_threads = atoi(optarg);
				if(global_context->config.all_threads <1) global_context->config.all_threads = 1;
				if(global_context->config.all_threads >32) global_context->config.all_threads = 32;

				break;
			case 'r':
				strncpy(global_context->config.first_read_file, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'R':
				global_context->input_reads.is_paired_end_reads = 1;
				strncpy(global_context->config.second_read_file, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'i':
				strncpy(global_context->config.index_prefix, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'o':
				strncpy(global_context->config.output_prefix, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'I':
				global_context->config.max_indel_length = atoi(optarg);

				if(global_context->config.max_indel_length <0)global_context->config.max_indel_length =0;
				if(global_context->config.max_indel_length > MAX_INSERTION_LENGTH)global_context->config.max_indel_length = MAX_INSERTION_LENGTH;
				if(global_context->config.max_indel_length > 16)
				{
					global_context->config.reassembly_subread_length = 12;
					global_context->config.reassembly_window_multiplex = 3;
					global_context->config.reassembly_start_read_number = 5;
					global_context->config.reassembly_tolerable_voting = 0;
					global_context->config.reassembly_window_alleles = 2;
					global_context->config.reassembly_key_length = 28;

					global_context->config.is_third_iteration_running = 1;
					global_context->config.max_mismatch_exonic_reads = 2;
					global_context->config.max_mismatch_junction_reads = 2;
					global_context->config.total_subreads = 28;
					global_context->config.do_big_margin_filtering_for_reads = 1;

					global_context->config.do_superlong_indel_detection = 0;
				}
				break;
			case 'P':
				if (optarg[0]=='3')
					global_context->config.phred_score_format = FASTQ_PHRED33;
				else
					global_context->config.phred_score_format = FASTQ_PHRED64;
				break;
			case 'p':
				global_context->config.minimum_subread_for_second_read = atoi(optarg);
				break;
			case 't':
				sprintf(global_context->config.temp_file_prefix, "%s/core-temp-sum-%06u-%05u", optarg, getpid(), myrand_rand()
				);
				break;
			case 'F':
				global_context->config.is_second_iteration_running = 0;
				global_context->config.report_sam_file = 0;
				break;
			case 'B':
				strcpy(global_context->config.exon_annotation_file, optarg);
				break;
			case 'c':
				global_context->config.space_type = GENE_SPACE_COLOR; 
				break;
				
			case 0:
				//if(strcmp("memory-optimisation",long_options[option_index].name)==0)
				//	memory_optimisation = atoi(optarg);
				break;
			case '?':
			default:
				return -1 ;
		}
	}


	return 0;
}

#ifdef MAKE_STANDALONE
int main_back(int argc , char ** argv)
{
	progress_report_callback = NULL;
#else
int subread_core_main(int argc , char ** argv)
{
	if(progress_report_callback) progress_report_callback(0,0,1);
#endif

	return core_main(argc, argv, parse_opts_core);
}


int check_configuration(global_context_t * global_context)
{
	int expected_type = FILE_TYPE_FAST_;
	if(global_context -> config.is_SAM_file_input && global_context -> config.is_BAM_input)
		expected_type = FILE_TYPE_BAM;
	else if(global_context -> config.is_SAM_file_input && !global_context -> config.is_BAM_input)
		expected_type = FILE_TYPE_SAM;
	else if(global_context -> config.is_gzip_fastq)
		expected_type = FILE_TYPE_GZIP_FAST_;
	
	if(global_context -> config.max_indel_length > 16)
		warning_file_limit();

	int wret = 0, wret2 = 0;
	if( global_context ->config.scRNA_input_mode == 0 )  warning_file_type(global_context -> config.first_read_file, expected_type);
	if(global_context -> config.second_read_file[0])
	{
		if(expected_type==FILE_TYPE_FAST_ || expected_type==FILE_TYPE_GZIP_FAST_)
			wret2 = warning_file_type(global_context -> config.second_read_file, expected_type);
		else
			print_in_box(80,0,0,"Only one input SAM or BAM file is needed. The second input file is ignored.");
	}

	if(wret == -1 || wret2 == -1){
		return -1;
	}

	if(4 == sizeof(void *) && global_context->config.sort_reads_by_coordinates){
		SUBREADputs("ERROR: the sort by coordinate function only works on 64-bit computers.");
		return -1;
	}

	if(global_context->config.is_input_read_order_required && global_context->config.sort_reads_by_coordinates){
		SUBREADputs("ERROR: you shouldn't specify keep input order and sort by coordinate at same time.");
		return -1;
	}

	if((global_context -> config.is_BAM_output == 0 || global_context->config.output_prefix[0]==0) && global_context->config.sort_reads_by_coordinates){
		if(global_context -> config.is_BAM_output==0)SUBREADputs("ERROR: SAM output doesn't support read sorting by coordinates.");
		else SUBREADputs("ERROR: STDOUT output doesn't support read sorting by coordinates.");
		return -1;
	}

	return 0;

}

unsigned long long myrand_seed = 0;
unsigned long long myrand_seed2 = 0;

void myrand_srand(unsigned long long seed){

	// when it is in R, the seed is NOT used but the R random numbers are used as the seed.
	myrand_seed = 0;
	myrand_seed2 = 0;

	#ifndef MAKE_STANDALONE
	GetRNGstate();
	seed = 0;
	seed = (seed << 16)+ (int)(unif_rand()*65521);
	seed = (seed << 16)+ (int)(unif_rand()*65521);
	seed = (seed << 16)+ (int)(unif_rand()*65521);
	seed = (seed << 16)+ (int)(unif_rand()*65521);

	PutRNGstate();	
	#endif
	myrand_seed = seed;
	myrand_rand();
	myrand_rand();
	myrand_rand();
	myrand_rand();
	myrand_rand();
	myrand_rand();
	myrand_rand();
	myrand_rand();
	myrand_rand();
}

#define THE_966666773323RD_PRIME 28962406029617llu
int myrand_rand(){
	//if(myrand_seed % 3133LLU == 0) myrand_srand(0);

	myrand_seed ^= (myrand_seed % 858173 + myrand_seed % 104729);
	myrand_seed2 += myrand_seed;
	myrand_seed ^= myrand_seed2 << 13;
	return (int)(((myrand_seed^myrand_seed2) %THE_966666773323RD_PRIME) % (1LLU+RAND_MAX));
}


int core_main(int argc , char ** argv, int (parse_opts (int , char **, global_context_t * )))
{
	struct timeval xtime; 
	gettimeofday(&xtime,NULL);
	myrand_srand(time(NULL)^xtime.tv_usec);

	global_context_t * global_context;
	global_context = (global_context_t*)malloc(sizeof(global_context_t));
	memset(global_context, 0, sizeof(global_context_t));
	init_global_context(global_context);


	int ret = parse_opts(argc , argv, global_context);
	init_core_temp_path(global_context);
	//global_context->config.reads_per_chunk = 200*1024;

	if(global_context->config.max_indel_length > 20 && !global_context->input_reads.is_paired_end_reads)
	{
		global_context->config.total_subreads = 28;
		global_context->config.reassembly_start_read_number = 3;
		global_context->config.do_superlong_indel_detection = 1;
	}

	if(global_context->config.fast_run){
		global_context -> config.top_scores = 1;
		global_context -> config.max_vote_combinations = 1;
		global_context -> config.max_vote_simples = 1;
		global_context -> config.multi_best_reads = 1;
	}

	ret = ret || print_configuration(global_context);
	ret = ret || check_configuration(global_context);
	ret = ret || load_global_context(global_context);
	ret = ret || init_modules(global_context);
	ret = ret || read_chunk_circles(global_context);
	ret = ret || write_final_results(global_context);
	ret = ret || destroy_modules(global_context);
	ret = ret || destroy_global_context(global_context);
	ret = ret || show_summary(global_context);

	free(global_context);

	return ret;
}

// the new file name is written into fname then.

int convert_BAM_to_SAM(global_context_t * global_context, char * fname, int is_bam)
{
	char temp_file_name[MAX_FILE_NAME_LENGTH+80], *fline=malloc(3000), tmp_readname[MAX_READ_NAME_LEN];
	short tmp_flags;
	SamBam_FILE * sambam_reader;

	char * env_no_sort = getenv("SUBREAD_DO_NOT_CHECK_INPUT");
	if(env_no_sort){
		global_context->input_reads.is_paired_end_reads = 1;
		SUBREADprintf("\nWARNING: The SAM input file was assumed to be sorted by name.\nENV: SUBREAD_DO_NOT_CHECK_INPUT\n");
		return 0;
	}

	int is_file_sorted = 1;
	unsigned long long int read_no = 0;
	sambam_reader = SamBam_fopen(fname, is_bam?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM);
	if(!sambam_reader)
	{
		SUBREADprintf("ERROR: Unable to open file '%s'. The file is not in the specified format.\n", fname);
		free(fline);
		return -1;
	}
	while(1)
	{
		char * is_ret = SamBam_fgets(sambam_reader, fline, 2999, 1);
		if(!is_ret) break;
		if(sambam_reader ->is_bam_broken ) {
			SUBREADputs("ERROR: the BAM format is broken.");
			return -1;
		}
		if(fline[0]=='@')continue;
		if(is_SAM_unsorted(fline, tmp_readname, &tmp_flags , read_no)){
			if(tmp_flags & 1) global_context->input_reads.is_paired_end_reads = 1;
			is_file_sorted = 0;
			read_no++;
			break;
		}
		read_no++;
		if(tmp_flags & 1) global_context->input_reads.is_paired_end_reads = 1;
		else global_context->input_reads.is_paired_end_reads = 0;

		if(! global_context->input_reads.is_paired_end_reads) break;
	}

	if(read_no<1)
	{
		SUBREADprintf("ERROR: unable to find any reads from file '%s'. The file is not in the specified format or is empty.\n", fname);
		return -1;
	}

	SamBam_fclose(sambam_reader);
	print_in_box(80,0,0,"The input %s file contains %s-end reads.", is_bam?"BAM":"SAM",  global_context->input_reads.is_paired_end_reads?"paired":"single");
	if(!is_file_sorted)
		print_in_box(80,0,0,"The input %s file is unsorted. Reorder it...", is_bam?"BAM":"SAM");
	else if(is_bam)
		print_in_box(80,0,0,"Convert the input BAM file...");

	int disk_is_full = 0;
	if(is_bam || (global_context->input_reads.is_paired_end_reads && !is_file_sorted))
	{
		sprintf(temp_file_name, "%s.sam", global_context->config.temp_file_prefix);
		sambam_reader = SamBam_fopen(fname, is_bam?SAMBAM_FILE_BAM: SAMBAM_FILE_SAM);
		if(!sambam_reader){
			SUBREADprintf("Unable to open %s.\n", fname);
			return -1;
		}
		FILE * sam_fp = NULL;
		SAM_sort_writer writer;
		int writer_opened = 0;

		if(is_file_sorted) sam_fp = f_subr_open(temp_file_name,"w");
		else writer_opened = sort_SAM_create(&writer, temp_file_name, NULL);

		if((is_file_sorted && !sam_fp) || (writer_opened && !is_file_sorted)){
			SUBREADprintf("Failed to write to the directory. You may not have permission to write to this directory or the disk is full.\n");
			return -1;
		}

		while(1)
		{
			char * is_ret = SamBam_fgets(sambam_reader, fline, 2999, 1);
			if(!is_ret) break;
			if(is_file_sorted){
				int wlen = fputs(fline, sam_fp);
				if(wlen <0){
					SUBREADprintf("ERROR: unable to write into the temporary SAM file. Please check the disk space in the output directory.\n");
					disk_is_full = 1;
					break;
				}
			}else{
				int ret = sort_SAM_add_line(&writer, fline, strlen(fline));
				if(ret== -1) {
					print_in_box(80,0,0,"ERROR: read name is too long; check input format.");
					break;
				}
				if(ret== -2) {
					SUBREADprintf("ERROR: unable to write into the temporary SAM file. Please check the disk space in the output directory.\n");
					disk_is_full = 1;
					break;
				}
			}
		}

		if(is_file_sorted) 
			fclose(sam_fp);
		else{
			int ret = sort_SAM_finalise(&writer);
			if(writer.unpaired_reads)
				#ifdef __MINGW32__
				print_in_box(80,0,0,"%I64d single-end mapped reads in reordering.", writer.unpaired_reads);
				#else
				print_in_box(80,0,0,"%lld single-end mapped reads in reordering.", writer.unpaired_reads);
				#endif
			if(ret) {
				disk_is_full = 1;
				SUBREADprintf("ERROR: unable to create the temporary file. Please check the disk space in the output directory.\n");
			}
		}

		SamBam_fclose(sambam_reader);
		strcpy(fname, temp_file_name);
		global_context -> will_remove_input_file = 1;
	}

	free(fline);
	return disk_is_full;
}

int convert_GZ_to_FQ(global_context_t * global_context, char * fname, int half_n)
{
	int is_OK = 0;
	char temp_file_name[MAX_FILE_NAME_LENGTH+30];
	char * linebuff=malloc(3001);
	gzFile rawfp = gzopen(fname, "r");
	
	if(rawfp)
	{
		print_in_box(80,0,0,"Decompress %s...", fname);
		sprintf(temp_file_name, "%s-%d.fq", global_context->config.temp_file_prefix, half_n);
		FILE * outfp = fopen(temp_file_name, "w");
		if(outfp)
		{
			while(1)
			{
				char * bufr =gzgets(rawfp, linebuff, 3000);
				if(!bufr) break;
				fputs(bufr, outfp);
			}
			is_OK = 1;
			fclose(outfp);
		}
		else{
			SUBREADprintf("Unable to create temporary file '%s'\nPlease run the program in a directory where you have the privilege to create files.\n", temp_file_name);
		}

		gzclose(rawfp);
	}

	strcpy(fname, temp_file_name);
	global_context -> will_remove_input_file |= (1<< (half_n-1));

	return !is_OK;
}

int core_geinput_open(global_context_t * global_context, gene_input_t * fp, int half_number)
{
	char *fname;
	if(global_context->config.is_SAM_file_input) {
		fname = global_context ->config.first_read_file;
		if(half_number == 1)
			if(convert_BAM_to_SAM(global_context, global_context ->config.first_read_file, global_context ->config.is_BAM_input)) return -1;
		if(!global_context->input_reads.is_paired_end_reads) half_number=0;
		return geinput_open_sam(fname, fp, half_number);
	} else {
		int rv = -1;
		//SUBREADprintf("SCRNA_MODE=%d\n", global_context->config.scRNA_input_mode );
		if(global_context -> config.is_gzip_fastq)
			if(convert_GZ_to_FQ(global_context, (half_number==2)? global_context ->config.second_read_file : global_context ->config.first_read_file, half_number)) return -1;
		fname = (half_number == 2)?global_context -> config.second_read_file:global_context -> config.first_read_file;
		if(global_context->config.scRNA_input_mode == GENE_INPUT_BCL)
			rv = geinput_open_bcl(fname , fp, global_context -> config.reads_per_chunk, global_context -> config.all_threads );
		else if(global_context->config.scRNA_input_mode == GENE_INPUT_SCRNA_FASTQ)
			rv = geinput_open_scRNA_fqs(fname , fp, global_context -> config.reads_per_chunk, global_context -> config.all_threads );
		else if(global_context->config.scRNA_input_mode == GENE_INPUT_SCRNA_BAM)
			rv = geinput_open_scRNA_BAM(fname , fp, global_context -> config.reads_per_chunk, global_context -> config.all_threads );
		else
			rv = geinput_open(fname, fp);

		if(global_context->input_reads.is_paired_end_reads && global_context->config. scRNA_input_mode){
			SUBREADprintf("ERROR: No paired-end input is allowed on scRNA mode.\n");
			return -1;
		}
		return rv;
	}
}


int fetch_next_read_pair(global_context_t * global_context, thread_context_t * thread_context, gene_input_t* ginp1, gene_input_t* ginp2, int *read_len_1, int *read_len_2, char * read_name_1, char * read_name_2, char * read_text_1, char * read_text_2, char * qual_text_1, char *qual_text_2, int remove_color_head, subread_read_number_t * read_no_in_chunk)
{
	int rl1=0, rl2=0;
	int is_second_R1, is_second_R2;
	subread_read_number_t this_number = -1;

	geinput_preload_buffer(ginp1, &global_context -> input_reads.input_lock);
	if(ginp2) geinput_preload_buffer(ginp2, &global_context -> input_reads.input_lock);

	subread_lock_occupy(&global_context -> input_reads.input_lock); 
	if(global_context -> running_processed_reads_in_chunk < global_context -> config.reads_per_chunk)
	{
		do
		{
			is_second_R1 = 0; is_second_R2 = 0;
			rl1 = geinput_next_read_trim(ginp1, read_name_1, read_text_1 , qual_text_1, global_context->config.read_trim_5, global_context->config.read_trim_3, &is_second_R1);
			//if(global_context -> running_processed_reads_in_chunk < 100)SUBREADprintf("GETRNAMES %s LEN=%d\n", read_name_1, rl1);
			if(global_context->config.space_type == GENE_SPACE_COLOR && remove_color_head)
			{
				if(isalpha(read_text_1[0]))
				{
					int xk1;
					for(xk1=2; read_text_1[xk1]; xk1++)
						read_text_1[xk1-2]=read_text_1[xk1];
					read_text_1[xk1-2]=0;
				}
			}

			if(ginp2)
			{
				rl2 = geinput_next_read_trim(ginp2, read_name_2, read_text_2 , qual_text_2, global_context->config.read_trim_5, global_context->config.read_trim_3, &is_second_R2);
				if(global_context->config.space_type == GENE_SPACE_COLOR && remove_color_head)
				{
					if(isalpha(read_text_2[0]))
					{
						int xk1;
						for(xk1=2; read_text_2[xk1]; xk1++)
							read_text_2[xk1-2]=read_text_2[xk1];
						read_text_2[xk1-2]=0;
					}
				}
			}
			if(rl1 <= 0 || (rl2 <= 0 && ginp2)) break;
		} while(is_second_R1||is_second_R2) ;


		//SUBREADprintf("IMBRL1=%d\tRL2=%d\n", rl1, rl2);
		//
		if(rl1 >0 || (rl2 > 0 && ginp2)){
			this_number = global_context -> running_processed_reads_in_chunk;
			global_context -> running_processed_reads_in_chunk ++;
		}
	}
	subread_lock_release(&global_context -> input_reads.input_lock); 

	if( global_context->config.space_type == GENE_SPACE_COLOR) {
		rl1-=1;rl2-=1;
	}

	if(ginp2 && rl1 * rl2 <=0 && (rl1>0 || rl2>0)){
		if(!global_context-> input_reads.is_internal_error)
			SUBREADprintf("\nERROR: two input files have different amounts of reads.\n\n");
		global_context-> input_reads.is_internal_error = 1;
		*read_no_in_chunk = -1;
		return 1;
	} else if(rl1>0 && (rl2>0 || !ginp2) && this_number>=0) {
		if(global_context->config.is_first_read_reversed)
		{
			reverse_read(read_text_1, rl1, global_context->config.space_type);
			if(qual_text_1)
				reverse_quality(qual_text_1, rl1);
		}

		if(ginp2 && global_context->config.is_second_read_reversed)
		{
			reverse_read(read_text_2, rl2, global_context->config.space_type);
			if(qual_text_2)
				reverse_quality(qual_text_2, rl2);
		}

		*read_no_in_chunk = this_number;
		*read_len_1 = rl1;
		if(ginp2)
			*read_len_2 = rl2;
		return 0;
	} else {
		// normally finished
		*read_no_in_chunk = -1;
		return 1;
	}
}

int write_final_results(global_context_t * context)
{

	if((context ->  config.do_fusion_detection || context ->  config.do_long_del_detection) && context -> config.do_structural_variance_detection)
		finalise_structural_variances(context);

	if(context -> config.output_prefix[0] && 0==context -> input_reads.is_internal_error && !(context -> config.is_BAM_output && context -> output_bam_writer -> is_internal_error) ) {
		write_indel_final_results(context);

		if(context -> config.entry_program_name == CORE_PROGRAM_SUBJUNC && (context -> config.prefer_donor_receptor_junctions||!(context ->  config.do_fusion_detection || context ->  config.do_long_del_detection)))
			write_junction_final_results(context);

		if((context ->  config.do_fusion_detection || context ->  config.do_long_del_detection))
			write_fusion_final_results(context);
	}
	
	return 0;
}

int get_soft_clipping_length(char* CIGAR)
{
	int nch;
	int cigar_cursor, tmp_int = 0;

	for(cigar_cursor = 0; (nch = CIGAR[cigar_cursor]) > 0;cigar_cursor++)
	{
		if(isdigit(nch)) tmp_int = tmp_int*10+(nch-'0');
		else
		{
			if(nch=='S') return tmp_int;
			return 0;
		}
	}
	return 0;
}


#define CIGAR_PERFECT_SECTIONS 12

typedef struct{
	char current_cigar_decompress[CORE_MAX_CIGAR_STR_LEN + 1];
	char cigar [CORE_MAX_CIGAR_STR_LEN];

	unsigned short chimeric_sections;
	unsigned int out_poses[CIGAR_PERFECT_SECTIONS];
	short out_lens[CIGAR_PERFECT_SECTIONS];
	char out_cigars[CIGAR_PERFECT_SECTIONS][60];
	char out_strands[CIGAR_PERFECT_SECTIONS];

	char additional_information[CORE_ADDITIONAL_INFO_LENGTH + 1];
	mapping_result_t * raw_result;

	unsigned int linear_position;
	short soft_clipping_movements;
	char * chro;
	int offset;
	int strand;
	int is_first_section_jumpped;

	int mapping_quality;
	int is_NM_appied;
}subread_output_tmp_t;

typedef struct{
	long long int *PE_distance;
	char * out_cigar_buffer[CIGAR_PERFECT_SECTIONS];
	subread_output_tmp_t *r1;
	subread_output_tmp_t *r2;
	subread_output_tmp_t ** out_pairs;
	mapping_result_t ** out_raws;

} subread_output_context_t;

void init_output_context(global_context_t * global_context, subread_output_context_t * out_context)
{
	int xk1;
	memset(out_context, 0, sizeof(subread_output_context_t));

	out_context -> r1 = malloc(sizeof(subread_output_tmp_t));
	for(xk1=0;xk1<CIGAR_PERFECT_SECTIONS;xk1++)
		out_context -> out_cigar_buffer[xk1] = malloc(60);

	out_context -> out_pairs = malloc(sizeof( subread_output_context_t *) * global_context->config.multi_best_reads * 2);
	out_context -> out_raws = malloc(sizeof( mapping_result_t * *) * global_context->config.multi_best_reads * 2);

	if(global_context -> input_reads.is_paired_end_reads)
	{
		out_context -> PE_distance = malloc(sizeof(long long) * global_context->config.multi_best_reads);
		out_context -> r2 = malloc(sizeof(subread_output_tmp_t));
	} else {
		out_context -> PE_distance = NULL;
		out_context -> r2 = NULL;
	}
}

void destroy_output_context(global_context_t * global_context , subread_output_context_t * out_context)
{
	int xk1;
	for(xk1=0;xk1<CIGAR_PERFECT_SECTIONS;xk1++)
		free(out_context -> out_cigar_buffer[xk1]);

	free(out_context -> out_raws);
	free(out_context -> out_pairs);
	free(out_context -> r1 );
	if(global_context -> input_reads.is_paired_end_reads)
	{
		free(out_context -> r2);
		free(out_context -> PE_distance);
	}
}

int locate_current_value_index(global_context_t * global_context, thread_context_t * thread_context, mapping_result_t * result, int rlen);
int calc_edit_dist(global_context_t * global_context, mapping_result_t * current_result, char * cigar , unsigned int pos ,  char * read_text, int all_mm)
{
	unsigned int cigar_cursor = 0, tmpi=0;

	while(1)
	{
		char nch = cigar[cigar_cursor++];
		if(!nch)break;

		if('0'<=nch && '9' >=nch)
		{
			tmpi = tmpi*10+nch-'0';
		}
		else{
			if(nch == 'I' || nch == 'D')
			{
				all_mm+=tmpi;
			}
	
			tmpi = 0;
		}
	}
	return all_mm;
}

unsigned int move_to_read_head(unsigned int tailpos, char * cigar){
	int cigar_i = 0, nch;
	unsigned int tmpi = 0;
	while(0<(nch = cigar[cigar_i++])){
		if(isdigit(nch)){
			tmpi = tmpi * 10 + nch - '0';
		}else{
			if(nch == 'N' || nch == 'M' || nch == 'D') tailpos -= tmpi;
			tmpi = 0;
		}
	}
	return tailpos;
} 


// This function returns 1 if the cut was added.
// It returns 0 if the head or tail cut is not able to be added (e.g., the cigar ends up like "50M4I10S" or "10S30N90M") 
int add_head_tail_cut_softclipping(global_context_t * global_context, unsigned int linear, char * cigar, int rlen, int head_cut, int tail_cut){
	//SUBREADprintf("ADD_SOFT: %s , %d, %d\n", cigar,head_cut,tail_cut);

	char cigar_added [CORE_MAX_CIGAR_STR_LEN];
	int cigar_cursor = 0, read_cursor = 0, next_read_cursor = 0;
	int tmpi = 0, nch, has_M = 0;
	unsigned int linear_cursor = linear;

	cigar_added[0]=0;

	while(1){
		nch = cigar[cigar_cursor++];
		if(nch == 0) break;

		if(isdigit(nch)){
			tmpi = 10 * tmpi + (nch - '0');
		}else{
			if('M' == nch || 'S' == nch || 'I' == nch)
				next_read_cursor = read_cursor + tmpi;

			if('M' == nch || 'D' == nch || 'S' == nch || 'N' == nch){
				linear_cursor += tmpi;
			}

			int head_S = 0, tail_S =0, remainder_tmpi = tmpi, is_skip = 0;
			
			if(next_read_cursor <= head_cut) is_skip = 1;
			if(read_cursor >= rlen - tail_cut) is_skip = 1;
			if(!is_skip){
				if(nch != 'S'){
					if(read_cursor <= head_cut)head_S = head_cut;
					if(next_read_cursor >= rlen - tail_cut) tail_S = tail_cut;
					remainder_tmpi = tmpi;
					if(head_S) remainder_tmpi-= (head_S - read_cursor);
					if(tail_S) remainder_tmpi-= next_read_cursor - (rlen - tail_S);
				}
				//SUBREADprintf("OPTR=%c, THREE-LENS=%d,%d,%d\n", nch, head_S, remainder_tmpi, tail_S);
				if(head_S >0 || tail_S > 0) if(nch !='M') return 0; 

				if(head_S > 0) sprintf(cigar_added + strlen(cigar_added), "%dS", head_S );
				if(remainder_tmpi > 0){
					sprintf(cigar_added + strlen(cigar_added), "%d%c", remainder_tmpi , nch );
					if( nch == 'M' ) has_M = 1;
				}
				if(tail_S > 0) sprintf(cigar_added + strlen(cigar_added), "%dS", tail_S );
			}

			read_cursor = next_read_cursor;
			tmpi = 0;
		}
	}
	//SUBREADprintf("RPCIGAR: %s => %s\n", cigar, cigar_added);
	strcpy(cigar, cigar_added);
	return has_M;
}


int convert_read_to_tmp(global_context_t * global_context , subread_output_context_t * output_context, int read_number, int is_second_read, int read_len, char * read_text, char * qual_text, realignment_result_t * current_result, subread_output_tmp_t * r, char * read_name)
{

	int is_r_OK;
	r -> raw_result = current_result -> mapping_result;
	r -> additional_information[0]=0;

	is_r_OK = (current_result -> mapping_result -> result_flags & CORE_IS_FULLY_EXPLAINED) > 0;

	if(0 && FIXLENstrcmp("V0112_0155:7:1101:5279:29143#ATCACG", read_name) == 0)
		SUBREADprintf("%s  R_%d CPOINT1 : is_OK=%d\n", read_name , is_second_read + 1 , is_r_OK);

	if(is_r_OK) {

		int current_strand = (current_result -> mapping_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;
		int is_first_section_jumped = current_result -> first_base_is_jumpped;
		if(is_first_section_jumped) assert((global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection));
		strcpy(r-> current_cigar_decompress, current_result -> cigar_string);

		int chimeric_sections = 0;

		r->is_first_section_jumpped = is_first_section_jumped;
		r->linear_position = current_result -> first_base_position;
		//r->mapping_quality = current_result -> final_quality;
		r->mapping_quality = 40;	// changed on 10-JUN-2016
		if(current_result -> realign_flags & CORE_IS_BREAKEVEN)
			r -> mapping_quality = 0;
		else
			r -> mapping_quality /= current_result->MAPQ_adjustment;
		//SUBREADprintf("REP=%d\n", current_result->MAPQ_adjustment);

		strcpy(r->cigar, r -> current_cigar_decompress);
		r->strand = (current_result -> mapping_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;

		//sprintf(r->additional_information, "\tSM:i:%d",current_result -> final_mismatched_bases);


		//printf("SM='%s'\n", r->additional_information);

		if((global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection))
		{			
			chimeric_sections = chimeric_cigar_parts(global_context, r->linear_position, is_first_section_jumped ^ current_strand, is_first_section_jumped, r->current_cigar_decompress, r->out_poses, output_context->out_cigar_buffer, r->out_strands, read_len, r->out_lens, read_name);

			if(chimeric_sections > 0){
				int xk1;
				r->chimeric_sections = chimeric_sections;
				strcpy(r->out_cigars[0], output_context->out_cigar_buffer[0]);
				for(xk1=1; xk1<chimeric_sections; xk1++)
				{
					int chimeric_pos;
					char * chimaric_chr;
					strcpy(r->out_cigars[xk1], output_context->out_cigar_buffer[xk1]);
					unsigned int vitual_head_pos = r->out_poses[xk1];
					char strand_xor = (r->out_strands[xk1] == '-') != is_second_read;//!= (r->out_strands[0]=='-') ;

					//if(( r->out_strands[xk1] == '-' ) != (r->out_strands[0]=='-' )) vitual_head_pos = move_to_read_head(vitual_head_pos, r->out_cigars[xk1]);

					if(0==locate_gene_position_max(vitual_head_pos ,& global_context -> chromosome_table, & chimaric_chr, & chimeric_pos, NULL, NULL, 0+r->out_lens[xk1]))
					{
						int soft_clipping_movement = 0;
						soft_clipping_movement = get_soft_clipping_length(r->out_cigars[xk1]);
						assert(chimaric_chr);
						sprintf(r->additional_information + strlen(r->additional_information), "\tCG:Z:%s\tCP:i:%u\tCT:Z:%c\tCC:Z:%s", r->out_cigars[xk1] , max(1,chimeric_pos + soft_clipping_movement + 1), strand_xor?'-':'+' , chimaric_chr );
					}
					else is_r_OK = 0;
				}
				r->linear_position = r->out_poses[0];
				r->strand = r->out_strands[0]=='-';

				strcpy(r->cigar , r->out_cigars[0]);
			}else is_r_OK = 0;
		}
		if(is_r_OK)
			r->soft_clipping_movements = get_soft_clipping_length(r->cigar);
	}

	if(is_r_OK)
	{
		int head_cut = 0 , tail_cut = 0;

		if(locate_gene_position_max(r->linear_position + r->soft_clipping_movements,& global_context -> chromosome_table, &r-> chro , &r -> offset, &head_cut, &tail_cut,(global_context-> config.do_fusion_detection || global_context-> config.do_long_del_detection)?read_len:(current_result->chromosomal_length - r->soft_clipping_movements))){
			is_r_OK = 0;
		} else {
			int is_added_OK = 1;
			if(head_cut!=0 || tail_cut!=0)
				is_added_OK = add_head_tail_cut_softclipping(global_context, r->linear_position, r->cigar , read_len, head_cut, tail_cut);

			if(is_added_OK){
				r -> offset++;
				assert(r-> chro);
			}else{
				is_r_OK = 0;
			}
		}

		if(global_context -> config.do_breakpoint_detection && !(current_result -> realign_flags & CORE_NOTFOUND_DONORS))
		{
			sprintf(r->additional_information + strlen(r->additional_information), "\tXS:A:%c", (current_result -> realign_flags & CORE_IS_GT_AG_DONORS)?'+':'-');
		}

		/*
		if(global_context -> config.more_accurate_fusions)
		{
			
			sprintf(r->additional_information + strlen(r->additional_information), "\tDF:i:%d", current_result -> best_second_diff_bases >0?current_result -> best_second_diff_bases:9999 );
		}*/


	}

	return is_r_OK;
}


int add_event_detected_from_cigar(global_context_t * global_context, unsigned int left_last_base, unsigned int right_first_base, int is_indel, int indel_len, int left_extend, int right_extend, int read_length, char * read_name)
{	
	HashTable * event_table = NULL;
	chromosome_event_t * event_space = NULL;
	event_table = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_entry_table;
	event_space = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;

	chromosome_event_t *site_events[MAX_EVENT_ENTRIES_PER_SITE];
	int search_types ;
	if(is_indel)	search_types = CHRO_EVENT_TYPE_INDEL;
	else		search_types = CHRO_EVENT_TYPE_JUNCTION | CHRO_EVENT_TYPE_FUSION;

	if(left_last_base > right_first_base) {
		unsigned int ttt;
		ttt = right_first_base;
		right_first_base = left_last_base;
		left_last_base = ttt;
	}

	int xk1, site_events_no = search_event(global_context, event_table , event_space , left_last_base, EVENT_SEARCH_BY_BOTH_SIDES , search_types , site_events);

	int is_found = 0;
	for(xk1 = 0; xk1 < site_events_no; xk1++){
		chromosome_event_t * this_event = site_events[xk1];
		if(is_indel && this_event -> event_type != CHRO_EVENT_TYPE_INDEL) continue;
		if((!is_indel) && this_event -> event_type == CHRO_EVENT_TYPE_INDEL) continue;

		#ifndef __MINGW32__
		if(0 && ( memcmp(read_name, "V0112_0155:7:1104:17648:117432",28)==0 || ( this_event -> event_small_side  > 2613701363 - 200 && this_event -> event_small_side  < 2613701363 + 100) ))
			SUBREADprintf("EVENT: L=%u, R=%u (TABLE L=%u , R=%u), LEXT=%d, REXT=%d -- %s\n", left_last_base, right_first_base, this_event -> event_small_side, this_event -> event_large_side, left_extend, right_extend, read_name);
		#endif

		if(this_event -> event_small_side == left_last_base &&  this_event -> event_large_side == right_first_base){
			if(is_indel && this_event -> indel_length == indel_len){
				this_event -> final_counted_reads ++;
				is_found = 1;
				break;
			}
			if(!is_indel){
				this_event -> final_counted_reads ++;

				//
				if(read_length >= 80)
				{
					if(left_extend > 30 && right_extend > 30) this_event -> critical_supporting_reads ++;
				}else{
					if(left_extend > 18 && right_extend > 18) this_event -> critical_supporting_reads ++;
				}
				this_event -> junction_flanking_left = max(this_event -> junction_flanking_left, left_extend);
				this_event -> junction_flanking_right = max(this_event -> junction_flanking_right, right_extend);
				is_found = 1;
				break;
			}
		}
	}

	if(!is_found){
		//SUBREADprintf("\nEVENT NOT FOUND.\n\n");
		return 1;
	}
	return 0;
}

int get_junction_right_extension(char * remainder_cigar){
	int ret = 0;
	unsigned int tmpi = 0;
	int cigar_cursor, nch;

	for(cigar_cursor = 0;;cigar_cursor++){
		nch = remainder_cigar[cigar_cursor];
		if(0 == nch)break;
		if(isdigit(nch))tmpi = tmpi*10 + nch - '0';
		else{
			if(nch == 'M' || nch == 'D')
				ret += tmpi;
			if(nch == 'N' || nch == 'n' || nch == 'B' || nch == 'b')
				break;

			tmpi = 0;
		}
	}
	return ret;
}

void remove_soft_clipping(char * dst, char * src){
	unsigned int tmpi = 0, head_clip = 0, tail_clip = 0, last_M = 0;
	int cigar_cursor, is_first_s = 1, nch, nch2;

	dst[0] = 0;

	for(cigar_cursor = 0;;cigar_cursor++){
		nch = src[cigar_cursor];
		nch2 = src[cigar_cursor+1];

		if(0 == nch)break;
		if(isdigit(nch))tmpi = tmpi*10 + nch - '0';
		else{
			if('S' == nch){
				if(is_first_s)	head_clip = tmpi;
				if(nch2 == 0)	tail_clip = tmpi;
			}
			else if('M' == nch){
				last_M = tmpi;
			}
			else{
				if(last_M){
					sprintf(dst + strlen(dst), "%uM", last_M + head_clip);
					last_M = 0;
					head_clip = 0;
				}
				sprintf(dst + strlen(dst), "%u%c", tmpi, nch);
			}
			tmpi=0;
			is_first_s = 0;
		}
	}
	if(last_M){
		sprintf(dst + strlen(dst), "%uM" , last_M + tail_clip + head_clip);
	}
	
}

int getFirstM(char * cig){
	int tmpi = 0, nch, x1;

	for(x1=0; 0!=(nch = cig[x1]) ; x1++){
		if(isdigit(nch)){
			tmpi=tmpi*10 + nch - '0';
		}else break;
	}
	return tmpi;
}

int calc_tlen(global_context_t * global_context, subread_output_tmp_t * rec1 , subread_output_tmp_t * rec2, int read_len_1, int  read_len_2);

int calc_flags(global_context_t * global_context, thread_context_t * thread_context, subread_output_tmp_t * rec1 , subread_output_tmp_t * rec2, int is_second_read, int read_len_1, int read_len_2, int current_location_no , int tlen, int this_OK, int mate_OK)
{
	int ret, is_TLEN_wrong = 0;

	if(global_context->input_reads.is_paired_end_reads)
	{
		ret  = SAM_FLAG_PAIRED_TASK;
		ret |= is_second_read?SAM_FLAG_SECOND_READ_IN_PAIR:SAM_FLAG_FIRST_READ_IN_PAIR;

		subread_output_tmp_t * this_rec = is_second_read?rec2:rec1;
		subread_output_tmp_t * mate_rec = is_second_read?rec1:rec2;

		if(this_OK == 0)  ret |= SAM_FLAG_UNMAPPED;
		else
			if(this_rec->strand + is_second_read == 1) ret |= SAM_FLAG_REVERSE_STRAND_MATCHED;

		if(mate_OK == 0)  ret |= SAM_FLAG_MATE_UNMATCHED;
		else
			if(mate_rec->strand + is_second_read != 1) ret |= SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;

		if(rec1 && rec2)
		{
			int TLEN = tlen;//calc_tlen(global_context , rec1, rec2, read_len_1, read_len_2);
			int is_PEM = 0;
			if( rec1 -> chro == rec2 -> chro /* two pointers can be directly compared. */ && TLEN >= global_context->config.minimum_pair_distance &&
			    TLEN <= global_context-> config.maximum_pair_distance &&  this_rec->strand == mate_rec->strand )
			{
				if(global_context -> config.is_first_read_reversed && !(global_context -> config.is_second_read_reversed))
				{
					if(this_rec -> strand == 0) {
						if((is_second_read + (mate_rec-> offset > this_rec -> offset) == 1) || mate_rec-> offset == this_rec -> offset)
							is_PEM = 1;
						else is_TLEN_wrong = 1;
					}
				}
				else
				{
					if(this_rec -> strand) {
						if((is_second_read + (mate_rec-> offset < this_rec -> offset) == 1) || mate_rec-> offset == this_rec -> offset) is_PEM = 1;
						else is_TLEN_wrong = 1;
					}else {
						if((is_second_read + (mate_rec-> offset > this_rec -> offset) == 1) || mate_rec-> offset == this_rec -> offset) is_PEM = 1;
						else is_TLEN_wrong = 1;
					}
				}
			}
			if(is_PEM) ret |= SAM_FLAG_MATCHED_IN_PAIR;
			else if(is_second_read){
				if(rec1 -> chro != rec2 -> chro){
					if(thread_context) thread_context -> not_properly_pairs_different_chro ++;
					else global_context -> not_properly_pairs_different_chro++;
				}else if(this_rec->strand != mate_rec->strand){
					if(thread_context) thread_context -> not_properly_different_strands ++;
					else global_context -> not_properly_different_strands ++;
				}else if(TLEN < global_context->config.minimum_pair_distance || TLEN > global_context-> config.maximum_pair_distance){
					if(thread_context) thread_context -> not_properly_pairs_TLEN_wrong ++;
					else global_context -> not_properly_pairs_TLEN_wrong ++;
				}else if(is_TLEN_wrong){
					if(thread_context) thread_context -> not_properly_pairs_wrong_arrangement ++;
					else global_context  -> not_properly_pairs_wrong_arrangement ++;
				}

			}
		}
	}
	else
	{
		ret = 0;

		if(this_OK == 0)  ret |= SAM_FLAG_UNMAPPED;
		else
			if(rec1->strand) ret |= SAM_FLAG_REVERSE_STRAND_MATCHED;
	}

	if(current_location_no>0){
		if((rec1 && !is_second_read)||(rec2 && is_second_read))
			ret |= SAM_FLAG_SECONDARY_MAPPING;
	}

	return ret;
}


int calc_tlen(global_context_t * global_context, subread_output_tmp_t * rec1 , subread_output_tmp_t * rec2, int read_len_1, int  read_len_2)
{
	int ret = -1;

	unsigned int r1_head_pos = rec1 -> offset;// - rec1 -> soft_clipping_movements;
	unsigned int r2_head_pos = rec2 -> offset;// - rec2 -> soft_clipping_movements;

	if(r1_head_pos == r2_head_pos)
		ret = max(read_len_1, read_len_2);	// the two reads are fully overlapping
	else{
		int is_r2_smaller = r2_head_pos < r1_head_pos;
		subread_output_tmp_t * smaller_rec = is_r2_smaller?rec2:rec1;
		unsigned int smaller_head_pos = is_r2_smaller?r2_head_pos:r1_head_pos;
		unsigned int larger_head_pos = is_r2_smaller?r1_head_pos:r2_head_pos;

		int len_larger_rec = is_r2_smaller?read_len_1:read_len_2;
		int len_smaller_rec = is_r2_smaller?read_len_2:read_len_1;

		unsigned int tmpi = 0;
		unsigned int chro_cursor = smaller_head_pos;
		unsigned int section_end = 0;
		unsigned int read_cursor = 0;
		int cigar_cursor, nch, nch2;
		for(cigar_cursor = 0; ; cigar_cursor++){
			nch = smaller_rec -> cigar[cigar_cursor];
			nch2 = smaller_rec -> cigar[cigar_cursor+1];

			if(nch > 0){
				if(isdigit(nch)) tmpi = tmpi * 10 + nch - '0';
				else{
					if(nch == 'M' || nch == 'S'){
						chro_cursor += tmpi;
						read_cursor += tmpi;
						section_end = chro_cursor;
					}
					if(nch == 'N' || nch == 'B' || nch == 'b' || nch == 'n' || nch == 'I' || nch == 'D' || nch2 == 0){
						if(nch == 'N' || nch == 'D')
							chro_cursor += tmpi;

						//SUBREADprintf("NCH=%c, NCH2=%d, SEC_END=%u, LARGE_HEAD=%u\n", nch, nch2, section_end, larger_head_pos);
						if(section_end >= larger_head_pos)
						{
							ret = read_cursor + larger_head_pos - section_end + len_larger_rec;
							break;
						}
					}

					if(nch == 'I'){
						read_cursor += tmpi;
					}

					if(nch == 'B' || nch == 'b' || nch == 'n') break;	// fusion block, unable to determine the TLEN after the fusion point. 

					tmpi = 0;
				}
			}else{
				break;
			}

		} 

		if(ret < 0)
			ret  = larger_head_pos - section_end + len_larger_rec + len_smaller_rec;
	}


	//SUBREADprintf("TLEN = %d :: \tREC1 : %s:%u (%s  - %d) ; REC2: %s:%u (%s  - %d )\n", ret, rec1 -> chro, rec1 -> offset, rec1 -> cigar, rec1 -> soft_clipping_movements, 
	//								       rec2 -> chro, rec2 -> offset, rec2 -> cigar, rec2 -> soft_clipping_movements);

	return ret;

}

int calc_should_reverse(global_context_t * global_context, subread_output_tmp_t * rec1 , subread_output_tmp_t * rec2, int is_second_read)
{
	subread_output_tmp_t * cur_rec = is_second_read?rec2:rec1;

	if(cur_rec)
		return cur_rec->strand;
	else
	{
		return is_second_read;
		//subread_output_tmp_t * mate_rec = is_second_read?rec2:rec1;
		//return mate_rec?mate_rec->strand:is_second_read;
	}
}

void remove_nm_i(char * st)
{
	char * st_nm = strstr(st, "\tNM:i:");
	if(st_nm)
	{
		char * write_cursor = st_nm;
		int started = 0;
		st_nm++;
		while(1)
		{
			char nch = *st_nm;
			if(!nch)
			{
				*write_cursor=0;
				break;
			}
			if(nch == '\t') started = 1;
			if(started)
			{
				*write_cursor=*st_nm;
				write_cursor++;
			}
			st_nm++;
		}
	}
}

double last_t1 = 0;
int for_other_thread = 0;
int for_one_threads = 0;


void add_buffered_fragment(global_context_t * global_context, thread_context_t * thread_context, subread_read_number_t pair_number ,
	char * read_name1, unsigned int flags1, char * chro_name1, unsigned int chro_position1, int mapping_quality1, char * cigar1,
	char * next_chro_name1, unsigned int next_chro_pos1, int temp_len1, int read_len1,
	char * read_text1, char * qual_text1, char * additional_columns1,
	char * read_name2, unsigned int flags2, char * chro_name2, unsigned int chro_position2, int mapping_quality2, char * cigar2,
	char * next_chro_name2, unsigned int next_chro_pos2, int temp_len2, int read_len2,
	char * read_text2, char * qual_text2, char * additional_columns2,
	int all_locations, int this_location
	){

	assert(global_context -> config.all_threads > 1);

	if( global_context -> config.is_BAM_output  && !global_context -> config.is_input_read_order_required){
		//SUBREADprintf("MTBAM : THR_%d CORE: %s\n", thread_context -> thread_id, read_name1);
		SamBam_writer_add_read(global_context -> output_bam_writer, thread_context -> thread_id, read_name1, flags1, chro_name1, chro_position1, mapping_quality1, cigar1, next_chro_name1, next_chro_pos1, temp_len1, read_len1, read_text1, qual_text1, additional_columns1, !global_context->input_reads.is_paired_end_reads);
		if(global_context->input_reads.is_paired_end_reads)
			SamBam_writer_add_read(global_context -> output_bam_writer, thread_context -> thread_id, read_name2, flags2, chro_name2, chro_position2, mapping_quality2, cigar2, next_chro_name2, next_chro_pos2, temp_len2, read_len2, read_text2, qual_text2, additional_columns2, 1);
	}else{
		while(1){
			int is_fin = 0;
			subread_lock_occupy(&global_context -> output_lock);
			//SUBREADprintf("THWT : TH_%d IS AFTER %lld ; LAST=%lld (%d/%d)\n", thread_context -> thread_id, pair_number -1, global_context -> last_written_fragment_number, this_location, all_locations);
			if(global_context -> last_written_fragment_number == pair_number-1){
				if(global_context -> config.is_BAM_output){
					SamBam_writer_add_read(global_context -> output_bam_writer, -1, read_name1, flags1,  chro_name1 , chro_position1, mapping_quality1, cigar1, next_chro_name1 , next_chro_pos1, temp_len1, read_len1, read_text1, qual_text1, additional_columns1, !global_context->input_reads.is_paired_end_reads);

					if(global_context->input_reads.is_paired_end_reads)
						SamBam_writer_add_read(global_context -> output_bam_writer, -2, read_name2, flags2,  chro_name2 , chro_position2, mapping_quality2, cigar2, next_chro_name2 , next_chro_pos2, temp_len2, read_len2, read_text2, qual_text2, additional_columns2, 1);

				}else{
					int write_len_2 = 100, write_len = sambamout_fprintf(global_context -> output_sam_fp , "%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%d\t%s\t%s%s%s\n", read_name1, flags1, chro_name1, chro_position1, mapping_quality1, cigar1, next_chro_name1, next_chro_pos1,  temp_len1, read_text1, qual_text1, additional_columns1[0]?"\t":"", additional_columns1);
					if(global_context->input_reads.is_paired_end_reads)
						write_len_2 = sambamout_fprintf(global_context -> output_sam_fp , "%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%d\t%s\t%s%s%s\n", read_name2, flags2, chro_name2, chro_position2, mapping_quality2, cigar2, next_chro_name2, next_chro_pos2,  temp_len2, read_text2, qual_text2, additional_columns2[0]?"\t":"", additional_columns2);
			
					if( write_len < 10 || write_len_2 < 10 ){
						global_context -> output_sam_is_full = 1;
					}
				}
				if(all_locations <= this_location +1)global_context -> last_written_fragment_number = pair_number;
				is_fin = 1;
			}
			subread_lock_release(&global_context -> output_lock);
			if(is_fin)break;
			usleep(2);
		}
	}
}

#define write_chunk_results_145 write_chunk_results

// rec1 or rec2 is OK if they are not NULL.
void write_single_fragment(global_context_t * global_context, thread_context_t * thread_context, subread_output_tmp_t * rec1, realignment_result_t * raw_r1, subread_output_tmp_t * rec2, realignment_result_t * raw_r2, int all_locations , int current_location , char * read_name_1, char * read_name_2, int read_len_1, int read_len_2, char * read_text_1, char * read_text_2, char * qual_text_1, char * qual_text_2, subread_read_number_t pair_number, int non_informative_subread_r1, int non_informative_subread_r2, int is_R1_OK, int is_R2_OK)
{

	//assert(all_locations <= global_context -> config.reported_multi_best_reads);
	int tlen = 0;

	if(is_R1_OK && is_R2_OK && rec1->chro == rec2->chro)tlen = calc_tlen(global_context , rec1, rec2,  read_len_1,  read_len_2);

	//#warning "SUBREAD_151 ============= REMOVE '+1' when not doing structure variance detection for the proposal  =========="
	if(0)
	if((global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection) && rec1 && rec2 && current_location == 0 && rec1 -> chimeric_sections == 1 && rec2 -> chimeric_sections == 1){
		// when current_location == 0, both r1 and r2 were on the same strand.
		int is_funky = is_funky_fragment(global_context, read_name_1, rec1 -> chro, rec1 -> offset, read_len_1, rec1 -> strand, rec1 -> out_cigars[0], read_text_1, read_name_2, rec2 -> chro, rec2 -> offset, read_len_2, rec2 -> strand, rec2 -> out_cigars[0] , read_text_2, tlen);
		if (0 &&is_funky)
			SUBREADprintf("RNAME %s is %d, CHRO = '%s' %p  '%s' %p, POS=%u, %u\n", read_name_1, is_funky, rec1 -> chro,rec1 -> chro,rec2 -> chro,rec2 -> chro, rec1 -> offset, rec2 -> offset);

		if(is_funky & FUNKY_FRAGMENT_A){
			fraglist_append(&global_context -> funky_list_A, pair_number);
		}
		if(1)if(is_funky & FUNKY_FRAGMENT_BC){
			//#warning "LOGIC WRONG: R1 AND R2 SHOULD BE DECIDED BY THEIR MAPPING POSITIONS"
			bktable_append(&global_context -> funky_table_BC, rec1 -> chro, rec1 -> offset , NULL + (2*pair_number));
			bktable_append(&global_context -> funky_table_BC, rec2 -> chro, rec2 -> offset , NULL + (2*pair_number+1));
		}
		if(1)if(is_funky & FUNKY_FRAGMENT_DE){
			fraglist_append(&global_context -> funky_list_DE, pair_number);
			bktable_append(&global_context -> funky_table_DE, rec1 -> chro, rec1 -> offset , NULL + (2*pair_number + (rec1 -> offset > rec2 -> offset ? 1:0)));
			bktable_append(&global_context -> funky_table_DE, rec2 -> chro, rec2 -> offset , NULL + (2*pair_number + (rec1 -> offset < rec2 -> offset ? 1:0)));
		}
	}

	int flag1 = calc_flags( global_context, thread_context , rec1, rec2, 0,  read_len_1,  read_len_2, current_location, tlen, is_R1_OK, is_R2_OK);

	int flag2 = -1;

	if(global_context->input_reads.is_paired_end_reads)
	{
		flag2 = calc_flags( global_context , thread_context, rec1, rec2, 1,  read_len_1,  read_len_2, current_location, tlen, is_R2_OK, is_R1_OK);
		if((0 == current_location) && (flag2 & SAM_FLAG_MATCHED_IN_PAIR)){
			if(thread_context)thread_context->all_correct_PE_reads  ++;
			else global_context->all_correct_PE_reads  ++;
		}
	}

	int applied_reverse_space;
	applied_reverse_space = global_context->config.space_type;
	if(global_context -> config.convert_color_to_base)
	{

		//SUBREADprintf("BEST_READ_NO = %d / %d\n", current_location, all_locations);
		//SUBREADprintf("ORGI: %s\n", read_text_1);
		colorread2base(read_text_1, read_len_1+1);
		//SUBREADprintf("CONV: %s\n\n", read_text_1);

		if(global_context->input_reads.is_paired_end_reads) 
			colorread2base(read_text_2, read_len_2+1);
		
		applied_reverse_space = GENE_SPACE_BASE;
	}

	int should_1_reverse = calc_should_reverse( global_context , rec1, rec2, 0);

	if(should_1_reverse)
	{
		reverse_read(read_text_1, read_len_1 + global_context->config.convert_color_to_base, applied_reverse_space);
		reverse_quality(qual_text_1, read_len_1);
	}

	if(global_context->input_reads.is_paired_end_reads)
	{
		int should_2_reverse = calc_should_reverse( global_context , rec1, rec2, 1);

		if(should_2_reverse)
		{
			reverse_read(read_text_2, read_len_2 + global_context->config.convert_color_to_base, applied_reverse_space);
			reverse_quality(qual_text_2, read_len_2);
		}
	}
	remove_backslash(read_name_1);
	if(global_context->input_reads.is_paired_end_reads) remove_backslash(read_name_2);

	int display_offset1 = 0, display_tailgate1 = 0;
	int display_offset2 = 0, display_tailgate2 = 0;

	if(global_context -> config.space_type == GENE_SPACE_COLOR)
	{

		if(rec1 && rec1 -> strand)
		{
			display_offset1 = 0;
			display_tailgate1 = 1;
		}
		else
		{
			if(!global_context -> config.convert_color_to_base)
			{
				// the first base was a fake prime base; the second base is the first meaningful base.
				char second_char = read_text_1[1];
				read_text_1[1] = color2char(second_char, read_text_1[0]);
			}
			display_offset1 = 1;
			display_tailgate1 = 0;
		}

		if(rec2 && rec2 -> strand)
		{
			if(!global_context -> config.convert_color_to_base)
			{
				// the first base was a fake prime base; the second base is the first meaningful base.
				char second_char = read_text_2[1];
				read_text_2[1] = color2char(second_char, read_text_2[0]);
			}
			display_offset2 = 1;
			display_tailgate2 = 0;
		}
		else
		{
			display_offset2 = 0;
			display_tailgate2 = 1;
		}
		assert(display_offset1 + display_tailgate1 == 1);
		assert(display_offset2 + display_tailgate2 == 1);
	}
	if(!qual_text_1[0])
	{
		int xi2;
		for(xi2=display_offset1;read_text_1[xi2 + display_tailgate1];xi2++) qual_text_1[xi2 - display_offset1] = 'I';
		qual_text_1[xi2 - display_offset1]=0;
		for(xi2=display_offset2;read_text_2[xi2 + display_tailgate2];xi2++) qual_text_2[xi2 - display_offset2] = 'I';
		qual_text_2[xi2 - display_offset2]=0;

	}



	if(display_tailgate1)
		read_text_1[strlen(read_text_1)-1]=0;
	if(display_tailgate2)
		read_text_2[strlen(read_text_2)-1]=0;


	if(global_context -> config.space_type == GENE_SPACE_BASE){ 
		if(is_R1_OK) {// && !strstr(rec1->additional_information, "\tNM:i:")){
			short rec1_edit = calc_edit_dist(global_context, rec1->raw_result, rec1->cigar , rec1->linear_position, read_text_1, raw_r1 -> final_mismatched_bases);
			sprintf(rec1->additional_information + strlen( rec1->additional_information), "\tNM:i:%d", rec1_edit );
		}
		if(global_context->input_reads.is_paired_end_reads && is_R2_OK) //&& !strstr(rec2->additional_information, "\tNM:i:"))
		{
			short rec2_edit = calc_edit_dist(global_context, rec2->raw_result, rec2->cigar , rec2->linear_position, read_text_2, raw_r2 -> final_mismatched_bases);
			sprintf(rec2->additional_information + strlen( rec2->additional_information), "\tNM:i:%d", rec2_edit );
		}
	}

	char extra_additional_1 [1000+CORE_ADDITIONAL_INFO_LENGTH], extra_additional_2[1000+CORE_ADDITIONAL_INFO_LENGTH];
	int extra_additional_1_ptr=0, extra_additional_2_ptr=0;

	extra_additional_1[0] = extra_additional_2[0]=0;

	if(is_R1_OK || is_R2_OK){	// no NH and HI only if both ends are unmapped.
		extra_additional_1_ptr +=snprintf(extra_additional_1, 310, "HI:i:%d\tNH:i:%d", current_location+1, all_locations);
		extra_additional_2_ptr +=snprintf(extra_additional_2, 310, "HI:i:%d\tNH:i:%d", current_location+1, all_locations);
	}

	if(global_context->config.read_group_id[0])
	{
		extra_additional_1_ptr += snprintf(extra_additional_1+ extra_additional_1_ptr, 310, "\tRG:Z:%s", global_context->config.read_group_id);
		extra_additional_2_ptr += snprintf(extra_additional_2+ extra_additional_2_ptr, 310, "\tRG:Z:%s", global_context->config.read_group_id);
	}

	char * out_chro1, * out_chro2, *out_cigar1, *out_cigar2;
	if(is_R1_OK)
	{
		strcpy(extra_additional_1 + extra_additional_1_ptr, rec1->additional_information);
		out_chro1 = rec1->chro;
		out_cigar1 = rec1->cigar;
	}
	else
	{
		out_chro1 = "*";
		out_cigar1 = "*";
	}

	if(is_R2_OK)
	{
		strcpy(extra_additional_2 + extra_additional_2_ptr, rec2->additional_information);
		out_chro2 = rec2->chro;
		out_cigar2 = rec2->cigar;
	}
	else
	{
		out_chro2 = "*";
		out_cigar2 = "*";
	}

	int out_offset1=0, out_offset2=0;
	long long int out_tlen1, out_tlen2;
	int  out_mapping_quality1 = 0, out_mapping_quality2 = 0;

	out_tlen1 = tlen;
	out_tlen2 = tlen;
	if(is_R1_OK && is_R2_OK)
	{
		if( rec1->offset >  rec2->offset) out_tlen1 = - out_tlen1;
		else if(rec2->offset > rec1->offset) out_tlen2 = -out_tlen2;
		else{
			if( rec1 -> strand )  out_tlen1 = - out_tlen1;
			else out_tlen2 = -out_tlen2;
		}
	}

	if(0==current_location)
	{
		if(global_context -> input_reads.is_paired_end_reads)
		{
			if(is_R1_OK || is_R2_OK){
				//SUBREADprintf("MAPTOLEN\n");
				if(thread_context)
					thread_context -> all_mapped_reads ++;
				else
					global_context -> all_mapped_reads ++;
			}
		}
		else if(is_R1_OK){
			if(thread_context)
				thread_context -> all_mapped_reads ++;
			else
				global_context -> all_mapped_reads ++;
		}
	}

	if(is_R1_OK){
		out_offset1 = max(1, rec1->offset);
		out_mapping_quality1 = rec1->mapping_quality;
	}
	if(is_R2_OK){
		out_offset2 = max(1, rec2->offset);
		out_mapping_quality2 = rec2->mapping_quality;
	}


	char * mate_chro_for_1  = out_chro2;
	char * mate_chro_for_2  = out_chro1;

	if(out_chro1 == out_chro2 && out_chro1 && out_chro1[0]!='*'){
		mate_chro_for_1="=";
		mate_chro_for_2="=";
	}

	//if(pair_number<10)SUBREADprintf("FIN_ADD #%d: %s in %p \n", pair_number, read_name_1, thread_context);

	if(thread_context)
		add_buffered_fragment(global_context, thread_context, pair_number,
			read_name_1, flag1,  out_chro1 , out_offset1, out_mapping_quality1, out_cigar1, mate_chro_for_1 , out_offset2, out_tlen1, read_len_1,
			read_text_1 + display_offset1, qual_text_1, extra_additional_1,
			read_name_2, flag2,  out_chro2 , out_offset2, out_mapping_quality2, out_cigar2, mate_chro_for_2 , out_offset1, out_tlen2, read_len_2,
			read_text_2 + display_offset2, qual_text_2, extra_additional_2, all_locations, current_location);
	else{
		//subread_lock_occupy(&global_context -> output_lock);
		if(global_context -> config.is_BAM_output){
			SamBam_writer_add_read(global_context -> output_bam_writer, -1, read_name_1, flag1,  out_chro1 , out_offset1, out_mapping_quality1, out_cigar1, mate_chro_for_1 , out_offset2, out_tlen1, read_len_1, read_text_1 + display_offset1, qual_text_1, extra_additional_1, !global_context->input_reads.is_paired_end_reads);

			if(global_context->input_reads.is_paired_end_reads)
				SamBam_writer_add_read(global_context -> output_bam_writer, -1, read_name_2, flag2,  out_chro2 , out_offset2, out_mapping_quality2, out_cigar2, mate_chro_for_2 , out_offset1, out_tlen2, read_len_2, read_text_2 + display_offset2, qual_text_2, extra_additional_2, 1);
		}
		else
		{
			int write_len_2 = 100, 
			#ifdef __MINGW32__
			write_len = sambamout_fprintf(global_context -> output_sam_fp , "%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%I64d\t%s\t%s%s%s\n", read_name_1, flag1, out_chro1, out_offset1, out_mapping_quality1, out_cigar1, mate_chro_for_1, out_offset2, out_tlen1, read_text_1 + display_offset1, qual_text_1, extra_additional_1[0]?"\t":"", extra_additional_1);
			#else
			write_len = sambamout_fprintf(global_context -> output_sam_fp , "%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%lld\t%s\t%s%s%s\n", read_name_1, flag1, out_chro1, out_offset1, out_mapping_quality1, out_cigar1, mate_chro_for_1, out_offset2, out_tlen1, read_text_1 + display_offset1, qual_text_1, extra_additional_1[0]?"\t":"", extra_additional_1);
			#endif
			if(global_context->input_reads.is_paired_end_reads)
			#ifdef __MINGW32__
				write_len_2 = sambamout_fprintf(global_context -> output_sam_fp , "%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%I64d\t%s\t%s%s%s\n", read_name_2, flag2, out_chro2, out_offset2, out_mapping_quality2, out_cigar2, mate_chro_for_2, out_offset1, out_tlen2, read_text_2 + display_offset2, qual_text_2, extra_additional_2[0]?"\t":"", extra_additional_2);
			#else
				write_len_2 = sambamout_fprintf(global_context -> output_sam_fp , "%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%lld\t%s\t%s%s%s\n", read_name_2, flag2, out_chro2, out_offset2, out_mapping_quality2, out_cigar2, mate_chro_for_2, out_offset1, out_tlen2, read_text_2 + display_offset2, qual_text_2, extra_additional_2[0]?"\t":"", extra_additional_2);
			#endif

			if( write_len < 10 || write_len_2 < 10 ){
				global_context -> output_sam_is_full = 1;
			}
		}
		//subread_lock_release(&global_context -> output_lock);
	}


	if(is_R1_OK) rec1->raw_result->selected_position += rec1->soft_clipping_movements;
	if(is_R2_OK) rec2->raw_result->selected_position += rec2->soft_clipping_movements;

}


void init_chunk_scanning_parameters(global_context_t * global_context, thread_context_t * thread_context, gene_input_t ** ginp1, gene_input_t ** ginp2)
{
	*ginp2 = NULL;
	*ginp1 = &global_context->input_reads.first_read_file;
	if(global_context->input_reads.is_paired_end_reads)
		*ginp2 = &global_context->input_reads.second_read_file;
}

gene_value_index_t * find_current_value_index(global_context_t * global_context, unsigned int pos, int len)
{
	int block_no;
	if(global_context->index_block_number<2)
	{
		unsigned index_begin = global_context -> all_value_indexes [0] . start_base_offset; 
		unsigned index_end = global_context -> all_value_indexes [0] . start_base_offset + global_context -> all_value_indexes [0] . length;

		if(pos>=index_begin && pos + len<index_end)
			return & global_context->all_value_indexes [0];
		else return NULL;
	}
	else
		for(block_no=0; block_no<global_context->index_block_number; block_no++)
		{
			unsigned index_begin = global_context -> all_value_indexes [block_no] . start_base_offset; 
			unsigned index_end = global_context -> all_value_indexes [block_no] . start_base_offset + global_context -> all_value_indexes [block_no] . length;
			if((block_no == 0 && pos >=  index_begin && pos < index_end - 1000000) ||
			   (block_no>0 && block_no<global_context->index_block_number-1 && pos >= index_begin+ 1000000 && pos < index_end - 1000000) ||
			   (block_no == global_context->index_block_number-1 && pos >= index_begin + 1000000 && pos < index_end ))
			{
				return & global_context -> all_value_indexes [block_no];
			}
		}
	return NULL;
}
//this function selects the correct all_value_indexes from global_context and put it to global_context or thread_context if thread_context is not NULL.
int locate_current_value_index(global_context_t * global_context, thread_context_t * thread_context, mapping_result_t * result, int rlen)
{
	int block_no;

	if(global_context->index_block_number<2)
	{
		unsigned index_begin = global_context -> all_value_indexes [0] . start_base_offset; 
		unsigned index_end = global_context -> all_value_indexes [0] . start_base_offset + global_context -> all_value_indexes [0] . length;

		//SUBREADprintf("RESET2 : %u should <= %u\n",  result->selected_position + rlen, index_end);

		if(result->selected_position>=index_begin && result->selected_position + rlen<=index_end)
		{
			if(thread_context)thread_context->current_value_index = & global_context->all_value_indexes [0];
			else global_context->current_value_index =& global_context->all_value_indexes [0];
			return 0;
		}
		else return 1;
	}
	for(block_no=0; block_no<global_context->index_block_number; block_no++)
	{
		unsigned index_begin = global_context -> all_value_indexes [block_no] . start_base_offset; 
		unsigned index_end = global_context -> all_value_indexes [block_no] . start_base_offset + global_context -> all_value_indexes [block_no] . length;
		if((block_no == 0 && result->selected_position >=  index_begin && result->selected_position < index_end - 1000000) ||
		   (block_no>0 && block_no<global_context->index_block_number-1 && result->selected_position >= index_begin+ 1000000 && result->selected_position < index_end - 1000000) ||
		   (block_no == global_context->index_block_number-1 && result->selected_position >= index_begin + 1000000 && result->selected_position < index_end ))
		{
			if(thread_context)thread_context->current_value_index =& global_context -> all_value_indexes [block_no];
			else global_context->current_value_index =& global_context -> all_value_indexes [block_no];
			return 0;
		}
	}
	return 1;
}




int do_iteration_one(global_context_t * global_context, thread_context_t * thread_context)
{
	assert(0); // not needed anymore
	return 0;
}



int finish_iteration_three(global_context_t * global_context, thread_context_t * thread_context)
{
	return 0;
}
int do_iteration_three(global_context_t * global_context, thread_context_t * thread_context)
{
	int ret;
	gene_input_t * ginp1 = NULL , * ginp2 = NULL;
	subread_read_number_t current_read_number=0;
	char read_text_1[MAX_READ_LENGTH+1], read_text_2[MAX_READ_LENGTH+1];
	char qual_text_1[MAX_READ_LENGTH+1], qual_text_2[MAX_READ_LENGTH+1];
	char read_name_1[MAX_READ_NAME_LEN+1], read_name_2[MAX_READ_NAME_LEN+1];
	int read_len_1, read_len_2=0, sqr_interval, sqr_read_number = 0;

	//unsigned int low_index_border = global_context -> current_value_index -> start_base_offset;
	//unsigned int high_index_border = global_context -> current_value_index -> start_base_offset + global_context -> current_value_index -> length; 

	print_in_box(80,0,0,"Prepare for long indel deleteion...");
	init_chunk_scanning_parameters(global_context,thread_context, & ginp1, & ginp2);
	sqr_interval = max(5000,global_context -> processed_reads_in_chunk/10/ global_context -> config.all_threads);

	while(1)
	{
		int is_second_read;

		sqr_read_number++;
		ret = fetch_next_read_pair(global_context, thread_context, ginp1, ginp2, &read_len_1, &read_len_2, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2, 1, &current_read_number);
		if(ret) break;

		int best_read_id = 0;
		for (is_second_read = 0; is_second_read < 1 + global_context -> input_reads.is_paired_end_reads; is_second_read ++)
		{
			int is_reversed_already = 0;
			for(best_read_id = 0; best_read_id < global_context -> config.multi_best_reads; best_read_id++)
			{
				mapping_result_t current_result_body;
				mapping_result_t mate_result_body;

				mapping_result_t * current_result = &current_result_body, *mate_result = NULL;
				bigtable_readonly_result(global_context, NULL, current_read_number, best_read_id, is_second_read, current_result, NULL);

				if(global_context -> input_reads.is_paired_end_reads){
					mate_result   = &mate_result_body;
					bigtable_readonly_result(global_context, NULL, current_read_number, best_read_id,!is_second_read, mate_result, NULL);
				}

				//if(best_read_id && (current_result->selected_votes <1)) break;

				char * current_read_name =  is_second_read?read_name_2 : read_name_1;
				char * current_read =  is_second_read?read_text_2 : read_text_1;
				char * current_qual =  is_second_read?qual_text_2 : qual_text_1;
				int current_rlen = is_second_read?read_len_2:read_len_1;
				int mate_rlen = is_second_read?read_len_1:read_len_2;

				if((current_result->result_flags & CORE_IS_FULLY_EXPLAINED) || (global_context -> input_reads.is_paired_end_reads &&(mate_result -> result_flags & CORE_IS_FULLY_EXPLAINED)))
				{
					//SUBREADprintf("BUILD BLOCK FOR READ_%d '%s' BEST=%d\n", is_second_read + 1, current_read_name , best_read_id);
					if(current_result->result_flags & CORE_IS_FULLY_EXPLAINED)
					{
						int is_negative_strand = ((current_result  -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0);

						if(is_negative_strand + is_reversed_already ==1)
						{
							reverse_read(current_read, current_rlen,  global_context->config.space_type);
							reverse_quality(current_qual, current_rlen);
							is_reversed_already=!is_reversed_already;
						}

						build_local_reassembly(global_context , thread_context , current_read_number, current_read_name , current_read , current_qual , current_rlen, 0 , is_second_read, best_read_id, 0, current_result, mate_result);

					}
					else if(global_context -> input_reads.is_paired_end_reads && (mate_result -> result_flags & CORE_IS_FULLY_EXPLAINED))
					{
						int is_negative_strand = ((mate_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0);
						if(is_negative_strand+is_reversed_already==1)
						{
							reverse_read(current_read, current_rlen,  global_context->config.space_type);
							reverse_quality(current_qual, current_rlen);
							is_reversed_already=!is_reversed_already;
						}

						build_local_reassembly(global_context , thread_context , current_read_number , current_read_name , current_read , current_qual , current_rlen, mate_rlen , is_second_read, best_read_id, 1, current_result, mate_result);
					}
				}	
			}
		}

		if(!thread_context || thread_context->thread_id == 0)
		{
			if(sqr_read_number>sqr_interval)
			{
				show_progress(global_context, thread_context, current_read_number, STEP_ITERATION_THREE);
				sqr_read_number = 0;
			}
		}
	}

	return 0;

}


void add_realignment_event_support(global_context_t * global_context , realignment_result_t * res){
	int xk1;
	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID];

	for(xk1 = 0; xk1 < MAX_EVENTS_IN_READ ; xk1++){
		chromosome_event_t *sup = res -> supporting_chromosome_events[xk1];
		if(!sup)break;

		int lock_hash = sup -> global_event_id % EVENT_BODY_LOCK_BUCKETS;
		subread_lock_occupy(indel_context -> event_body_locks+lock_hash);
		sup -> final_counted_reads ++;
		sup -> junction_flanking_left = max(sup -> junction_flanking_left, res -> flanking_size_left[xk1]);
		sup -> junction_flanking_right = max(sup -> junction_flanking_right, res -> flanking_size_right[xk1]);
		subread_lock_release(indel_context -> event_body_locks+lock_hash);
	}
}

unsigned int calc_end_pos(unsigned int p, char * cigar, unsigned int * all_skipped_len, int * is_exonic_regions, global_context_t * global_context);
void test_PE_and_same_chro_align(global_context_t * global_context , realignment_result_t * res1, realignment_result_t * res2, int * is_exonic_regions, int * is_PE_distance, int * is_same_chromosome, int read_len_1, int read_len_2, char * rname, int * res_tlen);
void write_realignments_for_fragment(global_context_t * global_context, thread_context_t * thread_context, subread_output_context_t * out_context, unsigned int read_number, realignment_result_t * res1, realignment_result_t * res2, char * read_name_1, char * read_name_2, char * read_text_1, char * read_text_2, char * qual_text_1, char * qual_text_2 , int rlen1 , int rlen2, int multi_mapping_number, int this_multi_mapping_i, int non_informative_subreads_r1, int non_informative_subreads_r2){

	int is_2_OK = 0, is_1_OK = 0;

	if(res1){
		is_1_OK = convert_read_to_tmp(global_context , out_context, read_number,  0, rlen1, read_text_1, qual_text_1, res1, out_context -> r1, read_name_1);
		if(is_1_OK) add_realignment_event_support(global_context, res1);
	}
	if(res2){
		is_2_OK = convert_read_to_tmp(global_context , out_context, read_number,  1, rlen2, read_text_2, qual_text_2, res2, out_context -> r2, read_name_2);
		if(is_2_OK) add_realignment_event_support(global_context, res2);
	}

	subread_output_tmp_t * r1_output = NULL;
	subread_output_tmp_t * r2_output = NULL;

	if(res1){
		r1_output = out_context -> r1;
	}

	if(res2){
		r2_output = out_context -> r2;
	}

	if(this_multi_mapping_i < 1){
		if( is_1_OK == 0 && is_2_OK == 0)
			if(thread_context) thread_context -> all_unmapped_reads ++;
			else global_context -> all_unmapped_reads ++;
		else if( is_1_OK == 0 || is_2_OK == 0){
			if(thread_context) thread_context -> not_properly_pairs_only_one_end_mapped ++;
			else global_context -> not_properly_pairs_only_one_end_mapped ++;

			if((is_1_OK && (res1 -> realign_flags & CORE_IS_BREAKEVEN )) ||
			   (is_2_OK && (res2 -> realign_flags & CORE_IS_BREAKEVEN))) {
				if(thread_context) thread_context -> all_multimapping_reads ++;
				else global_context -> all_multimapping_reads ++;
			}else{
				if(thread_context) thread_context -> all_uniquely_mapped_reads ++;
				else global_context -> all_uniquely_mapped_reads ++;
			}
		}else{
			if(res1 -> realign_flags & CORE_IS_BREAKEVEN){
				if(thread_context) thread_context -> all_multimapping_reads ++;
				else global_context -> all_multimapping_reads ++;
			}else{
				if(thread_context) thread_context -> all_uniquely_mapped_reads ++;
				else global_context -> all_uniquely_mapped_reads ++;
			}
		}
	}


	if((!global_context->config.ignore_unmapped_reads) || (is_2_OK || is_1_OK))
		write_single_fragment(global_context, thread_context, r1_output, res1, r2_output, res2, multi_mapping_number , this_multi_mapping_i , read_name_1, read_name_2, rlen1, rlen2, read_text_1, read_text_2, qual_text_1, qual_text_2, read_number, non_informative_subreads_r1, non_informative_subreads_r2, is_1_OK, is_2_OK);
}


void clear_repeated_buffer(global_context_t * global_context, unsigned int * repeated_buffer_position, char ** repeated_buffer_cigar, int * repeated_count){
	(*repeated_count) = 0;
}
int add_repeated_buffer(global_context_t * global_context, unsigned int * repeated_buffer_position, char ** repeated_buffer_cigar, int * repeated_count, realignment_result_t * res1, realignment_result_t * res2){
	int x1, is_repeated = 0;
	char * r1_cigar = "*";
	unsigned int r1_pos = 0;

	if(res1){
		r1_cigar = res1 -> cigar_string;
		r1_pos = res1 -> first_base_position;
	}

	char * r2_cigar = "*";
	unsigned int r2_pos = 0;

	if(res2){
		r2_cigar = res2 -> cigar_string;
		r2_pos = res2 -> first_base_position;
	}

	for(x1 = 0; x1 < (*repeated_count); x1 += 2){
		//if(strcmp(r1_cigar, "48M2769N53M")==0)
		//	SUBREADprintf("KYKY %u=%u  %s=%s\n" , r1_pos, repeated_buffer_position[x1], r1_cigar , repeated_buffer_cigar[x1]);
		if(repeated_buffer_position[x1] == r1_pos && repeated_buffer_position[x1+1] == r2_pos)
			if(strcmp(repeated_buffer_cigar[x1], r1_cigar) == 0 && strcmp(repeated_buffer_cigar[x1+1], r2_cigar) == 0)
			{
				is_repeated = 1;
				break;
			}
	}

	if( (!is_repeated) && (*repeated_count) < MAX_ALIGNMENT_PER_ANCHOR * 2 * global_context -> config.reported_multi_best_reads){


		//if(strcmp(r1_cigar, "48M2769N53M")==0)
		//	SUBREADprintf("CPCP %u=%u  %s=%s\n" , r1_pos, repeated_buffer_position[x1], r1_cigar , repeated_buffer_cigar[x1]);
		repeated_buffer_position[*repeated_count] = r1_pos;
		repeated_buffer_position[1 + *repeated_count] = r2_pos;
		strcpy(repeated_buffer_cigar[*repeated_count], r1_cigar);
		strcpy(repeated_buffer_cigar[1 + *repeated_count], r2_cigar);
		(*repeated_count) +=2;
	}

	return is_repeated;
}
int do_iteration_two(global_context_t * global_context, thread_context_t * thread_context)
{
	gene_input_t * ginp1 = NULL , * ginp2 = NULL;
	subread_read_number_t current_read_number=0;
	char read_text_1[MAX_READ_LENGTH+1], read_text_2[MAX_READ_LENGTH+1];
	char qual_text_1[MAX_READ_LENGTH+1], qual_text_2[MAX_READ_LENGTH+1];

	char raw_read_text_1[MAX_READ_LENGTH+1], raw_read_text_2[MAX_READ_LENGTH+1];
	char raw_qual_text_1[MAX_READ_LENGTH+1], raw_qual_text_2[MAX_READ_LENGTH+1];
	char read_name_1[MAX_READ_NAME_LEN+1], read_name_2[MAX_READ_NAME_LEN+1];
	char * repeated_buffer_cigars[MAX_ALIGNMENT_PER_ANCHOR *  2 * global_context -> config.reported_multi_best_reads];
	unsigned int repeated_buffer_pos[MAX_ALIGNMENT_PER_ANCHOR *  2 * global_context -> config.reported_multi_best_reads];
	int repeated_count;
	int read_len_1=0, read_len_2=0;
	int sqr_interval, sqr_read_number=0;

	int * final_MATCH_buffer1, *final_MISMATCH_buffer1, *final_PENALTY_buffer1;
	int * final_MATCH_buffer2, *final_MISMATCH_buffer2, *final_PENALTY_buffer2, non_informative_subreads_r1=0, non_informative_subreads_r2=0;
	int * final_realignment_index1, *final_realignment_index2;
	unsigned int *final_realignment_number;
	unsigned long long * final_SCORE_buffer;
	unsigned int * final_TLEN_buffer;

	realignment_result_t * final_realignments;
	subread_output_context_t out_context;

	//SUBREADprintf("THREAD_OPT2_start\n");
	init_output_context(global_context, &out_context);	

	for(repeated_count = 0;repeated_count < MAX_ALIGNMENT_PER_ANCHOR  *  2 * global_context -> config.reported_multi_best_reads ; repeated_count ++ ){
		repeated_buffer_cigars[repeated_count] = malloc(2*CORE_MAX_CIGAR_STR_LEN);
	}

	final_MATCH_buffer1 = malloc(sizeof(int) * 2 * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR);
	final_MISMATCH_buffer1 = malloc(sizeof(int) * 2 * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR);
	final_PENALTY_buffer1 = malloc(sizeof(int) * 2 * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR);
	final_realignment_index1 = malloc(sizeof(int) * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR);

	final_MATCH_buffer2 = malloc(sizeof(int) * 2 * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR);
	final_MISMATCH_buffer2 = malloc(sizeof(int) * 2 * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR);
	final_PENALTY_buffer2 = malloc(sizeof(int) * 2 * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR);
	final_realignment_index2 = malloc(sizeof(int) * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR);

	final_SCORE_buffer = malloc(sizeof(long long) * global_context -> config.multi_best_reads * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR*MAX_ALIGNMENT_PER_ANCHOR);
	final_TLEN_buffer = malloc(sizeof(int) * global_context -> config.multi_best_reads * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR*MAX_ALIGNMENT_PER_ANCHOR);

	final_realignments = malloc(sizeof(realignment_result_t) * global_context -> config.multi_best_reads * 2 * MAX_ALIGNMENT_PER_ANCHOR);
	final_realignment_number = malloc(sizeof(int) * global_context -> config.multi_best_reads * 2);

	mapping_result_t * r1_align_result_buffer, * r2_align_result_buffer;
	subjunc_result_t * r1_subjunc_result_buffer, * r2_subjunc_result_buffer;

	r1_align_result_buffer = malloc(sizeof(mapping_result_t) * global_context -> config.multi_best_reads);
	r2_align_result_buffer = malloc(sizeof(mapping_result_t) * global_context -> config.multi_best_reads);

	r1_subjunc_result_buffer = malloc(sizeof(subjunc_result_t) * global_context -> config.multi_best_reads);
	r2_subjunc_result_buffer = malloc(sizeof(subjunc_result_t) * global_context -> config.multi_best_reads);

	//unsigned int low_index_border = global_context -> current_value_index -> start_base_offset;
	//unsigned int high_index_border = global_context -> current_value_index -> start_base_offset + global_context -> current_value_index -> length; 

	init_chunk_scanning_parameters(global_context,thread_context, & ginp1, & ginp2);
	sqr_interval = max(5000,global_context -> processed_reads_in_chunk/10/ global_context -> config.all_threads);


	while(1)
	{
		int is_second_read;
		int max_votes;

		sqr_read_number++;
		fetch_next_read_pair(global_context, thread_context, ginp1, ginp2, &read_len_1, &read_len_2, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2, 0, &current_read_number);
		//if(current_read_number%500000==0)SUBREADprintf("THREAD_OPT2_fetch %s  %s ; RNO=%lld\n", read_name_1, read_text_1, current_read_number);
		strcpy(raw_read_text_1, read_text_1);
		strcpy(raw_qual_text_1, qual_text_1);
		//printf("OCT27-STEPB - %s\n", read_name_1);

		if(global_context -> input_reads.is_paired_end_reads){
			strcpy(raw_read_text_2, read_text_2);
			strcpy(raw_qual_text_2, qual_text_2);
		}

		if(global_context->config.space_type == GENE_SPACE_COLOR){
			if(isalpha(read_text_1[0])){
				//SUBREADprintf("FIRST : %s\n\n", read_text_1);
				int xk1;
				for(xk1=2; read_text_1[xk1]; xk1++){
					read_text_1[xk1-2] = read_text_1[xk1];
				}
				read_text_1[xk1-2] = 0;
			}

			if(global_context -> input_reads.is_paired_end_reads){
				//SUBREADprintf("SECOND: %s\n\n", read_text_2);
				if(isalpha(read_text_2[0])){
					int xk1;
					for(xk1=1; read_text_2[xk1]; xk1++){
						read_text_2[xk1-1] = read_text_2[xk1];
					}
					read_text_2[xk1-1] = 0;
				}
			}
		}

		if(raw_qual_text_1[0] && global_context -> config.phred_score_format == FASTQ_PHRED64){
			fastq_64_to_33(raw_qual_text_1);
			if(global_context->input_reads.is_paired_end_reads)
				fastq_64_to_33(raw_qual_text_2);
		}

		if((global_context -> config.is_BAM_output && global_context -> output_bam_writer -> is_internal_error) ||  global_context -> input_reads.is_internal_error || 
		   (global_context -> output_sam_is_full))break;
		if(current_read_number < 0) break;

		// if no more reads
		if( global_context -> input_reads.is_paired_end_reads)
			max_votes = max(_global_retrieve_alignment_ptr(global_context, current_read_number, 0, 0)->selected_votes, _global_retrieve_alignment_ptr(global_context, current_read_number, 1, 0)->selected_votes);
		else	max_votes = _global_retrieve_alignment_ptr(global_context, current_read_number, 0, 0)->selected_votes;

		//if(  current_read_number < 10 ) SUBREADprintf("BSSS2 : %s : %d\n", read_name_1, max_votes);

		int best_read_id=0;

		//memset(final_MATCH_buffer, 0, sizeof(int) * 2 * global_context -> config.multi_best_reads);
		//memset(final_MISMATCH_buffer, 0, sizeof(int) * 2 * global_context -> config.multi_best_reads);
		int r1_candidate_locations = 0, r2_candidate_locations = 0;
		int r1_step2_locations = 0, r2_step2_locations = 0;

		clear_repeated_buffer(global_context, repeated_buffer_pos, repeated_buffer_cigars, &repeated_count);
		for (is_second_read = 0; is_second_read < 1 + global_context -> input_reads.is_paired_end_reads; is_second_read ++)
		{
			char * current_read_name =  is_second_read?read_name_2 : read_name_1;
			int is_reversed_already = 0, realignment_i;
			int * current_candidate_locations = is_second_read?&r2_candidate_locations:&r1_candidate_locations;

			int * current_MATCH_buffer = is_second_read?final_MATCH_buffer2:final_MATCH_buffer1;
			int * current_MISMATCH_buffer = is_second_read?final_MISMATCH_buffer2:final_MISMATCH_buffer1;
			int * current_PENALTY_buffer = is_second_read?final_PENALTY_buffer2:final_PENALTY_buffer1;
			int * current_realignment_index = is_second_read?final_realignment_index2:final_realignment_index1;

	
			for(best_read_id = 0; best_read_id < global_context -> config.multi_best_reads; best_read_id++)
			{
				mapping_result_t *current_result = _global_retrieve_alignment_ptr(global_context, current_read_number, is_second_read, best_read_id);

				if(best_read_id == 0){
					if(is_second_read)
						non_informative_subreads_r2 = current_result -> noninformative_subreads_in_vote;
					else
						non_informative_subreads_r1 = current_result -> noninformative_subreads_in_vote;
				}

				char * current_read =  is_second_read?read_text_2 : read_text_1;
				char * current_qual =  is_second_read?qual_text_2 : qual_text_1;
				int current_rlen = is_second_read?read_len_2:read_len_1;

				if(current_result->selected_votes < global_context->config.minimum_subread_for_second_read || max_votes < global_context->config.minimum_subread_for_first_read)
				{
				//	if(current_read_number == 111 || current_read_number == 112)
					//SUBREADprintf("RESET0 [%lld] R_%d SEL=%d MAX=%d\n", current_read_number, is_second_read + 1, current_result->selected_votes, max_votes);
					current_result -> selected_votes = 0;
					continue;
				}
				if(locate_current_value_index(global_context, thread_context, current_result, current_rlen))
				{
					current_result -> selected_votes = 0;
					//SUBREADprintf("RESET1 Read position excesses index boundary.\n");
					continue;
				}

				//#warning ">>>>>>> COMMENT THIS <<<<<<<"
				//printf("OCT27-STEP21-%s:%d-ALN%02d-THRE %d in MM=%d\n", current_read_name, is_second_read + 1, best_read_id, thread_context -> thread_id, global_context -> config.multi_best_reads);

				int is_negative_strand = (current_result  -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;

				if(is_negative_strand + is_reversed_already == 1)
				{
					reverse_read(current_read, current_rlen,  global_context->config.space_type);
					reverse_quality(current_qual, current_rlen);
					is_reversed_already = !is_reversed_already;
				}

				current_result -> result_flags &= ~CORE_IS_FULLY_EXPLAINED;


				if(is_second_read) r2_step2_locations = best_read_id + 1;
				else r1_step2_locations = best_read_id + 1;

				unsigned int final_alignments = explain_read(global_context, thread_context , final_realignments + (is_second_read + 2 * best_read_id) * MAX_ALIGNMENT_PER_ANCHOR,
							   current_read_number, current_rlen, current_read_name, current_read, current_qual, is_second_read, best_read_id, is_negative_strand);

				final_realignment_number[ best_read_id * 2 + is_second_read ] = final_alignments;

				for(realignment_i = 0 ; realignment_i < final_alignments ; realignment_i ++){
					if((* current_candidate_locations) >= global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR) break;

					int realign_index = (is_second_read + 2 * best_read_id)* MAX_ALIGNMENT_PER_ANCHOR + realignment_i;
					int final_MATCH = final_realignments[realign_index].final_matched_bases;
					int final_PENALTY = final_realignments[realign_index].final_penalty;

					if((current_result -> result_flags & CORE_IS_FULLY_EXPLAINED) && final_MATCH>0) {
						int final_MISMATCH = final_realignments[realign_index].final_mismatched_bases;
						//#warning ">>>>>>>>>>>>>>>>>>>>>>>>>>>>> REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
						//printf("OCT27-STEPMSM-MMM %s M=%d MM=%d\n", read_name_1, final_MATCH ,final_MISMATCH);

						current_PENALTY_buffer[*current_candidate_locations] = final_PENALTY;
						current_MATCH_buffer[*current_candidate_locations] = final_MATCH;
						current_MISMATCH_buffer[*current_candidate_locations] = final_MISMATCH;
						current_realignment_index[*current_candidate_locations] = realign_index;
						(*current_candidate_locations) ++;
					}
				}
				//#warning ">>>>>>> COMMENT THIS <<<<<<<"
				//printf("OCT27-FI-%s:%d-ALN%02d-THRE %d\n", current_read_name, is_second_read + 1, best_read_id, thread_context -> thread_id);
			}
		}

		//if(161430 <= current_read_number) SUBREADprintf("LOC1=%d, LOC2=%d\n", r1_candidate_locations, r2_candidate_locations);

		int need_expect_TLEN = r2_candidate_locations && r1_candidate_locations && global_context -> config.reported_multi_best_reads<2 && global_context -> expected_TLEN_read_numbers < READPAIRS_FOR_CALC_EXPT_TLEN;
		int output_cursor = 0;
		if(r2_candidate_locations == 0 || r1_candidate_locations == 0) {
			int is_second_read, highest_score_occurence =0;
			for(is_second_read = 0; is_second_read < 1+ global_context -> input_reads.is_paired_end_reads; is_second_read++)
			{
				int current_candidate_locations = is_second_read?r2_candidate_locations:r1_candidate_locations;
				if(current_candidate_locations > 0){

					int * current_MATCH_buffer = is_second_read?final_MATCH_buffer2:final_MATCH_buffer1;
					int * current_PENALTY_buffer = is_second_read?final_PENALTY_buffer2:final_PENALTY_buffer1;
					int * current_MISMATCH_buffer = is_second_read?final_MISMATCH_buffer2:final_MISMATCH_buffer1;
					int * current_realignment_index = is_second_read?final_realignment_index2:final_realignment_index1;

					int read_record_i;
					unsigned long long int best_score_highest = 0;
					unsigned long long int scores_array [global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR];

					for(read_record_i = 0; read_record_i < current_candidate_locations; read_record_i++){
						realignment_result_t * current_realignment_result = final_realignments + current_realignment_index[read_record_i];

						unsigned int this_MATCH = current_MATCH_buffer[read_record_i];
						unsigned int this_MISMATCH = current_MISMATCH_buffer[read_record_i];
						unsigned int this_PENALTY = current_PENALTY_buffer[read_record_i];
						unsigned long long this_SCORE;

						if(global_context -> config.experiment_type == CORE_EXPERIMENT_DNASEQ){
							this_SCORE = this_MATCH * 100000llu + (10000 - this_MISMATCH);
						}else{
							unsigned int skip = 0; int is_exonic_regions = 1;
							if(global_context -> exonic_region_bitmap)calc_end_pos(current_realignment_result -> first_base_position, current_realignment_result -> cigar_string, &skip, &is_exonic_regions, global_context  );

							if(global_context -> config.scRNA_input_mode)
								this_SCORE =((100000llu * (10000 - this_MISMATCH + 2*is_exonic_regions) + this_MATCH)*50llu - this_PENALTY)*20llu+ current_realignment_result -> known_junction_supp;
							else this_SCORE =((100000llu * (10000 - this_MISMATCH) + this_MATCH)*50llu - this_PENALTY)*20llu+ current_realignment_result -> known_junction_supp;
						}

						best_score_highest = max(best_score_highest, this_SCORE);
						scores_array[read_record_i] = this_SCORE;
					}

					for(read_record_i = 0; read_record_i <  current_candidate_locations ; read_record_i++){
						realignment_result_t * current_realignment_result = final_realignments + current_realignment_index[read_record_i];
						if( scores_array[read_record_i] >= best_score_highest && (current_realignment_result -> realign_flags & CORE_TOO_MANY_MISMATCHES)==0) {
							int is_repeated = 0;

							is_repeated = add_repeated_buffer(global_context, repeated_buffer_pos, repeated_buffer_cigars, &repeated_count, is_second_read?NULL:current_realignment_result , is_second_read?current_realignment_result:NULL);
							if(is_repeated)
								scores_array[read_record_i] = 0;
							else
								highest_score_occurence ++;
						}
					}

					if(highest_score_occurence<2 ||  global_context -> config.report_multi_mapping_reads){
						int is_break_even = 0;
						if(highest_score_occurence>1)	is_break_even = 1;

						highest_score_occurence = min(highest_score_occurence, global_context -> config.reported_multi_best_reads);
						for(read_record_i = 0; read_record_i <  current_candidate_locations ; read_record_i++){
							realignment_result_t * current_realignment_result = final_realignments + current_realignment_index[read_record_i];

							if( scores_array[read_record_i] >= best_score_highest && (current_realignment_result -> realign_flags & CORE_TOO_MANY_MISMATCHES)==0
							  && output_cursor < global_context -> config.reported_multi_best_reads){
								strcpy(read_text_1, raw_read_text_1);
								strcpy(read_text_2, raw_read_text_2);
								strcpy(qual_text_1, raw_qual_text_1);
								strcpy(qual_text_2, raw_qual_text_2);

								if(is_break_even) current_realignment_result -> realign_flags |= CORE_IS_BREAKEVEN; 
								current_realignment_result -> MAPQ_adjustment = current_MISMATCH_buffer [read_record_i] + ( is_second_read?(r2_step2_locations): (r1_step2_locations));

//		SUBREADprintf("THREAD_OPT1_write %s\n", read_name_1);
								if(is_second_read)
									write_realignments_for_fragment(global_context, thread_context, &out_context, current_read_number, NULL, current_realignment_result, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2, read_len_1, read_len_2, highest_score_occurence, output_cursor,  non_informative_subreads_r1, non_informative_subreads_r2);
								else write_realignments_for_fragment(global_context, thread_context, &out_context , current_read_number, current_realignment_result, NULL, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2, read_len_1, read_len_2, highest_score_occurence, output_cursor,  non_informative_subreads_r1, non_informative_subreads_r2);
								output_cursor ++;
							}
						}

						assert(output_cursor == highest_score_occurence);
					}
				}
			}
		} else {
			int r1_best_id, r2_best_id, highest_score_occurence = 0;
			int expected_TLEN;
			if(global_context -> expected_TLEN_read_numbers >= READPAIRS_FOR_CALC_EXPT_TLEN)
				expected_TLEN = global_context -> expected_TLEN_sum / global_context -> expected_TLEN_read_numbers;
			else expected_TLEN = (global_context -> config.minimum_pair_distance + global_context -> config.maximum_pair_distance)/2;

			unsigned long long highest_score = 0;
			if(need_expect_TLEN) memset(final_TLEN_buffer, 0 , sizeof(int)  * global_context -> config.multi_best_reads * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR*MAX_ALIGNMENT_PER_ANCHOR );
			memset(final_SCORE_buffer, 0 , sizeof(long long)  * global_context -> config.multi_best_reads * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR*MAX_ALIGNMENT_PER_ANCHOR );
			for(r1_best_id = 0; r1_best_id < r1_candidate_locations; r1_best_id ++) {
				int r1_matched = final_MATCH_buffer1[r1_best_id];
				if(r1_matched < 1) continue;
				realignment_result_t * realignment_result_R1 = final_realignments + final_realignment_index1[r1_best_id];

				for(r2_best_id = 0; r2_best_id < r2_candidate_locations; r2_best_id ++) {
					int r2_matched = final_MATCH_buffer2[r2_best_id];
					if(r2_matched < 1) continue;

					realignment_result_t * realignment_result_R2 = final_realignments + final_realignment_index2[r2_best_id];
					int is_PE = 0, tlen=0;
					int is_same_chro = 0, is_exonic_regions = 0;
					unsigned long long final_SCORE = 0;

					test_PE_and_same_chro_align(global_context , realignment_result_R1 , realignment_result_R2, &is_exonic_regions, &is_PE, &is_same_chro , read_len_1, read_len_2, read_name_1, &tlen);

					unsigned long long TLEN_exp_score = 0;
					if(is_PE && global_context -> config.no_TLEN_preference == 0 && global_context -> config.reported_multi_best_reads<2){
						TLEN_exp_score = (tlen > expected_TLEN)?(tlen - expected_TLEN):(expected_TLEN - tlen);
						if(TLEN_exp_score>999) TLEN_exp_score=0;
						else TLEN_exp_score = 999-TLEN_exp_score;
					}
					if(global_context -> config.experiment_type == CORE_EXPERIMENT_DNASEQ){
						int weight;
						
						//#warning " ============ USE THE FIRST THREE WEIGHTS! ======== "

						if(is_PE)
							weight = 120;
							//weight = 300;
						else if(is_same_chro)
							weight = 100;
						else	weight = 80;
							//weight = 30;
						final_SCORE = weight * (final_MATCH_buffer1[r1_best_id] + final_MATCH_buffer2[r2_best_id]);
						//#warning "=========== ADD BY YANG LIAO FOR MORE MAPPED READS WITH '-u' OPTION ================"
						final_SCORE = final_SCORE * 1000llu - final_MISMATCH_buffer1[r1_best_id] - final_MISMATCH_buffer2[r2_best_id];
						final_SCORE = final_SCORE * 20llu - final_PENALTY_buffer1[r1_best_id] - final_PENALTY_buffer2[r2_best_id];
						final_SCORE = final_SCORE * 1000llu + TLEN_exp_score;
						if(is_PE && need_expect_TLEN) final_TLEN_buffer[r1_best_id * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR + r2_best_id] = tlen;
						//SUBREADprintf("%s FINAL_SCORE=%llu  TLEN_SCORE=%llu    DIST=%d    EXP_TLEN=%d\n", read_name_1, final_SCORE, TLEN_exp_score, tlen, expected_TLEN);
					} else if (global_context -> config.experiment_type == CORE_EXPERIMENT_RNASEQ) {
						int weight;


						if(1){

							if(is_exonic_regions && is_PE)
								weight = 5000;
							else if(is_PE || is_exonic_regions)
								weight = 3000;
							else if(is_same_chro)
								weight = 1000;
							else	weight = 300;
						}else{
							if(is_PE)
								weight = 3000;
							else if(is_same_chro)
								weight = 1000;
							else	weight = 300;

						}

						//#warning "=========== ADD BY YANG LIAO  ' + 2' ===================="
						final_SCORE = weight / (final_MISMATCH_buffer1[r1_best_id] + final_MISMATCH_buffer2[r2_best_id] + 1 + 2);

						//#warning "=========== ADD BY YANG LIAO FOR MORE MAPPED READS WITH '-u' OPTION ================"
						final_SCORE = final_SCORE * 3000llu + (final_MATCH_buffer1[r1_best_id] + final_MATCH_buffer2[r2_best_id]);
						final_SCORE = final_SCORE * 20 + realignment_result_R2 -> known_junction_supp + realignment_result_R1 -> known_junction_supp;
						final_SCORE = final_SCORE * 20 - final_PENALTY_buffer1[r1_best_id] - final_PENALTY_buffer2[r2_best_id];
						final_SCORE = final_SCORE * 1000 + TLEN_exp_score;

						//#warning ">>>>>>>>>>>>>>>>>>>>>>>>>>>>> REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
						//printf("OCT27-STEPMSM-FNSC %s M=%d,%d MM=%d,%d  CHROSAME=%d FSCR=%llu\n", read_name_1, final_MATCH_buffer1[r1_best_id], final_MATCH_buffer2[r2_best_id], final_MISMATCH_buffer1[r1_best_id] , final_MISMATCH_buffer2[r2_best_id], is_same_chro, final_SCORE );
						if(is_PE && need_expect_TLEN) final_TLEN_buffer[r1_best_id * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR + r2_best_id] = tlen;
					} else assert(0);

					assert(final_SCORE > 0);
					final_SCORE_buffer[r1_best_id * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR + r2_best_id] = final_SCORE;
	
					if(final_SCORE > highest_score) {
						//#warning ">>>>>>>>>>>> COMMENT THIS <<<<<<<<<<<<<<<<<<<"
						//printf("OCT27-STEPMSM-REPL %s : SCORE %llu -> %llu, pos=%u,%u\n", read_name_1, final_SCORE, highest_score , realignment_result_R1->mapping_result->selected_position, realignment_result_R2->mapping_result->selected_position );
						highest_score_occurence = 1;
						highest_score = final_SCORE;
						//SUBREADprintf("29NOV2016-FIRSTBEST %s [%d : score %llu] LOC = %d,%d\n", read_name_1, highest_score_occurence, final_SCORE, r1_best_id, r2_best_id);
						clear_repeated_buffer(global_context, repeated_buffer_pos, repeated_buffer_cigars, &repeated_count);
						add_repeated_buffer(global_context, repeated_buffer_pos, repeated_buffer_cigars, &repeated_count, realignment_result_R1, realignment_result_R2);
					} else if(final_SCORE == highest_score) {
						int is_repeat = 0;

						if(global_context -> config.reported_multi_best_reads)
							is_repeat = add_repeated_buffer(global_context, repeated_buffer_pos, repeated_buffer_cigars, &repeated_count, realignment_result_R1, realignment_result_R2);

						if(is_repeat)
							final_SCORE_buffer[r1_best_id * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR + r2_best_id] = 0;
						else{
							highest_score_occurence ++;
						//	SUBREADprintf("29NOV2016-BEST %s [%d : score %llu] LOC = %d,%d\n", read_name_1, highest_score_occurence, final_SCORE, r1_best_id, r2_best_id);
						}
					}
				}
			}

			//SUBREADprintf("Highest score = %llu, Occurance = %d\n", highest_score , highest_score_occurence);
			// Then, copy the (R1, R2) that have the highest score into the align_res buffer.
			
			//#warning ">>>>>>> COMMENT THIS <<<<<<<"
			//printf("OCT27-WRITE-BOTHMAOOED-%s-THRE %d ; occu=%d\n", read_name_1, thread_context -> thread_id, highest_score_occurence);

			if(highest_score_occurence <= 1 || global_context -> config.report_multi_mapping_reads){

				int is_break_even = 0;
				if(highest_score_occurence>1)is_break_even = 1;

				highest_score_occurence = min(highest_score_occurence, global_context -> config.reported_multi_best_reads);
				for(r1_best_id = 0; r1_best_id < r1_candidate_locations; r1_best_id ++)
				{
					int r1_matched = final_MATCH_buffer1[r1_best_id];
					if(r1_matched < 1) continue;

					for(r2_best_id = 0; r2_best_id < r2_candidate_locations; r2_best_id ++)
					{
						int r2_matched = final_MATCH_buffer2[r2_best_id];
						if(r2_matched < 1) continue;

						if(final_SCORE_buffer[r1_best_id * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR + r2_best_id] == highest_score && 
							output_cursor < global_context -> config.reported_multi_best_reads){
								if(need_expect_TLEN){
									int this_tlen = final_TLEN_buffer[r1_best_id * global_context -> config.multi_best_reads * MAX_ALIGNMENT_PER_ANCHOR + r2_best_id];
									if(this_tlen>0){
										subread_lock_occupy(&global_context -> output_lock);
										global_context -> expected_TLEN_read_numbers++;
										global_context -> expected_TLEN_sum += this_tlen;
										if(global_context -> expected_TLEN_read_numbers == READPAIRS_FOR_CALC_EXPT_TLEN)
											print_in_box(80,0,0,"  Estimated fragment length : %d bp\n",  (int)(global_context -> expected_TLEN_sum / global_context -> expected_TLEN_read_numbers ));
										subread_lock_release(&global_context -> output_lock);
									}
								}
								realignment_result_t * r1_realign = final_realignments + final_realignment_index1[r1_best_id];
								realignment_result_t * r2_realign = final_realignments + final_realignment_index2[r2_best_id];

								strcpy(read_text_1, raw_read_text_1);
								strcpy(read_text_2, raw_read_text_2);
								strcpy(qual_text_1, raw_qual_text_1);
								strcpy(qual_text_2, raw_qual_text_2);

								if(is_break_even){
									r1_realign -> realign_flags |= CORE_IS_BREAKEVEN;
									r2_realign -> realign_flags |= CORE_IS_BREAKEVEN;
								}

								r1_realign -> MAPQ_adjustment = r1_step2_locations + final_MISMATCH_buffer1[r1_best_id];
								r2_realign -> MAPQ_adjustment = r2_step2_locations + final_MISMATCH_buffer2[r2_best_id];

//		SUBREADprintf("THREAD_OPT2_write %s\n", read_name_1);
								write_realignments_for_fragment(global_context, thread_context, &out_context, current_read_number, r1_realign, r2_realign, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2, read_len_1, read_len_2, highest_score_occurence, output_cursor,  non_informative_subreads_r1, non_informative_subreads_r2);
								output_cursor ++;
						}
					}
				}
				assert(output_cursor == highest_score_occurence );
			}
		}

		//if(  current_read_number < 10 ) SUBREADprintf("BEEE2 : %s : %d\n", read_name_1, max_votes);
		//#warning ">>>>>>> COMMENT THIS <<<<<<<"
		//printf("OCT27-WRITE-UNMAP?-%s-THRE %d\n", read_name_1, thread_context -> thread_id);
		if(output_cursor<1) {
			strcpy(read_text_1, raw_read_text_1);
			strcpy(read_text_2, raw_read_text_2);
			strcpy(qual_text_1, raw_qual_text_1);
			strcpy(qual_text_2, raw_qual_text_2);
//		SUBREADprintf("THREAD_OPT3_write %s\n", read_name_1);
			write_realignments_for_fragment(global_context, thread_context, &out_context, current_read_number, NULL, NULL, read_name_1, read_name_2, read_text_1, read_text_2, raw_qual_text_1, raw_qual_text_2, read_len_1, read_len_2, 0, 0, non_informative_subreads_r1, non_informative_subreads_r2);
		}
		
		if(!thread_context || thread_context->thread_id == 0)
		{
			if(sqr_read_number>sqr_interval)
			{
				show_progress(global_context, thread_context, current_read_number, STEP_ITERATION_TWO);
				sqr_read_number=0;
			}
		}
		//#warning ">>>>>>> COMMENT THIS <<<<<<<"
		//printf("OCT27-FIN-%s-THRE %d\n", read_name_1, thread_context -> thread_id);
		//bigtable_release_result(global_context, thread_context, current_read_number, 1);
	}

	free(final_realignments);
	free(final_realignment_number);
	free(final_MATCH_buffer1);
	free(final_MISMATCH_buffer1);
	free(final_PENALTY_buffer1);
	free(final_realignment_index1);

	free(final_MATCH_buffer2);
	free(final_MISMATCH_buffer2);
	free(final_PENALTY_buffer2);
	free(final_realignment_index2);

	free(final_SCORE_buffer);
	free(final_TLEN_buffer);

	free(r1_align_result_buffer);
	free(r1_subjunc_result_buffer);
	free(r2_align_result_buffer);
	free(r2_subjunc_result_buffer);

	for(repeated_count = 0;repeated_count < MAX_ALIGNMENT_PER_ANCHOR  *  2 * global_context -> config.reported_multi_best_reads ; repeated_count ++ ){
		free(repeated_buffer_cigars[repeated_count]);
	}

	destroy_output_context(global_context, &out_context);	
	//SUBREADprintf("OPT2-FINISHED\n");
	if(global_context ->config.all_threads >1 && global_context ->config.is_BAM_output && !global_context ->config.is_input_read_order_required)SamBam_writer_finalise_thread(global_context -> output_bam_writer, thread_context -> thread_id);
	if(global_context ->config.all_threads <=1 && global_context ->config.is_BAM_output && global_context ->config.sort_reads_by_coordinates)SamBam_writer_finalise_one_thread(global_context -> output_bam_writer);
	return 0;
}

int core_get_subread_quality(global_context_t * global_context, thread_context_t * thread_context, char * qual, int qual_len)
{
	int x1, ret=1;

	if(!qual)return 1;
	if(!qual[0])return 1;

	int offset = (global_context->config.phred_score_format == FASTQ_PHRED33)?33:64; 

	for(x1=0; (x1 < qual_len) && qual[x1]; x1++)
		ret +=  (qual[x1] - offset);

	return  ret;
}

int has_better_mapping(global_context_t * global_context, thread_context_t * thread_context, subread_read_number_t current_read_number, int is_second_read, int this_aln_id){
	int better_read_id;

	mapping_result_t * this_r = _global_retrieve_alignment_ptr(global_context, current_read_number, is_second_read, this_aln_id);
	for(better_read_id = 0; better_read_id < this_aln_id; better_read_id++){
		mapping_result_t * better_r = _global_retrieve_alignment_ptr(global_context, current_read_number, is_second_read, better_read_id);
		if( this_r -> selected_position >= better_r -> selected_position - global_context -> config.max_indel_length - 1 
		 && this_r -> selected_position <=  better_r -> selected_position +  global_context -> config.max_indel_length +1 )
		if( this_r -> confident_coverage_start >= better_r -> confident_coverage_start
		 && this_r -> confident_coverage_end <= better_r -> confident_coverage_end) return 1;
	}
	return 0;
}

int do_voting(global_context_t * global_context, thread_context_t * thread_context)
{
	int xk1;
	gene_input_t * ginp1 = NULL , * ginp2 = NULL;
	subread_read_number_t current_read_number=0, last_shown_curr = 0;
	char * read_text_1, * read_text_2;
	char * qual_text_1, * qual_text_2;
	char read_name_1[MAX_READ_NAME_LEN+1], read_name_2[MAX_READ_NAME_LEN+1];
	int read_len_1, read_len_2=0;
	int min_first_read_votes = global_context -> config.minimum_subread_for_first_read; 
	int voting_max_indel_length = min(16, global_context->config.max_indel_length);
	int sqr_interval=10000, sqr_read_number = 0;

	read_text_1 = malloc(MAX_READ_LENGTH+1);
	read_text_2 = malloc(MAX_READ_LENGTH+1);
	qual_text_1 = malloc(MAX_READ_LENGTH+1);
	qual_text_2 = malloc(MAX_READ_LENGTH+1);

	gene_vote_t * vote_1, * vote_2, * vote_fg;

	vote_1 = (gene_vote_t *) malloc(sizeof(gene_vote_t));
	vote_2 = (gene_vote_t *) malloc(sizeof(gene_vote_t));
	vote_fg = (gene_vote_t *) malloc(sizeof(gene_vote_t));

	if(!(vote_1&&vote_2))
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_FATAL,"Cannot allocate voting memory.");
		return 1;
	}

	init_chunk_scanning_parameters(global_context,thread_context, & ginp1, & ginp2);

	unsigned int low_index_border = global_context -> current_value_index -> start_base_offset;
	unsigned int high_index_border = global_context -> current_value_index -> start_base_offset + global_context -> current_value_index -> length; 
	int has_second_read = 1 + global_context -> input_reads.is_paired_end_reads;

	if(thread_context)
		thread_context -> current_value_index = global_context -> current_value_index;

	int GENE_SLIDING_STEP = global_context -> current_index -> index_gap;
	int need_junction_step = global_context -> config.do_breakpoint_detection || global_context ->  config.do_fusion_detection ||global_context ->  config.do_long_del_detection;

	while(1)
	{
		int is_second_read;
		int subread_no;
		int is_reversed, applied_subreads = 0, v1_all_subreads=0, v2_all_subreads=0;

		fetch_next_read_pair(global_context, thread_context, ginp1, ginp2, &read_len_1, &read_len_2, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2,1, &current_read_number);
		//SUBREADprintf("VOTING READ # %lld: %s %s %d\n", current_read_number, read_name_1, read_text_1, read_len_1);

		if(current_read_number < 0){
		//	#warning ">>>>>>>>>>>>> COMMENT THIS <<<<<<<<<<<<<<<<<<<"
		//	printf("OCT27-STEPB-QUIT :T%d\n", thread_context -> thread_id);
			break;
		}

		for(is_reversed = 0; is_reversed<2; is_reversed++)
		{
			//printf("MAP_REA = %s / %s\n", read_text_1, read_text_2);

			for (is_second_read = 0; is_second_read < has_second_read; is_second_read ++)
			{
				gene_vote_t * current_vote = is_second_read?vote_2: vote_1;
				char * current_read =  is_second_read?read_text_2 : read_text_1;
				int current_rlen = is_second_read?read_len_2:read_len_1;
				int subread_step;
				if(current_rlen< 16) continue;
				int CR15GLS = (current_rlen - 15 - GENE_SLIDING_STEP)<<16;
				if(current_rlen<= EXON_LONG_READ_LENGTH)
				{
					subread_step =  CR15GLS /(global_context -> config.total_subreads -1);
					if(subread_step<(GENE_SLIDING_STEP<<16))subread_step = GENE_SLIDING_STEP<<16;
				}else{
					subread_step = 6<<16;
					if( CR15GLS / subread_step > 62)
						subread_step = CR15GLS/62;
				}

				int noninformative_subreads_for_each_gap[GENE_SLIDING_STEP];

				applied_subreads = 1 + CR15GLS / subread_step;
				if(is_second_read) v2_all_subreads = applied_subreads;
				else	 v1_all_subreads = applied_subreads;

				//SUBREADprintf("NSUBR=%d\tAPPLIED_SUBR=%d\tSTEP=%d\n", global_context -> config.total_subreads, applied_subreads, subread_step);

				unsigned int current_high_border = high_index_border -  current_rlen;

				if(global_context->config.do_breakpoint_detection && current_rlen > EXON_LONG_READ_LENGTH)//&& global_context->config.all_threads<2)
				{
					char * current_qual =  is_second_read?qual_text_2 : qual_text_1;
					core_fragile_junction_voting(global_context, thread_context, read_name_1, current_read, current_qual, current_rlen, is_reversed, global_context->config.space_type, low_index_border, current_high_border, vote_fg);

				}
				if(global_context->config.SAM_extra_columns)
				{
					for(xk1=0;xk1<GENE_SLIDING_STEP;xk1++)
						noninformative_subreads_for_each_gap[xk1]=0;
				}

				int allow_indel_i;
				unsigned int shift_indel_locs[ GENE_VOTE_TABLE_SIZE * GENE_VOTE_SPACE ];
				unsigned int shift_indel_NO = 0;
				for(allow_indel_i = 0; allow_indel_i<2; allow_indel_i++){
					if(is_reversed==0 || !(global_context-> config.do_fusion_detection || global_context-> config.do_long_del_detection)){
						if(is_second_read == 0)init_gene_vote(vote_1);
						if(is_second_read == 1)init_gene_vote(vote_2);
					}

					for(subread_no=0; subread_no < applied_subreads ; subread_no++)
					{
						for(xk1=0; xk1<GENE_SLIDING_STEP ; xk1++)
						{
	
							if(global_context->config.SAM_extra_columns)
							{
								current_vote -> noninformative_subreads = noninformative_subreads_for_each_gap[xk1];
							}
	
							int subread_offset = ((subread_step * subread_no) >> 16);
							if(GENE_SLIDING_STEP > 1)
								subread_offset -= subread_offset%(GENE_SLIDING_STEP) - xk1; 
	
							char * subread_string = current_read + subread_offset;
	
							gehash_key_t subread_integer = genekey2int(subread_string, global_context->config.space_type);
							if(global_context->config.is_methylation_reads)
								 gehash_go_q_CtoT(global_context->current_index, subread_integer , subread_offset, current_rlen, is_reversed, current_vote, 1, 0xffffff, voting_max_indel_length, subread_no, 1,  low_index_border, high_index_border - current_rlen);
							else
								 gehash_go_X(global_context->current_index, subread_integer , subread_offset, current_rlen, is_reversed, current_vote,  voting_max_indel_length, subread_no,  low_index_border, current_high_border, allow_indel_i, shift_indel_locs, &shift_indel_NO);
							if(global_context->config.SAM_extra_columns)
								noninformative_subreads_for_each_gap[xk1] = current_vote -> noninformative_subreads;
						}

					}
					if(shift_indel_NO == 0 || global_context-> config.do_fusion_detection || global_context-> config.do_long_del_detection)break;
				}

				if(global_context->config.SAM_extra_columns)
				{
					short max_noninformative_subreads = -1;

					for(xk1=0;xk1<GENE_SLIDING_STEP;xk1++)
						if(noninformative_subreads_for_each_gap[xk1] > max_noninformative_subreads)
							max_noninformative_subreads = noninformative_subreads_for_each_gap[xk1];

					current_vote -> noninformative_subreads = max_noninformative_subreads;
				}
			}



			if(is_reversed==1 || !(global_context-> config.do_fusion_detection || global_context-> config.do_long_del_detection))
			{
//#warning "====== CHECK PRINTING !!! =========="
				if(0&&current_read_number == 4000){
					SUBREADprintf(">>>%llu<<<\n%s [%d]  %s\n%s [%d]  %s  VOTE1_MAX=%d >= %d\n", current_read_number, read_name_1, read_len_1, read_text_1, read_name_2, read_len_2, read_text_2, vote_1->max_vote, min_first_read_votes);
					SUBREADprintf(" ======= PAIR %s = %llu ; NON_INFORMATIVE = %d, %d =======\n", read_name_1, current_read_number, vote_1 -> noninformative_subreads, vote_2 -> noninformative_subreads);
					print_votes(vote_1, global_context -> config.index_prefix);
					print_votes(vote_2, global_context -> config.index_prefix);
				}

				if(global_context -> input_reads.is_paired_end_reads)
					process_voting_junction(global_context, thread_context, current_read_number, vote_1, vote_2, read_name_1, read_name_2, read_text_1, read_text_2, read_len_1, read_len_2, is_reversed, v1_all_subreads, v2_all_subreads);
				else{
					if(vote_1->max_vote >= min_first_read_votes)
						process_voting_junction(global_context, thread_context, current_read_number, vote_1, vote_2, read_name_1, NULL ,  read_text_1, NULL, read_len_1, 0, is_reversed, v1_all_subreads, 0);
					else if(_global_retrieve_alignment_ptr(global_context, current_read_number, 0,0) -> selected_votes < 1)
					{
						mapping_result_t * allll = _global_retrieve_alignment_ptr(global_context, current_read_number, 0,0);
						allll -> noninformative_subreads_in_vote = 0;


						_global_retrieve_alignment_ptr(global_context, current_read_number, 0,0)->used_subreads_in_vote = max(_global_retrieve_alignment_ptr(global_context, current_read_number, 0,0)->used_subreads_in_vote, applied_subreads);
						_global_retrieve_alignment_ptr(global_context, current_read_number, 0,0)->noninformative_subreads_in_vote = max(_global_retrieve_alignment_ptr(global_context, current_read_number, 0,0)->noninformative_subreads_in_vote, vote_1 -> noninformative_subreads);
					}
				}
			}
		
			if(is_reversed == 0)
			{
				reverse_read(read_text_1, read_len_1,  global_context->config.space_type);
				if(global_context -> input_reads.is_paired_end_reads)
					reverse_read(read_text_2, read_len_2,  global_context->config.space_type);
			}
		}

		int read_1_reversed = 1;
		int read_2_reversed = 1;
		int best_read_id;

		if(global_context -> is_final_voting_run){
			for(is_second_read = 0; is_second_read < 1 + global_context -> input_reads.is_paired_end_reads; is_second_read ++)
			{
				int * this_read_has_reversed = is_second_read ? &read_2_reversed:&read_1_reversed;
				char * current_read_name =  is_second_read?read_name_2 : read_name_1;
				char * current_read =  is_second_read?read_text_2 : read_text_1;
				char * current_qual =  is_second_read?qual_text_2 : qual_text_1;
				int current_rlen = is_second_read?read_len_2:read_len_1;
	
				for(best_read_id = 0; best_read_id < global_context -> config.multi_best_reads; best_read_id++)
				{
					mapping_result_t * current_r = _global_retrieve_alignment_ptr(global_context, current_read_number, is_second_read,best_read_id);
					if(current_r -> selected_votes < 1) continue;
					
					int this_read_should_be_reversed = (current_r -> result_flags & CORE_IS_NEGATIVE_STRAND) ? 1:0;

					//SUBREADprintf("DETECT INDEL: should_reverse = %d, this_has_reversed = %d\n", this_read_should_be_reversed, *this_read_has_reversed);

					if(this_read_should_be_reversed != (*this_read_has_reversed))
					{
						(*this_read_has_reversed) = !(*this_read_has_reversed);
						reverse_read(current_read, current_rlen, global_context->config.space_type);
						if(current_qual)
							reverse_quality(current_qual , current_rlen);
					}

					gene_value_index_t * curr_val_index = thread_context? thread_context -> current_value_index: global_context -> current_value_index;
					locate_current_value_index(global_context, thread_context, current_r, current_rlen);

			//		#warning "==== UNCOMMENT THE NEXT LINE ===="

					if(0 && FIXLENstrcmp("R00000003493", current_read_name)==0)
						SUBREADprintf("TEST_BETTER: better_id = %d, votes = %d, has=%d\n", best_read_id, current_r -> selected_votes, has_better_mapping(global_context, thread_context, current_read_number,  is_second_read,best_read_id));

					if(!has_better_mapping(global_context, thread_context, current_read_number,  is_second_read,best_read_id))
						find_new_indels(global_context, thread_context, current_read_number, current_read_name, current_read, current_qual, current_rlen, is_second_read, best_read_id);
					if(need_junction_step)
						find_new_junctions(global_context, thread_context, current_read_number, current_read_name, current_read, current_qual, current_rlen, is_second_read, best_read_id);

			//		#warning "==== UNCOMMENT THE NEXT LINE ===="
					if(0 && FIXLENstrcmp("R00000003493", current_read_name)==0){
						SUBREADprintf("TEST_BETTER_END: better_id = %d, votes = %d, has=%d\n", best_read_id, current_r -> selected_votes, has_better_mapping(global_context, thread_context, current_read_number,  is_second_read,best_read_id));
						exit(1);
					}

					if(thread_context) thread_context -> current_value_index = curr_val_index;
					else global_context -> current_value_index = curr_val_index;
				}
			}
		}
		
		if(!thread_context || thread_context->thread_id == 0)
		{
			if(0 == global_context -> config.scRNA_input_mode  && sqr_read_number > sqr_interval)
			{
				show_progress(global_context, thread_context, current_read_number, STEP_VOTING);
				sqr_read_number = 0;
				unsigned long long total_file_size = global_context -> input_reads.first_read_file_size;
				unsigned long long guessed_all_reads = total_file_size / global_context -> input_reads . avg_read_length;// / (1+global_context -> config.is_SAM_file_input);
				sqr_interval = max(10000,guessed_all_reads / global_context -> config.all_threads/10);
			}
			if(global_context -> config.scRNA_input_mode && current_read_number - last_shown_curr >= 1000000){
				show_progress(global_context, thread_context, current_read_number + global_context -> all_processed_reads, STEP_VOTING);
				last_shown_curr = current_read_number;
			}
		}


		//bigtable_release_result(global_context, thread_context, current_read_number, 1);
		sqr_read_number++;

	}

	free(vote_1);
	free(vote_2);
	free(vote_fg);
	free(read_text_1);
	free(read_text_2);
	free(qual_text_1);
	free(qual_text_2);

	return 0;
}

void subread_init_topKbuff(global_context_t * global_context, topK_buffer_t * topKbuff){
	topKbuff -> vote_simple_1_buffer = malloc(global_context -> config.max_vote_simples * sizeof(simple_mapping_t));
	topKbuff -> vote_simple_2_buffer = malloc(global_context -> config.max_vote_simples * sizeof(simple_mapping_t));

	topKbuff -> junction_tmp_r1 = malloc(sizeof(subjunc_result_t) * global_context->config.multi_best_reads);
	topKbuff -> junction_tmp_r2 = malloc(sizeof(subjunc_result_t) * global_context->config.multi_best_reads);

	topKbuff -> alignment_tmp_r1 = malloc(sizeof(mapping_result_t) * global_context->config.multi_best_reads);
	topKbuff -> alignment_tmp_r2 = malloc(sizeof(mapping_result_t) * global_context->config.multi_best_reads);

	topKbuff -> comb_buffer = malloc(global_context -> config.max_vote_combinations * sizeof(vote_combination_t));
}

void subread_free_topKbuff(global_context_t * global_context, topK_buffer_t * topKbuff){
	
	free(topKbuff ->junction_tmp_r1);
	free(topKbuff ->junction_tmp_r2);
	free(topKbuff ->alignment_tmp_r1);
	free(topKbuff ->alignment_tmp_r2);
	free(topKbuff ->comb_buffer);
	free(topKbuff ->vote_simple_1_buffer);
	free(topKbuff ->vote_simple_2_buffer);
}



void * run_in_thread(void * pthread_param)
{
	void ** parameters = (void **)pthread_param;
	global_context_t * global_context = (global_context_t * ) parameters[0];
	thread_context_t * thread_context = (thread_context_t *) parameters[1];
	int task = *((int *)(parameters[2]));
	subread_lock_t * thread_init_lock = (subread_lock_t * ) parameters[3];
	int * ret_value_pointer = (int *)parameters[4];

	if(thread_init_lock)
		subread_lock_release(thread_init_lock);


	switch(task)
	{
		case STEP_VOTING:
			*ret_value_pointer = do_voting(global_context, thread_context);
		break;
		case STEP_ITERATION_TWO:
			*ret_value_pointer = do_iteration_two(global_context, thread_context);
		break;
	
	}

	//sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_DETAILS, "finished running %d", task);
	return NULL;
}

int run_maybe_threads(global_context_t *global_context, int task)
{
	void * thr_parameters [5];
	int ret_value =0;

	if(task == STEP_ITERATION_TWO){
		global_context -> last_written_fragment_number = -1;
	}


	if(global_context->config.all_threads<2) {
		thr_parameters[0] = global_context;
		thr_parameters[1] = NULL;
		thr_parameters[2] = &task;
		thr_parameters[3] = NULL;
		thr_parameters[4] = &ret_value;

		run_in_thread(thr_parameters);
		if(STEP_VOTING == task){
			// sort and merge events from all threads and the global event space.
			sort_global_event_table(global_context);
			// sort the event entry table at each location.
			//#warning "==== UNCOMMENT NEXT LINE ======"
			sort_junction_entry_table(global_context);
		}
	} else {
		int current_thread_no ;
		thread_context_t thread_contexts[64];
		int ret_values[64];

		memset(thread_contexts, 0, sizeof(thread_context_t)*64);
		global_context -> all_thread_contexts = thread_contexts;

		for(current_thread_no = 0 ; current_thread_no < global_context->config.all_threads ; current_thread_no ++) {
			thread_contexts[current_thread_no].thread_id = current_thread_no;
			init_indel_thread_contexts(global_context, thread_contexts+current_thread_no, task);
			if((global_context->config.do_breakpoint_detection || global_context-> config.do_fusion_detection || global_context->config.do_breakpoint_detection || global_context-> config.do_long_del_detection))
				init_junction_thread_contexts(global_context, thread_contexts+current_thread_no, task);

			if(STEP_VOTING == task) subread_init_topKbuff(global_context,&thread_contexts[current_thread_no].topKbuff);

			subread_lock_occupy(&global_context -> thread_initial_lock);
			thr_parameters[0] = global_context;
			thr_parameters[1] = thread_contexts+current_thread_no;
			thr_parameters[2] = &task;
			thr_parameters[3] = (void *)&global_context -> thread_initial_lock;
			thr_parameters[4] = ret_values + current_thread_no;

			pthread_create(&thread_contexts[current_thread_no].thread, NULL,  run_in_thread, &thr_parameters);
		}

		for(current_thread_no = 0 ; current_thread_no < global_context->config.all_threads ; current_thread_no ++) {
			pthread_join(thread_contexts[current_thread_no].thread, NULL);
			
			if(STEP_ITERATION_TWO == task){
				global_context -> all_mapped_reads += thread_contexts[current_thread_no].all_mapped_reads;
				global_context -> all_correct_PE_reads += thread_contexts[current_thread_no].all_correct_PE_reads;
				global_context -> not_properly_pairs_wrong_arrangement += thread_contexts[current_thread_no].not_properly_pairs_wrong_arrangement;
				global_context -> not_properly_pairs_different_chro += thread_contexts[current_thread_no].not_properly_pairs_different_chro;
				global_context -> not_properly_different_strands += thread_contexts[current_thread_no].not_properly_different_strands;
				global_context -> not_properly_pairs_TLEN_wrong += thread_contexts[current_thread_no].not_properly_pairs_TLEN_wrong;
				global_context -> all_unmapped_reads += thread_contexts[current_thread_no].all_unmapped_reads;
				global_context -> not_properly_pairs_only_one_end_mapped += thread_contexts[current_thread_no].not_properly_pairs_only_one_end_mapped;
				global_context -> all_multimapping_reads += thread_contexts[current_thread_no].all_multimapping_reads;
				global_context -> all_uniquely_mapped_reads += thread_contexts[current_thread_no].all_uniquely_mapped_reads;

			}
			if(STEP_VOTING == task) subread_free_topKbuff(global_context,&thread_contexts[current_thread_no].topKbuff);
			ret_value += *(ret_values + current_thread_no);
			if(ret_value)break;
		}

		// sort and merge events from all threads and the global event space.
		finalise_indel_and_junction_thread(global_context, thread_contexts, task);
		if(STEP_VOTING == task)
			// sort the event entry table at each location.
			sort_junction_entry_table(global_context);
	}

	if(CORE_SOFT_BR_CHAR == '\r')
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "");
	return ret_value;
}

void clean_context_after_chunk(global_context_t * context)
{
	context -> running_processed_reads_in_chunk = 0;
	context -> processed_reads_in_chunk = 0;
	init_bigtable_results(context, 1);

	indel_context_t * indel_context = (indel_context_t *)context -> module_contexts[MODULE_INDEL_ID];
	chromosome_event_t * event_space = indel_context -> event_space_dynamic;

	int event_id;
	//memset(context -> big_margin_record  , 0 , sizeof(*context -> big_margin_record) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.big_margin_record_size);
	for(event_id = 0; event_id < indel_context->total_events; event_id++){
		chromosome_event_t * event_body = event_space + event_id;
		event_body -> critical_read_id = 0xffffffffffffffffllu; 
	}
}

#define SKIP_CORE_NOEMPTY(fp_loc, buf_loc)	{ while(1){char *ret_loc = fgets(buf_loc, 3000, fp_loc); if(buf_loc[0]!='\n' || !ret_loc) break; } }

void locate_read_files(global_context_t * global_context, int type)
{
	// The BCL input module uses its own chunking algorithm.
	if(global_context -> input_reads.first_read_file.file_type == GENE_INPUT_BCL) return;

	if(type==SEEK_SET) {
		geinput_tell(&global_context -> input_reads.first_read_file, &global_context -> current_circle_start_position_file1);
		if(global_context ->input_reads.is_paired_end_reads)
			geinput_tell(&global_context -> input_reads.second_read_file, &global_context -> current_circle_start_position_file2);
	} else {
		geinput_tell(&global_context -> input_reads.first_read_file, &global_context -> current_circle_end_position_file1);
		if(global_context ->input_reads.is_paired_end_reads)
			geinput_tell(&global_context -> input_reads.second_read_file, &global_context -> current_circle_end_position_file2);
	}
	if(global_context -> input_reads.first_read_file.file_type == GENE_INPUT_SCRNA_FASTQ) return;
	if(global_context -> input_reads.first_read_file.file_type == GENE_INPUT_SCRNA_BAM) return;

	// This variable is only used to display the running progress.
	// That is why the scRNA-seq mode doesn't need to keep the start pos.
	if(type==SEEK_SET)
		global_context -> current_circle_start_abs_offset_file1 = geinput_file_offset(&(global_context -> input_reads.first_read_file));

}

void rewind_read_files(global_context_t * global_context, int type)
{
	if(type==SEEK_SET)
	{
		geinput_seek(&global_context -> input_reads.first_read_file, & global_context -> current_circle_start_position_file1);
		if(global_context ->input_reads.is_paired_end_reads)
			geinput_seek(&global_context -> input_reads.second_read_file, & global_context -> current_circle_start_position_file2);
	}
	else
	{
		geinput_seek(&global_context -> input_reads.first_read_file, & global_context -> current_circle_end_position_file1);
		if(global_context ->input_reads.is_paired_end_reads)
			geinput_seek(&global_context -> input_reads.second_read_file, & global_context -> current_circle_end_position_file2);
	
	}
}

void go_chunk_start( global_context_t * global_context ){
	if(global_context -> input_reads.first_read_file.file_type == GENE_INPUT_BCL)
		cacheBCL_go_chunk_start(&global_context -> input_reads.first_read_file.bcl_input);
	else 
		rewind_read_files(global_context, SEEK_SET);
	global_context -> running_processed_reads_in_chunk=0;
}

void go_chunk_nextchunk( global_context_t * global_context ){
	if(global_context -> input_reads.first_read_file.file_type == GENE_INPUT_BCL)
		cacheBCL_go_chunk_end(&global_context -> input_reads.first_read_file.bcl_input);
	else
		rewind_read_files(global_context, SEEK_END);
	global_context -> running_processed_reads_in_chunk=0;
}

int read_chunk_circles(global_context_t *global_context)
{
	int block_no;

//	printf("GINP1 AT %llu\n", ftello(global_context -> input_reads.first_read_file.input_fp));
	
	unsigned int chunk_no = 0;

	global_context -> current_index = (gehash_t*) malloc(sizeof(gehash_t));
	global_context -> current_value_index = global_context -> all_value_indexes;
	global_context -> running_processed_reads_in_chunk=0;
	global_context -> processed_reads_in_chunk=0;

	double time_load_index = miltime();
	for(block_no = 0; block_no< global_context->index_block_number; block_no++)
	{
		char tmp_fname[MAX_FILE_NAME_LENGTH+ 30];
		sprintf(tmp_fname, "%s.%02d.%c.array", global_context->config.index_prefix, block_no,  global_context->config.space_type == GENE_SPACE_COLOR?'c':'b');
		if(gvindex_load(&global_context -> all_value_indexes[block_no], tmp_fname)) return -1;
	}
	double period_load_index = miltime() - time_load_index;
	global_context -> timecost_load_index += period_load_index;

	while(1)
	{
		int ret;

		locate_read_files(global_context, SEEK_SET);
		for(global_context->current_index_block_number = 0; global_context->current_index_block_number < global_context->index_block_number; global_context->current_index_block_number++)
		{
			char tmp_fname[MAX_FILE_NAME_LENGTH+30];
			time_load_index = miltime();
			if(global_context->index_block_number > 1 || chunk_no == 0)
			{
				sprintf(tmp_fname, "%s.%02d.%c.tab", global_context->config.index_prefix, global_context->current_index_block_number,  global_context->config.space_type == GENE_SPACE_COLOR?'c':'b');
				print_in_box(80,0,0, "Load the %d-th index block...",1+ global_context->current_index_block_number);

				if(gehash_load(global_context -> current_index, tmp_fname)) return -1;
				print_in_box(80,0,0, "The index block has been loaded.");
				sprintf(tmp_fname, "%s.%02d.%c.array", global_context->config.index_prefix, global_context->current_index_block_number, global_context->config.space_type == GENE_SPACE_COLOR?'c':'b');
			}
			period_load_index = miltime() - time_load_index;
			global_context -> timecost_load_index += period_load_index;
			global_context -> current_value_index = global_context -> all_value_indexes + global_context->current_index_block_number;

			if(global_context->current_index_block_number ==0 && global_context -> all_processed_reads==0)
				global_context->align_start_time = miltime();

			if(global_context->index_block_number == global_context->current_index_block_number + 1)
				global_context -> is_final_voting_run = 1;
			else	global_context -> is_final_voting_run = 0;

			print_in_box(80,0,0, "Start read mapping in chunk.");
			double time_start_voting = miltime();
			ret = run_maybe_threads(global_context, STEP_VOTING);
			double period_voting = miltime() - time_start_voting;
			global_context -> timecost_voting += period_voting;

			if(0 == global_context->current_index_block_number){
				locate_read_files(global_context, SEEK_END);
				global_context -> processed_reads_in_chunk = global_context -> running_processed_reads_in_chunk;
			}

			if(global_context->current_index_block_number < global_context->index_block_number -1)
				go_chunk_start(global_context);
			//	rewind_read_files(global_context, SEEK_SET);

			int is_last_chunk = global_context -> processed_reads_in_chunk < global_context->config.reads_per_chunk;

			if(global_context->index_block_number > 1 || is_last_chunk)
				gehash_destory_fast(global_context -> current_index);

			if(ret) break;
			if(!global_context -> processed_reads_in_chunk) break;
		}


		// after the voting step, all subread index blocks are released and all base index blocks are loaded at once.
		double time_before_realign = miltime();

		if(0 == chunk_no && global_context->config.exon_annotation_file[0] && global_context->config.do_breakpoint_detection){
			ret = ret || load_known_junctions(global_context);
			if(!ret){
				// sort and merge events from all threads and the global event space.
				sort_global_event_table(global_context);
				// sort the event entry table at each location.
				sort_junction_entry_table(global_context);
			}
		}

		ret = ret || anti_supporting_read_scan(global_context);
		remove_neighbour(global_context);

		go_chunk_start(global_context);
	//	rewind_read_files(global_context, SEEK_SET);
		double period_before_realign = miltime() - time_before_realign;
		global_context -> timecost_before_realign += period_before_realign;

		double time_realign = miltime();
		ret = ret || run_maybe_threads(global_context, STEP_ITERATION_TWO);
		double period_realign = miltime() - time_realign;

		global_context -> timecost_for_realign += period_realign;

		if(global_context -> config.is_third_iteration_running)
		{
			go_chunk_start(global_context);
			//rewind_read_files(global_context, SEEK_SET);
			ret = ret || do_iteration_three(global_context, NULL);
		}

		go_chunk_nextchunk(global_context);
		//rewind_read_files(global_context, SEEK_END);

		global_context -> all_processed_reads+= global_context ->processed_reads_in_chunk;

		if(ret) return ret;

		if(global_context -> processed_reads_in_chunk < global_context->config.reads_per_chunk ||
		  (global_context -> config.is_BAM_output && global_context -> output_bam_writer -> is_internal_error) || global_context -> input_reads.is_internal_error  || 
		  (global_context -> output_sam_is_full))
			// base value indexes loaded in the last circle are not destroyed and are used in writting the indel VCF.
			// the indexes will be destroyed in destroy_global_context
			break;
		
		clean_context_after_chunk(global_context);
		chunk_no++;
	}

	free(global_context -> current_index);

	if(global_context -> config.is_third_iteration_running)
		finalise_long_insertions(global_context);

	return 0;
}

void char_strftime(char * tbuf){
	time_t rawtime;
	struct tm * timeinfo;

	time (&rawtime);
	timeinfo = localtime (&rawtime);
	strftime (tbuf,80,"%d-%b-%Y %H:%M:%S",timeinfo);

}

void print_subread_logo()
{
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m ========== %c[0m%c[36m    _____ _    _ ____  _____  ______          _____  ", CHAR_ESC, CHAR_ESC, CHAR_ESC);
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m =====      %c[0m%c[36m   / ____| |  | |  _ \\|  __ \\|  ____|   /\\   |  __ \\ ", CHAR_ESC, CHAR_ESC, CHAR_ESC);
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m   =====    %c[0m%c[36m  | (___ | |  | | |_) | |__) | |__     /  \\  | |  | |", CHAR_ESC, CHAR_ESC, CHAR_ESC);
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m     ====   %c[0m%c[36m   \\___ \\| |  | |  _ <|  _  /|  __|   / /\\ \\ | |  | |", CHAR_ESC, CHAR_ESC, CHAR_ESC);
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m       ==== %c[0m%c[36m   ____) | |__| | |_) | | \\ \\| |____ / ____ \\| |__| |", CHAR_ESC, CHAR_ESC, CHAR_ESC);
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m ========== %c[0m%c[36m  |_____/ \\____/|____/|_|  \\_\\______/_/    \\_\\_____/%c[0m", CHAR_ESC, CHAR_ESC, CHAR_ESC, CHAR_ESC);
	#ifdef MAKE_STANDALONE
	char * spaces = "";
	if(strlen(SUBREAD_VERSION) == 8) spaces = "";
	else if(strlen(SUBREAD_VERSION) == 5) spaces = "  ";
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"	%sv%s",spaces,SUBREAD_VERSION);
	#else
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %s",SUBREAD_VERSION);
	#endif
}



int print_configuration(global_context_t * context)
{
	setlocale(LC_NUMERIC, "");
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"");
	print_subread_logo();
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"");
	print_in_box(80, 1, 1, context->config.entry_program_name == CORE_PROGRAM_SUBJUNC?"setting":"setting");
	print_in_box(80, 0, 1, "");

	if(context->config.do_breakpoint_detection)
	{
	    if(context-> config.do_fusion_detection)
			print_in_box(80, 0, 0, "Function      : Read alignment + Junction/Fusion detection%s", context->config.experiment_type == CORE_EXPERIMENT_DNASEQ?" (DNA-Seq)":" (RNA-Seq)");
	    else if(context-> config.do_long_del_detection)
			print_in_box(80, 0, 0, "Function      : Read alignment + Long Deletion detection%s", context->config.experiment_type == CORE_EXPERIMENT_DNASEQ?" (DNA-Seq)":" (RNA-Seq)");
		else
			print_in_box(80, 0, 0, "Function      : Read alignment + Junction detection (%s)", context->config.experiment_type == CORE_EXPERIMENT_DNASEQ?"DNA-Seq":"RNA-Seq");
	}
	else
	        print_in_box(80, 0, 0, "Function      : Read alignment%s", context->config.experiment_type == CORE_EXPERIMENT_DNASEQ?" (DNA-Seq)":" (RNA-Seq)");
	if( context->config.second_read_file[0]){
	        print_in_box(80, 0, 0, "Input file 1  : %s", get_short_fname(context->config.first_read_file));
	        print_in_box(80, 0, 0, "Input file 2  : %s", get_short_fname(context->config.second_read_file));
	}else if(context->config.scRNA_input_mode == GENE_INPUT_SCRNA_FASTQ){
			const char *strmm_tmp = context->config.first_read_file;
			int sample_no = 1;
			while((strmm_tmp = strstr(strmm_tmp, SCRNA_FASTA_SPLIT1))!=NULL) {
				sample_no++;
				strmm_tmp++;
			}
	        print_in_box(80, 0, 0, "Input file    : %d samples from scRNA-seq", sample_no);
	}else if(context->config.scRNA_input_mode == GENE_INPUT_BCL){
	        print_in_box(80, 0, 0, "Input file    : %s%s", get_short_fname(context->config.first_read_file), " (scRNA)");
	}else if(context->config.scRNA_input_mode == GENE_INPUT_SCRNA_BAM){
	        print_in_box(80, 0, 0, "Input file    : %s%s", get_short_fname(context->config.first_read_file), " (10X BAM)");
	}else{
	        print_in_box(80, 0, 0, "Input file    : %s%s", get_short_fname(context->config.first_read_file), context->config.is_SAM_file_input?(context->config.is_BAM_input?" (BAM)":" (SAM)"):(""));
	}

	if(context->config.output_prefix [0])
	        print_in_box(80, 0, 0, "Output file   : %s (%s)%s", get_short_fname(context->config.output_prefix), context->config.is_BAM_output?"BAM":"SAM", context->config.is_input_read_order_required?", Keep Order":(context->config.sort_reads_by_coordinates?", Sorted":""));
	else
	        print_in_box(80, 0, 0, "Output method : STDOUT (%s)" , context->config.is_BAM_output?"BAM":"SAM");

	print_in_box(80, 0, 0,         "Index name    : %s", get_short_fname(context->config.index_prefix));
	print_in_box(80, 0, 0, "");
	print_in_box(80, 0, 1, "------------------------------------");
	print_in_box(80, 0, 0, "");
	print_in_box(80, 0, 0, "                              Threads : %d", context->config.all_threads);
	print_in_box(80, 0, 0, "                         Phred offset : %d", (context->config.phred_score_format == FASTQ_PHRED33)?33:64);
	if( context->config.second_read_file[0])
	{
	print_in_box(80, 0, 0, "              # of extracted subreads : %d", context->config.total_subreads);
	print_in_box(80, 0, 0, "                       Min read1 vote : %d", context->config.minimum_subread_for_first_read);
	print_in_box(80, 0, 0, "                       Min read2 vote : %d", context->config.minimum_subread_for_second_read);
	print_in_box(80, 0, 0, "                    Max fragment size : %d", context->config.maximum_pair_distance);
	print_in_box(80, 0, 0, "                    Min fragment size : %d", context->config.minimum_pair_distance);
	}
	else
	print_in_box(80, 0, 0, "                            Min votes : %d / %d", context->config.minimum_subread_for_first_read, context->config.total_subreads);

	print_in_box(80, 0, 0, "                       Max mismatches : %d", context->config.max_mismatch_exonic_reads);
	print_in_box(80, 0, 0, "                     Max indel length : %d", context->config.max_indel_length);
	print_in_box(80, 0, 0, "           Report multi-mapping reads : %s", context->config.report_multi_mapping_reads?"yes":"no");
	print_in_box(80, 0, 0, "Max alignments per multi-mapping read : %d", context->config.multi_best_reads);

	if(context->config.exon_annotation_file[0]){
	if(context->config.exon_annotation_file_screen_out[0])
	print_in_box(80, 0, 0, "                          Annotations : %s", context->config.exon_annotation_file_screen_out);
	else
	print_in_box(80, 0, 0, "                          Annotations : %s (%s)",get_short_fname( context->config.exon_annotation_file), context->config.exon_annotation_file_type==FILE_TYPE_GTF?"GTF":"SAF");
	}

	if(context->config.max_insertion_at_junctions)
	print_in_box(80, 0, 0, "                   Insertions at junc : %d", context->config.max_insertion_at_junctions);

	if(context->config.read_group_id[0])
	print_in_box(80, 0, 0, "                      Read group name : %s", context->config.read_group_id);

	print_in_box(80, 0, 1, "");
	print_in_box(80, 2, 1, "");
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"");


	if(!context->config.experiment_type){
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"You have to specify the experiment type by using the '-t' option.\n");
		return -1;
	}

	if(!context->config.first_read_file[0])
	{
	        sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"You have to specify at least one input file in the FASTQ/FASTA/PLAIN format using the '-r' option.\n");
	        return -1;
	}

	if(0 && !context->config.output_prefix[0])
	{
	        sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"You have to specify the path of output using the '-o' option.\n");
	        return -1;
	}

	if(!context->config.index_prefix[0])
	{
	        sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"You have to specify the prefix of the index files using the '-i' option.\n");
	        return -1;
	}
	char tbuf[90];
	char_strftime(tbuf);

	print_in_box(80,1,1,"Running (%s, pid=%d)", tbuf, getpid());
	print_in_box(80,0,1,"");


	return 0;
}



int init_paired_votes(global_context_t *context)
{
	init_bigtable_results(context, 0);
	return 0;
}

#define FETCH_SEQ_LEN(x1) (context->chromosome_table.read_offsets[xk1] - last_offset+16 - ( 2 * context->chromosome_table.padding ))

void write_sam_headers(global_context_t * context)
{
	char *  sorting_str = context -> config.sort_reads_by_coordinates?"coordinate":"unsorted";
	if(context -> config.is_BAM_output)
	{
		char header_buff[100];
		sprintf(header_buff, "@HD\tVN:1.0\tSO:%s", sorting_str);
		SamBam_writer_add_header(context -> output_bam_writer, header_buff, 0);
		int xk1;
		int last_offset = 0;
		char * obuf = malloc(10000+5000);
		for(xk1=0; xk1< context->chromosome_table.total_offsets; xk1++)
		{
			int seq_len = FETCH_SEQ_LEN(x1);
			//context->chromosome_table.read_offsets[xk1] - last_offset+16
			SamBam_writer_add_chromosome(context -> output_bam_writer, context->chromosome_table.read_names+ xk1 * MAX_CHROMOSOME_NAME_LEN, seq_len, 1);
			last_offset = context->chromosome_table.read_offsets[xk1];
		}


		if(context->config.read_group_id[0])
		{
			snprintf(obuf,10000, "@RG\tID:%s%s",context->config.read_group_id, context->config.read_group_txt);
			SamBam_writer_add_header(context -> output_bam_writer,obuf, 0);
		}
		snprintf(obuf,9999+4900, "@PG\tID:subread\tPN:subread\tVN:%s\tCL:%s", SUBREAD_VERSION, context->rebuilt_command_line);
		SamBam_writer_add_header(context -> output_bam_writer,obuf, 0);
		SamBam_writer_finish_header(context -> output_bam_writer);
		free(obuf);
	}
	else
	{
		sambamout_fprintf(context -> output_sam_fp, "@HD\tVN:1.0\tSO:%s\n", sorting_str);
		int xk1;
		int last_offset = 0;
		for(xk1=0; xk1< context->chromosome_table.total_offsets; xk1++)
		{
			int seq_len = FETCH_SEQ_LEN(x1);
			sambamout_fprintf(context -> output_sam_fp, "@SQ\tSN:%s\tLN:%u\n", context->chromosome_table.read_names+ xk1 * MAX_CHROMOSOME_NAME_LEN, seq_len);
			last_offset = context->chromosome_table.read_offsets[xk1];
		}

		if(context->config.read_group_id[0])
			sambamout_fprintf(context -> output_sam_fp, "@RG\tID:%s%s\n",context->config.read_group_id, context->config.read_group_txt);
		sambamout_fprintf(context -> output_sam_fp, "@PG\tID:subread\tPN:subread\tVN:%s\tCL:%s\n", SUBREAD_VERSION, context->rebuilt_command_line);
		
	}
}

#define EXONIC_REGION_RESOLUTION 16

int is_pos_in_annotated_exon_regions(global_context_t * global_context, unsigned int pos){ 
	int exonic_map_byte = pos / EXONIC_REGION_RESOLUTION / 8;
	int exonic_map_bit = (pos / EXONIC_REGION_RESOLUTION) % 8;

	return (global_context ->exonic_region_bitmap [exonic_map_byte] & (1<<exonic_map_bit))?1:0;

}

char * get_sam_chro_name_from_alias(HashTable * tab, char * anno_chro){
	KeyValuePair * cursor = NULL;
	long x1;
	for(x1 = 0; x1 < tab -> numOfBuckets; x1 ++){
		cursor = tab -> bucketArray[x1];
		while(1){
			if(NULL == cursor)break;
			char * tab_anno_chro = cursor -> value;
			if(strcmp(tab_anno_chro, anno_chro) == 0) return (char *)cursor -> key;
			cursor = cursor -> next;
		}
	}
	return NULL;
}

int do_anno_bitmap_add_feature(char * gene_name, char * tracnscript_id, char * chro_name, unsigned int feature_start, unsigned int feature_end, int is_negative_strand, void * context){
	global_context_t * global_context = context;
	if(global_context -> sam_chro_to_anno_chr_alias){
		char * sam_chro = get_sam_chro_name_from_alias(global_context -> sam_chro_to_anno_chr_alias, chro_name);
		if(sam_chro!=NULL) chro_name = sam_chro;
	}

	if(!HashTableGet(global_context -> annotation_chro_table,chro_name))
		HashTablePut(global_context -> annotation_chro_table, memstrcpy( chro_name ), NULL+1);

	char tmp_chro_name[MAX_CHROMOSOME_NAME_LEN];
	int access_n = HashTableGet( global_context -> chromosome_table.read_name_to_index, chro_name ) - NULL;
	if(access_n < 1){
		if(chro_name[0]=='c' && chro_name[1]=='h' && chro_name[2]=='r'){
			chro_name += 3;
		}else{
			strcpy(tmp_chro_name, "chr");
			strcat(tmp_chro_name, chro_name);
			chro_name = tmp_chro_name;
		}
	}

	unsigned int exonic_map_start = linear_gene_position(&global_context->chromosome_table , chro_name, feature_start);
	unsigned int exonic_map_stop = linear_gene_position(&global_context->chromosome_table , chro_name, feature_end);
	if(exonic_map_start > 0xffffff00 || exonic_map_stop > 0xffffff00){
		return -1;
	}

	exonic_map_start -= exonic_map_start%EXONIC_REGION_RESOLUTION;
	exonic_map_stop -= exonic_map_stop%EXONIC_REGION_RESOLUTION;

	for(; exonic_map_start <= exonic_map_stop; exonic_map_start+=EXONIC_REGION_RESOLUTION){
		int exonic_map_byte = exonic_map_start / EXONIC_REGION_RESOLUTION / 8;
		int exonic_map_bit = (exonic_map_start / EXONIC_REGION_RESOLUTION) % 8;

		global_context ->exonic_region_bitmap[exonic_map_byte] |= (1<<exonic_map_bit);
	}
	return 0;
}

void convert_table_key_to_alias(void * ky ,void * va, HashTable * tab){
	HashTable * alias_index_to_annot = tab -> appendix2;
	HashTable * result_annot_chros = tab -> appendix1;

	char * index_name = ky;
	char * annot_name = HashTableGet(alias_index_to_annot, index_name);
	if(annot_name==NULL) annot_name = index_name;
	HashTablePut( result_annot_chros, annot_name, va ); 
}

void warning_anno_vs_index(HashTable * anno_chros_tab, gene_offset_t * index_chros_offset, HashTable * alias_index_to_annot){

	HashTable * index_chros_tab = index_chros_offset -> read_name_to_index;
	HashTable * chros_in_annot_using_index_names = anno_chros_tab ;
	HashTable * chros_in_index_using_annot_names = index_chros_tab ;
	HashTable * alias_annot_to_index = NULL;

	if(alias_index_to_annot){
		chros_in_annot_using_index_names = StringTableCreate(1000);
		chros_in_index_using_annot_names = StringTableCreate(1000);
	
		alias_annot_to_index = StringTableReverse(alias_index_to_annot);
		anno_chros_tab -> appendix1 = chros_in_annot_using_index_names;
		anno_chros_tab -> appendix2 = alias_annot_to_index;
		HashTableIteration(anno_chros_tab, convert_table_key_to_alias);
	
		index_chros_tab -> appendix1 = chros_in_index_using_annot_names;
		index_chros_tab -> appendix2 = alias_index_to_annot;
		HashTableIteration(index_chros_tab, convert_table_key_to_alias);
	}

	warning_hash_hash(anno_chros_tab, chros_in_index_using_annot_names, "Chromosomes/contigs in annotation but not in index :");
	warning_hash_hash(index_chros_tab, chros_in_annot_using_index_names, "Chromosomes/contigs in index but not in annotation :");

	if(alias_index_to_annot){
		HashTableDestroy(alias_annot_to_index);
		HashTableDestroy(chros_in_annot_using_index_names);
		HashTableDestroy(chros_in_index_using_annot_names);
	}
}

int load_annotated_exon_regions(global_context_t * global_context){
	int bitmap_size = (4096 / EXONIC_REGION_RESOLUTION / 8)*1024*1024;
	global_context ->exonic_region_bitmap = malloc(bitmap_size);
	memset( global_context ->exonic_region_bitmap , 0, bitmap_size );
	global_context -> annotation_chro_table = HashTableCreate(1003);
	HashTableSetDeallocationFunctions( global_context -> annotation_chro_table, free, NULL);
	HashTableSetKeyComparisonFunction( global_context -> annotation_chro_table, my_strcmp);
	HashTableSetHashFunction( global_context -> annotation_chro_table, fc_chro_hash);

	int loaded_features = load_features_annotation(global_context->config.exon_annotation_file, global_context->config.exon_annotation_file_type, global_context->config.exon_annotation_gene_id_column, NULL, global_context->config.exon_annotation_feature_name_column, global_context, do_anno_bitmap_add_feature);
	if(loaded_features < 0)return -1;
	else print_in_box(80,0,0,"%d annotation records were loaded.\n", loaded_features);
	return 0;
}

int load_global_context(global_context_t * context)
{
	char tmp_fname [MAX_FILE_NAME_LENGTH + 50];
	int min_phred_score = -1 , max_phred_score = -1;
	context -> is_phred_warning = 0; 
	
	if(context->config.multi_best_reads>1 && ! context->config.report_multi_mapping_reads){
		print_in_box(80,0,0,"WARNING: Multi-mapping reads are reported.");
		context->config.report_multi_mapping_reads = 1;
	}

	if(context->config.scRNA_input_mode){
		context -> config.multi_best_reads = 3;
		context -> config.do_remove_neighbour_for_scRNA = 1;
		context -> config.reads_per_chunk = 30000000*context -> config.multi_best_reads;
		context -> config.multi_best_reads = max(context -> config.multi_best_reads , context -> config.reported_multi_best_reads);
		// opening a BCL input needs the exact chunk size. 
		if(context->config.multi_best_reads>1) context -> config.reads_per_chunk /= context->config.multi_best_reads;
		char * reads_per_chunk = getenv("CC_READS_PER_CHUNK");
		char * remove_beighbour = getenv("CC_REMOVE_NEIGHBOUR");
		if(remove_beighbour) context -> config.do_remove_neighbour_for_scRNA = remove_beighbour[0]-'0';
		if(reads_per_chunk) context -> config.reads_per_chunk = atoi(reads_per_chunk);
	}
	print_in_box(80,0,0,"Check the input reads.");
	subread_init_lock(&context->input_reads.input_lock);
	if(core_geinput_open(context, &context->input_reads.first_read_file, 1)) {
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Unable to open the input file '%s'\n", context->config.first_read_file);
		return -1;
	}

	context->config.space_type = context->input_reads.first_read_file.space_type;
	print_in_box(89,0,0,"The input file contains %c[36m%s%c[0m space reads.", CHAR_ESC, context->config.space_type == GENE_SPACE_COLOR?"color":"base", CHAR_ESC);
	if(context->config.space_type == GENE_SPACE_COLOR && context->config.is_BAM_output && !context->config.convert_color_to_base)
	{
		print_in_box(80,0,0,"The color-space bases will be converted into base space in the BAM output.");
		context->config.convert_color_to_base=1;
	}
	else if(context->config.space_type == GENE_SPACE_BASE && context->config.convert_color_to_base)
	{
		print_in_box(80,0,0,"The reads will not be converted into base space.");
		context->config.convert_color_to_base=0;
	}



	if(context->input_reads.is_paired_end_reads)
	{
		if(core_geinput_open(context, &context->input_reads.second_read_file, 2))
		{
			//sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Unable to open '%s' as input. Please check if it exists, you have the permission to read it, and it is in the correct format.\n", context->config.second_read_file);
			return -1;
		}

		context -> config.max_vote_combinations = 3;
		context -> config.multi_best_reads = 3;
		context -> config.max_vote_simples = 64;
		context -> config.max_vote_number_cutoff = 2;
	}else{
		context -> config.max_vote_combinations = 3;
		context -> config.multi_best_reads = 3;
		context -> config.max_vote_simples = 3;
		context -> config.max_vote_number_cutoff = 2;
	}

	print_in_box(80,0,0,"Initialise the memory objects.");
	context -> config.max_vote_simples = max(context -> config.max_vote_simples ,  context -> config.reported_multi_best_reads);
	context -> config.max_vote_combinations = max(context -> config.max_vote_combinations ,  context -> config.reported_multi_best_reads);

	if(!context->config.scRNA_input_mode){
		// if it is BCL input then these two parameters have been decided .
		if(context->config.multi_best_reads>1) context->config.reads_per_chunk /= context->config.multi_best_reads;
		if(context->input_reads.is_paired_end_reads) context->config.reads_per_chunk /= 2;
	}

	struct stat ginp1_stat;
	int guess_tested_reads = 0;
	stat(context->config.first_read_file , &ginp1_stat);
	context->input_reads.first_read_file_size = ginp1_stat.st_size;

	print_in_box(80,0,0,"Estimate the mean read length.");
	context -> input_reads.avg_read_length = 100;
	if(context->config.scRNA_input_mode != GENE_INPUT_SCRNA_BAM && context->config.scRNA_input_mode != GENE_INPUT_SCRNA_FASTQ) context -> input_reads.avg_read_length = guess_reads_density_format(context->config.first_read_file , context->config.is_SAM_file_input?1:0, &min_phred_score, &max_phred_score , &guess_tested_reads);
	if(context -> input_reads.avg_read_length<0 )context -> input_reads.avg_read_length = 250;
//	SUBREADprintf("QR=[%d,%d]; ALEN=%f\n",  min_phred_score, max_phred_score, context -> input_reads.avg_read_length);
	if(max_phred_score>=0)
	{
		int inferred_offset;

		if(abs(min_phred_score - 33) < abs(min_phred_score - 64))  inferred_offset = 33;
		else	inferred_offset = 64;

		if(context -> input_reads.first_read_file.file_type != GENE_INPUT_GZIP_FASTA){
			if((context->config.phred_score_format == FASTQ_PHRED64 && inferred_offset == 33) ||
			   (context->config.phred_score_format == FASTQ_PHRED33 && inferred_offset == 64))
			{
				print_in_box(80,0,0, "WARNING  - The specified Phred score offset (%d) seems incorrect.", context->config.phred_score_format == FASTQ_PHRED33?33:64);
				print_in_box(80,0,0, "           ASCII values of the quality scores of read bases included in");
				print_in_box(80,0,0, "           the first %d reads were found to be within the range of",guess_tested_reads);
				print_in_box(80,0,0, "           [%d,%d].", min_phred_score, max_phred_score);
				print_in_box(80,0,0, "");
				context -> is_phred_warning = 1;
			}
			else{
				print_in_box(80,0,0, "The range of Phred scores observed in the data is [%d,%d]", min_phred_score - inferred_offset, max_phred_score - inferred_offset);
			}
		}
	}

	subread_init_lock(&context -> output_lock);

	if(context->config.report_sam_file && context -> config.output_prefix[0])
	{
		// ====== open output files ======
		// Only the sam file is opened here; other files like bed, indel and etc are opened in init_modules()
		sprintf(tmp_fname,"%s", context->config.output_prefix);

		print_in_box(80,0,0,"Create the output %s file.", context -> config.is_BAM_output?"BAM":"SAM");
		if(context -> config.is_BAM_output)
		{
			context -> output_bam_writer = malloc(sizeof(SamBam_Writer));
			SamBam_writer_create(context -> output_bam_writer , tmp_fname, context -> config.is_input_read_order_required?1:context -> config.all_threads,  context -> config.sort_reads_by_coordinates, context->config.scRNA_input_mode, context -> config.temp_file_prefix);
			context -> output_sam_fp = NULL;
		}
		else
		{
			context -> output_sam_fp = f_subr_open(tmp_fname,"wb");
			//context -> output_sam_inner_buffer = malloc(OUTPUT_BUFFER_SIZE);
			//setvbuf (context -> output_sam_fp, context -> output_sam_inner_buffer, _IOFBF, OUTPUT_BUFFER_SIZE);
			context -> output_bam_writer = NULL;
		}
		if((!context -> output_bam_writer) && (!context->output_sam_fp))
		{
			sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Unable to open '%s' to write. Please check if the path exists and you have the permission to create/write this file", tmp_fname);
			return -1;
		}
	}
	else
	{
		if(context -> config.is_BAM_output)
		{
			context -> output_bam_writer = malloc(sizeof(SamBam_Writer));
			SamBam_writer_create(context -> output_bam_writer ,NULL, context -> config.is_input_read_order_required?1:context -> config.all_threads,  context -> config.sort_reads_by_coordinates, 0,  context -> config.temp_file_prefix);
		}
		context->output_sam_fp = NULL;
	}
	
	// ====== check index files, count blocks and load chro table ======
	print_in_box(80,0,0,"Check the index.");
	sprintf(tmp_fname, "%s.reads", context->config.index_prefix);
	if(!does_file_exist(tmp_fname))
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Unable top open index '%s'. Please make sure that the correct prefix is specified and you have the permission to read these files. For example, if there are files '/opt/my_index.reads', '/opt/my_index.files' and etc, the index prefix should be specified as '/opt/my_index' without any suffix. \n", context->config.index_prefix);
		return -1;
	}
	
	context->current_index_block_number = 0;
	if(load_offsets(&context->chromosome_table, context->config.index_prefix)){
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"\nThe index was built by using an old version of Subread; its format is no longer supported. Please use the current version of the index builder to rebuild it.\n");
		return 1;
	}

	if(context->config.space_type == GENE_SPACE_COLOR)
		sprintf(tmp_fname, "%s.00.c.tab", context->config.index_prefix);
	else
		sprintf(tmp_fname, "%s.00.b.tab", context->config.index_prefix);
	if(!does_file_exist(tmp_fname))
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Your reads are in the %s space but the index was not built in the same space. Unable to precess the reads.\n", context->config.space_type == GENE_SPACE_COLOR?"color":"base");
		return -1;
	}

	context->index_block_number = 0; 
	while(1)
	{
		sprintf(tmp_fname, "%s.%02d.%c.tab", context->config.index_prefix, context->index_block_number, context->config.space_type == GENE_SPACE_COLOR?'c':'b');
		if(!does_file_exist(tmp_fname))break;
		context->index_block_number ++;
		if(context->index_block_number>=2 && context->config.max_indel_length > 16)
		{
			print_in_box(80,0,0,"ERROR You cannot use multi-block index for very-long indel detection.");
			print_in_box(80,0,0,"Please set the maximum indel length <= 16.");
			return -1;
		}
	}

	if(context->config.report_sam_file)
		write_sam_headers(context);

	print_in_box(80,0,0,"Init the voting space.");
	// ====== init other variables ======
	if(context -> config.all_threads>1)
		subread_init_lock(&context -> thread_initial_lock);
	else subread_init_topKbuff(context, &context -> topKbuff);

	if(init_paired_votes(context))
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Cannot initialise the voting space. You need at least 2GB of empty physical memory to run this program.\n");
		return 1;
	}
	sublog_printf(SUBLOG_STAGE_DEV1, SUBLOG_LEVEL_DEBUG, "load_global_context: finished");

	memset( context->all_value_indexes , 0 , 100 * sizeof(gene_value_index_t));

	context -> sam_chro_to_anno_chr_alias = NULL;
	if(context->config.exon_annotation_file[0]){ 
		print_in_box(80,0,0,"Load the annotation file.");
		if( load_annotated_exon_regions( context ) ) return -1;
		if(context->config.exon_annotation_alias_file[0])
			 context -> sam_chro_to_anno_chr_alias = load_alias_table(context->config.exon_annotation_alias_file);
		warning_anno_vs_index(context -> annotation_chro_table, &context -> chromosome_table, context -> sam_chro_to_anno_chr_alias);
		HashTableDestroy(context -> annotation_chro_table);
	} else context -> exonic_region_bitmap = NULL;
	print_in_box(80,0,0,"Global environment is initialised.");
	return 0;
}


int init_modules(global_context_t * context)
{
	sublog_printf(SUBLOG_STAGE_DEV1, SUBLOG_LEVEL_DEBUG, "init_modules: begin");
	int ret = init_indel_tables(context);
	if((context->config.do_breakpoint_detection || context-> config.do_fusion_detection || context->config.do_breakpoint_detection || context-> config.do_long_del_detection))
		ret = ret || init_junction_tables(context);

	sublog_printf(SUBLOG_STAGE_DEV1, SUBLOG_LEVEL_DEBUG, "init_modules: finished: %d",ret);
	return ret;
}

int destroy_modules(global_context_t * context)
{
	destroy_indel_module(context);
	if((context->config.do_breakpoint_detection || context-> config.do_fusion_detection || context->config.do_breakpoint_detection || context-> config.do_long_del_detection))
		destroy_junction_tables(context);
	return 0;
}

int destroy_global_context(global_context_t * context)
{
	int xk1, block_no, ret = 0;

	if(context -> exonic_region_bitmap) free(context -> exonic_region_bitmap);

	for(block_no = 0; block_no< context->index_block_number; block_no++)
		gvindex_destory(&context -> all_value_indexes[block_no]);

	if(context->config.all_threads<2){
		subread_free_topKbuff( context,&context -> topKbuff );
	}
	if(context->output_sam_fp) {
		if(context -> output_sam_is_full){
			unlink(context->config.output_prefix);
			SUBREADprintf("\nERROR: cannot finish the SAM file. Please check the disk space in the output directory.\nNo output file was generated.\n");
			ret = 1;
		}
		fclose(context -> output_sam_fp);
	}

	if(context->input_reads.is_internal_error){
		unlink(context->config.output_prefix);
		return 1;
	}

	if(context->output_bam_writer) {
		SamBam_writer_close(context->output_bam_writer);
		if(context->output_bam_writer -> is_internal_error){
			unlink(context->config.output_prefix);
			SUBREADprintf("\nERROR: cannot finish the BAM file. Please check the disk space in the output directory.\nNo output file was generated.\n");
			ret = 1;
		}
		free(context->output_bam_writer);
		context->output_bam_writer=NULL;
	}
	//free(context->big_margin_record);
	
	for(xk1=0; xk1<5; xk1++)
		if(context->module_contexts[xk1])free(context->module_contexts[xk1]);
	geinput_close(&context -> input_reads.first_read_file);
	if(context->input_reads.is_paired_end_reads) geinput_close(&context -> input_reads.second_read_file);
	destroy_offsets(&context->chromosome_table);
	finalise_bigtable_results(context);

	if((context -> will_remove_input_file & 1) && strstr(context ->config.first_read_file, "/core-temp")) unlink(context ->config.first_read_file);
	if((context -> will_remove_input_file & 2) && strstr(context ->config.second_read_file, "/core-temp")) unlink(context ->config.second_read_file);
	free(context -> rebuilt_command_line);
	
	return ret;
}

void subread_rebuild_cmd(int argc, char ** argv, global_context_t * global_context){
	global_context -> rebuilt_command_line_size = rebuild_command_line(&global_context -> rebuilt_command_line , argc, argv);
}


int write_bincigar_part(char * bincigar, int chropt, unsigned int optlen, int bincigar_len)
{
	int binopt, binbytes, x1;

	if(optlen<1) return -1;

	if(optlen < 4) binbytes=1;
	else if(optlen < 1024) binbytes=2;
	else if(optlen < 262144) binbytes=3;
	else if(optlen < 67108864) binbytes=4;
	else binbytes=5;

	if(bincigar_len<binbytes) return -1; 

	switch(chropt)
	{
		case 'S':
			binopt = CORE_CIGAR_OPT_S;
			break;
		case 'M':
			binopt = CORE_CIGAR_OPT_M;
			break;
		case 'I':
			binopt = CORE_CIGAR_OPT_I;
			break;
		case 'D':
			binopt = CORE_CIGAR_OPT_D;
			break;
		case 'N':
			binopt = CORE_CIGAR_OPT_N;
			break;
		case 'n':
			binopt = CORE_CIGAR_OPT_NN;
			break;
		case 'B':
			binopt = CORE_CIGAR_OPT_B;
			break;
		case 'b':
			binopt = CORE_CIGAR_OPT_BB;
			break;
		default:
			return -1;
	}

	bincigar[0]=binopt | (binbytes << 3) | ((optlen & 3)<< 6);
	optlen >>= 2;
	for(x1=1;x1<binbytes; x1++)
	{
		bincigar[x1] = optlen&0xff;
		optlen>>=8;
	}

	return binbytes;
}

// function returns the actual length of bincigar, or -1 if anything is wrong, e.g., bincigar_len is too short or unrecognized operations.
int OLD_cigar2bincigar(char *cigar, char *bincigar, int bincigar_len)
{
	int xk1=0;
	unsigned int tmpv=0, bincigar_cursor=0;
	while(1)
	{
		int nch = cigar[xk1];
		if(!nch) break;
		xk1++;

		if(isdigit(nch)) tmpv=tmpv*10+(nch-'0');
		else
		{
			int bincigar_sec_len = write_bincigar_part(bincigar+bincigar_cursor, nch, tmpv, bincigar_len-bincigar_cursor);
			if(bincigar_sec_len<0){
				bincigar[0]=0;
				return -1;
			}
			bincigar_cursor += bincigar_sec_len;
			tmpv=0;
		}
	}

	if(bincigar_cursor<bincigar_len) bincigar[bincigar_cursor] = 0;

	//printf("%s : BL=%d\n", cigar, bincigar_cursor);

	return bincigar_cursor;
}


int write_cigar_part(char *bincigar, char *cigar, int cigar_len , int * bincigar_move)
{
	int binbytes, x1, binopt, charopt;
	unsigned int tmpv = 0;
	char sec_buf[13];

	binbytes = 7& (bincigar[0] >> 3);
	binopt = 7 & bincigar[0];

	switch(binopt)
	{
		case CORE_CIGAR_OPT_D:
			charopt='D'; 
			break;
		case CORE_CIGAR_OPT_I:
			charopt='I'; 
			break;
		case CORE_CIGAR_OPT_M:
			charopt='M'; 
			break;
		case CORE_CIGAR_OPT_S:
			charopt='S'; 
			break;
		case CORE_CIGAR_OPT_B:
			charopt='B'; 
			break;
		case CORE_CIGAR_OPT_BB:
			charopt='b'; 
			break;
		case CORE_CIGAR_OPT_N:
			charopt='N'; 
			break;
		default:
			charopt='n'; 
			break;
	}

	tmpv = (bincigar[0]>>6) & 3;
	for(x1 = 1; x1 < binbytes; x1++)
	{
		unsigned int dtmpv = 0xff & bincigar[x1];
		dtmpv <<= (x1*8 - 6);
		tmpv += dtmpv; 
	}

	int added_len = sprintf(sec_buf, "%u%c", tmpv, charopt);
	if(added_len > cigar_len)
		return -1;
	memcpy(cigar, sec_buf, added_len);
	(*bincigar_move) = binbytes;

	return added_len;
}

int OLD_bincigar2cigar(char * cigar, int cigar_len, char * bincigar, int bincigar_max_len, int read_len)
{

	int cigar_cursor = 0, bincigar_cursor = 0;

	while(1)
	{
		int bincigar_move = 0;
		int cigar_sec_len = write_cigar_part(bincigar + bincigar_cursor, cigar+cigar_cursor, cigar_len-cigar_cursor-1, &bincigar_move);
		if(cigar_sec_len<0){
			sprintf(cigar,"%dM", read_len);
			return -1;
		}
		//printf("NPC=%s\n", cigar);
		cigar_cursor += cigar_sec_len;
		bincigar_cursor += bincigar_move;
		if(bincigar_cursor>=bincigar_max_len) break;
		if(bincigar[bincigar_cursor] == 0) break;
	}
	cigar[cigar_cursor] = 0;

	return cigar_cursor;
}

int term_strncpy(char * dst, char * src, int max_dst_mem)
{
	int i;

	for(i=0; i<max_dst_mem; i++)
	{
		if(!src[i]) break;
		dst[i]=src[i];
		if(i == max_dst_mem-1)
			SUBREADprintf("String out of memory limit: '%s'\n", src);
	}
	if(i >= max_dst_mem) i = max_dst_mem-1;
	dst[i] = 0;

	return 0;
}
void absoffset_to_posstr(global_context_t * global_context, unsigned int pos, char * res){
	char * ch;
	int off;

	locate_gene_position(pos, &global_context -> chromosome_table, &  ch, &off);

	sprintf(res, "%s:%u", ch, off);
}

int chimeric_cigar_parts(global_context_t * global_context, unsigned int sel_pos, int is_first_section_negative_strand, int is_first_section_reversed, char * in_cigar, unsigned int * out_poses, char ** out_cigars, char * out_strands, int read_len, short * perfect_lens, char * read_name)
{
	unsigned int current_perfect_map_start = sel_pos;
	int current_perfect_section_no = 0;
	unsigned int current_perfect_cursor = sel_pos;
	int is_reversed = is_first_section_reversed;
	int is_negative = is_first_section_negative_strand;
	int read_cursor = 0;
	int out_cigar_writer_ptr = 0;
	unsigned int tmpi = 0;

	short perfect_len = 0;

	int cigar_cursor;

	out_poses[0] = current_perfect_map_start - (is_reversed?1:0);
	out_strands[0] = is_negative?'-':'+';
	char main_piece_strand = (is_first_section_negative_strand == is_first_section_reversed)?'+':'-';

	for(cigar_cursor=0;;cigar_cursor++)
	{
		char ncch = in_cigar[cigar_cursor];
		int is_chimeric_section_end = 0;

		if(!ncch){
			perfect_lens [current_perfect_section_no] = perfect_len ;
			current_perfect_section_no++;
			break;
		}

		if(toupper(ncch)=='N'||toupper(ncch)=='B')
		{

			unsigned int jummped_location;
			int is_chro_jump = 0, is_long_jump = 0;

			if(is_reversed)
			{
				if(toupper(ncch)=='N')
					jummped_location = current_perfect_map_start - 1 + tmpi;
				else
					jummped_location = current_perfect_map_start - 1 - tmpi;
			}
			else
			{
				if(toupper(ncch)=='N')
					jummped_location = current_perfect_cursor + tmpi;
				else
					jummped_location = current_perfect_cursor - tmpi;

			}

			if(ncch == 'N')
			{
				char * curr_chr, * new_chr;
				int curr_offset, new_offset;
				locate_gene_position_max(current_perfect_cursor, &global_context -> chromosome_table, & curr_chr, & curr_offset, NULL, NULL, 1);
				locate_gene_position_max(jummped_location      , &global_context -> chromosome_table, &  new_chr, &  new_offset, NULL, NULL, 1);
				if( curr_chr == NULL || new_chr == NULL ){
					/*
					char outpos[100];
					absoffset_to_posstr(global_context, sel_pos + 1, outpos);
					SUBREADprintf("Wrong CIGAR: mapped to %s, CIGAR=%s\n", outpos , in_cigar);
					*/
					return -1;
				}
				assert(curr_chr);
				assert(new_chr);
				is_chro_jump = (curr_chr != new_chr);

				long long int dist = (long long int)current_perfect_cursor;
				dist -= (long long int)jummped_location;
				if(abs(dist) >= global_context -> config.maximum_intron_length)
					is_long_jump = 1;

			//#warning ">>>>>> COMMENT WHEN RELEASE <<<<<<"
			if(0 && FIXLENstrcmp("R000000007", read_name)==0)
				SUBREADprintf("dist=%lld, abs=%lld, %u, %u\n",dist,abs(dist), current_perfect_cursor, jummped_location);
				// A long jump is the jump longer than 2^27.
				// Picard does not like it!!
			}

			// is_long_jump is true only if the two sections are on different chromosomes.
			//#warning ">>>>>> COMMENT WHEN RELEASE <<<<<<"
			if(0 && FIXLENstrcmp("R000000007", read_name)==0)
				SUBREADprintf("CHR_JMP=%d, NCCH=%c, LONG_JMP=%d\n", is_chro_jump, ncch, is_long_jump);
			if(is_chro_jump || islower(ncch) || ncch == 'B' || is_long_jump)
			{
				current_perfect_cursor = jummped_location;

				if(islower(ncch)){
					is_reversed = !is_reversed;
					is_negative = !is_negative;

					if(is_reversed)
						current_perfect_cursor --;
				}

				current_perfect_map_start = current_perfect_cursor;
				tmpi = 0;
				if(read_cursor<read_len)
					sprintf(out_cigars[current_perfect_section_no] + out_cigar_writer_ptr,"%dS", read_len - read_cursor);

				perfect_lens [current_perfect_section_no] = perfect_len ;
				perfect_len = 0;

				current_perfect_section_no++;
				if(current_perfect_section_no>=CIGAR_PERFECT_SECTIONS)break;

				out_poses[current_perfect_section_no] = current_perfect_map_start - read_cursor;
				out_strands[current_perfect_section_no] = is_negative?'-':'+';
				out_cigar_writer_ptr = sprintf(out_cigars[current_perfect_section_no],"%dS", read_cursor);
				is_chimeric_section_end  = 1;
			}
		}

		if(!is_chimeric_section_end)
		{
			if(isalpha(ncch))
			{
				out_cigar_writer_ptr+=sprintf(out_cigars[current_perfect_section_no]+out_cigar_writer_ptr, "%u%c", tmpi, ncch);
			}
			if(ncch == 'M'|| ncch == 'S')
			{
				read_cursor += tmpi;
				if(ncch == 'M')
					perfect_len += tmpi;
				if(is_reversed)
					out_poses[current_perfect_section_no] += tmpi;
				else
					current_perfect_cursor += tmpi;
				tmpi = 0;
			}
			else if(ncch == 'D' || ncch == 'N')
			{
				if(is_reversed)
					out_poses[current_perfect_section_no] += tmpi;
				else
					current_perfect_cursor += tmpi;
				tmpi = 0;
			}
			else if(ncch == 'I')
			{
				read_cursor += tmpi;
				tmpi = 0;
			}
			else if(isdigit(ncch))
				tmpi = tmpi*10+(ncch-'0');
		}
	}

	int xk1 = 0, best_match = -9999;

	for(xk1=0; xk1<current_perfect_section_no;xk1++)
	{
		if(best_match < 0 || (main_piece_strand == out_strands[xk1] && perfect_lens[xk1]>perfect_lens[best_match]))
			best_match = xk1;
	}

	if(best_match>0)
	{
		unsigned int tmpv;
		char cigar_sec[100];
		tmpv = out_poses[0];
		out_poses[0]=out_poses[best_match];
		out_poses[best_match] = tmpv;

		tmpv = out_strands[0];
		out_strands[0] = out_strands[best_match];
		out_strands[best_match] = tmpv;

		strcpy(cigar_sec, out_cigars[0]);
		strcpy(out_cigars[0], out_cigars[best_match]);
		strcpy(out_cigars[best_match] , cigar_sec);

		tmpv = perfect_lens[0];
		perfect_lens[0] = perfect_lens[best_match];
		perfect_lens[best_match] = tmpv;
	}

	return current_perfect_section_no;
}

void quick_sort_run(void * arr, int spot_low,int spot_high, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r));

void quick_sort(void * arr, int arr_size, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r))
{
	quick_sort_run(arr, 0, arr_size-1, compare, exchange);
}
 
 
void quick_sort_run(void * arr, int spot_low,int spot_high, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r))
{
	// https://en.wikipedia.org/wiki/Quicksort
	// Lomuto partition scheme

	int pivot,j,i;

	if(spot_high <= spot_low) return;
	pivot = spot_high;
	i = spot_low;

	for(j = spot_low+1; j < spot_high; j++)
		if(compare(arr, j, pivot)<=0) {
			exchange(arr,i,j);
			i++;
		}

	exchange(arr, i, spot_high);

	quick_sort_run(arr, spot_low, i-1, compare, exchange);
	quick_sort_run(arr, i+1, spot_high, compare, exchange);
	
}
void basic_sort_run(void * arr, int start, int items, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r)){
	int i, j;
	for(i=start; i< start + items - 1; i++)
	{
		int min_j = i;
		for(j=i + 1; j< start + items; j++)
		{
			if(compare(arr, min_j, j) > 0)
				min_j = j;
		}
		if(i!=min_j)
			exchange(arr, i, min_j);
	}
}

void basic_sort(void * arr, int items, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r)){
	basic_sort_run(arr, 0, items, compare, exchange);
}


void merge_sort_run(void * arr, int start, int items, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2))
{
	if(items > 11)
	{
		int half_point = items/2;
		merge_sort_run(arr, start, half_point, compare, exchange, merge);
		merge_sort_run(arr, start + half_point, items - half_point, compare, exchange, merge);
		merge(arr, start, half_point, items - half_point);
	}
	else
	{
		basic_sort_run(arr, start, items, compare, exchange);
	}
}
void merge_sort(void * arr, int arr_size, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2))
{
	merge_sort_run(arr, 0, arr_size, compare, exchange, merge);
}

unsigned int calc_end_pos(unsigned int p, char * cigar, unsigned int * all_skipped_len, int * is_exonic_regions, global_context_t * global_context){
	unsigned int cursor = p, tmpi=0;
	int nch, cigar_cursor;
	for(cigar_cursor = 0; 0!=(nch = cigar[cigar_cursor]); cigar_cursor++){
		if(isdigit(nch)){
			tmpi = tmpi * 10 + (nch - '0');
		}else{
			if((nch == 'S' && cursor == p) || nch == 'M' || nch == 'N' || nch == 'D'){
				if(nch == 'M' && global_context -> exonic_region_bitmap){
					if(global_context -> config.do_breakpoint_detection){
						if(is_pos_in_annotated_exon_regions(global_context, cursor) == 0 || is_pos_in_annotated_exon_regions(global_context, cursor + tmpi - 1) == 0)  ( * is_exonic_regions) = 0;
					} else {
						if(is_pos_in_annotated_exon_regions(global_context, cursor + tmpi/2) == 0) ( * is_exonic_regions) = 0;
					}
				}

				cursor += tmpi;
				if(nch == 'N' || nch == 'D')(*all_skipped_len) += tmpi;
			}
			tmpi = 0;
		}
	}
	return cursor;

}

void test_PE_and_same_chro_cigars(global_context_t * global_context , unsigned int pos1, unsigned int pos2, int * is_exonic_regions, int * is_PE_distance, int * is_same_chromosome, int read_len_1, int read_len_2, char * cigar1, char * cigar2, char *read_name, int * res_tlen){
 	char * r1_chr = NULL, * r2_chr = NULL;
	int r1_pos, r2_pos;

	(*is_same_chromosome) = 0;
	(*is_PE_distance) = 0;
	(*is_exonic_regions) = 1;

	locate_gene_position(pos1, &global_context -> chromosome_table, & r1_chr, & r1_pos);
	locate_gene_position(pos2, &global_context -> chromosome_table, & r2_chr, & r2_pos);

	if(r1_chr == r2_chr){
		unsigned int skip_1 = 0;
		unsigned int skip_2 = 0;
		unsigned int r1_end_pos = calc_end_pos(pos1, cigar1, &skip_1, is_exonic_regions, global_context  );
		unsigned int r2_end_pos = calc_end_pos(pos2, cigar2, &skip_2, is_exonic_regions, global_context  );

		unsigned int tlen = max(r1_end_pos, r2_end_pos) - min(pos1, pos2);
		if(tlen > skip_1) tlen -= skip_1;
		if(tlen > skip_2) tlen -= skip_2;

		(*is_same_chromosome) = 1;

		if(tlen >= global_context -> config.minimum_pair_distance && tlen <= global_context -> config.maximum_pair_distance)
			(* is_PE_distance) = 1;
		*res_tlen = tlen;
	}else *res_tlen = 0x7fffffff;
}

void test_PE_and_same_chro_align(global_context_t * global_context , realignment_result_t * res1, realignment_result_t * res2, int * is_exonic_regions, int * is_PE_distance, int * is_same_chromosome, int read_len_1, int read_len_2, char * read_name, int *res_tlen){
	return test_PE_and_same_chro_cigars(global_context, res1 -> first_base_position, res2 -> first_base_position, is_exonic_regions, is_PE_distance, is_same_chromosome , read_len_1 , read_len_2, res1 -> cigar_string, res2 -> cigar_string, read_name, res_tlen);
}



void test_PE_and_same_chro(global_context_t * global_context , unsigned int pos1, unsigned int pos2, int * is_PE_distance, int * is_same_chromosome, int read_len_1, int read_len_2)
{
 	char * r1_chr, * r2_chr;
	int r1_pos, r2_pos;

	(*is_same_chromosome) = 0;
	(*is_PE_distance) = 0;

	int re1 = locate_gene_position(pos1, &global_context -> chromosome_table, & r1_chr, & r1_pos);
	int re2 = locate_gene_position(pos2, &global_context -> chromosome_table, & r2_chr, & r2_pos);

	if(re1 ==0 && 0 ==re2) {
		long long tlen = r1_pos;
		tlen  -= r2_pos;
		tlen = abs(tlen);
		tlen += (r1_pos > r2_pos)?read_len_1:read_len_2;
		unsigned int tlenI = (unsigned int) tlen;

		//SUBREADprintf("TEST PE: %p == %p , TLEN=%u\n", r1_chr, r2_chr, tlenI);

		if(r1_chr == r2_chr){
			(*is_same_chromosome) = 1;
			if(tlenI >= global_context -> config.minimum_pair_distance && tlenI <= global_context -> config.maximum_pair_distance)
				(* is_PE_distance) = 1;
		}
	}
}

int FIXLENstrcmp(char * flen, char * rname){
	char * tt=NULL;
	char * fixed_len = malloc(strlen(flen)+1);
	strcpy(fixed_len, flen);
	char * s1 = strtok_r(fixed_len, "\n", &tt);

	int found = 0;
	while(1){
		int this_ok = 1,x=0;
		for(; s1[x]; x++){
			if(rname[x]!=s1[x]) {
				this_ok = 0;
				break;
			}
		}
		if(this_ok){
			found = 1;
			break;
		}
		s1 = strtok_r(NULL, "\n", &tt);
		if(!s1)break;
	}
	free(fixed_len);
	return !found;
}
