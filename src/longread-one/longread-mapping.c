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
#include <stdarg.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/time.h>
#include <sys/types.h>
#ifndef __MINGW32__
#include <sys/resource.h>
#endif
#include <sys/timeb.h>
#include <sys/stat.h>
#include <locale.h>
#include <ctype.h>
#include <unistd.h>
#include <getopt.h>

#include "LRMconfig.h"
#include "LRMhelper.h"
#include "LRMsorted-hashtable.h"
#include "LRMbase-index.h"
#include "LRMchro-event.h"
#include "LRMfile-io.h"

int LRMvalidate_and_init_context(LRMcontext_t ** context, int argc, char ** argv);

#ifdef MAKE_STANDALONE
int main(int argc, char ** argv){
#else
int longread_mapping_R(int argc, char ** argv){
#endif
	int retv=0;
	
	LRMcontext_t *context = NULL;
	retv = retv || LRMvalidate_and_init_context(&context, argc, argv);
	retv = retv || LRMshow_conf(context);
	retv = retv || LRMrun_task(context);
	retv = retv || LRMfinalise(context);
	retv = retv || LRMprint_mapping_summary(context);
	retv = retv || LRMdestroy_context(context);
	context = NULL;
	
	return retv;
}

int LRMprint_mapping_summary(LRMcontext_t * context){
	LRMprintf("\n\nAll finished.\n\nTotal processed reads : %d\n", context -> all_processed_reads);
	LRMprintf("Mapped reads: %u (%.1f%%)\n", context->mapped_reads, context->mapped_reads*100./context->all_processed_reads);
	LRMprintf("Time: %.1f minutes\n\n" , (LRMmiltime() - context->start_running_time)/60);
	return 0;
}


void LRMprint_version(){
	LRMprintf("\nSublong v%s\n", SUBREAD_VERSION);
	LRMputs("");
}
void LRMprint_usage(){
	LRMprint_version();
	LRMputs("Usage:");
	LRMputs("");
	LRMputs("./sublong [options] -i <index_name> -r <input> -o <output>");
	LRMputs("");
	LRMputs("## Mandatory arguments:");
	LRMputs("");
	LRMputs(" -i <string>   Base name of the index. The index must be built as a full index");
	LRMputs("               and has only one block.");
	LRMputs("");
	LRMputs(" -r <string>   Name of an input read file. Acceptable formats include gzipped");
    LRMputs("               FASTQ and FASTQ (automatically identified).");
	LRMputs("");
	LRMputs(" -o <string>   Name of an output file in BAM format.");
	LRMputs("");
	LRMputs("## Optional arguments:");
	LRMputs("# input reads and output");
	LRMputs("");
	LRMputs(" --SAMoutput   Save mapping results in SAM format.");
	LRMputs("");
	LRMputs("# number of CPU threads");
	LRMputs("");
	LRMputs(" -T <int>      Number of CPU threads used. 1 by default.");
	LRMputs("");
	LRMputs("# others");
	LRMputs("");
	LRMputs(" -v            Output version of the program.");
	LRMputs("");
	LRMputs(" -X            Turn on the RNA-seq mode.");
	LRMputs("");
	LRMputs("Refer to Users Manual for detailed description to the arguments.");
	LRMputs("");
}


static struct option long_options[] ={
	{"SAMoutput", no_argument, 0, 's'},
	{0, 0, 0, 0}
};


int LRMvalidate_and_init_context(LRMcontext_t ** context, int argc, char ** argv){
	int c;
	
	(*context) = malloc(sizeof(LRMcontext_t));
	memset((*context), 0, sizeof(LRMcontext_t));
	LRMset_default_values_context(*context);

	(*context) -> input_file_name[0] = 0;
	(*context) -> output_file_name[0] = 0;	
	(*context) -> index_prefix[0] = 0;
	
	optind = 0;
	opterr = 1;
	optopt = 63;
	int option_index = 0;
	while ((c = getopt_long (argc, argv, "Xr:i:o:B:T:v", long_options, &option_index))!=-1){
		switch(c){
			case 'B':
				(*context)->multi_best_read_alignments = min(max(1,atoi(optarg)),20);
				break;
			case 'X':
				(*context) -> is_RNAseq_mode = 1;
				(*context)->dynamic_programming_score_match = 6;
				(*context)->dynamic_programming_score_mismatch = 0;
				(*context)->dynamic_programming_score_create_gap = -6;
				(*context)->dynamic_programming_score_extend_gap = -1;

				break;
			case 'r':
				strcpy((*context) -> input_file_name, optarg);
				break;
			case 'i':
				strcpy((*context) -> index_prefix, optarg);
				break;
			case 's':
				(*context) -> max_cigars_in_read = 99900;
				(*context) -> is_SAM_output = 1;
				break;
			case 'o':
				strcpy((*context) -> output_file_name, optarg);
				break;
			case 'T':
				(*context) -> threads = min(max(1,atoi(optarg)),LRMMAX_THREADS);
				break;
			case 'm':
				(*context) -> min_voting_number = atoi(optarg);
				break;
			case 'v':
				LRMprint_version();
				return 1;
			case '?':
			default:
				LRMprint_usage();
				return 1;
		}
	}

	if((*context) -> input_file_name[0] ==0 || (*context) -> output_file_name[0]==0 || (*context) -> index_prefix[0]==0){
		LRMprintf("Please specify the input, output files and the index.\n");
		LRMprint_usage();
		return 1;
	}

	(*context) -> user_command_line[0]=0;
	for(c = 0; c<argc;c++)
		sprintf((*context) -> user_command_line+strlen( (*context) -> user_command_line), "\"%s\" ", argv[c]);
	

	LRMthread_lockinit(&(*context) -> input_lock);
	LRMthread_lockinit(&(*context) -> sam_bam_file_lock);
	
	(*context)-> sam_bam_chromosome_table = LRMHashTableCreate(199);
	LRMHashTableSetKeyComparisonFunction((*context)-> sam_bam_chromosome_table, LRMhash_strcmp);
	LRMHashTableSetHashFunction((*context)-> sam_bam_chromosome_table, LRMhash_strhash);
	LRMHashTableSetDeallocationFunctions((*context)-> sam_bam_chromosome_table, NULL, NULL);
	

	(*context)-> chromosome_size_list = LRMArrayListCreate(29);

	(*context)-> chromosome_size_table = LRMHashTableCreate(199);
	LRMHashTableSetKeyComparisonFunction((*context)-> chromosome_size_table, LRMhash_strcmp);
	LRMHashTableSetHashFunction((*context)-> chromosome_size_table, LRMhash_strhash);
	LRMHashTableSetDeallocationFunctions((*context)-> chromosome_size_table, free, NULL);

	(*context) -> sam_bam_chromosome_list = LRMArrayListCreate(29);
	
	LRMload_offsets(*context);

	int retv = LRMgeinput_open((*context)->input_file_name,&(*context) -> input_file);

	(*context)->sam_bam_file = fopen( (*context) -> output_file_name, "w");
	if(NULL == (*context)->sam_bam_file) retv = 1;
	
	(*context)->event_space = malloc(sizeof(LRMevent_t)*20000);
	(*context)->event_space_size = 20000;
	LRMthread_lockinit(&(*context) -> event_space_lock);
	(*context)->events_realignment  = LRMHashTableCreate(320000);

	return retv;
}

#ifndef MAKE_STANDALONE
#define LRM_CLOCK_USE_GETTIME
#endif
double LRMmiltime(){
	double ret;
	#ifdef LRM_CLOCK_USE_GETTIME
	struct timespec tsc;
	clock_gettime(CLOCK_REALTIME, &tsc);
	ret = tsc.tv_sec*1. + tsc.tv_nsec*1./1000/1000/1000;
	#else
	struct timeb trp;
	ftime(&trp);
	ret = trp.time*1.0+(trp.millitm*1.0/1000.0);
	#endif

	return ret;
}


void LRMset_default_values_context(LRMcontext_t * context){
	context->threads = 1;
	context->start_running_time = LRMmiltime();
	context->multi_best_read_alignments = 1;

	context->max_dynamic_indel_length = 5000;
	context->max_read_indel_length = 0;
	context->max_junction_distance = 100000;
	context->min_matched_bases_in_alignment =40;
	context->max_mismatched_bases_in_subread = 0;
	context->min_voting_number = 1;
	context->result_merge_tolerance = LRMDYNAMIC_MAXIMUM_GAP_LENGTH;
	context->max_cigars_in_read = 65000;

	context->dynamic_programming_score_match = 6;
	context->dynamic_programming_score_mismatch = 0;
	context->dynamic_programming_score_create_gap = -6;
	context->dynamic_programming_score_extend_gap = -1;
}

int LRMshow_conf(LRMcontext_t * context){
	LRMprintf("\n ====== Subread long read mapping %s======\n\n", context -> is_RNAseq_mode? "(RNA-seq) ": "");
	LRMprintf("Threads: %d\n" , context->threads);
	LRMprintf("Input file: %s\n" , context->input_file_name);
	LRMprintf("Output file: %s (%s)\n" , context->output_file_name,  context->is_SAM_output?"SAM":"BAM");
	LRMprintf("Index: %s\n\n" , context->index_prefix);
	
	return 0;
}

int LRMinit_chunk(LRMcontext_t * context){
	return 0;
}

int LRMrun_task(LRMcontext_t * context){
	int retv = 0;
	retv = LRMload_index( context );
	LRMprintf("Index was loaded; the gap bewteen subreads is %d bases\n", context -> current_index.index_gap );
	while(!(retv ||LRMinput_has_finished( context ))){
		retv=retv || LRMinit_chunk(context);
		retv=retv || LRMsave_input_pos(context);
		retv=retv || LRMiterate_reads( context , LRMRUNNING_STEP_VOTING);
		retv=retv || LRMfinalise_chunk_reads(context);
	}
	return retv;
}

int LRMfinalise(LRMcontext_t * context){
	return 0;
}

int LRMdestroy_context(LRMcontext_t * context){

	LRMgehash_destory(&(context -> current_index));
	LRMgvindex_destory(&(context -> current_base_index));
	
	LRMHashTableDestroy(context -> chromosome_size_table);
	LRMArrayListDestroy(context -> chromosome_size_list);

	LRMHashTableDestroy(context -> sam_bam_chromosome_table);
	LRMArrayListDestroy(context -> sam_bam_chromosome_list);
	
	LRMHashTableSetDeallocationFunctions(context -> events_realignment, NULL, free);
	LRMHashTableDestroy(context -> events_realignment);

	free(context -> event_space);

	if(!context -> is_SAM_output){
		fwrite(context -> bam_file_tail_binary,1, context -> bam_file_tail_length, context->sam_bam_file);
	}

	LRMgeinput_close(&context->input_file);
	fclose(context->sam_bam_file);
	//free(context->user_command_line);
	free(context);
	return 0;
}


int LRMinput_has_finished( LRMcontext_t * context ){
	return context -> input_exhausted ;
}

int LRMload_index(LRMcontext_t * context){
	int retv = 0;
	char indextab_fname[LRMMAX_FILENAME_LENGTH + 20];

	sprintf(indextab_fname, "%s.00.b.tab", context -> index_prefix);
	retv = retv || LRMgehash_load(&(context -> current_index), indextab_fname);

	sprintf(indextab_fname, "%s.00.b.array", context -> index_prefix);
	retv = retv || LRMgvindex_load(&(context -> current_base_index), indextab_fname);
	
	return retv;
} 


int LRMiterate_reads( LRMcontext_t * context, int task ){
	int retv = 0;
	retv = retv || LRMstart_thread( context , task );
	retv = retv || LRMwait_threads( context );
	retv = retv || LRMmerge_threads( context, task );
	return retv;
}

void * LRM_thread_runner (void * args){
	void ** argv = args;
	LRMcontext_t * context = argv[0];
	int thid = argv[1]-NULL;
	int task = argv[2]-NULL;	
	free(args);

	LRMchunk_read_iteration(context, thid, task);	
	
	return NULL;
}

int LRMstart_thread_init_context(LRMcontext_t * context, int thread_id, int step){
	LRMthread_context_t * thread_context = context -> thread_contexts+thread_id;
	memset(thread_context, 0, sizeof(LRMthread_context_t));
	thread_context->thread_id = thread_id;
	
	if( thread_context -> thread_id == 0 )LRMsambam_write_header(context, thread_context);


	thread_context -> mapped_reads = 0;
	thread_context -> out_SAMBAM_buffer = malloc(3000000); 
	if(thread_context -> out_SAMBAM_buffer == NULL) return 1;
	thread_context -> out_buff_used = 0;		     
	thread_context -> out_buff_capacity = 3000000;

	thread_context -> dynamic_programming_movement_buffer = malloc(( 2* LRMINDEL_DYNAMIC_CHANNEL_TOLERANCE + 1) * LRMDYNAMIC_MAXIMUM_GAP_LENGTH);
	thread_context -> dynamic_programming_score_buffer = malloc(sizeof(int) * ( 2 * LRMINDEL_DYNAMIC_CHANNEL_TOLERANCE + 1) *( LRMDYNAMIC_MAXIMUM_GAP_LENGTH+1));
	thread_context -> dynamic_programming_indel_movement_buf = malloc( LRMMAX_READ_LENGTH + 1 );
	thread_context -> final_cigar_string = malloc( LRMMAX_READ_LENGTH + 1 );
	return 0;
}

int LRMstart_thread(LRMcontext_t * context, int task ){
	int th_id, retv=0;
		
	for(th_id=0; th_id<context -> threads; th_id++){
		
		retv = retv || LRMstart_thread_init_context(context,th_id,task);
		if(retv)
			break;
		else {
			void ** th_args=malloc(sizeof(void *)*3); // to be freed in the thread.
			th_args[0] = context;
			th_args[1] = NULL + th_id;
			th_args[2] = NULL + task;
			LRMpthread_create(context -> running_threads+th_id, NULL, LRM_thread_runner, th_args);
		}
	}
	
	return retv;
}

int LRMwait_threads( LRMcontext_t * context ){
	int th_id;
	for(th_id=0; th_id<context -> threads; th_id++)
		LRMpthread_join(context -> running_threads[th_id], NULL);
	return 0;
}

void LRMmerge_threads_destroy_context(LRMcontext_t * context, LRMthread_context_t * thread_context, int task){
	free(thread_context -> dynamic_programming_movement_buffer);
	free(thread_context -> dynamic_programming_score_buffer);
	free(thread_context -> dynamic_programming_indel_movement_buf);
	free(thread_context -> final_cigar_string);
}

int LRMmerge_threads( LRMcontext_t * context , int step){
	int retv = 0;
	int th_id;

	for(th_id=0; th_id<context -> threads; th_id++){
		LRMwrite_chunk_check_buffer_write(context,  context -> thread_contexts+th_id, 1);
		LRMmerge_threads_destroy_context(context,  context -> thread_contexts+th_id, step);
		if(th_id == context -> threads-1)LRMbam_generate_tail_binary(context,  context -> thread_contexts+th_id);
		context -> mapped_reads += context -> thread_contexts[th_id].mapped_reads;
	}
	
	return retv;
}

int LRMrewind_input_pos(LRMcontext_t * context){
	context -> processed_reads_in_chunk = 0;
	if(context->input_file.file_type == LRMGENE_INPUT_GZIP_FASTQ)
		LRMseekgz_seek(context->input_file.input_fp, &context->last_saved_zlib_pos);
	else
		fseeko(context->input_file.input_fp, context->last_saved_raw_pos, SEEK_SET);
	return 0;
}

int LRMsave_input_pos( LRMcontext_t * context){
	context -> processed_reads_in_chunk = 0;
	if(context->input_file.file_type == LRMGENE_INPUT_GZIP_FASTQ)
		LRMseekgz_tell(context->input_file.input_fp, &context->last_saved_zlib_pos);
	else
		context -> last_saved_raw_pos = ftello(context->input_file.input_fp);
		
	return 0;
}

int LRMdo_one_voting_read_process_samechro(LRMcontext_t * context, unsigned int p1, unsigned int p2){
	char * chro_name1, *chro_name2;
	int chro_pos1, chro_pos2;
	LRMlocate_gene_position(context, p1, &chro_name1, & chro_pos1);
	LRMlocate_gene_position(context, p2, &chro_name2, & chro_pos2);
	
	return chro_name1 == chro_name2; // they can be compared in this way because they are pointers in the sam_bam_chromosome_list.
}



void LRMreverse_read_and_qual(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context){
	LRMreverse_read(iteration_context -> read_text, iteration_context -> read_length);
	LRMreverse_quality(iteration_context -> qual_text, iteration_context -> read_length);
}


void LRMcalc_total_subreads(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context){
	if(iteration_context->read_length < 16){
		iteration_context -> subread_total = 0;
		return;
	}
	int maximum_subreads = iteration_context->read_length - 16 + 1;
	maximum_subreads = maximum_subreads/3;
	//LRMprintf("TOTAL SUBREADS=%d\n",  min(maximum_subreads, LRMMAX_SUBREAD_PER_READ_HARDLIMIT));
	iteration_context -> subread_total = min(maximum_subreads, LRMMAX_SUBREAD_PER_READ_HARDLIMIT);
	int subread_last =  iteration_context -> subread_total -1;
	iteration_context -> subread_extract_gap = (iteration_context->read_length - 16)*1.0 / ( subread_last ) + 0.000001;
}

int LRMcalc_subread_start(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, int subread_no){
	if(iteration_context->read_length < 16)return 0;
	int subread_last =  iteration_context -> subread_total -1;
	if(subread_no < subread_last){
		int ret = (int)(iteration_context -> subread_extract_gap * subread_no);
		assert(ret <= iteration_context->read_length - 16);
		return ret;
	}
	return iteration_context->read_length - 16;
}

void LRMdo_one_voting_read_single(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context){
	int subread_i;
	LRMcalc_total_subreads(context,thread_context,iteration_context);
	for(subread_i=0; subread_i< iteration_context -> subread_total; subread_i++){
		int subread_start = LRMcalc_subread_start(context,thread_context,iteration_context,subread_i);
		char * subread_text = iteration_context->read_text + subread_start;
		unsigned int subread_int = LRMgenekey2int(subread_text);

		LRMgehash_go_QQ(context, thread_context, iteration_context,& context->current_index, subread_int , subread_start, iteration_context->read_length, iteration_context -> is_reversed, & iteration_context-> vote_table, context -> max_read_indel_length, subread_i);
	}
}

void LRMdo_one_voting_read(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context){
	LRMinit_gene_vote((& iteration_context-> vote_table));
	
	for(iteration_context->is_reversed = 0; iteration_context->is_reversed<2; iteration_context->is_reversed++){
		LRMdo_one_voting_read_single(context, thread_context, iteration_context);
		if(0 == iteration_context->is_reversed) LRMreverse_read_and_qual(context, thread_context, iteration_context);
	}
	iteration_context->is_reversed=1;
}

int LRM_longvote_location_compare(void * arr, int l, int r){
	LRMread_iteration_context_t * iteration_context = arr;
	if(iteration_context -> sorting_vote_locations[l] > iteration_context -> sorting_vote_locations[r])return 1;
	if(iteration_context -> sorting_vote_locations[l] < iteration_context -> sorting_vote_locations[r])return -1;
	return 0;
}

void LRM_longvote_location_exchange(void * arr, int l, int r){
	LRMread_iteration_context_t * iteration_context = arr;
	unsigned int t;
	unsigned short ts;
	t = iteration_context -> sorting_vote_locations[l];
	iteration_context -> sorting_vote_locations[l] = iteration_context -> sorting_vote_locations[r];
	iteration_context -> sorting_vote_locations[r] = t;

	t = iteration_context -> sorting_subread_nos[l];
	iteration_context -> sorting_subread_nos[l] = iteration_context -> sorting_subread_nos[r];
	iteration_context -> sorting_subread_nos[r] = t;

	ts = iteration_context -> sorting_subread_votes[l];
	iteration_context -> sorting_subread_votes[l] = iteration_context -> sorting_subread_votes[r];
	iteration_context -> sorting_subread_votes[r] = ts;

	ts = iteration_context -> sorting_is_negative_strand[l];
	iteration_context -> sorting_is_negative_strand[l] = iteration_context -> sorting_is_negative_strand[r];
	iteration_context -> sorting_is_negative_strand[r] = ts;
}

void LRM_longvote_location_merge(void * arr,  int start, int items, int items2){
	LRMread_iteration_context_t * iteration_context = arr;
	unsigned int * tmp_locs = malloc(sizeof(unsigned int) * (items + items2));
	unsigned int * tmp_subread_nos = malloc(sizeof(unsigned int) * (items + items2));
	unsigned int * tmp_negative = malloc(sizeof(unsigned int) * (items + items2));
	unsigned short * tmp_votes = malloc(sizeof(unsigned short) * (items + items2));

	int rcur_1 = start, rcur_2 = start + items, wcur = 0;
	for(wcur=0; wcur < items + items2 ; wcur++){
		int select_from_1 = rcur_1 < start + items;
		if(select_from_1 && rcur_2 < start + items + items2)select_from_1 = (LRM_longvote_location_compare(arr, rcur_1, rcur_2)<0);
		if(select_from_1){
			tmp_locs[wcur] = iteration_context -> sorting_vote_locations[rcur_1];
			tmp_subread_nos[wcur] = iteration_context -> sorting_subread_nos[rcur_1];
			tmp_votes[wcur] = iteration_context -> sorting_subread_votes[rcur_1];
			tmp_negative[wcur] = iteration_context -> sorting_is_negative_strand[rcur_1++];
		}else{
			tmp_locs[wcur] = iteration_context -> sorting_vote_locations[rcur_2];
			tmp_subread_nos[wcur] = iteration_context -> sorting_subread_nos[rcur_2];
			tmp_votes[wcur] = iteration_context -> sorting_subread_votes[rcur_2];
			tmp_negative[wcur] = iteration_context -> sorting_is_negative_strand[rcur_2++];
		}
	}
	memcpy(iteration_context -> sorting_vote_locations +  start , tmp_locs, sizeof(unsigned int) * (items + items2));
	memcpy(iteration_context -> sorting_subread_nos +  start , tmp_subread_nos, sizeof(unsigned int) * (items + items2));
	memcpy(iteration_context -> sorting_subread_votes +  start , tmp_votes, sizeof(unsigned short) * (items + items2));
	memcpy(iteration_context -> sorting_is_negative_strand +  start , tmp_negative, sizeof(unsigned int) * (items + items2));

	free(tmp_locs);
	free(tmp_subread_nos);
	free(tmp_votes);
	free(tmp_negative);
}

void LRMprint_longvote( LRMcontext_t * context, LRMthread_context_t * thread_context,LRMread_iteration_context_t * iteration_context){
	int ii;
	LRMprintf("\n+++ READ %s ++++++++++++++++++++\n", iteration_context -> read_name);

	if(0)
	for(ii = 0; ii< iteration_context -> sorting_total_votes; ii++){
		char postxt[55];
		LRMgene_vote_t * vtab = &(iteration_context->vote_table);
		if( iteration_context -> sorting_subread_votes[ii] >0){
			int bb = iteration_context -> sorting_subread_nos[ii] >> 16;
			int tt = iteration_context -> sorting_subread_nos[ii] & 0xffff;
			LRMpos2txt(context,  iteration_context -> sorting_vote_locations[ii], postxt);
			LRMprintf("  %c %d ~ %d  %u (%s) : %d votes    Array %d %d\n", iteration_context -> sorting_is_negative_strand[ii]?'-':'+', vtab -> coverage_start[bb][tt], vtab -> coverage_end[bb][tt], iteration_context -> sorting_vote_locations[ii], postxt, iteration_context -> sorting_subread_votes[ii], bb, tt);
		}
	}
	LRMprintf("===========================================================\n");
	for(ii=0; ii<LRMMAX_SORTED_WINDOW_NO_HARD_LIMIT; ii++){
		if(iteration_context -> sorted_window_vote_stop[ii]<1)break;

		int sorted_no = iteration_context -> sorted_window_vote_start[ii];
		char postxt[55];
		LRMpos2txt(context,  iteration_context -> sorting_vote_locations[sorted_no], postxt);
		LRMprintf("  Window %s subread no: %d ~ %d : Votes = %d ; Pos = %u (%s)\n", iteration_context -> sorted_window_is_negative_strand[ii]?"NEG":"POS", iteration_context -> sorted_window_vote_start[ii], iteration_context -> sorted_window_vote_stop[ii], iteration_context -> sorted_window_total_votes[ii], iteration_context -> sorting_vote_locations[sorted_no], postxt);
	}

	LRMprintf(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n");
	unsigned int prev_end_chro = 0;
	int prev_end_read = 0;
	for(ii = 0; ii < iteration_context -> chain_total_items ; ii++){
		char postxt[55];
		long long to_prev_delta = prev_end_chro + iteration_context -> chain_cov_start[ii] - prev_end_read;
		to_prev_delta -=  iteration_context -> chain_chro_at_cov_start[ii];
		if(ii < 1) to_prev_delta = 0;
		LRMpos2txt(context,  iteration_context -> chain_chro_at_cov_start[ii], postxt);
		#ifndef __MINGW32__
		LRMprintf("  >> COV = %d ~ %d  : Pos = %u (%s)  DELTA_PREV=%lld\n", iteration_context -> chain_cov_start[ii], iteration_context -> chain_cov_end[ii], iteration_context -> chain_chro_at_cov_start[ii], postxt, to_prev_delta);
		#endif
		prev_end_chro = iteration_context -> chain_chro_at_cov_start[ii] + (iteration_context -> chain_cov_end[ii] - iteration_context -> chain_cov_start[ii]);
		prev_end_read = iteration_context -> chain_cov_end[ii];
	}

	LRMprintf("\n");
}

void LRMcopy_longvotes_to_itr( LRMcontext_t * context, LRMthread_context_t * thread_context,LRMread_iteration_context_t * iteration_context){
	int bb, ii;
	iteration_context -> sorting_total_votes = 0;
	LRMgene_vote_t * vtab = &(iteration_context->vote_table);
	for(bb =0; bb< LRMGENE_VOTE_TABLE_SIZE; bb++){
		for(ii = 0; ii< vtab -> items[bb]; ii++){
		//	if( vtab -> votes[bb][ii]<2)continue;
			iteration_context -> sorting_vote_locations[iteration_context -> sorting_total_votes] = vtab -> pos[bb][ii] + vtab -> coverage_start[bb][ii];
			iteration_context -> sorting_is_negative_strand[iteration_context -> sorting_total_votes] = (vtab ->masks[bb][ii] & LRMIS_NEGATIVE_STRAND)? 1:0;
			iteration_context -> sorting_subread_nos[iteration_context -> sorting_total_votes] = ii | (bb<<16);
			iteration_context -> sorting_subread_votes[iteration_context -> sorting_total_votes] = vtab -> votes[bb][ii];
			iteration_context -> sorting_total_votes ++;
		}
	}
}

int LRMis_gap_in_used_gap(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, unsigned int vote_start_pos, unsigned int vote_cov_len){
	unsigned int vote_end_pos = vote_start_pos + vote_cov_len;
	int xxi;
	for(xxi = 0; xxi < iteration_context -> chain_used_gaps -> numOfElements; xxi+=2){
		unsigned int used_start = LRMArrayListGet(iteration_context -> chain_used_gaps, xxi) - NULL;
		unsigned int used_end = LRMArrayListGet(iteration_context -> chain_used_gaps, xxi+1) - NULL;

		//LRMprintf("ADD_VOTE? OLD: %u ~ %u ; NEXT: %u ~ %u\n", used_start, used_end, vote_start_pos, vote_end_pos);

		if(used_start <= vote_start_pos && used_end >= vote_start_pos) return 1;
		if(vote_start_pos <= used_start && vote_end_pos >= used_start) return 1;
	}
	return 0;
}



#define LRM_SORTED_VOTE_WINDOW_SIZE 1500000

void LRMfind_top_windows( LRMcontext_t * context, LRMthread_context_t * thread_context,LRMread_iteration_context_t * iteration_context){
	int ii, max_window_no = LRMMAX_SORTED_WINDOW_NO_HARD_LIMIT, is_reversed; 
	int testing_window_width = min(LRM_SORTED_VOTE_WINDOW_SIZE, iteration_context -> read_length *4/3);
	if(iteration_context -> read_length < 10000) testing_window_width = 20000;
	if(context -> is_RNAseq_mode ) testing_window_width = min(1000000, testing_window_width);
	if(iteration_context -> sorting_total_votes<1){
		iteration_context -> sorted_window_total_votes[0]=0;
		return;
	}

	memset(iteration_context -> sorted_window_total_votes, 0, sizeof(int) * LRMMAX_SORTED_WINDOW_NO_HARD_LIMIT);

	for(is_reversed = 0; is_reversed < 2; is_reversed++){
		unsigned int votes_in_win = 0, win_start_ii = 0;
		for(win_start_ii = 0; win_start_ii < iteration_context -> sorting_total_votes && iteration_context -> sorting_is_negative_strand[win_start_ii] != is_reversed; win_start_ii++)
			win_start_ii ++;

		if(win_start_ii == iteration_context -> sorting_total_votes) continue;

		for(ii = win_start_ii; ii< iteration_context -> sorting_total_votes; ii++){
			if( iteration_context -> sorting_is_negative_strand[ii] != is_reversed )
				continue;


			if(0 && LRMFIXLENstrcmp("R4", iteration_context->read_name) == 0)
					LRMprintf("    PASSING %d += %d : Vote in win = %d\n", ii, iteration_context -> sorting_subread_votes[ii], votes_in_win);
			votes_in_win += iteration_context -> sorting_subread_votes[ii];

			while(win_start_ii < iteration_context -> sorting_total_votes && win_start_ii < ii){
				while(iteration_context -> sorting_is_negative_strand[win_start_ii] != is_reversed) win_start_ii ++;

				if(iteration_context -> sorting_vote_locations[ii] - iteration_context -> sorting_vote_locations[win_start_ii] < testing_window_width && !LRMis_gap_in_used_gap(context, thread_context, iteration_context,iteration_context -> sorting_vote_locations[win_start_ii], 15))break;

				votes_in_win -= iteration_context -> sorting_subread_votes[win_start_ii];
				win_start_ii ++;
				while(iteration_context -> sorting_is_negative_strand[win_start_ii] != is_reversed){
					if(win_start_ii >= iteration_context -> sorting_total_votes) break;
					win_start_ii ++;
				}
			}

			if(0 && LRMFIXLENstrcmp("R4", iteration_context->read_name) == 0)
					LRMprintf("    PASSING2 %d += %d : Vote in win = %d\n", ii, iteration_context -> sorting_subread_votes[ii], votes_in_win);
	
			if(0 && LRMFIXLENstrcmp("R4", iteration_context->read_name) == 0 && votes_in_win > 8200){
				LRMprintf("    EXCHANGE: %d > %d ; window = %d ~ %d\n=================================\n", votes_in_win , iteration_context -> sorted_window_total_votes[max_window_no-1], win_start_ii, ii);
			}

			//LRMprintf("LRMtest_window_in : %u ~ %u\n", iteration_context->sorting_vote_locations[win_start_ii], iteration_context->sorting_vote_locations[ii] - iteration_context -> sorting_vote_locations[win_start_ii] + 16);
			if(votes_in_win > iteration_context -> sorted_window_total_votes[max_window_no-1] && ! LRMis_gap_in_used_gap(context, thread_context, iteration_context, iteration_context->sorting_vote_locations[win_start_ii], iteration_context->sorting_vote_locations[ii] - iteration_context -> sorting_vote_locations[win_start_ii] + 15)){
				int tins = -1;
				int tt;
				int has_in_window = 0;
	
				for(tt = max_window_no - 1; tt >= 0 ; tt --){
					if( is_reversed == iteration_context -> sorted_window_is_negative_strand[tt] &&
						( (iteration_context -> sorted_window_vote_start[tt] <=win_start_ii && iteration_context -> sorted_window_vote_stop[tt] > win_start_ii)||
						( win_start_ii <= iteration_context -> sorted_window_vote_start[tt] && ii > iteration_context -> sorted_window_vote_start[tt] ))){
						if(votes_in_win >iteration_context -> sorted_window_total_votes[tt]){
							int mm;
							for(mm = tt; mm < max_window_no - 1; mm++){
								iteration_context -> sorted_window_total_votes[mm] = iteration_context -> sorted_window_total_votes[mm+1];
								iteration_context -> sorted_window_vote_start[mm] = iteration_context -> sorted_window_vote_start[mm+1];
								iteration_context -> sorted_window_vote_stop[mm] = iteration_context -> sorted_window_vote_stop[mm+1];
								iteration_context -> sorted_window_is_negative_strand[mm] = iteration_context -> sorted_window_is_negative_strand[mm+1];
							}
							iteration_context -> sorted_window_total_votes[max_window_no - 1] = 0;
						}else{
							has_in_window = 1;
							break;
						}
					}
				}
	
				if(!has_in_window){
					for(tt = max_window_no - 1; tt >= 0 ; tt --){
						if(votes_in_win > iteration_context -> sorted_window_total_votes[tt]) tins = tt;
					}
	
					for(tt = max_window_no-1; tt > tins; tt--){
						iteration_context -> sorted_window_total_votes[tt] = iteration_context -> sorted_window_total_votes[tt -1] ;
						iteration_context -> sorted_window_vote_start[tt] = iteration_context -> sorted_window_vote_start[tt -1] ;
						iteration_context -> sorted_window_vote_stop[tt] = iteration_context -> sorted_window_vote_stop[tt -1] ;
						iteration_context -> sorted_window_is_negative_strand[tt] = iteration_context -> sorted_window_is_negative_strand[tt -1];
					}
	
					iteration_context -> sorted_window_total_votes[tins] = votes_in_win;
					iteration_context -> sorted_window_vote_start[tins] = win_start_ii;
					iteration_context -> sorted_window_vote_stop[tins] = ii;
					iteration_context -> sorted_window_is_negative_strand[tins] =  iteration_context -> sorting_is_negative_strand[ii];
				}
			}
		}
	}
}

#define MAX_INTRON_LENGTH 700000

int LRM_test_chain_extension(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, int window_no, int go_larger, long long N2W_read_dist, long long N2W_chro_dist, unsigned int last_chro_pos, unsigned int this_chro_pos, unsigned int this_cov_len){
	//LRMprintf("TESTING CHAIN: N2W_READ: %lld ; N2W_CHRO: %lld  %s\n", N2W_read_dist, N2W_chro_dist, go_larger?"GO_LARGE":"GO_SMALL");
	if(LRMis_gap_in_used_gap( context, thread_context, iteration_context, this_chro_pos, this_cov_len )) return 0;
	
	if(N2W_read_dist == 0 || N2W_chro_dist == 0) return 0;

	if(go_larger){
		if(N2W_read_dist < 0 || N2W_chro_dist < 0) return 0;
	}else{
		if(N2W_read_dist > 0 || N2W_chro_dist > 0) return 0;
	}

	if(!context -> is_RNAseq_mode)if( abs(N2W_chro_dist) > 50 || abs(N2W_read_dist) > 50 )if(abs(N2W_chro_dist) >= abs(N2W_read_dist)* 14/10 || abs(N2W_chro_dist) <= abs(N2W_read_dist) * 6 / 10) return 0;
	if(context -> is_RNAseq_mode)if( abs(N2W_chro_dist) > 50 || abs(N2W_read_dist) > 50 )if( abs(N2W_chro_dist) >= abs(N2W_read_dist) + MAX_INTRON_LENGTH || abs(N2W_chro_dist) <= abs(N2W_read_dist) * 6 / 10) return 0;

	if(abs(N2W_read_dist) > context -> result_merge_tolerance) return 0;

	if(!LRMdo_one_voting_read_process_samechro(context , last_chro_pos, this_chro_pos) ) return 0;

	return 1;
}

void LRMfix_extension_overlapping(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, int window_no ){
	int ii;
	int last_end_read = iteration_context -> chain_cov_end[0];
	unsigned int last_chro_end = iteration_context -> chain_chro_at_cov_start[0] + iteration_context -> chain_cov_end[0] - iteration_context -> chain_cov_start[0];
	if(0){
		for(ii = 1; ii < iteration_context -> chain_total_items ; ii++){
			LRMprintf("BEFORE_FIX READ COV: %d ~ %d\n", iteration_context -> chain_cov_start[ii], iteration_context -> chain_cov_end[ii]);
		}
	}
	for(ii = 1; ii < iteration_context -> chain_total_items ; ii++){
		int this_start_read = iteration_context -> chain_cov_start[ii];
		int add_to_this = last_end_read - this_start_read;
		if(add_to_this >= 0){
			iteration_context -> chain_cov_start[ii] += add_to_this+1;
			iteration_context -> chain_chro_at_cov_start[ii] += add_to_this+1;
			//LRMprintf("FIXING OVERLAPPING MOVE : %d\n", add_to_this);
		}

		//LRMprintf("ITEM %d / %d ; LAST_END %u <= THIS_START %u  ; READ START END = %d, %d ; LAST_READ_END = %d\n", ii, iteration_context -> chain_total_items , last_chro_end,  iteration_context -> chain_chro_at_cov_start[ii], iteration_context -> chain_cov_start[ii], iteration_context -> chain_cov_end[ii] , last_end_read );
		if(iteration_context -> chain_cov_start[ii]>=iteration_context -> chain_cov_end[ii] || iteration_context -> chain_cov_start[ii] < last_end_read || iteration_context -> chain_chro_at_cov_start[ii] < last_chro_end){
		//	LRMprintf("DELETED COV\n");
			//assert( iteration_context -> chain_cov_start[ii] < iteration_context -> chain_cov_end[ii] );
			int ti;
			for(ti=ii; ti<iteration_context -> chain_total_items-1; ti++){
				iteration_context -> chain_cov_start[ti] = iteration_context -> chain_cov_start[ti+1];
				iteration_context -> chain_cov_end[ti] = iteration_context -> chain_cov_end[ti+1];
				iteration_context -> chain_chro_at_cov_start [ti] = iteration_context -> chain_chro_at_cov_start [ti+1];
			}
			iteration_context -> chain_total_items--;
			last_end_read = iteration_context -> chain_cov_end[ii-1];
			last_chro_end = iteration_context -> chain_chro_at_cov_start [ii-1] + iteration_context -> chain_cov_end[ii-1] - iteration_context -> chain_cov_start[ii-1];
			if(ii < iteration_context -> chain_total_items) ii--;
		}else{
			last_end_read = iteration_context -> chain_cov_end[ii];
			last_chro_end = iteration_context -> chain_chro_at_cov_start [ii] + iteration_context -> chain_cov_end[ii] - iteration_context -> chain_cov_start[ii];
		}
	}
}

void LRMbuild_chains(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, int window_no ){
	int testing_window_width = min(LRM_SORTED_VOTE_WINDOW_SIZE, iteration_context -> read_length *4/3);
	if(iteration_context -> read_length < 10000) testing_window_width = 20000;
	if(iteration_context -> sorted_window_total_votes[window_no] < 1){
		iteration_context -> chain_total_items = iteration_context -> chain_tosmall_items = iteration_context -> chain_tolarge_items = 0;
		return;
	}
	int seed_subread_no = iteration_context -> sorted_window_vote_start[window_no];

	int test_subread_no, go_large, has_overlapping = 0;
	for(test_subread_no = iteration_context -> sorted_window_vote_start[window_no] ; test_subread_no < iteration_context -> sorted_window_vote_stop[window_no] ; test_subread_no++){
		if(iteration_context -> sorted_window_is_negative_strand[window_no] != iteration_context -> sorting_is_negative_strand[test_subread_no])
			continue;
		if(iteration_context -> sorting_subread_votes[test_subread_no] > iteration_context -> sorting_subread_votes[seed_subread_no]) seed_subread_no = test_subread_no;
	}

	LRMgene_vote_t * vtab = &(iteration_context->vote_table);
	
	for(go_large = 0; go_large < 2; go_large ++){
		int last_added_tt, last_added_ii;
		last_added_tt = iteration_context -> sorting_subread_nos[seed_subread_no] >> 16;
		last_added_ii = iteration_context -> sorting_subread_nos[seed_subread_no] & 0xffff;
		unsigned int last_added_read_pos = vtab -> coverage_start[last_added_tt][last_added_ii];
		long long last_added_chro_pos = vtab -> pos[last_added_tt][last_added_ii] + vtab -> coverage_start[last_added_tt][last_added_ii];
		//LRMprintf("SEED %s : [%d:%d] ; pos=%u OR %u\n", iteration_context -> sorted_window_is_negative_strand[window_no]?"NEG":"POS", last_added_tt, last_added_ii, last_added_read_pos, iteration_context -> sorting_vote_locations[seed_subread_no]);

		int last_added_read_edge = go_large?vtab -> coverage_end[last_added_tt][last_added_ii]: vtab -> coverage_start[last_added_tt][last_added_ii];

		if(!go_large){
			iteration_context ->chain_tosmall_items = 1;
			iteration_context ->chain_tolarge_items = 0;
			iteration_context -> chain_cov_start[0] = last_added_read_pos;
			iteration_context -> chain_cov_end[0] = vtab -> coverage_end[last_added_tt][last_added_ii];
			iteration_context -> chain_chro_at_cov_start[0] = last_added_chro_pos;
		}

		int delta = go_large?1:-1;
		int next_test_subread_no = -1;
		for(test_subread_no = seed_subread_no ; test_subread_no >= 0 && test_subread_no < iteration_context -> sorting_total_votes; test_subread_no+= delta){
			int test_tt, test_ii;

			test_tt = iteration_context -> sorting_subread_nos[test_subread_no] >> 16;
			test_ii = iteration_context -> sorting_subread_nos[test_subread_no] &0xffff;

			if( test_ii >= LRMGENE_VOTE_SPACE || test_tt >= LRMGENE_VOTE_TABLE_SIZE) LRMprintf("Error: Table oversize %s , Subr: %d/%d : %d %d\n", iteration_context -> read_name, test_subread_no, iteration_context -> sorting_total_votes, test_tt, test_ii);
			assert(test_ii < LRMGENE_VOTE_SPACE);
			assert(test_tt < LRMGENE_VOTE_TABLE_SIZE);
			int this_read_pos = vtab -> coverage_start[test_tt][test_ii];
			int this_read_end = vtab -> coverage_end[test_tt][test_ii];
			long long this_chro_pos = iteration_context -> sorting_vote_locations[test_subread_no];

			int same_min_subread_test_no=-1;
			long long min_same_DELTA = 99999999999;
			for(same_min_subread_test_no = test_subread_no; same_min_subread_test_no <iteration_context -> sorting_total_votes && same_min_subread_test_no >=0; same_min_subread_test_no+= delta){
				int same_min_tt = iteration_context -> sorting_subread_nos[same_min_subread_test_no] >> 16;
				int same_min_ii = iteration_context -> sorting_subread_nos[same_min_subread_test_no] & 0xffff;
				int same_min_read_pos = vtab -> coverage_start[same_min_tt][same_min_ii];
				if(same_min_read_pos!= this_read_pos){
					next_test_subread_no = same_min_subread_test_no - delta;
					break;
				}
				long long same_min_chro_pos = iteration_context -> sorting_vote_locations[same_min_subread_test_no];
				long long int same_read_dist = same_min_read_pos;
				same_read_dist -= last_added_read_pos;
				long long same_chro_dist = same_min_chro_pos;
				same_chro_dist -= last_added_chro_pos;

				long long same_DELTA = same_chro_dist - same_read_dist;

				if(abs(same_DELTA) < min_same_DELTA){
					min_same_DELTA = abs(same_DELTA);
					test_subread_no = same_min_subread_test_no;
				}
				if(0){
					//int Asame_min_tt = iteration_context -> sorting_subread_nos[same_min_subread_test_no+delta] >> 16;
					//int Asame_min_ii = iteration_context -> sorting_subread_nos[same_min_subread_test_no+delta] & 0xffff;
					//int Asame_min_read_pos = vtab -> coverage_start[Asame_min_tt][Asame_min_ii];
					//if(Asame_min_read_pos == same_min_read_pos || same_min_subread_test_no > test_subread_no)
					//	LRMprintf("SAME_MIN %s TO %d - D=%lld ; LOC # %lld = %u\n", go_large?"LARGE":"SMALL" ,  same_min_read_pos, same_DELTA, same_min_subread_test_no, same_min_chro_pos);
				}

				if(same_min_subread_test_no == iteration_context -> sorting_total_votes -1 || same_min_subread_test_no == 0)
					next_test_subread_no = same_min_subread_test_no;
			}


			test_tt = iteration_context -> sorting_subread_nos[test_subread_no] >> 16;
			test_ii = iteration_context -> sorting_subread_nos[test_subread_no] &0xffff;
			this_read_pos = vtab -> coverage_start[test_tt][test_ii];
			this_read_end = vtab -> coverage_end[test_tt][test_ii];
			this_chro_pos = iteration_context -> sorting_vote_locations[test_subread_no];

			if( test_ii >= LRMGENE_VOTE_SPACE || test_tt >= LRMGENE_VOTE_TABLE_SIZE) LRMprintf("Error: Table oversize %s , Subr: %d/%d : %d %d\n", iteration_context -> read_name, test_subread_no, iteration_context -> sorting_total_votes, test_tt, test_ii);

			int N2W_read_dist = this_read_pos;
			N2W_read_dist -= last_added_read_pos;

			long long N2W_chro_dist = this_chro_pos;
			N2W_chro_dist -= last_added_chro_pos;
			if(abs(N2W_chro_dist) > testing_window_width) break;

			//LRMprintf("CALC N2W_chro_dist : %lld = %lld - %lld \n", N2W_chro_dist, last_added_chro_pos, this_chro_pos);

			int can_add_this = iteration_context -> sorted_window_is_negative_strand[window_no] == iteration_context -> sorting_is_negative_strand[test_subread_no] &&
                               LRM_test_chain_extension(context, thread_context, iteration_context, window_no, go_large, N2W_read_dist, N2W_chro_dist, last_added_chro_pos , this_chro_pos, this_read_end - this_read_pos);
			if(can_add_this){
				int idx = iteration_context ->chain_tosmall_items;
				if(go_large) idx += iteration_context ->chain_tolarge_items;

				int this_added_read_inner = go_large? this_read_pos:this_read_end;
				if( go_large && this_added_read_inner <= last_added_read_edge ){
					//LRMprintf("OVERLAPPING_LAGRE %s : %d <= %d\n", iteration_context -> read_name, this_added_read_inner, last_added_read_edge);
					has_overlapping = 1;
				}else if( (!go_large) && this_added_read_inner >= last_added_read_edge ){
					//LRMprintf("OVERLAPPING_SMALL %s : %d >= %d\n", iteration_context -> read_name, this_added_read_inner, last_added_read_edge);
					has_overlapping = 1;
				}

				iteration_context -> chain_cov_start[idx] = this_read_pos;
				iteration_context -> chain_cov_end[idx] = this_read_end;
				iteration_context -> chain_chro_at_cov_start[idx] = this_chro_pos;

				if(go_large) iteration_context ->chain_tolarge_items++; else iteration_context ->chain_tosmall_items++;

				last_added_read_pos = this_read_pos;
				last_added_chro_pos = this_chro_pos;
				last_added_read_edge = go_large? this_read_end : this_read_pos ;
			}

			//assert(next_test_subread_no >= 0);
			test_subread_no = next_test_subread_no;
		}
	}

	int ii;
	for(ii = 0; ii < iteration_context -> chain_tosmall_items/2 ; ii++){
		unsigned int tv;

		tv = iteration_context -> chain_cov_start[ii];
		iteration_context -> chain_cov_start[ii] = iteration_context -> chain_cov_start[iteration_context -> chain_tosmall_items - 1 - ii];
		iteration_context -> chain_cov_start[iteration_context -> chain_tosmall_items - 1 - ii] = tv;

		tv = iteration_context -> chain_cov_end[ii];
		iteration_context -> chain_cov_end[ii] = iteration_context -> chain_cov_end[iteration_context -> chain_tosmall_items - 1 - ii];
		iteration_context -> chain_cov_end[iteration_context -> chain_tosmall_items - 1 - ii] = tv;

		tv = iteration_context -> chain_chro_at_cov_start[ii];
		iteration_context -> chain_chro_at_cov_start[ii] = iteration_context -> chain_chro_at_cov_start[iteration_context -> chain_tosmall_items - 1 - ii];
		iteration_context -> chain_chro_at_cov_start[iteration_context -> chain_tosmall_items - 1 - ii] = tv;
	}
	iteration_context -> chain_total_items = iteration_context -> chain_tosmall_items + iteration_context -> chain_tolarge_items;

	if(has_overlapping) LRMfix_extension_overlapping(context, thread_context, iteration_context, window_no);

	unsigned int this_gap_start_corr = iteration_context -> chain_chro_at_cov_start[0];
	unsigned int this_gap_end_corr = iteration_context -> chain_chro_at_cov_start[iteration_context -> chain_total_items -1] + iteration_context -> chain_cov_end[iteration_context -> chain_total_items -1] - iteration_context -> chain_cov_start[ iteration_context -> chain_total_items -1 ];

	if(iteration_context -> chain_total_items>0){
		LRMArrayListPush(iteration_context -> chain_used_gaps, NULL+this_gap_start_corr);
		LRMArrayListPush(iteration_context -> chain_used_gaps, NULL+this_gap_end_corr);
		//LRMprintf("USED_CHAIN_RANGE ITEMS=%d  : %u ~ %u\n", iteration_context -> chain_total_items, this_gap_start_corr , this_gap_end_corr);
	}else{
		//LRMprintf("NO_CHAIN_IS_BUILT\n");
	}
}

void LRMfill_gaps_addNM(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, int window_no, int subread_no){
	int middle_delta = subread_no?1:1;
	thread_context -> dynamic_programming_indel_movement_start += sprintf( thread_context -> dynamic_programming_indel_movement_buf + thread_context -> dynamic_programming_indel_movement_start, "%dM/",  iteration_context -> chain_cov_end[subread_no] - iteration_context -> chain_cov_start[subread_no] - middle_delta);
	//LRMprintf("NORMAL MS: %dM\n",  iteration_context -> chain_cov_end[subread_no] - iteration_context -> chain_cov_start[subread_no] - middle_delta);
}

void LRMsave_mapping_result(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, int window_no){
	int flags = 4;
	//if(iteration_context -> total_matched_bases >= context -> min_matched_bases_in_alignment)
	char * chro_name;
	int chro_pos, map_quality = 10, mis_matched = 0, corr_report = 0;

	if( iteration_context -> chain_total_items > 0 ){
		flags = iteration_context -> sorted_window_is_negative_strand[window_no]?16:0;
		LRMlocate_gene_position(context, iteration_context -> final_mapping_location, &chro_name, &chro_pos);
		if(NULL != chro_name){
			corr_report = 1;
			thread_context -> mapped_reads++;
		}
	}

	if(!corr_report){
		strcpy(thread_context -> final_cigar_string, "*");
		chro_name = "*";
		chro_pos = -1;
		map_quality = 0;
	}
	//LRMprintf("Writing [%d] : %u => %s , %d\n", loc_ret, iteration_context -> final_mapping_location,  chro_name , chro_pos);
	LRMwrite_chunk_add_buffered_output(context, thread_context, iteration_context, flags, chro_name, chro_pos +1, map_quality, thread_context -> final_cigar_string, mis_matched);
}

void LRMfill_gaps_reduce_Cigar(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, int * PTRmapped_length, int window_no){
	//LRMtest_move_buff( context, thread_context, iteration_context, thread_context -> dynamic_programming_indel_movement_buf, strlen(thread_context -> dynamic_programming_indel_movement_buf), iteration_context -> read_length);
	int tmpi = -1;
	int ci, nch, repeat_i = 0, old_opt = 0, wcur=0, r_rebuilt_len=0, mapped_length = 0;
 	for(ci = 0; ; ci++){
		nch = thread_context -> dynamic_programming_indel_movement_buf[ci];
		if(!nch) break;
		if(nch == '/') continue;
		if(nch == '.') continue;
		if(nch == 'X') nch='M';

		if(isdigit(nch)){
			if(tmpi<0) tmpi = 0;
			tmpi = tmpi*10 + (nch-'0');
		}else{
			if(tmpi<0) tmpi = 1;
			if(old_opt != nch && repeat_i>0){
				wcur += sprintf( thread_context -> final_cigar_string + wcur, "%d%c", repeat_i, old_opt );
				if( old_opt == 'M' || old_opt == 'I' || old_opt == 'S' ) r_rebuilt_len += repeat_i;
				if( old_opt == 'M' ) mapped_length += repeat_i;
				repeat_i = 0;
			}
			repeat_i += tmpi; 
			tmpi = -1;
			old_opt = nch;
		}
	}
	if(repeat_i>0){
		if( old_opt == 'M' || old_opt == 'I' || old_opt == 'S' ) r_rebuilt_len += repeat_i;
		if( old_opt == 'M' ) mapped_length += repeat_i;
		sprintf( thread_context -> final_cigar_string + wcur, "%d%c", repeat_i, old_opt );
	}
	//LRMprintf("Rebuild Rlen of %s = %d\n", iteration_context -> read_name, r_rebuilt_len);
	if(r_rebuilt_len != iteration_context -> read_length)
		LRMprintf("WRONG_REBUILD : %s : %d != %d ; %s\n", iteration_context -> read_name, r_rebuilt_len , iteration_context -> read_length, thread_context -> dynamic_programming_indel_movement_buf);
	*PTRmapped_length = mapped_length;
}

unsigned int LRMfill_gaps_find_final_mapping_loc(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, int window_no){
	int ii;
	int rollback = 0;
	for(ii = thread_context -> dynamic_programming_indel_movement_start ; ii >=0; ii--){
		int nch = thread_context -> dynamic_programming_indel_movement_buf[ii];
		if( nch == 'I' || nch == '/' || nch == '.' ) continue;
		if( nch == 'S') break;
		rollback ++;
	}
	return iteration_context -> chain_chro_at_cov_start[0] - rollback;
}
void LRMfill_gaps(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, int window_no){
	int ii;	
	int total_read_len = 0;

	//LRMprintf("READ_REVERSING: %d == %d ?\n", iteration_context -> is_reversed ,  iteration_context -> sorted_window_is_negative_strand[window_no]);
	if(iteration_context -> is_reversed != iteration_context -> sorted_window_is_negative_strand[window_no] ){
		LRMreverse_read(iteration_context -> read_text, iteration_context -> read_length);
		iteration_context -> is_reversed=!iteration_context -> is_reversed;
	}

	int first_mapped_subread_start = iteration_context -> chain_cov_start[0];
	int new_moves = LRMdynamic_to_ends(context, thread_context, iteration_context, first_mapped_subread_start, iteration_context -> chain_chro_at_cov_start[0], 0);
	//LRMprintf("MOVES_START %s: moves=%d ; read_len = %d ; %s\n", iteration_context -> read_name, new_moves, first_mapped_subread_start, thread_context ->  dynamic_programming_indel_movement_buf);
	if(new_moves<0){
		LRMprintf("MINUS MOVES : %d\n", new_moves);
		new_moves = 0;
	}
	total_read_len += first_mapped_subread_start;
	thread_context -> dynamic_programming_indel_movement_start += new_moves;
	iteration_context -> final_mapping_location = LRMfill_gaps_find_final_mapping_loc(context, thread_context, iteration_context, window_no);

	LRMfill_gaps_addNM(context, thread_context, iteration_context, window_no, 0);
	total_read_len += iteration_context -> chain_cov_end[0] - iteration_context -> chain_cov_start[0];

	int last_subread_end = iteration_context -> chain_cov_end[0];
	unsigned int last_subread_chro_end = iteration_context -> chain_chro_at_cov_start[0] + (last_subread_end - first_mapped_subread_start);

	for(ii = 1; ii < iteration_context -> chain_total_items; ii++){
		unsigned int this_subread_chro_start = iteration_context -> chain_chro_at_cov_start[ii];
		int this_subread_start = iteration_context -> chain_cov_start[ii];
		int read_minus_chro_delta = (this_subread_start - last_subread_end) - ( this_subread_chro_start - last_subread_chro_end );
		//LRMprintf("read_minus_chro_delta== %d = (%d - %d) - (%d - %d)\n", read_minus_chro_delta, this_subread_start, last_subread_end, this_subread_chro_start, last_subread_chro_end);
		if(this_subread_start <=last_subread_end -1) LRMprintf("Error: Reversed order of %s : %d < %d\n", iteration_context -> read_name, this_subread_start, last_subread_end -1 );
		assert( this_subread_start > last_subread_end -1 );
		if( this_subread_start - last_subread_end >=  LRMDYNAMIC_MAXIMUM_GAP_LENGTH ){
			int indel_after_M = -read_minus_chro_delta;

			int gap_read_M = this_subread_start - last_subread_end +1 - max( 0, -indel_after_M );
			int gap_read_M_L = gap_read_M/2;
			int gap_read_M_R = gap_read_M - gap_read_M_L;
			int indel_move = indel_after_M<0?'I':'D';
			thread_context -> dynamic_programming_indel_movement_start += sprintf( thread_context -> dynamic_programming_indel_movement_buf + thread_context -> dynamic_programming_indel_movement_start, "%dM%d%c%dM/", gap_read_M_L,abs(indel_after_M), indel_move, gap_read_M_R);
			LRMprintf("LONG GAP %s : %d\n", iteration_context -> read_name, this_subread_start - last_subread_end );
			LRMprintf("LONG GAP CIGAR : %dM%d%c\n", gap_read_M,abs(indel_after_M), indel_move);
		}
		else{
			new_moves = LRMdynamic_in_middle(context, thread_context, iteration_context, last_subread_end -1, this_subread_start, last_subread_chro_end - 1, read_minus_chro_delta);

			if(new_moves<0){
				LRMprintf("MINUS MOVES : %d\n", new_moves);
				new_moves = 0;
			}
			thread_context -> dynamic_programming_indel_movement_start += new_moves;
		}

		total_read_len += iteration_context -> chain_cov_end[ii] - last_subread_end ;
		LRMfill_gaps_addNM(context, thread_context, iteration_context, window_no, ii);
		last_subread_end = iteration_context -> chain_cov_end[ii];
		last_subread_chro_end = this_subread_chro_start + (last_subread_end - this_subread_start);
	}

	unsigned int last_end_chro_pos =  iteration_context -> chain_chro_at_cov_start[ii-1]; 
	last_end_chro_pos += ( iteration_context -> chain_cov_end[ii-1] - iteration_context -> chain_cov_start[ii-1] );
	new_moves = LRMdynamic_to_ends(context, thread_context, iteration_context, last_subread_end -1, last_end_chro_pos -1, 1 );

	if(new_moves<0){
		LRMprintf("MINUS MOVES : %d\n", new_moves);
		new_moves = 0;
	}

	total_read_len += iteration_context -> read_length - last_subread_end;
	if( total_read_len!=iteration_context -> read_length )LRMprintf("WRONG LENGTH %s: %d (mapping) != %d (fastq length) \n", iteration_context -> read_name, total_read_len,  iteration_context -> read_length);

	thread_context -> dynamic_programming_indel_movement_start += new_moves;
	int mapped_length = 0;
	LRMfill_gaps_reduce_Cigar(context, thread_context, iteration_context, &mapped_length, window_no);
	if(mapped_length < context->min_matched_bases_in_alignment) iteration_context -> chain_total_items=0;
}

void LRMreset_iteration_context_before_read_one_alignment(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context){
	iteration_context -> total_matched_bases = 0;
	iteration_context -> chain_tolarge_items = 0;
	iteration_context -> chain_tosmall_items = 0;
	iteration_context -> chain_total_items = 0;
	memset(iteration_context -> sorted_window_total_votes, 0, sizeof(int) * LRMMAX_SORTED_WINDOW_NO_HARD_LIMIT);
	thread_context -> dynamic_programming_indel_movement_start = 1;
	thread_context -> dynamic_programming_indel_movement_buf[0]='/';
}

void LRMreset_iteration_context_before_read( LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context ){
	if(iteration_context -> chain_used_gaps) LRMArrayListDestroy(iteration_context -> chain_used_gaps);
	iteration_context -> chain_used_gaps = LRMArrayListCreate(40);
}

void LRMdo_dynamic_programming_read( LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context ){
	LRMreset_iteration_context_before_read(context , thread_context , iteration_context);
	LRMcopy_longvotes_to_itr(context, thread_context, iteration_context);
	LRMmerge_sort(iteration_context, iteration_context -> sorting_total_votes, LRM_longvote_location_compare, LRM_longvote_location_exchange, LRM_longvote_location_merge);

	for(iteration_context -> current_alignment_no = 0; iteration_context -> current_alignment_no < context -> multi_best_read_alignments; iteration_context -> current_alignment_no++){
		LRMreset_iteration_context_before_read_one_alignment(context , thread_context , iteration_context);
		LRMfind_top_windows(context, thread_context, iteration_context);
	
		int ww = 0;
		if(iteration_context -> sorted_window_total_votes [ww] >0){
			LRMbuild_chains(context, thread_context, iteration_context, ww);
	
			if(0&& LRMFIXLENstrcmp("R4", iteration_context->read_name) == 0 )
				LRMprint_longvote(context, thread_context, iteration_context);
	
			LRMfill_gaps(context, thread_context, iteration_context, ww);
		}else iteration_context -> chain_total_items = iteration_context -> chain_tosmall_items = iteration_context -> chain_tolarge_items = 0;
		LRMsave_mapping_result(context, thread_context, iteration_context, ww);
	}
}

int LRMchunk_read_iteration(LRMcontext_t * context, int thread_id, int task){
	LRMthread_context_t * thread_context = context -> thread_contexts+ thread_id;

	LRMread_iteration_context_t * iteration_context;
	iteration_context = malloc(sizeof(LRMread_iteration_context_t));
	memset(iteration_context, 0, sizeof(LRMread_iteration_context_t));

	while(1){
		int retv = LRMfetch_next_read(context, thread_context, &iteration_context-> read_length, iteration_context->read_name, iteration_context->read_text, iteration_context->qual_text, &iteration_context -> read_no_in_chunk);
		if(retv) break;

		LRMdo_one_voting_read(context, thread_context, iteration_context);
		LRMdo_dynamic_programming_read(context, thread_context, iteration_context);
		
		if(iteration_context -> read_no_in_chunk % 2000 == 0)
			LRMprintf("Processing %d-th read for task %d; used %.1f minutes\n", context -> all_processed_reads + iteration_context -> read_no_in_chunk, task, (LRMmiltime() - context -> start_running_time)/60);
	}
	if(iteration_context -> chain_used_gaps) LRMArrayListDestroy(iteration_context -> chain_used_gaps);
	iteration_context -> chain_used_gaps = NULL;
	free(iteration_context);
	return 0;
}

int LRMfinalise_chunk_reads(LRMcontext_t* context){
	context ->  all_processed_reads += context -> processed_reads_in_chunk;
	return 0;
}

int LRMFIXLENstrcmp(char * fixed_len, char * rname){
        int x=0;
        for(; fixed_len[x]; x++){
                if(rname[x]!=fixed_len[x]) return 1;
        }
        return 0;
}

