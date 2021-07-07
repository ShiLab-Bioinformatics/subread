#include <getopt.h>

#include "subread.h"
#include "gene-algorithms.h"
#include "input-files.h"
#include "input-blc.h"


typedef struct
{
	unsigned int selected_position;
	short result_flags;
	short read_length;
	// 4 bytes
	gene_vote_number_t selected_votes;
	gene_vote_number_t used_subreads_in_vote;
	char indels_in_confident_coverage;
	char is_fully_covered;
	gene_vote_number_t selected_indel_record [MAX_INDEL_SECTIONS*3 + 1];
	unsigned short confident_coverage_start;
	unsigned short confident_coverage_end;
} voting_location_t;

typedef struct{
	int thread_no;
	pthread_t thread;

	int event_space_capacity;
	int total_events;
	HashTable * event_entry_table;
	
	short ** dynamic_align_table;
	char ** dynamic_align_table_mask;
	chromosome_event_t * event_space_dynamic;

	topK_buffer_t topKbuff;
	short * final_reads_mismatches_array;
	short * final_counted_reads_array;
	gene_value_index_t * current_value_index;
} cellcounts_align_thread_t;


typedef struct{
	int total_threads;
	cellcounts_align_thread_t * all_thread_contexts;
	int reads_per_chunk;
	int max_voting_simples;
	int max_voting_locations;
	int max_best_alignments;
	int max_indel_length;
	int max_top_vote_simples;
	int max_vote_number_cutoff;
	int min_votes_per_mapped_read;
	int total_subreads_per_read;

	int processed_reads_in_chunk;
	int running_processed_reads_in_chunk;
	int is_final_voting_run;
	int output_binfiles_are_full;
	gene_inputfile_position_t current_circle_start_position, current_circle_end_position;
	int last_written_fragment_number;
	voting_location_t * voting_location_table;

	int total_events;
	int event_space_capacity;
	chromosome_event_t * event_space_dynamic;
	chromosome_event_t * event_space;
	HashTable * event_entry_table;

	char index_prefix[MAX_FILE_NAME_LENGTH];
	char output_prefix[MAX_FILE_NAME_LENGTH];
	char input_dataset_name[MAX_FILE_NAME_LENGTH * MAX_SCRNA_FASTQ_FILES * 3];
	int input_mode;

	int total_index_blocks;
	int current_index_block_number;
	gehash_t * current_index;
	gene_value_index_t * current_value_index;
	gene_value_index_t all_value_indexes[100];
	gene_input_t input_dataset;
	subread_lock_t input_dataset_lock;
	subread_lock_t output_lock;
	
	char features_annotation_file[MAX_FILE_NAME_LENGTH];
	char features_annotation_alias_file[MAX_FILE_NAME_LENGTH];
	int  features_annotation_file_type;
	char features_annotation_gene_id_column[MAX_READ_NAME_LEN];
	char features_annotation_feature_type[MAX_READ_NAME_LEN];
	gene_offset_t chromosome_table;

	unsigned char 	* features_sorted_strand;
	srInt_64 	* features_sorted_start, * features_sorted_stop;
	int 		* features_sorted_geneid;
	char 		** features_sorted_chr;
	HashTable 	* sam_chro_to_anno_chr_alias;
	
	char		* exonic_region_bitmap;
} cellcounts_global_t;

typedef struct{
	int thread_no;
} cellcounts_final_thread_t;

int features_load_one_line(char * gene_name, char * transcript_name, char * chro_name, unsigned int start, unsigned int end, int is_negative_strand, void * context){
	cellcounts_global_t * cct_context = context;
	if(cct_context -> sam_chro_to_anno_chr_alias){
		char * sam_chro = get_sam_chro_name_from_alias(cct_context -> sam_chro_to_anno_chr_alias, chro_name);
		if(sam_chro!=NULL) chro_name = sam_chro;
	}

	char tmp_chro_name[MAX_CHROMOSOME_NAME_LEN];

	int access_n = HashTableGet(cct_context -> chromosome_table.read_name_to_index, chro_name ) - NULL;
	if(access_n < 1){
		if(chro_name[0]=='c' && chro_name[1]=='h' && chro_name[2]=='r'){
			chro_name += 3;
		}else{
			strcpy(tmp_chro_name, "chr");
			strcat(tmp_chro_name, chro_name);
			chro_name = tmp_chro_name;
		}
	}

	return 0;
}

static struct option cellCounts_long_options[]={
	{"dataset", required_argument ,0,0},
	{"index", required_argument ,0,0},
	{"inputMode", required_argument ,0,0},
	{"output", required_argument ,0,0},
	{"threads", required_argument ,0,0},

	{"annotation", required_argument ,0,0},
	{"isGTFannotation", no_argument ,0,0},
	{"geneIdColumn", required_argument ,0,0},
	{"annotationType", required_argument ,0,0},
	{"annotationChroAlias", required_argument ,0,0},

	{0,0,0,0}
};


int cellCounts_args_context(cellcounts_global_t * cct_context, int argc, char** argv){
	int c , option_index=0;
	optind = 0;
	opterr = 1;
	optopt = 63;

	cct_context -> input_mode = GENE_INPUT_BCL;
	cct_context -> total_threads = 10;
	cct_context -> features_annotation_file_type = FILE_TYPE_RSUBREAD;
	cct_context -> reads_per_chunk = 30000000;
	cct_context -> max_best_alignments = 1;
	cct_context -> max_voting_simples = 5;
	cct_context -> max_voting_locations = 3;
	cct_context -> max_indel_length = 5;
	cct_context -> max_top_vote_simples = 3;
	cct_context -> max_vote_number_cutoff = 2;
	cct_context -> min_votes_per_mapped_read = 3;
	cct_context -> total_subreads_per_read = 10;

	while (1){
		c = getopt_long(argc, argv, "", cellCounts_long_options, &option_index);
		if(c<0 || c==255)break;

		if(strcmp("dataset", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> input_dataset_name, optarg, MAX_FILE_NAME_LENGTH * MAX_SCRNA_FASTQ_FILES * 3 -1);
		}
		if(strcmp("index", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> index_prefix, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("inputMode", cellCounts_long_options[option_index].name)==0){
			if(strcmp("FASTQ", optarg)==0) cct_context -> input_mode = GENE_INPUT_SCRNA_FASTQ;
		}
		if(strcmp("output", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> output_prefix, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("threads", cellCounts_long_options[option_index].name)==0){
			cct_context -> total_threads = min(64, max(1, atoi(optarg)));
		}
		if(strcmp("annotation", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> features_annotation_file, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("annotationChroAlias", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> features_annotation_alias_file, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("annotationType", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> features_annotation_feature_type, optarg, MAX_READ_NAME_LEN-1);
		}
		if(strcmp("geneIdColumn", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> features_annotation_gene_id_column, optarg, MAX_READ_NAME_LEN-1);
		}
		if(strcmp("isGTFannotation", cellCounts_long_options[option_index].name)==0){
			cct_context -> features_annotation_file_type = FILE_TYPE_GTF;
		}
	}
	return 0;
}


int determine_total_index_blocks(cellcounts_global_t * cct_context){
	char tmp_fname[MAX_FILE_NAME_LENGTH+ 30];
	cct_context-> total_index_blocks = 0;
	while(1){
		sprintf(tmp_fname, "%s.%02d.b.tab", cct_context->index_prefix, cct_context->total_index_blocks);
		if(!does_file_exist(tmp_fname))break;
		cct_context->total_index_blocks ++;
	}
	return 0;
}

#define EXONIC_REGION_RESOLUTION 16

int cellCounts_load_context(cellcounts_global_t * cct_context){
	int rv = 0, block_no;

	if(cct_context -> input_mode == GENE_INPUT_BCL)
		rv = rv || geinput_open_bcl(cct_context -> input_dataset_name , & cct_context -> input_dataset , cct_context -> reads_per_chunk, cct_context -> total_threads);
	else if(cct_context -> input_mode == GENE_INPUT_SCRNA_FASTQ)
		rv = rv || geinput_open_scRNA_fqs(cct_context -> input_dataset_name , & cct_context -> input_dataset , cct_context -> reads_per_chunk, cct_context -> total_threads );
	else 	rv = -1;

	rv = rv || load_offsets(& cct_context -> chromosome_table, cct_context -> index_prefix);
	rv = rv || determine_total_index_blocks(cct_context);

	int bitmap_size = (4096 / EXONIC_REGION_RESOLUTION / 8)*1024*1024;
	rv = rv || ((cct_context -> exonic_region_bitmap = calloc(bitmap_size, 1))==NULL);
	int voting_location_table_items = cct_context -> reads_per_chunk * max(cct_context -> max_voting_locations, cct_context -> max_best_alignments);
	rv = rv || ((cct_context->voting_location_table = calloc(sizeof(voting_location_t), voting_location_table_items))==NULL);

	if(!rv)for(block_no = 0; block_no< cct_context->total_index_blocks; block_no++) {
		char tmp_fname[MAX_FILE_NAME_LENGTH+ 30];
		sprintf(tmp_fname, "%s.%02d.b.array", cct_context ->index_prefix, block_no);
		rv = rv || gvindex_load(&cct_context -> all_value_indexes[block_no], tmp_fname);
	}

	return rv;
}

int cellCounts_destroy_context(cellcounts_global_t * cct_context){
	geinput_close(&cct_context -> input_dataset);
	destroy_offsets(&cct_context->chromosome_table);
	free(cct_context->voting_location_table);
	free(cct_context->exonic_region_bitmap);
	return 0;
}

void cellCounts_locate_read_files(cellcounts_global_t * cct_context, int type) {
	// The BCL input module uses its own chunking algorithm.
	if(cct_context -> input_mode == GENE_INPUT_BCL) return;

	if(type==SEEK_SET) geinput_tell(&cct_context -> input_dataset, &cct_context -> current_circle_start_position);
	else geinput_tell(&cct_context -> input_dataset, &cct_context -> current_circle_end_position);

	return;
}

void cellCounts_rewind_read_files(cellcounts_global_t * cct_context, int type)
{
	if(type==SEEK_SET) geinput_seek(&cct_context -> input_dataset, &cct_context -> current_circle_start_position);
	else geinput_seek(&cct_context -> input_dataset, &cct_context -> current_circle_end_position);
}

void cellCounts_go_chunk_start(cellcounts_global_t * cct_context){
	if(cct_context -> input_mode == GENE_INPUT_BCL)
		cacheBCL_go_chunk_start(&cct_context -> input_dataset.bcl_input);
	else
		cellCounts_rewind_read_files(cct_context, SEEK_SET);
	cct_context -> running_processed_reads_in_chunk=0;
}

void cellCounts_go_chunk_nextchunk(cellcounts_global_t * cct_context){
	if(cct_context -> input_mode == GENE_INPUT_BCL)
		cacheBCL_go_chunk_end(&cct_context -> input_dataset.bcl_input);
	else
		cellCounts_rewind_read_files(cct_context, SEEK_END);
	cct_context -> running_processed_reads_in_chunk=0;
}

void cellCounts_clean_context_after_chunk(cellcounts_global_t * cct_context) {
	cct_context -> running_processed_reads_in_chunk = 0;
	cct_context -> processed_reads_in_chunk = 0;

	int event_id;
	//memset(context -> big_margin_record  , 0 , sizeof(*context -> big_margin_record) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.big_margin_record_size);
	for(event_id = 0; event_id < cct_context -> total_events; event_id++){
		chromosome_event_t * event_body = cct_context -> event_space + event_id;
		event_body -> critical_read_id = 0xffffffffffffffffllu;
	}
}

void cellCounts_init_topKbuff(cellcounts_global_t * cct_context, int thread_no){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	topK_buffer_t * topKbuff = &thread_context -> topKbuff;
	topKbuff -> vote_simple_1_buffer = malloc(cct_context -> max_voting_simples * sizeof(simple_mapping_t));
	topKbuff -> alignment_tmp_r1 = malloc(sizeof(voting_location_t) * cct_context -> max_voting_locations);
}

void cellCounts_free_topKbuff(cellcounts_global_t * cct_context, int thread_no){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	topK_buffer_t * topKbuff = &thread_context -> topKbuff;
	free(topKbuff -> vote_simple_1_buffer);
	free(topKbuff -> alignment_tmp_r1);
}

int cellCounts_anti_supporting_read_scan(cellcounts_global_t * cct_context){
	SUBREADprintf("SHOULD COPY anti_supporting_read_scan from align()\n");
	return 0;
}

int cellCounts_remove_neighbour(cellcounts_global_t * cct_context){
	SUBREADprintf("SHOULD COPY remove_neighbour from align()\n");
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
}cellCounts_output_tmp_t;


typedef struct{
	long long int *PE_distance;
	char * out_cigar_buffer[CIGAR_PERFECT_SECTIONS];
	cellCounts_output_tmp_t *r1;
	cellCounts_output_tmp_t ** out_pairs;
	mapping_result_t ** out_raws;
} cellCounts_output_context_t;

voting_location_t * cellCounts_global_retrieve_alignment_ptr(cellcounts_global_t * cct_context, subread_read_number_t pair_number, int best_voting_loc_id){
	int item_idx = pair_number * cct_context -> max_voting_locations + best_voting_loc_id; 
	voting_location_t * ret = cct_context -> voting_location_table+item_idx;
	return ret;
}

void cellCounts_init_output_context(cellcounts_global_t * cct_context, cellCounts_output_context_t * out_context) {
	int xk1;
	memset(out_context, 0, sizeof(cellCounts_output_context_t));

	out_context -> r1 = malloc(sizeof(cellCounts_output_tmp_t));
	for(xk1=0;xk1<CIGAR_PERFECT_SECTIONS;xk1++)
		out_context -> out_cigar_buffer[xk1] = malloc(60);

	out_context -> out_pairs = malloc(sizeof( cellCounts_output_context_t *) * cct_context -> max_voting_locations);
	out_context -> out_raws = malloc(sizeof( voting_location_t **) * cct_context -> max_voting_locations);
}

void cellCounts_destroy_output_context(cellcounts_global_t * cct_context, cellCounts_output_context_t * out_context) {
	int xk1;
	for(xk1=0;xk1<CIGAR_PERFECT_SECTIONS;xk1++)
		free(out_context -> out_cigar_buffer[xk1]);

	free(out_context -> out_raws);
	free(out_context -> out_pairs);
	free(out_context -> r1 );
}

unsigned int cellCounts_calc_end_pos(unsigned int p, char * cigar, unsigned int * all_skipped_len, int * is_exonic_regions, cellcounts_global_t * cct_context){
	unsigned int cursor = p, tmpi=0;
	int nch, cigar_cursor;
	for(cigar_cursor = 0; 0!=(nch = cigar[cigar_cursor]); cigar_cursor++){
		if(isdigit(nch)){
			tmpi = tmpi * 10 + (nch - '0');
		}else{
			if((nch == 'S' && cursor == p) || nch == 'M' || nch == 'N' || nch == 'D'){
				if(nch == 'M' && cellCounts_is_pos_in_annotated_exon_regions(cct_context, cursor + tmpi/2) == 0) ( * is_exonic_regions) = 0;
				cursor += tmpi;
				if(nch == 'N' || nch == 'D')(*all_skipped_len) += tmpi;
			}
			tmpi = 0;
		}
	}
	return cursor;

}

void cellCounts_write_realignments_for_fragment(cellcounts_global_t * cct_context, int thread_no, cellCounts_output_context_t * out_context, unsigned int read_number, realignment_result_t * res, char * read_name, char * read_text, char * qual_text, int rlen, int multi_mapping_number, int this_multi_mapping_i){
}

unsigned int cellCounts_explain_read(cellcounts_global_t * cct_context, int thread_no, realignment_result_t * realigns, subread_read_number_t pair_number,int read_len, char * read_name , char *read_text, char *qual,  int voting_loc_no, int is_negative_strand){
	return 0;
}

int cellCounts_do_iteration_two(cellcounts_global_t * cct_context, int thread_no){
	subread_read_number_t current_read_number=0;
	char read_text[MAX_READ_LENGTH+1], qual_text[MAX_READ_LENGTH+1];

	char raw_read_text[MAX_READ_LENGTH+1], raw_qual_text[MAX_READ_LENGTH+1];
	char read_name[MAX_READ_NAME_LEN+1];
	char * repeated_buffer_cigars[MAX_ALIGNMENT_PER_ANCHOR *  2 * cct_context -> max_best_alignments];
	unsigned int repeated_buffer_pos[MAX_ALIGNMENT_PER_ANCHOR *  2 * cct_context -> max_best_alignments];
	int repeated_count, read_len=0;

	int * final_MATCH_buffer, *final_MISMATCH_buffer, *final_PENALTY_buffer;
	int * final_realignment_index;
	unsigned int *final_realignment_number;
	unsigned long long * final_SCORE_buffer;

	realignment_result_t * final_realignments;
	cellCounts_output_context_t out_context;

	cellCounts_init_output_context(cct_context, &out_context);	

	for(repeated_count = 0;repeated_count < MAX_ALIGNMENT_PER_ANCHOR  *  2 * cct_context -> max_best_alignments ; repeated_count ++)
		repeated_buffer_cigars[repeated_count] = malloc(2*CORE_MAX_CIGAR_STR_LEN);

	final_MATCH_buffer = malloc(sizeof(int) * 2 * cct_context -> max_voting_locations * MAX_ALIGNMENT_PER_ANCHOR);
	final_MISMATCH_buffer = malloc(sizeof(int) * 2 * cct_context -> max_voting_locations * MAX_ALIGNMENT_PER_ANCHOR);
	final_PENALTY_buffer = malloc(sizeof(int) * 2 * cct_context -> max_voting_locations * MAX_ALIGNMENT_PER_ANCHOR);
	final_realignment_index = malloc(sizeof(int) * cct_context -> max_voting_locations * MAX_ALIGNMENT_PER_ANCHOR);

	final_SCORE_buffer = malloc(sizeof(long long) * cct_context -> max_voting_locations * cct_context -> max_voting_locations * MAX_ALIGNMENT_PER_ANCHOR*MAX_ALIGNMENT_PER_ANCHOR);

	final_realignments = malloc(sizeof(realignment_result_t) * cct_context -> max_voting_locations * MAX_ALIGNMENT_PER_ANCHOR);
	final_realignment_number = malloc(sizeof(int) * cct_context -> max_voting_locations );

	mapping_result_t * align_result_buffer = malloc(sizeof(mapping_result_t) * cct_context -> max_voting_locations);

	while(1) {
		int max_votes;

		cellCounts_fetch_next_read_pair(cct_context, thread_no,  &read_len, read_name, read_text, qual_text, &current_read_number);
		// if no more reads
		if(current_read_number < 0) break;
		strcpy(raw_read_text, read_text);
		strcpy(raw_qual_text, qual_text);

		max_votes = cellCounts_global_retrieve_alignment_ptr(cct_context, current_read_number, 0)->selected_votes;

		int voting_loc_no=0;
		int step2_locations = 0;
		repeated_count=0;

		int is_reversed_already = 0, realignment_i;
		int candidate_locations = 0;

		for(voting_loc_no = 0; voting_loc_no < cct_context -> max_voting_locations; voting_loc_no++) {
			voting_location_t *current_result = cellCounts_global_retrieve_alignment_ptr(cct_context, current_read_number, voting_loc_no);

			if(max_votes < cct_context -> min_votes_per_mapped_read) {
				current_result -> selected_votes = 0;
				continue;
			}
			if(cellCounts_locate_current_value_index(cct_context, thread_no, current_result, read_len)) {
				current_result -> selected_votes = 0;
				continue;
			}

			int is_negative_strand = (current_result  -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;

			if(is_negative_strand + is_reversed_already == 1) {
				reverse_read(read_text, read_len, GENE_SPACE_BASE);
				reverse_quality(qual_text, read_len);
				is_reversed_already = !is_reversed_already;
			}

			current_result -> result_flags &= ~CORE_IS_FULLY_EXPLAINED;
			step2_locations = voting_loc_no + 1;

			unsigned int final_alignments = cellCounts_explain_read(cct_context, thread_no , final_realignments + voting_loc_no * MAX_ALIGNMENT_PER_ANCHOR,
						   current_read_number, read_len, read_name, read_text, qual_text, voting_loc_no, is_negative_strand);

			final_realignment_number[ voting_loc_no ] = final_alignments;

			for(realignment_i = 0 ; realignment_i < final_alignments ; realignment_i ++){
				if(candidate_locations >= cct_context -> max_voting_locations * MAX_ALIGNMENT_PER_ANCHOR) break;

				int realign_index = voting_loc_no* MAX_ALIGNMENT_PER_ANCHOR + realignment_i;
				int final_MATCH = final_realignments[realign_index].final_matched_bases;
				int final_PENALTY = final_realignments[realign_index].final_penalty;

				if((current_result -> result_flags & CORE_IS_FULLY_EXPLAINED) && final_MATCH>0) {
					int final_MISMATCH = final_realignments[realign_index].final_mismatched_bases;

					final_PENALTY_buffer[candidate_locations] = final_PENALTY;
					final_MATCH_buffer[candidate_locations] = final_MATCH;
					final_MISMATCH_buffer[candidate_locations] = final_MISMATCH;
					final_realignment_index[candidate_locations] = realign_index;
					candidate_locations ++;
				}
			}
		}

		int output_cursor = 0;
		int highest_score_occurence =0;
		if(candidate_locations > 0){
			int read_record_i;
			unsigned long long int best_score_highest = 0;
			unsigned long long int scores_array [cct_context -> max_voting_locations * MAX_ALIGNMENT_PER_ANCHOR];

			for(read_record_i = 0; read_record_i < candidate_locations; read_record_i++){
				realignment_result_t * final_realignment_result = final_realignments + final_realignment_index[read_record_i];

				unsigned int this_MATCH = final_MATCH_buffer[read_record_i];
				unsigned int this_MISMATCH = final_MISMATCH_buffer[read_record_i];
				unsigned int this_PENALTY = final_PENALTY_buffer[read_record_i];
				unsigned long long this_SCORE;

				unsigned int skip = 0; int is_exonic_regions = 1;
				cellCounts_calc_end_pos(final_realignment_result -> first_base_position, final_realignment_result -> cigar_string, &skip, &is_exonic_regions, cct_context);
				this_SCORE =((100000llu * (10000 - this_MISMATCH + 2*is_exonic_regions) + this_MATCH)*50llu - this_PENALTY)*20llu+ final_realignment_result -> known_junction_supp;

				best_score_highest = max(best_score_highest, this_SCORE);
				scores_array[read_record_i] = this_SCORE;
			}

			for(read_record_i = 0; read_record_i <  candidate_locations ; read_record_i++){
				realignment_result_t * final_realignment_result = final_realignments + final_realignment_index[read_record_i];
				if( scores_array[read_record_i] >= best_score_highest && (final_realignment_result -> realign_flags & CORE_TOO_MANY_MISMATCHES)==0) {
					int is_repeated = 0;
					is_repeated = cellCounts_add_repeated_buffer(cct_context, repeated_buffer_pos, repeated_buffer_cigars, &repeated_count, final_realignment_result);
					if(is_repeated)
						scores_array[read_record_i] = 0;
					else
						highest_score_occurence ++;
				}
			}

			if(highest_score_occurence<2 || cct_context -> max_best_alignments){
				int is_break_even = 0;
				if(highest_score_occurence>1)	is_break_even = 1;

				highest_score_occurence = min(highest_score_occurence, cct_context -> max_best_alignments);
				for(read_record_i = 0; read_record_i <  candidate_locations ; read_record_i++){
					realignment_result_t * final_realignment_result = final_realignments + final_realignment_index[read_record_i];

					if( scores_array[read_record_i] >= best_score_highest && (final_realignment_result -> realign_flags & CORE_TOO_MANY_MISMATCHES)==0
					  && output_cursor < cct_context -> max_best_alignments){
						strcpy(read_text, raw_read_text);
						strcpy(qual_text, raw_qual_text);

						if(is_break_even) final_realignment_result -> realign_flags |= CORE_IS_BREAKEVEN; 
						final_realignment_result -> MAPQ_adjustment = final_MISMATCH_buffer[read_record_i] + step2_locations;

						cellCounts_write_realignments_for_fragment(cct_context, thread_no, &out_context , current_read_number, final_realignment_result, read_name, read_text, qual_text, read_len, highest_score_occurence, output_cursor);
						output_cursor ++;
					}
				}

				assert(output_cursor == highest_score_occurence);
			}
		}

		if(output_cursor<1) {
			strcpy(read_text, raw_read_text);
			strcpy(qual_text, raw_qual_text);
			cellCounts_write_realignments_for_fragment(cct_context, thread_no, &out_context, current_read_number, NULL, read_name, read_text, qual_text, read_len, 0, 0);
		}
		if(current_read_number % 50000 == 0) SUBREADprintf("realign step: %d\n", current_read_number);
	}

	free(final_realignments);
	free(final_realignment_number);
	free(final_MATCH_buffer);
	free(final_MISMATCH_buffer);
	free(final_PENALTY_buffer);
	free(final_realignment_index);

	free(final_SCORE_buffer);
	free(align_result_buffer);

	for(repeated_count = 0;repeated_count < MAX_ALIGNMENT_PER_ANCHOR * cct_context -> max_best_alignments ; repeated_count ++ )
		free(repeated_buffer_cigars[repeated_count]);

	cellCounts_destroy_output_context(cct_context, &out_context);	
	return 0;
}

void * cellCounts_run_in_thread(void * params){
	void ** parameters = (void **)params;
	cellcounts_global_t * cct_context = (cellcounts_global_t * ) parameters[0];
	int thread_no = parameters[1]-NULL;
	int task = parameters[2]-NULL;
	int * ret_value_pointer = (int *)parameters[3];
	free(parameters);

	switch(task)
	{
		case STEP_VOTING:
			*ret_value_pointer = cellCounts_do_voting(cct_context, thread_no);
		break;
		case STEP_ITERATION_TWO:
			*ret_value_pointer = cellCounts_do_iteration_two(cct_context, thread_no);
		break;

	}
}

#define MAX_EVENT_NUMBER_INIT 200000

int cellCounts_prepare_context_for_align(cellcounts_global_t * cct_context, int thread_no, int task) {
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	if(task == STEP_VOTING) {
		thread_context -> event_entry_table = HashTableCreate(1599997);
		HashTableSetKeyComparisonFunction(thread_context->event_entry_table, localPointerCmp_forEventEntry);
		HashTableSetHashFunction(thread_context->event_entry_table, localPointerHashFunction_forEventEntry);

		thread_context -> total_events = 0;
		thread_context -> event_space_capacity = MAX_EVENT_NUMBER_INIT;
		thread_context -> event_space_dynamic = malloc(sizeof(chromosome_event_t)*thread_context -> event_space_capacity);

		thread_context -> dynamic_align_table = malloc(sizeof(short*)*MAX_READ_LENGTH);
		thread_context -> dynamic_align_table_mask = malloc(sizeof(char *)*MAX_READ_LENGTH);

		int xk1;
		for(xk1=0;xk1<MAX_READ_LENGTH; xk1++){
			thread_context -> dynamic_align_table[xk1] = malloc(sizeof(short)*MAX_READ_LENGTH);
			thread_context -> dynamic_align_table_mask[xk1] = malloc(sizeof(char)*MAX_READ_LENGTH);
		}
	} else if(task == STEP_ITERATION_TWO) {
		thread_context -> event_entry_table = cct_context -> event_entry_table;
		thread_context -> total_events = cct_context -> total_events;
		thread_context -> event_space_dynamic = cct_context -> event_space_dynamic;
		thread_context -> final_reads_mismatches_array = (unsigned short*)malloc(sizeof(unsigned short)*cct_context -> total_events);
		thread_context -> final_counted_reads_array = (unsigned short*)malloc(sizeof(unsigned short)*cct_context -> total_events);

		memset(thread_context -> final_reads_mismatches_array , 0, sizeof(unsigned short)*cct_context -> total_events);
		memset(thread_context -> final_counted_reads_array , 0, sizeof(unsigned short)*cct_context -> total_events);

		subread_init_lock(&cct_context -> output_lock);
	}
	return 0;
}


int cellCounts_run_maybe_threads(cellcounts_global_t * cct_context, int task){
	int ret_value =0;

	if(task == STEP_ITERATION_TWO){
		cct_context -> last_written_fragment_number = -1;
	}

	int current_thread_no ;
	cellcounts_align_thread_t thread_contexts[64];
	int ret_values[64];

	memset(thread_contexts, 0, sizeof(cellcounts_align_thread_t)*64);
	cct_context -> all_thread_contexts = thread_contexts;

	for(current_thread_no = 0 ; current_thread_no < cct_context->total_threads ; current_thread_no ++) {
		thread_contexts[current_thread_no].thread_no = current_thread_no;
		cellCounts_prepare_context_for_align(cct_context, current_thread_no, task);

		if(STEP_VOTING == task) cellCounts_init_topKbuff(cct_context, current_thread_no);
		void ** thr_parameters = malloc(sizeof(void*)*4);
		thr_parameters[0] = cct_context;
		thr_parameters[1] = NULL+current_thread_no;
		thr_parameters[2] = NULL+task;
		thr_parameters[3] = ret_values + current_thread_no;

		pthread_create(&thread_contexts[current_thread_no].thread, NULL, cellCounts_run_in_thread, thr_parameters);
	}

	for(current_thread_no = 0 ; current_thread_no < cct_context->total_threads ; current_thread_no ++) {
		pthread_join(thread_contexts[current_thread_no].thread, NULL);

		if(STEP_VOTING == task) cellCounts_free_topKbuff(cct_context, current_thread_no);
		ret_value += *(ret_values + current_thread_no);
		if(ret_value)break;
	}

	// sort and merge events from all threads and the global event space.
	cellCounts_finalise_indel_and_junction_thread(cct_context, thread_contexts, task);

	return ret_value;
}

typedef struct {
	unsigned int thread_bodytable_number;
	short thread_no;
} concatinating_events_record_t;

int cellCounts_locate_current_value_index(cellcounts_global_t * cct_context, int thread_no, voting_location_t * result, int rlen) {
	int block_no;
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	if(cct_context->total_index_blocks == 1){
		unsigned index_begin = cct_context -> all_value_indexes[0].start_base_offset; 
		unsigned index_end = cct_context -> all_value_indexes[0].start_base_offset + cct_context -> all_value_indexes [block_no] . length;
		thread_context->current_value_index = cct_context -> all_value_indexes;
		if(index_begin <= result->selected_position && index_end > result->selected_position) return 0;
		return 1;
	}

	for(block_no=0; block_no< cct_context->total_index_blocks; block_no++) {
		unsigned index_begin = cct_context -> all_value_indexes [block_no] . start_base_offset; 
		unsigned index_end = cct_context -> all_value_indexes [block_no] . start_base_offset + cct_context -> all_value_indexes [block_no] . length;
		SUBREADprintf("   TRYING BLOCK %d of %u - %u\n", block_no, index_begin, index_end);
		if((block_no == 0 && result->selected_position >=  index_begin && result->selected_position < index_end - 1000000) ||
		   (block_no>0 && block_no<cct_context->total_index_blocks && result->selected_position >= index_begin+ 1000000 && result->selected_position < index_end - 1000000) ||
		   (block_no == cct_context->total_index_blocks-1 && result->selected_position >= index_begin + 1000000 && result->selected_position < index_end )) {
			thread_context->current_value_index = cct_context -> all_value_indexes + block_no;
			return 0;
		}
	}
	return 1;
}

int cellCounts_add_repeated_buffer(cellcounts_global_t * cct_context, unsigned int * repeated_buffer_position, char ** repeated_buffer_cigar, int * repeated_count, realignment_result_t * res1){
	int x1, is_repeated = 0;

	char * r1_cigar = res1 -> cigar_string;
	unsigned int r1_pos = res1 -> first_base_position;

	for(x1 = 0; x1 < (*repeated_count); x1 ++)
		if(repeated_buffer_position[x1] == r1_pos && strcmp(repeated_buffer_cigar[x1], r1_cigar) == 0) {
			is_repeated = 1;
			break;
		}

	if((!is_repeated) && (*repeated_count) < MAX_ALIGNMENT_PER_ANCHOR * cct_context -> max_best_alignments){
		repeated_buffer_position[*repeated_count] = r1_pos;
		strcpy(repeated_buffer_cigar[*repeated_count], r1_cigar);
		(*repeated_count) +=1;
	}

	return is_repeated;
}


int cellCounts_fetch_next_read_pair(cellcounts_global_t * cct_context, int thread_no,int *read_len_1, char * read_name_1, char * read_text_1, char * qual_text_1, subread_read_number_t * read_no_in_chunk) {
	int rl1=0;
	int is_second_R1;
	subread_read_number_t this_number = -1;
	gene_input_t * ginp1 = &cct_context -> input_dataset;

	geinput_preload_buffer(ginp1, &cct_context -> input_dataset_lock);
	subread_lock_occupy(&cct_context -> input_dataset_lock); 
	if(cct_context -> running_processed_reads_in_chunk < cct_context -> reads_per_chunk) {
		do{
			is_second_R1 = 0;
			rl1 = geinput_next_read_trim(ginp1, read_name_1, read_text_1 , qual_text_1, 0, 0, &is_second_R1);
			if(rl1 <= 0) break;
		}while(is_second_R1) ;

		if(rl1 >0){
			this_number = cct_context -> running_processed_reads_in_chunk;
			cct_context -> running_processed_reads_in_chunk ++;
		}
	}
	subread_lock_release(&cct_context -> input_dataset_lock); 

	if(rl1>0 && this_number>=0) {
		*read_no_in_chunk = this_number;
		*read_len_1 = rl1;
		return 0;
	} else {
		*read_no_in_chunk = -1;
		return 1;
	}
}

int cellCounts_has_better_mapping(cellcounts_global_t * cct_context, int thread_no, subread_read_number_t current_read_number, int this_aln_id){
	int better_voting_loc_id;

	voting_location_t * this_r = cellCounts_global_retrieve_alignment_ptr(cct_context, current_read_number, this_aln_id);
	for(better_voting_loc_id = 0; better_voting_loc_id < this_aln_id; better_voting_loc_id++){
		voting_location_t * better_r = cellCounts_global_retrieve_alignment_ptr(cct_context, current_read_number, better_voting_loc_id);
		if( this_r -> selected_position >= better_r -> selected_position - cct_context -> max_indel_length - 1
		 && this_r -> selected_position <=  better_r -> selected_position +  cct_context -> max_indel_length +1 )
		if( this_r -> confident_coverage_start >= better_r -> confident_coverage_start
		 && this_r -> confident_coverage_end <= better_r -> confident_coverage_end) return 1;
	}
	return 0;
}


void cellCounts_update_top_three(cellcounts_global_t * cct_context, int * top_buffer_3i, int new_value){
	if(new_value > top_buffer_3i[cct_context -> max_top_vote_simples - 1]){
		int x1;
		for(x1 = 0;x1 < cct_context -> max_top_vote_simples ; x1++){
			if(new_value > top_buffer_3i[x1]){
				int x2;
				for(x2 = cct_context -> max_top_vote_simples - 1 ; x2 > x1 ; x2 --){
					top_buffer_3i[x2] = top_buffer_3i[x2-1];
				}
				top_buffer_3i[x1] = new_value;
				break;
			}else if(new_value == top_buffer_3i[x1]) break;
		}
	}
}

void cellCounts_set_insertion_sequence(cellcounts_global_t * cct_context, int thread_no , char ** binary_bases , char * read_text , int insertions) {
	int xk1;

	(*binary_bases) = malloc((1+insertions)/4+2);
	//SUBREADprintf("ALLOC PTR=%p\n", (*binary_bases) );

	assert(insertions <= MAX_INSERTION_LENGTH);
	memset((*binary_bases),0, (1+insertions)/4+2);

	for(xk1=0; xk1<insertions;xk1++)
	{
		int byte_no = xk1/4;
		int bit_no = 2*(xk1%4);

		*((*binary_bases)+byte_no) |= (base2int(read_text[xk1]))<<bit_no;
	}
}

int cellCounts_search_event(cellcounts_global_t * cct_context, HashTable * event_table, chromosome_event_t * event_space, unsigned int pos, int search_type, unsigned char event_type, chromosome_event_t ** return_buffer) {
	int ret = 0;

	if(pos < 1 || pos > 0xffff0000)return 0;

	if(event_table->appendix1) {
		if(search_type == EVENT_SEARCH_BY_SMALL_SIDE)
			if(! check_event_bitmap(event_table -> appendix1, pos))return 0;

		if(search_type == EVENT_SEARCH_BY_LARGE_SIDE)
			if(! check_event_bitmap(event_table -> appendix2, pos))return 0;

		if(search_type == EVENT_SEARCH_BY_BOTH_SIDES)
			if(!(check_event_bitmap(event_table -> appendix1, pos) || check_event_bitmap(event_table -> appendix2, pos)))
				return 0;
	}

	unsigned int * res = HashTableGet(event_table, NULL+pos); 
	if(res) {
		int xk2;
		int current_size = res[0]&0x0fffffff;
		for(xk2=1; xk2< current_size+1 ; xk2++) {
			if(!res[xk2])break;
			chromosome_event_t * event_body = &event_space[res[xk2]-1]; 

			if((event_body -> event_type & event_type) == 0)continue;
			if(search_type == EVENT_SEARCH_BY_SMALL_SIDE && event_body -> event_small_side != pos)continue;
			if(search_type == EVENT_SEARCH_BY_LARGE_SIDE && event_body -> event_large_side != pos)continue;
			if(search_type == EVENT_SEARCH_BY_BOTH_SIDES && event_body -> event_small_side != pos && event_body -> event_large_side != pos)continue;

			return_buffer[ret++] = event_body;
		}
	}

	return ret;
}



chromosome_event_t * cellCounts_local_add_indel_event(cellcounts_global_t * cct_context, int thread_no, HashTable * event_table, char * read_text, unsigned int left_edge, int indels, int score_supporting_read_added, int is_ambiguous, int mismatched_bases,int * old_event_id)
{
	chromosome_event_t * found = NULL;
	chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	int found_events = cellCounts_search_event(cct_context, event_table, thread_context ->event_space_dynamic, left_edge, EVENT_SEARCH_BY_SMALL_SIDE, CHRO_EVENT_TYPE_INDEL|CHRO_EVENT_TYPE_LONG_INDEL, search_return);

	if(found_events) {
		int kx1; 
		for(kx1 = 0; kx1 < found_events ; kx1++) {
			if(search_return[kx1] -> indel_length == indels) {
				found = search_return[kx1];	
				break;
			}
		}
	}

	if(found){
		if(old_event_id) (*old_event_id)=found->global_event_id;
		found -> supporting_reads += score_supporting_read_added;
		return NULL;
	} else {
		int event_no;

		event_no = thread_context -> total_events ++;
		if(event_no >=thread_context -> event_space_capacity){
			thread_context -> event_space_capacity = thread_context -> event_space_capacity * 1.6;
			thread_context -> event_space_dynamic = realloc(thread_context -> event_space_dynamic , sizeof(chromosome_event_t) * thread_context -> event_space_capacity );
		}

		chromosome_event_t * new_event = thread_context ->event_space_dynamic + event_no; 
		memset(new_event , 0 , sizeof(chromosome_event_t));
		if(indels <0)	// if it is an insertion:
			cellCounts_set_insertion_sequence(cct_context, thread_no , &new_event -> inserted_bases , read_text, -indels);

		new_event -> event_small_side = left_edge;
		new_event -> event_large_side = left_edge + 1 + max(0, indels);
		new_event -> event_type = CHRO_EVENT_TYPE_INDEL;
		new_event -> indel_length = indels;
		new_event -> supporting_reads = 1;
		new_event -> event_quality = 1;//pow(0.5 , 3*mismatched_bases);

		put_new_event(event_table, new_event , event_no);
		return new_event;
	}
}

#define INDEL_MASK_BY_INSERTION 1
#define INDEL_MASK_BY_DELETION 2
#define INDEL_MASK_BY_MATCH 0
#define INDEL_MASK_BY_MISMATCH 3

int CELLCTS_DPALIGN_CREATEGAP_PENALTY = -1;
int CELLCTS_DPALIGN_EXTENDGAP_PENALTY = 0;
int CELLCTS_DPALIGN_MATCH_SCORE = 2;
int CELLCTS_DPALIGN_MISMATCH_PENALTY = 0;


int cellCounts_dynamic_align(cellcounts_global_t * cct_context, int thread_no, char * read, int read_len, unsigned int begin_position, char * movement_buffer, int expected_offset, char * read_name){
// read must be converted to the positive strand.
// movement buffer: 0:match, 1: read-insert, 2: gene-insert, 3:mismatch
// the size of the movement buffer must be equal to the length of the read plus max_indel * 3.
	int max_indel = min(16 , cct_context->max_indel_length);
	int i,j;
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;

	if(read_len < 3 || abs(expected_offset) > max_indel)
		return 0;
	if(expected_offset < 0 && read_len < (3-expected_offset))
		return 0;

	gene_value_index_t * current_value_index = thread_context->current_value_index;

	//unsigned long long table_ptr = (unsigned long long) indel_context -> dynamic_align_table;

	short ** table = thread_context -> dynamic_align_table;
	char ** table_mask = thread_context -> dynamic_align_table_mask;

	// vertical move: deletion (1)
	// horizontal move: insertion (2)
	// cross move: match (0) or mismatch (3)
	// i: vertical move; j: horizontal move

	for (i=0; i<read_len +  expected_offset; i++) {
		for(j=0; j<read_len; j++) {
			table_mask[i][j]=0;
			if (j < i - max_indel || j > max_indel + i) {
				table[i][j]=-9999;
				continue;
			}

			short from_upper;

			if (i>0) from_upper = table[i-1][j] + (table_mask[i-1][j] == INDEL_MASK_BY_DELETION?CELLCTS_DPALIGN_EXTENDGAP_PENALTY:CELLCTS_DPALIGN_CREATEGAP_PENALTY);
			else     from_upper = -9999;

			short from_left;

			if (j>0) from_left = table[i][j-1] + (table_mask[i][j-1] == INDEL_MASK_BY_INSERTION?CELLCTS_DPALIGN_EXTENDGAP_PENALTY:CELLCTS_DPALIGN_CREATEGAP_PENALTY);
			else     from_left = -9999;

			char chromo_ch = gvindex_get(current_value_index, begin_position + i);
			char is_matched_ij = (chromo_ch == read[j])?CELLCTS_DPALIGN_MATCH_SCORE:CELLCTS_DPALIGN_MISMATCH_PENALTY;
			
			short from_upperleft;

			if (i>0 && j>0) from_upperleft = table[i-1][j-1] + is_matched_ij;
			else if(i==0 && j==0) from_upperleft = is_matched_ij;
			else	    from_upperleft = -9999; 

			if (from_upperleft == from_upper && from_upperleft > from_left) {
				table_mask[i][j]= INDEL_MASK_BY_DELETION;
				table[i][j] = from_upper;
			} else if(from_upperleft == from_left && from_upperleft > from_upper) {
				table_mask[i][j]= INDEL_MASK_BY_INSERTION;
				table[i][j] = from_left;
			} else if(from_upperleft > from_left && from_upperleft > from_upper) {
				table_mask[i][j]= (chromo_ch == read[j])?INDEL_MASK_BY_MATCH:INDEL_MASK_BY_MISMATCH;
				table[i][j] = from_upperleft;
			} else if(from_upperleft == from_left && from_upperleft == from_upper) {
				table_mask[i][j]= (chromo_ch == read[j])?INDEL_MASK_BY_MATCH:INDEL_MASK_BY_MISMATCH;
				table[i][j] = from_upperleft;
			} else if(from_left > from_upper) {
				table_mask[i][j]= INDEL_MASK_BY_INSERTION;
				table[i][j] = from_left;
			} else if(from_left <= from_upper) {
				table_mask[i][j]= INDEL_MASK_BY_DELETION;
				table[i][j] = from_upper;
			}

		}
	}
	short path_i = read_len + expected_offset - 1;
	int out_pos = 0, delta=0;
	j = read_len - 1;

	while(1) {
		if(table_mask[path_i][j] == INDEL_MASK_BY_INSERTION) {
			j--;
			delta --;
			movement_buffer[out_pos++] = 2;
		} else if(table_mask[path_i][j] == INDEL_MASK_BY_DELETION) {
			path_i--;
			delta ++;
			movement_buffer[out_pos++] = 1;
		} else if(table_mask[path_i][j] == INDEL_MASK_BY_MATCH || table_mask[path_i][j] == INDEL_MASK_BY_MISMATCH) {
			movement_buffer[out_pos++] = table_mask[path_i][j] == INDEL_MASK_BY_MATCH?0:3;
			path_i--;
			j--;
		}

		if(path_i == -1 && j == -1) break;
		if(j<0 || path_i<0) return 0;
	}

	if(expected_offset!=delta)return 0;
	for(i=0; i<out_pos/2; i++) {
		char tmp;
		tmp = movement_buffer[out_pos-1-i];
		movement_buffer[out_pos-1-i] = movement_buffer[i];
		movement_buffer[i] = tmp;
	}
	return out_pos;
}



int cellCounts_find_new_indels(cellcounts_global_t * cct_context, int thread_no, int pair_number, char * read_name, char * read_text, char * qual_text, int read_len, int voting_loc_no) {
	int i,j;

	int last_correct_subread = 0;
	int last_indel = 0;

	gene_vote_number_t * indel_recorder;
	unsigned int voting_position;

	HashTable * event_table = NULL;
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;

	event_table = thread_context -> event_entry_table; 
	voting_location_t *current_result = cellCounts_global_retrieve_alignment_ptr(cct_context, pair_number, voting_loc_no);
	indel_recorder = current_result -> selected_indel_record; 
	voting_position = current_result -> selected_position; 

	if(!indel_recorder[0])
		return 0;

	gene_value_index_t * current_value_index = thread_context->current_value_index;

	for(i=0; indel_recorder[i]  && (i<MAX_INDEL_SECTIONS); i+=3) {
		int indels = indel_recorder[i+2] - last_indel;
			// -1 : 1 insert; 1: 1 delete
		int next_correct_subread = indel_recorder[i] -1;
		//int next_correct_subread_last = indel_recorder[i+1] -1;

		if(indels) {
			int last_correct_base = find_subread_end(read_len, cct_context->total_subreads_per_read , last_correct_subread) - 9;
			int first_correct_base = find_subread_end(read_len, cct_context->total_subreads_per_read , next_correct_subread) - 16 + 9;

			last_correct_base = max(0, last_correct_base);
			last_correct_base = min(read_len-1, last_correct_base);
			first_correct_base = min(first_correct_base, read_len-1);
			first_correct_base = max(0, first_correct_base);
			first_correct_base = max(first_correct_base, last_correct_base);

			char movement_buffer[MAX_READ_LENGTH * 10 / 7];
			//chromosome_event_t * last_event = NULL;
			int last_event_id = -1;

			first_correct_base = min(first_correct_base+10, read_len);
			int x1, dyna_steps;

			dyna_steps = cellCounts_dynamic_align(cct_context, thread_no, read_text + last_correct_base, first_correct_base - last_correct_base, voting_position + last_correct_base + last_indel, movement_buffer, indels, read_name);
			movement_buffer[dyna_steps]=0;

			unsigned int cursor_on_chromosome = voting_position + last_correct_base + last_indel, cursor_on_read = last_correct_base;
			int last_mv = 0;
			unsigned int indel_left_boundary = 0;
			int is_in_indel = 0, current_indel_len = 0, total_mismatch = 0;

			for(x1=0; x1<dyna_steps;x1++) {
				int mv=movement_buffer[x1];
				if(mv==3) total_mismatch++;
			}

			if(total_mismatch<=2){
				for(x1=0; x1<dyna_steps;x1++)
				{
					int mv=movement_buffer[x1];

					if(last_mv != mv)
					{
						if( ( mv==1 || mv==2 ) && ! is_in_indel)
						{
							indel_left_boundary = cursor_on_chromosome;
							is_in_indel = 1;
							current_indel_len = 0;
						}
						else if ( is_in_indel && (mv == 0 || mv == 3)  )
						{
							gene_value_index_t * current_value_index = thread_context->current_value_index;

							int ambiguous_count=0;
							if(0){
								int ambiguous_i, best_matched_bases = match_chro(read_text + cursor_on_read - 6, current_value_index, indel_left_boundary - 6, 6, 0, GENE_SPACE_BASE)  +
											 match_chro(read_text + cursor_on_read - min(current_indel_len,0), current_value_index, indel_left_boundary + max(0, current_indel_len), 6, 0, GENE_SPACE_BASE);
								for(ambiguous_i=-5; ambiguous_i<=5; ambiguous_i++)
								{
									int left_match = match_chro(read_text + cursor_on_read - 6, current_value_index, indel_left_boundary - 6, 6+ambiguous_i, 0, GENE_SPACE_BASE);
									int right_match = match_chro(read_text + cursor_on_read + ambiguous_i - min(current_indel_len,0), current_value_index, indel_left_boundary + ambiguous_i + max(0, current_indel_len), 6-ambiguous_i, 0, GENE_SPACE_BASE);
									if(left_match+right_match == best_matched_bases) ambiguous_count ++;
								}
							}

							if(abs(current_indel_len)<=cct_context -> max_indel_length)
							{
								int  old_event_id = -1;

								chromosome_event_t * new_event = cellCounts_local_add_indel_event(cct_context, thread_no, event_table, read_text + cursor_on_read + min(0,current_indel_len), indel_left_boundary - 1, current_indel_len, 1, ambiguous_count, 0, &old_event_id);
								mark_gapped_read(current_result);

								chromosome_event_t * event_space = NULL;
								event_space = thread_context -> event_space_dynamic;

								if(last_event_id>=0){
									chromosome_event_t * last_event = event_space + last_event_id;
									int dist = indel_left_boundary - last_event -> event_large_side ;

									if(last_event -> connected_next_event_distance<1)
										last_event -> connected_next_event_distance = dist;
									else    last_event -> connected_next_event_distance = min(dist, last_event -> connected_next_event_distance);

									chromosome_event_t * old_new_event = new_event;
									if(old_event_id>=0){
										assert(!new_event);
										old_new_event = event_space + old_event_id;
									}

									if(old_new_event -> connected_previous_event_distance < 1)
										old_new_event -> connected_previous_event_distance = dist;
									else    old_new_event -> connected_previous_event_distance = min(dist, old_new_event -> connected_previous_event_distance);
								}

								if (new_event)
									last_event_id = new_event -> global_event_id;
								else	last_event_id = old_event_id;
							}
						}
						

						if(mv == 0 || mv == 3)
							is_in_indel = 0;
					}

					if(is_in_indel && mv == 1)
						current_indel_len += 1;
					if(is_in_indel && mv == 2)
						current_indel_len -= 1;

					if(mv == 1 || mv == 3 || mv == 0) cursor_on_chromosome++;
					if(mv == 2 || mv == 3 || mv == 0) cursor_on_read++;

					last_mv = mv;
				}

			}
		}

		last_correct_subread = indel_recorder[i+1]-1;
		last_indel = indel_recorder[i+2];
	}
	return 0;
}

int cellCounts_is_pos_in_annotated_exon_regions(cellcounts_global_t * cct_context, unsigned int pos){ 
	int range_start = (pos+5) / EXONIC_REGION_RESOLUTION;
	int range_stop = (pos-5) / EXONIC_REGION_RESOLUTION;
	
	for(; range_start <=range_stop; range_start++){
		int exonic_map_byte= range_start / 8;
		int exonic_map_bit = range_start % 8;
		if (cct_context ->exonic_region_bitmap [exonic_map_byte] & (1<<exonic_map_bit))return 1;
	}
	return 0;
}

#define SE_READ_IN_KNOWN_EXON_REWARD 1

void cellCounts_copy_vote_to_alignment_res(cellcounts_global_t * cct_context, int thread_no, voting_location_t * align_res, gene_vote_t * current_vote, int vote_i, int vote_j, int curr_read_len, char * read_name, char * curr_read_text, int used_subreads_in_vote, subread_read_number_t pair_number){
	int vv = current_vote -> votes[vote_i][vote_j];
	vv += SE_READ_IN_KNOWN_EXON_REWARD *cellCounts_is_pos_in_annotated_exon_regions(cct_context, current_vote -> pos[vote_i][vote_j]);
	align_res -> selected_position = current_vote -> pos[vote_i][vote_j];
	align_res -> selected_votes = vv;
	align_res -> indels_in_confident_coverage = indel_recorder_copy(align_res -> selected_indel_record, current_vote -> indel_recorder[vote_i][vote_j]);
	align_res -> confident_coverage_end = current_vote -> coverage_end[vote_i][vote_j];
	align_res -> confident_coverage_start = current_vote -> coverage_start[vote_i][vote_j];
	align_res -> result_flags = (current_vote -> masks[vote_i][vote_j] & IS_NEGATIVE_STRAND)?(CORE_IS_NEGATIVE_STRAND):0;
	align_res -> used_subreads_in_vote = used_subreads_in_vote;
}

int cellCounts_process_voting_junction(cellcounts_global_t * cct_context, int thread_no, subread_read_number_t pair_number, gene_vote_t * votetab, char * read_name, char * read_text, int read_len, int is_negative_strand, gene_vote_number_t all_subreads) {
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	topK_buffer_t * topbuf = &thread_context->topKbuff;

	simple_mapping_t * vote_simple_buffer = (simple_mapping_t *) topbuf -> vote_simple_1_buffer;

	int i,j;
	int top_three_buff[9];
	memset(top_three_buff, 0 , cct_context -> max_top_vote_simples * sizeof(int));

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++){
		for (j=0; j< votetab->items[i]; j++){
			int vv = votetab -> votes[i][j];
			vv += SE_READ_IN_KNOWN_EXON_REWARD*cellCounts_is_pos_in_annotated_exon_regions(cct_context, votetab -> pos[i][j]);
			cellCounts_update_top_three(cct_context, top_three_buff, vv);
		}
	}

	for(i = 0; i < cct_context -> max_voting_locations; i++){
		voting_location_t * old_result = cellCounts_global_retrieve_alignment_ptr(cct_context, pair_number, i);
		if(old_result -> selected_votes>0)
			cellCounts_update_top_three(cct_context, top_three_buff, old_result -> selected_votes);
	}

	int simple_record_numbers= 0, third_k;

	// populate the simple stub list
	for(third_k = 0 ; third_k < cct_context -> max_top_vote_simples; third_k ++){
		if(simple_record_numbers >= cct_context -> max_voting_simples)break;
		int this_vote_N = top_three_buff[third_k];
		if(this_vote_N < 1 || (top_three_buff[0] - this_vote_N > cct_context -> max_vote_number_cutoff )) break;

		for (i=0; i<GENE_VOTE_TABLE_SIZE; i++){
			if(simple_record_numbers >= cct_context -> max_voting_simples)break;
			for (j=0; j< votetab->items[i]; j++){
				if(simple_record_numbers >= cct_context -> max_voting_simples)break;

				int vv = votetab->votes[i][j];
				vv += SE_READ_IN_KNOWN_EXON_REWARD*cellCounts_is_pos_in_annotated_exon_regions(cct_context,  votetab -> pos[i][j]);
				if(vv == this_vote_N && votetab->votes[i][j] >= cct_context -> min_votes_per_mapped_read){
					vote_simple_buffer[simple_record_numbers].is_vote_t_item = 1;
					vote_simple_buffer[simple_record_numbers].item_index_i = i;
					vote_simple_buffer[simple_record_numbers].item_index_j = j;
					vote_simple_buffer[simple_record_numbers].read_start_base = votetab -> coverage_start[i][j];
					vote_simple_buffer[simple_record_numbers].mapping_position = votetab -> pos[i][j];
					vote_simple_buffer[simple_record_numbers].major_half_votes = vv;
					simple_record_numbers ++;
				}
			}
		}

		for(i = 0; i < cct_context -> max_voting_locations; i++){
			voting_location_t * old_result = cellCounts_global_retrieve_alignment_ptr(cct_context, pair_number, i);
			if(simple_record_numbers >= cct_context -> max_voting_simples)break;
			if(old_result -> selected_votes == this_vote_N) {
				vote_simple_buffer[simple_record_numbers].is_vote_t_item = 0;
				vote_simple_buffer[simple_record_numbers].item_index_i = i;
				vote_simple_buffer[simple_record_numbers].mapping_position = old_result -> selected_position;
				vote_simple_buffer[simple_record_numbers].major_half_votes = old_result -> selected_votes;
				vote_simple_buffer[simple_record_numbers].read_start_base = old_result -> confident_coverage_start;

				simple_record_numbers ++;
			}
		}
	}

	voting_location_t * alignment_tmp = (voting_location_t *) topbuf -> alignment_tmp_r1;
	memset(alignment_tmp, 0, sizeof(voting_location_t) * cct_context->max_voting_locations);

	int alignment_res_cursor = 0;

	if(simple_record_numbers>0) {
		for(i = 0; i < simple_record_numbers; i++){
			if(alignment_res_cursor >= cct_context -> max_voting_locations)break;
			simple_mapping_t * current_loc = vote_simple_buffer+i;

			if(current_loc -> major_half_votes < cct_context -> min_votes_per_mapped_read) continue;
			unsigned int current_pos = current_loc -> mapping_position;

			int is_exist = 0;
			for(j = 0; j < alignment_res_cursor; j++) {
				if(alignment_tmp[j].selected_position == current_pos){
					is_exist = 1;
					break;
				}
			}
			if(!is_exist){
				if(current_loc -> is_vote_t_item)
					cellCounts_copy_vote_to_alignment_res(cct_context, thread_no, alignment_tmp + alignment_res_cursor, votetab, current_loc -> item_index_i, current_loc -> item_index_j, read_len, read_name, read_text, all_subreads , pair_number);
				else{
					//checked:boundary
					memcpy(alignment_tmp + alignment_res_cursor, cellCounts_global_retrieve_alignment_ptr(cct_context, pair_number, current_loc -> item_index_i), sizeof(voting_location_t));
				}

				alignment_res_cursor++;
			}
		}
	}

	for(i = 0; i < cct_context->max_voting_locations ; i++){
		voting_location_t * cur_res = cellCounts_global_retrieve_alignment_ptr(cct_context, pair_number, i);
		if( i < alignment_res_cursor)
			memcpy(cur_res, alignment_tmp + i, sizeof(voting_location_t));
		else	cur_res -> selected_votes = 0;
	}

	return 0;
}

int cellCounts_do_voting(cellcounts_global_t * cct_context, int thread_no) {
	int xk1;
	subread_read_number_t current_read_number=0;
	char * read_text_1, * qual_text_1;
	char read_name_1[MAX_READ_NAME_LEN+1];
	int read_len_1=0;
	int min_first_read_votes = cct_context -> min_votes_per_mapped_read; 
	int voting_max_indel_length = min(16, cct_context->max_indel_length);
	int sqr_interval=10000;

	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	read_text_1 = malloc(MAX_READ_LENGTH+1);
	qual_text_1 = malloc(MAX_READ_LENGTH+1);
	gene_vote_t * vote_1 = malloc(sizeof(gene_vote_t));

	if(vote_1==NULL) {
		SUBREADprintf("Cannot allocate voting memory.\n");
		return -1;
	}

	unsigned int low_index_border = cct_context -> current_value_index -> start_base_offset;
	unsigned int high_index_border = cct_context -> current_value_index -> start_base_offset + cct_context -> current_value_index -> length; 

	int index_gap_width = cct_context -> current_index -> index_gap;

	while(1) {
		int subread_no;
		int is_reversed, applied_subreads = 0, v1_all_subreads=0;

		cellCounts_fetch_next_read_pair(cct_context, thread_no,  &read_len_1, read_name_1, read_text_1, qual_text_1, &current_read_number);
		if(current_read_number < 0) break;

		for(is_reversed = 0; is_reversed<2; is_reversed++) {
			gene_vote_t * current_vote = vote_1;
			char * current_read = read_text_1;
			int current_rlen = read_len_1;
			int subread_step;
			if(current_rlen< 16) continue;
			int CR15GLS = (current_rlen - 15 - index_gap_width)<<16;
			if(current_rlen<= EXON_LONG_READ_LENGTH)
			{
				subread_step =  CR15GLS /(cct_context -> total_subreads_per_read -1);
				if(subread_step<(index_gap_width<<16))subread_step = index_gap_width<<16;
			}else{
				subread_step = 6<<16;
				if( CR15GLS / subread_step > 62)
					subread_step = CR15GLS/62;
			}

			applied_subreads = 1 + CR15GLS / subread_step;
			v1_all_subreads = applied_subreads;

			unsigned int current_high_border = high_index_border -  current_rlen;

			int allow_indel_i;
			unsigned int shift_indel_locs[ GENE_VOTE_TABLE_SIZE * GENE_VOTE_SPACE ];
			unsigned int shift_indel_NO = 0;
			for(allow_indel_i = 0; allow_indel_i<2; allow_indel_i++){
				if(is_reversed==0) init_gene_vote(vote_1);

				for(subread_no=0; subread_no < applied_subreads ; subread_no++) {
					for(xk1=0; xk1<index_gap_width ; xk1++) {

						int subread_offset = ((subread_step * subread_no) >> 16);
						if(index_gap_width > 1) subread_offset -= subread_offset%(index_gap_width) - xk1; 

						char * subread_string = current_read + subread_offset;

						gehash_key_t subread_integer = genekey2int(subread_string, GENE_SPACE_BASE);
						gehash_go_X(cct_context->current_index, subread_integer , subread_offset, current_rlen, is_reversed, current_vote,  voting_max_indel_length, subread_no,  low_index_border, current_high_border, allow_indel_i, shift_indel_locs, &shift_indel_NO);
					}
				}
				if(shift_indel_NO == 0)break;
			}



			if(is_reversed) {
				if(current_read_number % 50000 == 0) SUBREADprintf("voting step: %d\n", current_read_number);
				if(0&&current_read_number == 4000){
					SUBREADprintf(">>>%llu<<<\n%s [%d]  %s VOTE1_MAX=%d >= %d\n", current_read_number, read_name_1, read_len_1, read_text_1, vote_1->max_vote, min_first_read_votes);
					SUBREADprintf(" ======= PAIR %s = %llu =======\n", read_name_1, current_read_number);
					print_votes(vote_1, cct_context -> index_prefix);
				}

				if(vote_1->max_vote >= min_first_read_votes)
					cellCounts_process_voting_junction(cct_context, thread_no, current_read_number, vote_1, read_name_1, read_text_1, read_len_1, is_reversed, v1_all_subreads);
			} else reverse_read(read_text_1, read_len_1, GENE_SPACE_BASE);
		}

		int read_1_reversed = 1;

		if(cct_context -> is_final_voting_run) {
			int * this_read_has_reversed = &read_1_reversed;
			char * current_read_name =  read_name_1;
			char * current_read =  read_text_1;
			char * current_qual =  qual_text_1;
			int current_rlen = read_len_1;
			int voting_loc_id;

			for(voting_loc_id = 0; voting_loc_id < cct_context -> max_voting_locations; voting_loc_id++) {
				voting_location_t * current_r = cellCounts_global_retrieve_alignment_ptr(cct_context, current_read_number, voting_loc_id);
				if(current_r -> selected_votes < 1) continue;
				
				int this_read_should_be_reversed = (current_r -> result_flags & CORE_IS_NEGATIVE_STRAND) ? 1:0;

				if(this_read_should_be_reversed != (*this_read_has_reversed)) {
					(*this_read_has_reversed) = !(*this_read_has_reversed);
					reverse_read(current_read, current_rlen, GENE_SPACE_BASE);
					if(current_qual)
						reverse_quality(current_qual , current_rlen);
				}

				int rv = cellCounts_locate_current_value_index(cct_context, thread_no, current_r, current_rlen);
				//SUBREADprintf("LOCATE in THR %d for pos=%u has %d\n", thread_no, current_r->selected_position, rv);
				assert(rv==0);
				if(!cellCounts_has_better_mapping(cct_context, thread_no, current_read_number, voting_loc_id))
					cellCounts_find_new_indels(cct_context, thread_no, current_read_number, current_read_name, current_read, current_qual, current_rlen, voting_loc_id);
			}
		}
	}

	free(vote_1);
	free(read_text_1);
	free(qual_text_1);

	return 0;
}




#define _test_conc_size if(conc_rec_items >= conc_rec_size - 1){\
			conc_rec_size *= 1.5;\
			records = realloc(records, conc_rec_size * sizeof(concatinating_events_record_t));\
			if(NULL == records) return -1;\
		}

#define _get_event_body(iii) (thread_contexts[records[iii].thread_no].event_space_dynamic + records[iii].thread_bodytable_number)

int cellCounts_conc_sort_compare(void * arr, int l, int r){
	void ** arrr = (void **)arr;
	concatinating_events_record_t * records = arrr[0];
	cellcounts_global_t * cct_context = arrr[1];
	cellcounts_align_thread_t * thread_contexts = arrr[2];

	chromosome_event_t * body_l = _get_event_body(l);
	chromosome_event_t * body_r = _get_event_body(r);


	if(body_l -> event_small_side > body_r -> event_small_side)return 3;
	if(body_l -> event_small_side < body_r -> event_small_side)return -3;
	if(body_l -> event_large_side > body_r -> event_large_side)return 3;
	if(body_l -> event_large_side < body_r -> event_large_side)return -3;

	if(abs(body_l -> indel_length) < abs(body_r -> indel_length)) return 2;
	if(abs(body_l -> indel_length) > abs(body_r -> indel_length)) return -2;
	if(body_l -> indel_length > body_r -> indel_length) return -2; // same length, but L is del and R is ins -- prefer del than ins
	if(body_l -> indel_length < body_r -> indel_length) return 2;

	if((body_l -> is_donor_found_or_annotation & 64)!=0 && (body_r -> is_donor_found_or_annotation & 64)==0) return 1; // prefer known truth.
	if((body_l -> is_donor_found_or_annotation & 64)==0 && (body_r -> is_donor_found_or_annotation & 64)!=0) return -1;
	if(body_l -> supporting_reads > body_r -> supporting_reads) return -1;
	if(body_l -> supporting_reads < body_r -> supporting_reads) return 1;
	return 0;
}

void cellCounts_conc_sort_exchange(void * arr, int l, int r){
	void ** arrr = (void **)arr;
	concatinating_events_record_t * records = arrr[0];

	unsigned int tmpi;
	tmpi = records[l].thread_bodytable_number;
	records[l].thread_bodytable_number =  records[r].thread_bodytable_number;
	records[r].thread_bodytable_number = tmpi;

	tmpi = records[l].thread_no;
	records[l].thread_no =  records[r].thread_no;
	records[r].thread_no = tmpi;
}

void cellCounts_conc_sort_merge(void * arr,  int start, int items, int items2){
	void ** arrr = (void **) arr;
	concatinating_events_record_t * records = arrr[0];

	int read_1_ptr = start, read_2_ptr = start+items, write_ptr;
	concatinating_events_record_t * merged_records = malloc(sizeof(concatinating_events_record_t) * (items+items2));

	for(write_ptr=0; write_ptr<items+items2; write_ptr++){
		if((read_1_ptr >= start+items)||(read_2_ptr < start+items+items2 && cellCounts_conc_sort_compare(arr, read_1_ptr, read_2_ptr) > 0))
			memcpy(merged_records+write_ptr, records+(read_2_ptr++), sizeof(concatinating_events_record_t));
		else
			memcpy(merged_records+write_ptr, records+(read_1_ptr++), sizeof(concatinating_events_record_t));
	}
	memcpy(records + start, merged_records, sizeof(concatinating_events_record_t) * (items+items2));
	free(merged_records);
}

int cellCounts_finalise_indel_and_junction_thread(cellcounts_global_t * cct_context, int task) {
	cellcounts_align_thread_t * thread_contexts = cct_context -> all_thread_contexts;
	
	if(task == STEP_VOTING) {
		concatinating_events_record_t * records;
		int conc_rec_size = 10000, conc_rec_items = 0;
		records = malloc(sizeof(concatinating_events_record_t) * conc_rec_size);

		int xk1, thn;
		for(thn =0; thn < cct_context->total_threads; thn++){
			cellcounts_align_thread_t * thread_context = thread_contexts+thn;

			for(xk1 = 0; xk1 < thread_context -> total_events; xk1++){
			chromosome_event_t * old_body = thread_context  ->  event_space_dynamic + xk1;
			if(old_body -> event_type == CHRO_EVENT_TYPE_REMOVED) continue;

			_test_conc_size;
			records[conc_rec_items].thread_no=thn;
			records[conc_rec_items++].thread_bodytable_number=xk1;
			}
		}

		for(xk1 = 0; xk1 < cct_context -> total_events; xk1++){
			chromosome_event_t * old_body = cct_context ->  event_space_dynamic + xk1;
			if(old_body -> event_type == CHRO_EVENT_TYPE_REMOVED) continue;
			_test_conc_size;
			records[conc_rec_items].thread_no=-1;
			records[conc_rec_items++].thread_bodytable_number=xk1;
		}

		void * sort_arr[3];
		sort_arr[0] = records;
		sort_arr[1] = cct_context;
		sort_arr[2] = thread_contexts;

		// many repeated elements -- do not use quick sort.
		merge_sort(sort_arr, conc_rec_items, cellCounts_conc_sort_compare, cellCounts_conc_sort_exchange, cellCounts_conc_sort_merge);

		chromosome_event_t * prev_env = NULL;
		int merge_target_size = 10000;
		int merge_target_items = 0;
		chromosome_event_t * merge_target = malloc(sizeof(chromosome_event_t) * merge_target_size);

		int merge_start = 0;
		for(xk1 = 0; xk1 <= conc_rec_items; xk1++){
			chromosome_event_t * this_event = (xk1 == conc_rec_items)?NULL:_get_event_body(xk1);
			if(xk1 > 0){
				int compret = 0;
				if(xk1 < conc_rec_items) compret = cellCounts_conc_sort_compare(sort_arr, xk1-1, xk1);
				if(abs(compret)>1 || xk1 == conc_rec_items){// different events or last one -- merge [ prev_event_record_no , xk1 - 1 ]
					int xk_merge;
					// find a new slot in the target space
					if(merge_target_items >= merge_target_size - 1){
						merge_target_size *= 1.5;
						merge_target = realloc(merge_target, sizeof(chromosome_event_t) * merge_target_size);
					}

					chromosome_event_t * merged_body = merge_target + merge_target_items;
					memcpy( merged_body, prev_env, sizeof(chromosome_event_t) );
					merged_body -> global_event_id = merge_target_items++;
					for(xk_merge = merge_start; xk_merge < xk1 - 1; xk_merge++){
						chromosome_event_t * old_body = _get_event_body(xk_merge);

						assert(merged_body -> event_type == old_body -> event_type);
						merged_body -> supporting_reads += old_body -> supporting_reads;
						merged_body -> anti_supporting_reads += old_body -> anti_supporting_reads;
						merged_body -> final_counted_reads += old_body -> final_counted_reads;
						merged_body -> final_reads_mismatches += old_body -> final_reads_mismatches;
						merged_body -> critical_supporting_reads += old_body -> critical_supporting_reads;
						merged_body -> junction_flanking_left = max(merged_body -> junction_flanking_left, old_body -> junction_flanking_left);
						merged_body -> junction_flanking_right = max(merged_body -> junction_flanking_right, old_body -> junction_flanking_right);
						merged_body -> is_donor_found_or_annotation |= old_body -> is_donor_found_or_annotation;

						if(merged_body -> connected_next_event_distance > 0 && old_body -> connected_next_event_distance > 0)
							merged_body -> connected_next_event_distance = min(merged_body -> connected_next_event_distance, old_body -> connected_next_event_distance);
						else	merged_body -> connected_next_event_distance = max(merged_body -> connected_next_event_distance, old_body -> connected_next_event_distance);

						if(merged_body -> connected_previous_event_distance > 0 && old_body -> connected_previous_event_distance > 0)
							merged_body -> connected_previous_event_distance = min(merged_body -> connected_previous_event_distance, old_body -> connected_previous_event_distance);
						else	merged_body -> connected_previous_event_distance = max(merged_body -> connected_previous_event_distance, old_body -> connected_previous_event_distance);

						if(old_body->inserted_bases && old_body -> event_type == CHRO_EVENT_TYPE_INDEL  && old_body -> indel_length<0){
							if(merged_body -> inserted_bases && merged_body -> inserted_bases != old_body -> inserted_bases){
								free(old_body -> inserted_bases);
							}
							else merged_body -> inserted_bases = old_body -> inserted_bases;
							old_body -> inserted_bases = NULL;
						}
					}
					merge_start = xk1;
				}
			}
			if(xk1 == conc_rec_items) break;
			prev_env = this_event;
		}

		free(records);

		for(thn =0; thn < cct_context->total_threads; thn++){
			cellcounts_align_thread_t * thread_context = thread_contexts+thn;
			
			destory_event_entry_table(thread_context -> event_entry_table);
			HashTableDestroy(thread_context -> event_entry_table);

			free(thread_context -> event_space_dynamic);

			for(xk1=0; xk1<MAX_READ_LENGTH; xk1++) {
				free(thread_context -> dynamic_align_table[xk1]);
				free(thread_context -> dynamic_align_table_mask[xk1]);
			}

			free(thread_context->dynamic_align_table);
			free(thread_context->dynamic_align_table_mask);
			free(thread_context);
		}
		if(cct_context -> event_space_dynamic ) free(cct_context -> event_space_dynamic);

		cct_context -> event_space_dynamic = merge_target;
		cct_context -> event_space_capacity = merge_target_size;
		cct_context -> total_events = merge_target_items;
	} else if(task == STEP_ITERATION_TWO) {
		int xk1, thn;
		for(thn =0; thn < cct_context->total_threads; thn++){
			cellcounts_align_thread_t * thread_context = thread_contexts+thn;

			for(xk1 = 0; xk1 < thread_context -> total_events; xk1++) {
				cct_context -> event_space_dynamic [xk1].final_counted_reads += thread_context -> final_counted_reads_array[xk1];
				cct_context -> event_space_dynamic [xk1].final_reads_mismatches += thread_context -> final_reads_mismatches_array[xk1];
			}
			free(thread_context -> final_counted_reads_array);
			free(thread_context -> final_reads_mismatches_array);
			free(thread_context);
		}
	}
	return 0;
}

int cellCounts_run_mapping(cellcounts_global_t * cct_context){
	int chunk_no = 0;

	cct_context -> current_index = (gehash_t*) malloc(sizeof(gehash_t));
	cct_context -> current_value_index = cct_context -> all_value_indexes;
	cct_context -> running_processed_reads_in_chunk=0;
	cct_context -> processed_reads_in_chunk=0;

	while(1) {
		int ret = 0;

		cellCounts_locate_read_files(cct_context, SEEK_SET);
		for(cct_context->current_index_block_number = 0; cct_context->current_index_block_number < cct_context->total_index_blocks; cct_context->current_index_block_number++) {	   
			char tmp_fname[MAX_FILE_NAME_LENGTH+30];

			if(cct_context->total_index_blocks > 1 || chunk_no == 0) {	   
				sprintf(tmp_fname, "%s.%02d.b.tab", cct_context->index_prefix, cct_context->current_index_block_number);
				print_in_box(80,0,0, "Load the %d-th index block...",1+ cct_context->current_index_block_number);
				
				if(gehash_load(cct_context -> current_index, tmp_fname)) return -1;
				print_in_box(80,0,0, "The index block has been loaded.");
				sprintf(tmp_fname, "%s.%02d.b.array", cct_context->index_prefix, cct_context->current_index_block_number);
			}
			cct_context -> current_value_index = cct_context -> all_value_indexes + cct_context->current_index_block_number;
			
			if(cct_context->total_index_blocks == cct_context->current_index_block_number + 1)
				cct_context -> is_final_voting_run = 1;
			else	cct_context -> is_final_voting_run = 0;
			
			ret = cellCounts_run_maybe_threads(cct_context, STEP_VOTING);
			
			if(0 == cct_context->current_index_block_number) {
				cellCounts_locate_read_files(cct_context, SEEK_END); 
				cct_context -> processed_reads_in_chunk = cct_context -> running_processed_reads_in_chunk;
			}
			
			if(cct_context->current_index_block_number < cct_context->total_index_blocks -1)
				cellCounts_go_chunk_start(cct_context);
			
			int is_last_chunk = cct_context -> processed_reads_in_chunk < cct_context-> reads_per_chunk;
			
			if(cct_context->total_index_blocks > 1 || is_last_chunk)
				gehash_destory_fast(cct_context -> current_index);
			
			if(ret) break;
			if(!cct_context -> processed_reads_in_chunk) break;
		}

		ret = ret || cellCounts_anti_supporting_read_scan(cct_context);
		cellCounts_remove_neighbour(cct_context);

		cellCounts_go_chunk_start(cct_context);
		ret = ret || cellCounts_run_maybe_threads(cct_context, STEP_ITERATION_TWO);

		cellCounts_go_chunk_nextchunk(cct_context);

		if(ret) return ret;

		if(cct_context -> processed_reads_in_chunk < cct_context -> reads_per_chunk ||
		  (cct_context -> output_binfiles_are_full))
			// base value indexes loaded in the last circle are not destroyed and are used in writting the indel VCF.
			// the indexes will be destroyed in destroy_global_context
			break;

		cellCounts_clean_context_after_chunk(cct_context);
		chunk_no++;
	}

	free(cct_context -> current_index);
	return 0;
}

int cellCounts_run_counting(cellcounts_global_t * cct_context){
	return 0;
}

#ifdef MAKE_STANDALONE
	#define cellCounts_main main
#endif
int cellCounts_main(int argc, char** argv){
	cellcounts_global_t * cct_context = calloc(sizeof(cellcounts_global_t),1);

	int ret = 0;

	ret = ret || cellCounts_args_context(cct_context, argc, argv);
	ret = ret || cellCounts_load_context(cct_context);
	ret = ret || cellCounts_run_mapping(cct_context);
	ret = ret || cellCounts_run_counting(cct_context);
	ret = ret || cellCounts_destroy_context(cct_context);

	free(cct_context);
	if(ret) SUBREADprintf("cellCounts terminates with errors.\n");
	return ret;
}
