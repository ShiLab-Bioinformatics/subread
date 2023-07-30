#define MAKE_CELLCOUNTS
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <locale.h>
#include <getopt.h>
#include <math.h>
#include <pthread.h>
#include <sys/stat.h>
#include "subread.h"
#include "gene-algorithms.h"
#include "input-files.h"
#include "interval_merge.h"
#include "input-blc.h"
#include "core-indel.h"
#include "core-junction.h"

//#define DO_STARSOLO_THING

#define CELLCOUNTS_REALIGNMENT_TRIES 15
#define MAX_FC_READ_LENGTH 10001
#define READ_BIN_BUF_SIZE 1000 // sufficient for a <=150bp read.
#define CELLBC_BATCH_NUMBER 149
#define CIGAR_PERFECT_SECTIONS 12 
#define MAX_UMI_LEN 14
#define MAX_SCRNA_SAMPLE_NUMBER 40 
#define MAX_SUBREADS_PER_READ 30
#define SCRNA_SUBREADS_HARD_LIMIT 20 


// the configurations used in our cellCounts paper
#define GENE_SCRNA_VOTE_SPACE 3
#define GENE_SCRNA_VOTE_TABLE_SIZE 17 


typedef struct{
	gene_vote_number_t max_vote;
	int max_vote_IJ;
	gehash_data_t max_position;
	gene_quality_score_t max_quality;
	gene_vote_number_t max_indel_recorder[MAX_INDEL_TOLERANCE*3];
	gene_vote_number_t * max_tmp_indel_recorder;
	int max_mask;
	gene_vote_number_t noninformative_subreads;

	unsigned short items[GENE_SCRNA_VOTE_TABLE_SIZE];
	unsigned int pos [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	int masks [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	int marked_shift_indel[GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	gene_vote_number_t votes [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	gene_quality_score_t quality [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	gene_vote_number_t last_subread_cluster [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	gene_vote_number_t indel_recorder [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE][MAX_INDEL_TOLERANCE*3];
	char current_indel_cursor[GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	char toli[GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	int topK_votes[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	int topK_IJ[SCRNA_HIGHEST_REPORTED_ALIGNMENTS ];

	short coverage_start [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	short coverage_end [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	short max_coverage_start;
	short max_coverage_end;
} gene_sc_vote_t;

typedef struct {
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

typedef struct {
	int votes[MAX_SUBREADS_PER_READ *2];
	int offsets[MAX_SUBREADS_PER_READ *2];
	unsigned int * start_location_in_index [MAX_SUBREADS_PER_READ *2]; // * 2 because positive and negative strand votes are evaluated together.
} temp_votes_per_read_t;

typedef struct{
	int thread_no;
	pthread_t thread;

	int event_space_capacity;
	int total_events;
	HashTable * event_entry_table;

	topK_buffer_t topKbuff;
	short * final_reads_mismatches_array;
	short * final_counted_reads_array;

	int hits_number_capacity;
	int * hits_start_pos;
	int * hits_length;
	char ** hits_chro;
	srInt_64 * hits_indices;

	srInt_64 mapped_reads_per_sample[MAX_SCRNA_SAMPLE_NUMBER];
	srInt_64 assigned_reads_per_sample[MAX_SCRNA_SAMPLE_NUMBER];
	srInt_64 reads_per_sample[MAX_SCRNA_SAMPLE_NUMBER+1];
	srInt_64 bcl_input_local_start_no;
	int bcl_input_local_filled, bcl_input_local_cached;
	char bcl_input_local_readbin[BCL_READBIN_ITEMS_LOCAL][BCL_READBIN_SIZE];
	int bcl_input_local_readlane[BCL_READBIN_ITEMS_LOCAL];

	srInt_64 hiconf_map;
	srInt_64 loconf_map;
	int reporting_count;
	int reporting_multi_alignment_no, reporting_this_alignment_no;
	srInt_64 reporting_scores[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	srInt_64 reporting_flags[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	unsigned int reporting_positions[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	int reporting_mapq[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	char reporting_cigars[SCRNA_HIGHEST_REPORTED_ALIGNMENTS][MAX_SCRNA_READ_LENGTH+20];
	int reporting_editing_distance[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];

} cellcounts_align_thread_t;

typedef struct{
	int total_threads;
	cellcounts_align_thread_t * all_thread_contexts;
	int reads_per_chunk;
	int allow_multi_overlapping_reads;
	int max_voting_simples;
	int max_voting_locations;
	int max_reported_alignments_per_read;
	int max_indel_length;
	int max_distinct_top_vote_numbers;
	int max_differential_from_top_vote_number;
	int max_mismatching_bases_in_reads;
	int min_mapped_length_for_mapped_read;
	int min_votes_per_mapped_read;
	int total_subreads_per_read;
	int report_multi_mapping_reads;
	int is_BAM_and_FQ_out_generated;
	int current_dataset_no;

	int processed_reads_in_chunk;
	int running_processed_reads_in_chunk;

	srInt_64 mapped_reads_per_sample[MAX_SCRNA_SAMPLE_NUMBER];
	srInt_64 assigned_reads_per_sample[MAX_SCRNA_SAMPLE_NUMBER];
	srInt_64 reads_per_sample[MAX_SCRNA_SAMPLE_NUMBER];

	srInt_64 hiconf_map;
	srInt_64 loconf_map;
	srInt_64 all_processed_reads_before_chunk;
	double program_start_time;
	int is_final_voting_run;
	int longest_chro_name;
	int output_binfiles_are_full;
	gene_inputfile_position_t current_circle_start_position, current_circle_end_position;
	int last_written_fragment_number;

	char index_prefix[MAX_FILE_NAME_LENGTH];
	char output_prefix[MAX_FILE_NAME_LENGTH];
	char temp_file_dir[MAX_FILE_NAME_LENGTH];
	char input_dataset_name[MAX_FILE_NAME_LENGTH * MAX_SCRNA_FASTQ_FILES * 3];
	int input_mode;

	int total_index_blocks;
	int current_index_block_number;
	gene_value_index_t * value_index;
	gehash_t * current_index;
	gene_input_t input_dataset;
	cellCounts_lock_t input_dataset_lock;

	char cell_barcode_list_file[MAX_FILE_NAME_LENGTH];
	char bcl_sample_sheet_file[MAX_FILE_NAME_LENGTH];
	int known_cell_barcode_length;
	int is_dual_index;
	HashTable * cell_barcode_head_tail_table;
	ArrayList * cell_barcodes_array;
	HashTable * sample_sheet_table;
	ArrayList * sample_barcode_list;
	ArrayList * sample_id_to_name;
	HashTable * lineno1B_to_sampleno1B_tab;
	FILE * batch_files[CELLBC_BATCH_NUMBER+2];
	cellCounts_lock_t batch_file_locks[CELLBC_BATCH_NUMBER+2];
	HashTable * sample_BAM_writers;

	parallel_gzip_writer_t fastq_unassigned_writer[4];
	cellCounts_lock_t fastq_unassigned_lock;
	pthread_t thread_delete_files;

	int UMI_length;
	int barcode_batched_max_genes;
	int barcode_batched_max_Rbin_len;
	float umi_cutoff;
	int applied_umi_cut[MAX_SCRNA_SAMPLE_NUMBER];
	int do_one_batch_runner_current;
	int report_excluded_barcodes;
	int has_error;
	
	char features_annotation_file[MAX_FILE_NAME_LENGTH];
	char features_annotation_alias_file[MAX_FILE_NAME_LENGTH];
	int  features_annotation_file_type;
	char features_annotation_gene_id_column[MAX_READ_NAME_LEN];
	char features_annotation_feature_type[MAX_READ_NAME_LEN];
	srInt_64 * block_min_start, *block_max_end, *block_end_index;
	gene_offset_t chromosome_table;
	ArrayList * all_features_array;

	HashTable * chromosome_exons_table;
	unsigned char ** gene_name_array;
	HashTable * gene_name_table; 
	char * unistr_buffer_space;
	srInt_64 unistr_buffer_size;
	srInt_64 unistr_buffer_used;
	char * cmd_rebuilt;


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

#define CHROMOSOME_NAME_LENGTH 256 
#define REVERSE_TABLE_BUCKET_LENGTH 131072
#define IMPOSSIBLE_MEMORY_SPACE 0x5CAFEBABE0000000llu



int cellCounts_reduce_Cigar(char * cigar, char * cigarout){
	//LRMtest_move_buff( context, thread_context, iteration_context, thread_context -> dynamic_programming_indel_movement_buf, strlen(thread_context -> dynamic_programming_indel_movement_buf), iteration_context -> read_length);
	int tmpi = -1;
	int ci, nch, repeat_i = 0, old_opt = 0, wcur=0, rlen=0;
	for(ci = 0; ; ci++){
		nch = cigar[ci];
		if(!nch) break;
		if(isdigit(nch)){
			if(tmpi<0) tmpi = 0;
			tmpi = tmpi*10 + (nch-'0');
		}else{
			if(tmpi<0) tmpi = 1;
			if(old_opt != nch && repeat_i>0){
				if(old_opt=='M' || old_opt=='S' || old_opt=='I') rlen+=repeat_i;
				wcur += SUBreadSprintf( cigarout + wcur, 11, "%d%c", repeat_i, old_opt );
				repeat_i = 0;
			}
			repeat_i += tmpi;
			tmpi = -1;
			old_opt = nch;
		}
	}
	if(repeat_i>0){
		SUBreadSprintf(cigarout + wcur, 11, "%d%c", repeat_i, old_opt);
		if(old_opt=='M' || old_opt=='S' || old_opt=='I') rlen+=repeat_i;
	}
	return rlen;
}

void cellCounts_cell_barcode_tabel_destroy(void *a){
	if(((a-NULL) & 0xfffffffff0000000llu ) ==IMPOSSIBLE_MEMORY_SPACE )return;
	ArrayListDestroy((ArrayList*)a);
}

int cellCounts_make_barcode_HT_table(cellcounts_global_t * cct_context){
	int xx1,xx2;
	cct_context -> cell_barcode_head_tail_table = StringTableCreate(600000);
	HashTableSetDeallocationFunctions(cct_context -> cell_barcode_head_tail_table, free, cellCounts_cell_barcode_tabel_destroy);

	for(xx1=0;xx1 < cct_context-> cell_barcodes_array -> numOfElements; xx1++){
		char * bc = ArrayListGet(cct_context-> cell_barcodes_array, xx1);
		int bcl =strlen(bc);
		if(cct_context -> known_cell_barcode_length==0) cct_context -> known_cell_barcode_length=bcl;
		if(bcl!=cct_context -> known_cell_barcode_length){
			SUBREADprintf("ERROR: the cell barcode list must contain equal-length strings!\n");
			return 1;
		}

		char bctmp[20];
		HashTablePut(cct_context -> cell_barcode_head_tail_table, strdup(bc), NULL+xx1+IMPOSSIBLE_MEMORY_SPACE);
		for(xx2=0; xx2<2; xx2++){
			bctmp[0] = xx2?'S':'F';
			int xx3;
			for(xx3 = 0; xx3< cct_context -> known_cell_barcode_length/2; xx3++)
				bctmp[xx3+1] = bc[ xx3*2+xx2 ];
			bctmp[bcl/2+1]=0;

			ArrayList * array_of_codes = HashTableGet(cct_context -> cell_barcode_head_tail_table, bctmp);
			if(!array_of_codes){
				array_of_codes = ArrayListCreate(4);
				HashTablePut(cct_context -> cell_barcode_head_tail_table, strdup(bctmp), array_of_codes);
			}
			ArrayListPush(array_of_codes, NULL+xx1);
		}
	}
	return 0;
}


void cellCounts_absoffset_to_posstr(cellcounts_global_t * cct_context, unsigned int pos, char * res);

#define SOFT_CLIPPING_WINDOW_SIZE 5
#define SOFT_CLIPPING_MAX_ERROR   1
#define gvindex_baseno2offset_m(base_number, index, offset_byte, offset_bit)    {offset_byte =  ((base_number) - index -> start_base_offset) >>2; offset_bit = (base_number) % 4 * 2;}
inline int cellCounts_get_index_int(gene_value_index_t * value_index, unsigned int pos){
	int offset_byte , offset_bit;
	gvindex_baseno2offset_m(pos, value_index,  offset_byte , offset_bit);
	return (value_index ->values [offset_byte] >> offset_bit)&3;
}

inline int cellCounts_get_read_int(char * read_bin, int read_offset){
	int read_byte , read_bit ;
	read_byte = read_offset/4;
	read_bit = read_offset%4 *2;
	return (read_bin[read_byte]>>read_bit)&3;
}

// it returns the number of bases to be clipped off.
int cellCounts_find_soft_clipping(cellcounts_global_t * cct_context, int thread_no, char * read_bin, int read_offset, unsigned int mapped_pos, int test_len,  int search_to_tail, int search_center) {
	int base_in_window = 0;
	int added_base_index = 0, removed_base_index = 0;
	int search_start = 0;
	int matched_in_window = SOFT_CLIPPING_WINDOW_SIZE;
	int last_matched_base_index = -1, delta;
	gene_value_index_t * current_value_index = cct_context->value_index;

	if(search_to_tail) {
		if(search_center < 0)
			search_start = 0;
		else if(search_center >= test_len)
			// SHOULD NOT HAPPEN!!!
			search_start = test_len - 1;
		else	search_start = search_center - 1;

		delta = 1;
	}else{
		if(search_center < 0)
			// SHOULD NOT HAPPEN!!!
			search_start = 0;
		else if(search_center >= test_len)
			search_start = test_len - 1;
		else	search_start = search_center + 1;

		delta = -1;
	}

//if(FIXLENstrcmp("GACTGACACATGAGCTGAGAATTTATTTTTTTAAGCAAGATGAAAGGGGCAGCTCTAAGACAGACAGGTCACGGGCTCCCTGATAAGTTTCTGAGCTC", read_text)==0)
//SUBREADprintf("SEARCH_SCLIP to_tail %d ; test_len %d ; mappos %u\n",search_to_tail , test_len, mapped_pos);
	for(added_base_index = search_start; added_base_index >= 0 && added_base_index < test_len; added_base_index += delta) {
		// add the new base
		char reference_base = cellCounts_get_index_int(current_value_index , added_base_index+mapped_pos);
		int added_is_matched = reference_base == cellCounts_get_read_int(read_bin, read_offset+added_base_index);

		matched_in_window += added_is_matched;
		if(added_is_matched)
			last_matched_base_index = added_base_index;

		base_in_window ++;

		if(base_in_window > SOFT_CLIPPING_WINDOW_SIZE){
			removed_base_index = added_base_index - delta * SOFT_CLIPPING_WINDOW_SIZE;
			char removing_ref_base = cellCounts_get_index_int(current_value_index, removed_base_index + mapped_pos );
			matched_in_window -= removing_ref_base == cellCounts_get_read_int(read_bin, removed_base_index + read_offset);
		}else{
			matched_in_window --;
		}

		if(matched_in_window < SOFT_CLIPPING_WINDOW_SIZE - SOFT_CLIPPING_MAX_ERROR){
			// clip, bondary is the last matched base.
			if(search_to_tail){
				if(last_matched_base_index < 0) return test_len - search_start;
				else return test_len - last_matched_base_index - 1;
			}else{
				if(last_matched_base_index >= 0) return last_matched_base_index;
				else return search_start - 1;
			}
		}
	}

	if(last_matched_base_index < 0) return test_len;

	if(search_to_tail){
		if(last_matched_base_index < 0) return test_len - search_start;
		else return test_len - last_matched_base_index - 1;
	}else{
		if(last_matched_base_index >= 0) return last_matched_base_index;
		else return search_start - 1;
	}
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
	{"reportExcludedBarcodes", required_argument ,0,0},

	{"cellBarcodeFile",required_argument, 0,0},
	{"sampleSheetFile",required_argument, 0,0},
	{"reportMultiMappingReads", no_argument ,0,0},

	{"maxDiffToTopVotes", required_argument ,0,0},
	{"maxMismatch", required_argument ,0,0},
	{"subreadsPerRead",required_argument,0,0},
	{"minVotesPerRead",required_argument,0,0},
	{"minMappedLength",required_argument,0,0},
	{"umiCutoff",required_argument,0,0},
	{"reportedAlignmentsPerRead",required_argument,0,0},

	{0,0,0,0}
};


void cellCounts_print_config(cellcounts_global_t * cct_context);
#define _gehash_hash(k) ((unsigned int)(k))
#define _gehash_get_bucket(tab, key)  ( (tab) -> buckets + _gehash_hash( key ) % (tab) -> buckets_number )

void prefill_votes(gehash_t * the_table, temp_votes_per_read_t * pnts, int applied_subreads, unsigned int subread, int offset, int subread_no, int is_negative_strand){
	struct gehash_bucket * current_bucket;
	current_bucket = _gehash_get_bucket (the_table, subread);
	int items = current_bucket -> current_items;
	int my_no = subread_no + applied_subreads * is_negative_strand;
	pnts -> votes[my_no] = 0;
	if(!items) return;

	short *current_keys = current_bucket -> new_item_keys;
	int imin=0, imax=items - 1;
	int last_accepted_index;
	short key = subread / the_table->buckets_number;
	while(1){
		last_accepted_index=(imin+imax)/2;
		short current_key = current_keys[last_accepted_index];
		if(current_key>key) imax = last_accepted_index - 1;
		else if(current_key<key) imin = last_accepted_index + 1;
		else break;

		if(imax<imin) return;
	}

	imax -= imin;
	int start_scan_idx = last_accepted_index, tl_last_accepted_index, prefill_big_search_step = imax/4;
	int stoploc;
	
	for(; prefill_big_search_step > 1 ; prefill_big_search_step /= 3)while(1){
		tl_last_accepted_index = last_accepted_index + prefill_big_search_step;
		if(tl_last_accepted_index >= items || current_keys[tl_last_accepted_index]!=key)break;
		last_accepted_index = tl_last_accepted_index;
	}

	while(1){
		last_accepted_index ++;
		if(last_accepted_index == items|| current_keys[last_accepted_index]!= key){
			stoploc = last_accepted_index ; 
			last_accepted_index = start_scan_idx;
			break;
		}
	}

	prefill_big_search_step = imax /4;
	for(; prefill_big_search_step > 1 ; prefill_big_search_step /= 3)while(1){
		tl_last_accepted_index = last_accepted_index - prefill_big_search_step;
		if(tl_last_accepted_index < imin|| current_keys[tl_last_accepted_index]!=key)break;
		last_accepted_index = tl_last_accepted_index;
	}

	while(1){
		if(last_accepted_index == imin || current_keys[last_accepted_index -1]!= key){
			pnts -> start_location_in_index[my_no] = current_bucket->item_values + last_accepted_index;
			pnts -> votes[my_no] = stoploc - last_accepted_index;
			pnts -> offsets[my_no] = offset;
			break;
		}
		last_accepted_index --;
	}
}

int cellCounts_args_context(cellcounts_global_t * cct_context, int argc, char** argv){
	int c , option_index=0;

	optind = 0;
	opterr = 1;
	optopt = 63;

	int cmd_rebuilt_size = 2000;
	char * cmd_rebuilt = malloc(cmd_rebuilt_size);

	cmd_rebuilt[0]=0;
	for(c = 0; c<argc;c++)
	{
		int needed_buff_len = strlen(cmd_rebuilt) + 100+strlen(argv[c]);
		if(needed_buff_len > cmd_rebuilt_size)
		{
			cmd_rebuilt_size = max(cmd_rebuilt_size *2, needed_buff_len);
			cmd_rebuilt = realloc(cmd_rebuilt, cmd_rebuilt_size);
		}
		SUBreadSprintf(cmd_rebuilt+strlen(cmd_rebuilt), cmd_rebuilt_size - strlen(cmd_rebuilt), "\"%s\" ", argv[c]);
	}

	cct_context -> input_mode = GENE_INPUT_BCL;
	cct_context -> total_threads = 10;
	cct_context -> features_annotation_file_type = FILE_TYPE_RSUBREAD;
	cct_context -> reads_per_chunk = 30000000;
	cct_context -> max_reported_alignments_per_read = 1;
	cct_context -> max_voting_simples = 3;
	cct_context -> max_mismatching_bases_in_reads = 3;
	cct_context -> max_voting_locations = 3;
	cct_context -> max_indel_length = 5;
	cct_context -> max_distinct_top_vote_numbers = 3;
	cct_context -> umi_cutoff = -1;
	cct_context -> max_differential_from_top_vote_number = 2;
	cct_context -> min_votes_per_mapped_read = 3;
	cct_context -> total_subreads_per_read = 10;
	cct_context -> is_BAM_and_FQ_out_generated = 1;
	cct_context -> current_dataset_no = 1;
	cct_context -> min_mapped_length_for_mapped_read = 40;
	cct_context -> cmd_rebuilt = cmd_rebuilt;
	strcpy(cct_context -> temp_file_dir, "./");

	while (1){
		c = getopt_long(argc, argv, "", cellCounts_long_options, &option_index);
		if(c<0 || c==255)break;

		if(strcmp("maxMismatch", cellCounts_long_options[option_index].name)==0){
			cct_context -> max_mismatching_bases_in_reads = min(100, max(0, atoi(optarg)));
		}
		if(strcmp("minMappedLength", cellCounts_long_options[option_index].name)==0){
			cct_context -> min_mapped_length_for_mapped_read = min(MAX_SCRNA_READ_LENGTH, max(-1, atoi(optarg)));
		}
		if(strcmp("minVotesPerRead", cellCounts_long_options[option_index].name)==0){
			cct_context -> min_votes_per_mapped_read = min(64, max(1, atoi(optarg)));
		}
		if(strcmp("subreadsPerRead", cellCounts_long_options[option_index].name)==0){
			cct_context -> total_subreads_per_read = min(SCRNA_SUBREADS_HARD_LIMIT, max(7, atoi(optarg)));
		}
		if(strcmp("reportExcludedBarcodes", cellCounts_long_options[option_index].name)==0){
			cct_context -> report_excluded_barcodes = atoi(optarg);
		}
		if(strcmp("dataset", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> input_dataset_name, optarg, MAX_FILE_NAME_LENGTH * MAX_SCRNA_FASTQ_FILES * 3 -1);
		}
		if(strcmp("maxDiffToTopVotes", cellCounts_long_options[option_index].name)==0){
			cct_context -> max_differential_from_top_vote_number = min(30, max(1, atoi(optarg)));
		}
		if(strcmp("index", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> index_prefix, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("inputMode", cellCounts_long_options[option_index].name)==0){
			if(strcmp("FASTQ", optarg)==0) cct_context -> input_mode = GENE_INPUT_SCRNA_FASTQ;
			if(strcmp("BAM", optarg)==0) cct_context -> input_mode = GENE_INPUT_SCRNA_BAM;
		}
		if(strcmp("output", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> output_prefix, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("reportedAlignmentsPerRead", cellCounts_long_options[option_index].name)==0){
			cct_context -> max_reported_alignments_per_read = min(SCRNA_HIGHEST_REPORTED_ALIGNMENTS, max(1, atoi(optarg)));
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
		if(strcmp("reportMultiMappingReads", cellCounts_long_options[option_index].name)==0){
			cct_context -> report_multi_mapping_reads = 1;
		}
		if(strcmp("geneIdColumn", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> features_annotation_gene_id_column, optarg, MAX_READ_NAME_LEN-1);
		}
		if(strcmp("isGTFannotation", cellCounts_long_options[option_index].name)==0){
			cct_context -> features_annotation_file_type = FILE_TYPE_GTF;
		}
		if(strcmp("cellBarcodeFile", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> cell_barcode_list_file, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("sampleSheetFile", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> bcl_sample_sheet_file, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("umiCutoff", cellCounts_long_options[option_index].name)==0){
			cct_context -> umi_cutoff = atof(optarg);
//			SUBREADprintf("UMI_CUT=%.2f\n", cct_context -> umi_cutoff);
		}
	}

	cct_context -> max_voting_simples = max(cct_context -> max_voting_simples, cct_context -> max_reported_alignments_per_read);
	cct_context -> max_voting_locations = max(cct_context -> max_voting_locations, cct_context -> max_reported_alignments_per_read);

	return 0;
}

void cellCounts_print_config(cellcounts_global_t * cct_context){

	SUBREADputs("        ==========     _____ _    _ ____  _____  ______          _____  ");
	SUBREADputs("        =====         / ____| |  | |  _ \\|  __ \\|  ____|   /\\   |  __ \\ ");
	SUBREADputs("          =====      | (___ | |  | | |_) | |__) | |__     /  \\  | |  | |");
	SUBREADputs("            ====      \\___ \\| |  | |  _ <|  _  /|  __|   / /\\ \\ | |  | |");
	SUBREADputs("              ====    ____) | |__| | |_) | | \\ \\| |____ / ____ \\| |__| |");
	SUBREADputs("        ==========   |_____/ \\____/|____/|_|  \\_\\______/_/    \\_\\_____/");
	SUBREADprintf("       %s\n",SUBREAD_VERSION);
	SUBREADputs("");


	print_in_box(80,1,PRINT_BOX_CENTER,"cellCounts settings");
	print_in_box(80,0,0,"");
	print_in_box(80,0,0,"         Index : %s", cct_context -> index_prefix);
	print_in_box(80,0,0,"    Input mode : %s", cct_context -> input_mode == GENE_INPUT_SCRNA_FASTQ?"FASTQ files":(
	    cct_context -> input_mode == GENE_INPUT_SCRNA_BAM?"BAM files":"Raw BCL files"));
	print_in_box(80,0,0,"");
	print_in_box(80,2,PRINT_BOX_CENTER,"");
	SUBREADputs(  "");
}


int determine_total_index_blocks(cellcounts_global_t * cct_context){
	char tmp_fname[MAX_FILE_NAME_LENGTH+ 30];
	cct_context-> total_index_blocks = 0;
	while(1){
		SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+ 30, "%s.%02d.b.tab", cct_context->index_prefix, cct_context->total_index_blocks);
		if(!does_file_exist(tmp_fname))break;
		cct_context->total_index_blocks ++;
	}
	if(cct_context->total_index_blocks> 1){
		SUBREADprintf("ERROR: cellCounts can only run with one-block index. Please build the index with indexSplit=FALSE.\n");
		return 1;
	}
	return 0;
}

void sheet_convert_ss_to_arr( void * key, void * hashed_obj, HashTable * tab ){
	ArrayList * hashed_arr = hashed_obj ;
	cellcounts_global_t * cct_context = tab->appendix1;
	ArrayListPush(cct_context -> sample_id_to_name, key);
	hashed_arr -> appendix1 = NULL+ cct_context -> sample_id_to_name -> numOfElements; // One-based
					
	srInt_64 xx1;		   
	for(xx1 =0; xx1< hashed_arr -> numOfElements; xx1++){
		char ** push_arr = malloc(sizeof(char*)*4); 
		char ** sbc_lane_sample = ArrayListGet(hashed_arr, xx1);
		srInt_64 lane_sample_int = sbc_lane_sample[0]-(char*)NULL;
			
		ArrayListPush(cct_context -> sample_barcode_list, push_arr);
		push_arr[0] = NULL + lane_sample_int;
		push_arr[1] = NULL + cct_context -> sample_id_to_name -> numOfElements;
		push_arr[2] = sbc_lane_sample[1]; // Sample Barcode
		push_arr[3] = NULL + (sbc_lane_sample[1]!=NULL && strlen(sbc_lane_sample[1])>12);
		int line_no_in_sheet = sbc_lane_sample[2] - (char*)NULL;
		HashTablePut(cct_context -> lineno1B_to_sampleno1B_tab , NULL+line_no_in_sheet, NULL + cct_context -> sample_id_to_name -> numOfElements);
	}
}       


void cellCounts_close_sample_SamBam_writers(void *v){
	void ** vv = v;
	simple_bam_writer * wtr = vv[0];
	simple_bam_close(wtr);

	if(vv[1]){
		parallel_gzip_writer_t* gzfp = vv[1];
		parallel_gzip_writer_close(gzfp);

		gzfp = vv[2];
		parallel_gzip_writer_close(gzfp);

		gzfp = vv[3];
		if(gzfp)parallel_gzip_writer_close(gzfp);

		gzfp = vv[4];
		parallel_gzip_writer_close(gzfp);
	}

	cellCounts_lock_t * gz_lock = vv[5];
	cellCounts_destroy_lock(gz_lock);
	free(gz_lock);
	
	free(vv);
}

void cellCounts_add_simple_writer_header(cellcounts_global_t * cct_context, simple_bam_writer * wtr){
	int x1, outcapa, outsize=0;
	unsigned int last_end = 0;
	int rebuilt_len = strlen(cct_context->cmd_rebuilt) + 200;
	outcapa = rebuilt_len + MAX_CHROMOSOME_NAME_LEN *2 + 200;
	char * outbuff = malloc(outcapa);

	outsize += SUBreadSprintf(outbuff, outcapa, "@HD\tVN:1.0\tSO:%s\n", "coordinate");
	for(x1=0;x1<cct_context->chromosome_table.total_offsets; x1++){
		if(outsize + rebuilt_len >= outcapa -MAX_CHROMOSOME_NAME_LEN -40){
			outcapa *=2;
			outbuff = realloc(outbuff, outcapa);
		}
		char * chro_name = cct_context->chromosome_table.read_names+x1*MAX_CHROMOSOME_NAME_LEN;
		unsigned int this_end = cct_context->chromosome_table.read_offsets[x1];
		int this_size = this_end - last_end;
		last_end = this_end;
		outsize += SUBreadSprintf(outbuff+outsize, outcapa-outsize,"@SQ\tSN:%s\tLN:%d\n",chro_name,this_size);
	}
	// if the command line is longer than 16K the remaining part is omitted.
	outsize += snprintf(outbuff+outsize,1024*32, "@PG\tID:cellCounts\tPN:cellCounts\tVN:%s\tCL:%s", SUBREAD_VERSION, cct_context->cmd_rebuilt);
	outbuff[outsize++]='\n';
	outbuff[outsize++]='\0';
	simple_bam_write(&outsize,4,wtr,0);
	simple_bam_write(outbuff, outsize, wtr, 1);
	free(outbuff);

	simple_bam_write(&cct_context->chromosome_table.total_offsets,4,wtr,0);
	for(x1=0;x1<cct_context->chromosome_table.total_offsets; x1++){
		char * chro_name = cct_context->chromosome_table.read_names+x1*MAX_CHROMOSOME_NAME_LEN;
		unsigned int this_end = cct_context->chromosome_table.read_offsets[x1];
		int this_size = this_end - last_end;
		last_end = this_end;
		int sqname_len = strlen(chro_name)+1;
		simple_bam_write(&sqname_len, 4, wtr, 0);
		simple_bam_write(chro_name, sqname_len, wtr, 0);
		simple_bam_write(&this_size, 4, wtr, 0);
	}
	simple_bam_write("",0,wtr,1);
	wtr -> total_chromosomes = cct_context->chromosome_table.total_offsets;
}

void cellCounts_sample_SamBam_writers_new_files(void *k, void *v, HashTable * tab){
	HashTable * fp_tab = tab -> appendix1;
	cellcounts_global_t * cct_context = tab -> appendix2;
	ArrayList * scRNA_sample_id_to_name = tab -> appendix3;

	char * samplename = k;
	char fname [MAX_FILE_NAME_LENGTH+20], fnamet[MAX_FILE_NAME_LENGTH+20];
	SUBreadSprintf(fname, MAX_FILE_NAME_LENGTH+20, "%s.bam", samplename);
	SUBreadSprintf(fnamet,MAX_FILE_NAME_LENGTH+20, "del4-cC-tmp0-%s.del", samplename);
	simple_bam_writer * wtr = simple_bam_create(fname);
	cellCounts_add_simple_writer_header(cct_context, wtr);
	parallel_gzip_writer_t * gzipR1fq=NULL, * gzipI1fq=NULL, * gzipR2fq=NULL, * gzipI2fq=NULL;

	if(cct_context -> input_mode == GENE_INPUT_BCL || cct_context -> input_mode == GENE_INPUT_SCRNA_BAM){
		gzipR1fq = calloc(sizeof(parallel_gzip_writer_t),1);
		gzipI1fq = calloc(sizeof(parallel_gzip_writer_t),1);
		if(cct_context -> is_dual_index)gzipI2fq = calloc(sizeof(parallel_gzip_writer_t),1);
		gzipR2fq = calloc(sizeof(parallel_gzip_writer_t),1);
		SUBreadSprintf(fname, MAX_FILE_NAME_LENGTH+20, "%s_R1.fastq.gz", samplename);
		parallel_gzip_writer_init(gzipR1fq, fname, cct_context -> total_threads);
		SUBreadSprintf(fname, MAX_FILE_NAME_LENGTH+20, "%s_I1.fastq.gz", samplename);
		parallel_gzip_writer_init(gzipI1fq, fname, cct_context -> total_threads);
		SUBreadSprintf(fname, MAX_FILE_NAME_LENGTH+20, "%s_I2.fastq.gz", samplename);
		if(gzipI2fq)parallel_gzip_writer_init(gzipI2fq, fname, cct_context -> total_threads);
		SUBreadSprintf(fname, MAX_FILE_NAME_LENGTH+20, "%s_R2.fastq.gz", samplename);
		parallel_gzip_writer_init(gzipR2fq, fname, cct_context -> total_threads);
	}

	cellCounts_lock_t * gzfp_lock = malloc(sizeof(cellCounts_lock_t));
	cellCounts_init_lock(gzfp_lock, 0);
	int x1;
	for(x1=0; x1<scRNA_sample_id_to_name -> numOfElements; x1++){
		char * sample_name = ArrayListGet( scRNA_sample_id_to_name, x1 );
		if(strcmp(sample_name, samplename)==0){
			void ** wtrptr = malloc(sizeof(void*)*6);
			wtrptr[0]=wtr;
			wtrptr[1]=gzipR1fq;
			wtrptr[2]=gzipI1fq;
			wtrptr[3]=gzipI2fq;
			wtrptr[4]=gzipR2fq;
			wtrptr[5]=gzfp_lock;
			HashTablePut(fp_tab, NULL+x1+1 , wtrptr);
			break;
		}
	}
}

int cellCounts_load_scRNA_tables(cellcounts_global_t * cct_context){
	int rv = 0;

	cct_context-> cell_barcodes_array = input_BLC_parse_CellBarcodes( cct_context-> cell_barcode_list_file );
	if(NULL == cct_context-> cell_barcodes_array){
		SUBREADprintf("ERROR: cannot find valid cell barcodes from the cell barcode list. Please check the content and the accessibility of the file.\n");
		rv = 1;
	}
	if(!rv){
		rv = cellCounts_make_barcode_HT_table( cct_context );
		if(!rv){
			cct_context-> sample_sheet_table = input_BLC_parse_SampleSheet( cct_context -> bcl_sample_sheet_file);
			if(NULL == cct_context-> sample_sheet_table) rv = 1;
			if(rv==0 && cct_context-> sample_sheet_table -> numOfElements > MAX_SCRNA_SAMPLE_NUMBER){
				SUBREADprintf("ERROR: too many samples in the sample sheet.\n");
				rv = 1;
			}
			if(!rv){
				cct_context -> sample_id_to_name = ArrayListCreate(64);
				cct_context -> lineno1B_to_sampleno1B_tab = HashTableCreate(40);

				cct_context -> sample_sheet_table -> appendix1 = cct_context;
				cct_context -> sample_barcode_list = ArrayListCreate(64);

				ArrayListSetDeallocationFunction(cct_context -> sample_barcode_list, free);
				HashTableIteration(cct_context-> sample_sheet_table, sheet_convert_ss_to_arr);

				if(cct_context -> is_BAM_and_FQ_out_generated){
					cct_context -> sample_BAM_writers = HashTableCreate(cct_context -> sample_sheet_table -> numOfElements);
					HashTableSetDeallocationFunctions(cct_context -> sample_BAM_writers, NULL, cellCounts_close_sample_SamBam_writers);
					cct_context -> sample_sheet_table ->appendix1 = cct_context -> sample_BAM_writers;
					cct_context -> sample_sheet_table ->appendix2 = cct_context;
					cct_context -> sample_sheet_table ->appendix3 = cct_context -> sample_id_to_name;
					HashTableIteration( cct_context -> sample_sheet_table, cellCounts_sample_SamBam_writers_new_files);
				}

			}
		}
	}
	return rv;
}

int cellCounts_load_base_value_indexes(cellcounts_global_t * cct_context){
	int rv=0;
	char tmp_fname[MAX_FILE_NAME_LENGTH+ 30];
	SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+30, "%s.%02d.b.array", cct_context ->index_prefix, 0);
	cct_context -> value_index = calloc(sizeof(gene_value_index_t),1);
	rv = rv || gvindex_load(cct_context -> value_index, tmp_fname);
	return rv;
}

typedef struct {
	srInt_64 feature_name_pos;
	unsigned int start;
	unsigned int end;
	unsigned int sorted_order;

	unsigned short chro_name_pos_delta;
	char is_negative_strand;
	char * extra_columns;
} fc_feature_info_t;

typedef struct {
	unsigned int chro_number;
	unsigned int chro_features;
	unsigned int chro_feature_table_start;
	unsigned int chro_block_table_end;
	unsigned int chro_possible_length;

	unsigned short chro_reverse_table_current_size;
	unsigned int * reverse_table_start_index;
	int reverse_table_start_index_size;
	//unsigned int * reverse_table_end_index;
} fc_chromosome_index_info;

srInt_64 cellCounts_unistr_cpy(cellcounts_global_t * cct_context, char * str, int strl)
{
	srInt_64 ret;
	if(cct_context->unistr_buffer_used + strl >= cct_context->unistr_buffer_size-1)
	{
		if( cct_context->unistr_buffer_size < (1000llu*1000u*1000u*32)) // 32GB
		{
			cct_context -> unistr_buffer_size = cct_context->unistr_buffer_size /2 *3;
			cct_context -> unistr_buffer_space = realloc(cct_context -> unistr_buffer_space, cct_context->unistr_buffer_size);
		}
		else
		{
			SUBREADprintf("Error: exceed memory limit (32GB) for storing feature names.\n");
			return 0xffffffffu;
		}
	}

	strcpy(cct_context -> unistr_buffer_space + cct_context->unistr_buffer_used, str);
	ret = cct_context->unistr_buffer_used;

	cct_context->unistr_buffer_used += strl +1;

	return ret;
}

int features_load_one_line(char * gene_name, char * transcript_name, char * chro_name, unsigned int start, unsigned int end, int is_negative_strand, void * context){
	cellcounts_global_t * cct_context = context;
	ArrayList * the_features = cct_context -> all_features_array;
	fc_feature_info_t * new_added = calloc(sizeof(fc_feature_info_t), 1);

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

	// for the featureCounts part.
	cct_context -> longest_chro_name = max(cct_context -> longest_chro_name, strlen(chro_name));
	new_added -> feature_name_pos = cellCounts_unistr_cpy(cct_context, gene_name, strlen(gene_name));
	new_added -> chro_name_pos_delta = cellCounts_unistr_cpy(cct_context, chro_name, strlen(chro_name)) - new_added -> feature_name_pos;
	new_added -> start = start;
	new_added -> end = end;
	new_added -> sorted_order = the_features -> numOfElements;
	new_added -> is_negative_strand = is_negative_strand;
	ArrayListPush(the_features, new_added);

	fc_chromosome_index_info * chro_stub = HashTableGet(cct_context -> chromosome_exons_table, chro_name);
	if(chro_stub){
		if(chro_stub -> chro_possible_length < end+1) chro_stub -> chro_possible_length = end+1;
	}else{
		chro_stub = calloc(sizeof(fc_chromosome_index_info),1);
		char * tmp_chro_name = malloc(CHROMOSOME_NAME_LENGTH);
		term_strncpy(tmp_chro_name, chro_name, CHROMOSOME_NAME_LENGTH);
		chro_stub -> chro_number = cct_context -> chromosome_exons_table -> numOfElements;
		chro_stub -> chro_possible_length = end+1;
		chro_stub -> reverse_table_start_index_size = 0;
		chro_stub -> reverse_table_start_index = NULL;
		HashTablePut(cct_context -> chromosome_exons_table, tmp_chro_name, chro_stub);
	}
	chro_stub -> chro_features ++;

	if(chro_stub -> reverse_table_start_index){
		int bin_location = start / REVERSE_TABLE_BUCKET_LENGTH;
		chro_stub -> reverse_table_start_index[bin_location]++;
	}

	// for the align part.
	unsigned int exonic_map_start = linear_gene_position(&cct_context->chromosome_table , chro_name, start);
	unsigned int exonic_map_stop = linear_gene_position(&cct_context->chromosome_table , chro_name, end), exonpos_i;
	if(exonic_map_start > 0xffffff00 || exonic_map_stop > 0xffffff00){
		return -1;
	}

	for(exonpos_i = exonic_map_start; exonpos_i <= exonic_map_stop; exonpos_i++){
		int exonic_map_byte = exonpos_i/ 8;
		int exonic_map_bit = exonpos_i % 8;

		//if(strcmp("IGKV3-31", gene_name)==0)SUBREADprintf("LINE1ADD %s of %d: %d > %d   (%s:%d ~ %d)\n", gene_name, exonpos_i, exonic_map_byte , exonic_map_bit, chro_name, start, end );
		cct_context ->exonic_region_bitmap[exonic_map_byte] |= (1<<exonic_map_bit);
	}

	for(exonpos_i = exonic_map_start -100; exonpos_i <= exonic_map_stop +100; exonpos_i++){
		int exonic_map_byte = exonpos_i / 8 + (4096 / 8)*1024*1024;
		int exonic_map_bit =  exonpos_i % 8;

		//if(strcmp("ENSMUSG00000024608", gene_name)==0)SUBREADprintf("LINE1ADD %s: %d > %d\n", gene_name, exonic_map_byte , exonic_map_bit );
		cct_context ->exonic_region_bitmap[exonic_map_byte] |= (1<<exonic_map_bit);
	}

	return 0;
}

int cellCounts_find_or_insert_gene_name(cellcounts_global_t * cct_context, unsigned char * feature_name) {
	HashTable * genetable = cct_context -> gene_name_table;

	srInt_64 gene_number = HashTableGet(genetable, feature_name) - NULL;
	if(gene_number>0)
		return gene_number-1;
	else {
		gene_number = genetable -> numOfElements; 
		HashTablePut(genetable, feature_name, NULL+gene_number+1);
		cct_context -> gene_name_array[gene_number] = feature_name;
			// real memory space of feature_name is in the "loaded_features" data structure.
			// now we only save its pointer.

		return gene_number;
	}
}

void cellCounts_register_reverse_table(int block_no, srInt_64 this_block_min_start, srInt_64 this_block_max_end, fc_chromosome_index_info * chro_inf) {
	unsigned int reversed_bucket_start = this_block_min_start /  REVERSE_TABLE_BUCKET_LENGTH;
	unsigned int reversed_bucket_end = this_block_max_end / REVERSE_TABLE_BUCKET_LENGTH;
	assert(this_block_min_start <= this_block_max_end);
	assert(reversed_bucket_end < chro_inf -> chro_possible_length);
	if(!chro_inf->reverse_table_start_index)return;
	int x1;
	for(x1 = reversed_bucket_start; x1 <= reversed_bucket_end; x1++)
		chro_inf->reverse_table_start_index[x1] = min(chro_inf->reverse_table_start_index[x1], block_no);
}

void cellCounts_feature_merge(void * arrv, int start, int items, int items2)
{

	void ** arr = (void **) arrv;

	srInt_64 * ret_start = (srInt_64 *) arr[0];
	srInt_64 * ret_end = (srInt_64 *) arr[1];
	unsigned char * ret_strand = (unsigned char *) arr[2];
	int * ret_entyrez = (int *) arr[3];
	fc_feature_info_t ** old_info_ptr = (fc_feature_info_t **) arr[4];

	int total_items = items+items2;
	srInt_64 * tmp_start = malloc(sizeof(srInt_64) * total_items);
	srInt_64 * tmp_end = malloc(sizeof(srInt_64) * total_items);
	unsigned char * tmp_strand = malloc(sizeof(char) * total_items);
	int * tmp_entyrez = malloc(sizeof(int) * total_items);
	fc_feature_info_t ** tmp_info_ptr = malloc(sizeof(fc_feature_info_t*) * total_items);

	int read_1_ptr = start;
	int read_2_ptr = start+items;
	int write_ptr;

	for(write_ptr=0; write_ptr<total_items; write_ptr++)
	{
		if((read_1_ptr >= start+items)||(read_2_ptr < start+total_items && ret_start[read_1_ptr] >= ret_start[read_2_ptr]))
		{
			tmp_start[write_ptr] = ret_start[read_2_ptr];
			tmp_end[write_ptr] = ret_end[read_2_ptr];
			tmp_strand[write_ptr] = ret_strand[read_2_ptr];
			tmp_entyrez[write_ptr] = ret_entyrez[read_2_ptr];
			tmp_info_ptr[write_ptr] = old_info_ptr[read_2_ptr];
			read_2_ptr++;
		}
		else
		{
			tmp_start[write_ptr] = ret_start[read_1_ptr];
			tmp_end[write_ptr] = ret_end[read_1_ptr];
			tmp_strand[write_ptr] = ret_strand[read_1_ptr];
			tmp_entyrez[write_ptr] = ret_entyrez[read_1_ptr];
			tmp_info_ptr[write_ptr] = old_info_ptr[read_1_ptr];
			read_1_ptr++;
		}
	}

	memcpy(ret_start+ start, tmp_start, sizeof(srInt_64) * total_items);
	memcpy(ret_end+ start, tmp_end, sizeof(srInt_64) * total_items);
	memcpy(ret_strand+ start, tmp_strand, sizeof(char) * total_items);
	memcpy(ret_entyrez+ start, tmp_entyrez, sizeof(int) * total_items);
	memcpy(old_info_ptr+ start, tmp_info_ptr, sizeof(fc_feature_info_t*) * total_items);

	free(tmp_start);
	free(tmp_end);
	free(tmp_strand);
	free(tmp_entyrez);
	free(tmp_info_ptr);
}


int cellCounts_feature_sort_compare(void * arrv, int l, int r)
{
	void ** arr = (void **) arrv;
	srInt_64 * ret_start = (srInt_64 *)arr[0];
	srInt_64 ll = ret_start[l];
	srInt_64 rl = ret_start[r];

	if(ll==rl) return 0;
	else if(ll>rl) return 1;
	else return -1;
}

void cellCounts_feature_sort_exchange(void * arrv, int l, int r)
{
	void ** arr = (void **) arrv;
	srInt_64 tmp;
	fc_feature_info_t * tmpptr;

	srInt_64 * ret_start = (srInt_64 *) arr[0];
	srInt_64 * ret_end = (srInt_64 *) arr[1];
	unsigned char * ret_strand = (unsigned char *) arr[2];
	int * ret_entyrez = (int *) arr[3];
	fc_feature_info_t ** old_info_ptr = (fc_feature_info_t **) arr[4];

	
	tmp = ret_start[r];
	ret_start[r]=ret_start[l];
	ret_start[l]=tmp;

	tmp = ret_end[r];
	ret_end[r]=ret_end[l];
	ret_end[l]=tmp;

	tmp = ret_strand[r];
	ret_strand[r]=ret_strand[l];
	ret_strand[l]=tmp;

	tmp = ret_entyrez[r];
	ret_entyrez[r]=ret_entyrez[l];
	ret_entyrez[l]=tmp;

	tmpptr = old_info_ptr[r];
	old_info_ptr[r]=old_info_ptr[l];
	old_info_ptr[l]=tmpptr;

}



void cellCounts_sort_feature_info(cellcounts_global_t * cct_context, unsigned int features, ArrayList * loaded_features, char *** sorted_chr_names, int ** sorted_entrezid, srInt_64 ** sorted_start, srInt_64 ** sorted_end, unsigned char ** sorted_strand, srInt_64 ** block_end_index, srInt_64 ** block_min_start_pos, srInt_64 ** block_max_end_pos) {
	unsigned int chro_pnt;
	unsigned int xk1,xk2;
	int * ret_entrez = malloc(sizeof(int) * features);
	srInt_64 * ret_start = malloc(sizeof(srInt_64) * features);
	srInt_64 * ret_end = malloc(sizeof(srInt_64) * features);
	int current_block_buffer_size = 2000;

	srInt_64 * ret_block_end_index = malloc(sizeof(srInt_64) * current_block_buffer_size);
	srInt_64 * ret_block_min_start = malloc(sizeof(srInt_64) * current_block_buffer_size);
	srInt_64 * ret_block_max_end = malloc(sizeof(srInt_64) * current_block_buffer_size);
	unsigned char * ret_strand = malloc(features);
	char ** ret_char_name = malloc(sizeof(void *) * features);
	fc_feature_info_t ** old_info_ptr = malloc(sizeof(void *) * features);
	unsigned int * chro_feature_ptr = calloc(sizeof(int) , cct_context -> chromosome_exons_table -> numOfElements);
	fc_chromosome_index_info ** tmp_chro_info_ptrs = calloc(sizeof(fc_chromosome_index_info *), cct_context -> chromosome_exons_table -> numOfElements);

	cct_context -> gene_name_array = malloc(sizeof(char *) * features);	// there should be much less identical names.
	cct_context -> gene_name_table = HashTableCreate(5000);
	HashTableSetHashFunction(cct_context -> gene_name_table, HashTableStringHashFunction);
	HashTableSetKeyComparisonFunction(cct_context -> gene_name_table, my_strcmp);

	// init start positions of each chromosome block.
	if(1) {
		KeyValuePair * cursor;
		int bucket;
		unsigned int sum_ptr = 0;
		for(bucket=0; bucket < cct_context -> chromosome_exons_table -> numOfBuckets; bucket++) {
			cursor = cct_context -> chromosome_exons_table -> bucketArray[bucket];
			while (1) {
				if (!cursor) break;
				fc_chromosome_index_info * tmp_chro_inf = cursor -> value;
				cursor = cursor->next;
				assert(tmp_chro_inf -> chro_number <  cct_context -> chromosome_exons_table -> numOfElements);
				chro_feature_ptr [tmp_chro_inf -> chro_number] = tmp_chro_inf -> chro_features;
				tmp_chro_info_ptrs[tmp_chro_inf -> chro_number] = tmp_chro_inf;
			}
		}

		for(xk1 = 0; xk1 < cct_context -> chromosome_exons_table -> numOfElements; xk1++) {
			unsigned int tmpv = chro_feature_ptr[xk1];
			chro_feature_ptr[xk1] = sum_ptr;
			tmp_chro_info_ptrs[xk1] -> chro_feature_table_start = sum_ptr;
			sum_ptr += tmpv;
		}
	}
	int current_block_id = 0, sort_i = 0;

	(*sorted_chr_names) = ret_char_name;
	(*sorted_entrezid) = ret_entrez;
	(*sorted_start) = ret_start;
	(*sorted_end) = ret_end;
	(*sorted_strand) = ret_strand;

	for(chro_pnt=0; chro_pnt < features; chro_pnt++) {
		fc_feature_info_t * cur_feature = ArrayListGet(loaded_features, chro_pnt);
		char * this_chro_name = cct_context -> unistr_buffer_space + cur_feature -> feature_name_pos + cur_feature -> chro_name_pos_delta;
		
		fc_chromosome_index_info * chro_stub = HashTableGet(cct_context -> chromosome_exons_table, this_chro_name);
		int this_chro_number = chro_stub -> chro_number;
		assert( this_chro_number < cct_context -> chromosome_exons_table -> numOfElements );
		unsigned int this_chro_table_ptr = chro_feature_ptr[this_chro_number];

		ret_char_name[this_chro_table_ptr] = this_chro_name;// (char *)cur_feature -> chro;
		ret_entrez[this_chro_table_ptr] = cellCounts_find_or_insert_gene_name(cct_context, (unsigned char *)(cct_context -> unistr_buffer_space + cur_feature -> feature_name_pos));
		ret_start[this_chro_table_ptr] = cur_feature -> start;
		ret_end[this_chro_table_ptr] = cur_feature -> end;
		ret_strand[this_chro_table_ptr] = cur_feature -> is_negative_strand;
		old_info_ptr[this_chro_table_ptr] = cur_feature;

		chro_feature_ptr[this_chro_number]++;
	}

	print_in_box(80,0,0,"Sort the %d genes...", cct_context -> gene_name_table -> numOfElements);
	for(xk1 = 0; xk1 < cct_context -> chromosome_exons_table -> numOfElements; xk1++) {
		fc_chromosome_index_info * tmp_chro_inf = tmp_chro_info_ptrs[xk1];
		int bins_in_chr = ( tmp_chro_inf->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2);
		short * features_per_block_bins = malloc(sizeof(short)*bins_in_chr);

		if(tmp_chro_inf -> reverse_table_start_index){
			for(xk2=0; xk2<bins_in_chr; xk2++)
				features_per_block_bins[xk2] = max(1,min(1000,(int)(0.9999999+sqrt(tmp_chro_inf -> reverse_table_start_index[xk2]))));
			memset(tmp_chro_inf -> reverse_table_start_index, 0xff, sizeof(int) *bins_in_chr);
		}else continue;

		unsigned int this_block_items = 0;
		srInt_64 this_block_min_start = 0x7fffffff, this_block_max_end = 0;
		unsigned int this_chro_tab_end =  tmp_chro_inf -> chro_features + tmp_chro_inf -> chro_feature_table_start;

		void * in_array[5];
		in_array[0] = ret_start + tmp_chro_inf -> chro_feature_table_start; 
		in_array[1] = ret_end + tmp_chro_inf -> chro_feature_table_start; 
		in_array[2] = ret_strand + tmp_chro_inf -> chro_feature_table_start; 
		in_array[3] = ret_entrez + tmp_chro_inf -> chro_feature_table_start; 
		in_array[4] = old_info_ptr + tmp_chro_inf -> chro_feature_table_start; 

		merge_sort(in_array, this_chro_tab_end - tmp_chro_inf -> chro_feature_table_start, cellCounts_feature_sort_compare, cellCounts_feature_sort_exchange, cellCounts_feature_merge);

		for(sort_i = tmp_chro_inf -> chro_feature_table_start; sort_i< this_chro_tab_end ; sort_i++) {
			old_info_ptr[sort_i]->sorted_order = sort_i;
			int feature_bin_location = ret_start[sort_i] / REVERSE_TABLE_BUCKET_LENGTH;
			int block_bin_location = this_block_min_start / REVERSE_TABLE_BUCKET_LENGTH;

			if(this_block_items && (this_block_items > features_per_block_bins[block_bin_location] || feature_bin_location != block_bin_location)){
				if(current_block_id >= current_block_buffer_size - 1) {
					current_block_buffer_size *= 1.3;
					ret_block_min_start = realloc(ret_block_min_start, sizeof(srInt_64)*current_block_buffer_size);
					ret_block_max_end = realloc(ret_block_max_end, sizeof(srInt_64)*current_block_buffer_size);
					ret_block_end_index = realloc(ret_block_end_index, sizeof(srInt_64)*current_block_buffer_size);
				}
				ret_block_end_index[current_block_id] = sort_i;	// FIRST UNWANTED ID
				ret_block_min_start[current_block_id] = this_block_min_start;
				ret_block_max_end[current_block_id] = this_block_max_end;
				cellCounts_register_reverse_table(current_block_id, this_block_min_start, this_block_max_end, tmp_chro_inf);
				current_block_id++;
				this_block_max_end = 0;
				this_block_items = 0;
				this_block_min_start = 0x7fffffff;
			}

			this_block_max_end = max(this_block_max_end, ret_end[sort_i]);
			this_block_min_start = min(this_block_min_start, ret_start[sort_i]);
			this_block_items ++;
		
		}

		if(this_block_items) {
			if(current_block_id >= current_block_buffer_size) {
				current_block_buffer_size *= 1.3;
				ret_block_min_start = realloc(ret_block_min_start, sizeof(srInt_64)*current_block_buffer_size);
				ret_block_max_end = realloc(ret_block_max_end, sizeof(srInt_64)*current_block_buffer_size);
				ret_block_end_index = realloc(ret_block_end_index, sizeof(srInt_64)*current_block_buffer_size);
			}

			ret_block_end_index[current_block_id] = this_chro_tab_end;	// FIRST UNWANTED ID
			ret_block_min_start[current_block_id] = this_block_min_start;
			ret_block_max_end[current_block_id] = this_block_max_end;
			cellCounts_register_reverse_table(current_block_id, this_block_min_start, this_block_max_end, tmp_chro_inf);
			current_block_id++;
		}

		tmp_chro_inf -> chro_block_table_end = current_block_id; 
		free(features_per_block_bins);
	}

	(*block_end_index) = ret_block_end_index;
	(*block_min_start_pos) = ret_block_min_start;
	(*block_max_end_pos) = ret_block_max_end;

	free(old_info_ptr);
	free(tmp_chro_info_ptrs);
	free(chro_feature_ptr);
}



int cellCounts_load_annotations(cellcounts_global_t * cct_context){
	int rv = 0;
	unsigned int last_loc = 0;
			
	if(cct_context -> features_annotation_alias_file[0])
		rv = (NULL != (cct_context -> sam_chro_to_anno_chr_alias = load_alias_table(cct_context->features_annotation_alias_file)));
	if(!rv){
		int x1;
		cct_context -> unistr_buffer_size = 1024*1024*2;
		cct_context -> unistr_buffer_space = malloc(cct_context -> unistr_buffer_size);
		cct_context -> chromosome_exons_table = HashTableCreate(163);
		HashTableSetHashFunction(cct_context -> chromosome_exons_table, HashTableStringHashFunction);
		HashTableSetKeyComparisonFunction(cct_context -> chromosome_exons_table, my_strcmp);

		for(x1 = 0; x1 < cct_context -> chromosome_table.total_offsets; x1++){
			fc_chromosome_index_info * chro_stub = calloc(sizeof(fc_chromosome_index_info),1);
			char * tmp_chro_name = malloc(CHROMOSOME_NAME_LENGTH);
			char * seq_name = cct_context -> chromosome_table.read_names + MAX_CHROMOSOME_NAME_LEN*x1;
			term_strncpy(tmp_chro_name, seq_name, CHROMOSOME_NAME_LENGTH);
			chro_stub -> chro_number = HashTableGet(cct_context -> chromosome_table.read_name_to_index, seq_name ) - NULL -1;
			chro_stub -> chro_possible_length = cct_context -> chromosome_table.read_offsets[x1] - last_loc;
			last_loc = cct_context -> chromosome_table.read_offsets[x1];
			chro_stub -> reverse_table_start_index_size = chro_stub -> chro_possible_length + 1024*1024;
			chro_stub -> reverse_table_start_index = calloc( chro_stub -> reverse_table_start_index_size  / REVERSE_TABLE_BUCKET_LENGTH +2, sizeof(int));
			HashTablePut(cct_context -> chromosome_exons_table, tmp_chro_name, chro_stub);
		}
		
		cct_context -> all_features_array = ArrayListCreate(350000);
		ArrayListSetDeallocationFunction(cct_context -> all_features_array, free);
		int loaded_features = load_features_annotation(cct_context->features_annotation_file, cct_context->features_annotation_file_type, cct_context->features_annotation_gene_id_column, NULL, cct_context-> features_annotation_feature_type, cct_context, features_load_one_line);
		if(loaded_features<1) rv = 1;

		if(!rv){
			int anno_index_matched=0;
			ArrayList * annot_chros = HashTableKeys(cct_context -> chromosome_exons_table);
			for(x1=0; x1<annot_chros -> numOfElements; x1++){
				char * t1chro = (char*)ArrayListGet(annot_chros,x1);
				fc_chromosome_index_info * chro_stub = HashTableGet(cct_context -> chromosome_exons_table, t1chro);
				if(chro_stub -> chro_features<1) ArrayListSet(annot_chros,x1,NULL);
			}
			int all_chro_unmatched = warning_array_hash_numbers(annot_chros, cct_context-> chromosome_table.read_name_to_index, & anno_index_matched);
			rv=all_chro_unmatched;
			ArrayListDestroy(annot_chros);

			if(all_chro_unmatched) SUBREADprintf("ERROR: no matched chromosomes/contigs found between reference sequences and gene annotation.\n"); else{
				char tbuf[90];
				char_strftime(tbuf);
				SUBREADprintf("Number of chromosomes/contigs matched between reference sequences and gene annotation is %d.\n\n", anno_index_matched);
				cellCounts_print_config(cct_context);
				print_in_box(80,1,1,"Running (%s, pid=%d)", tbuf, getpid());
				print_in_box(80,0,0,"");
			}
			if(!rv) cellCounts_sort_feature_info(cct_context, loaded_features, cct_context -> all_features_array, &cct_context -> features_sorted_chr, &cct_context -> features_sorted_geneid, &cct_context -> features_sorted_start, &cct_context -> features_sorted_stop, &cct_context -> features_sorted_strand, &cct_context -> block_end_index, &cct_context -> block_min_start, &cct_context -> block_max_end);
		}
	}
	return rv;
}

int cellCounts_open_cellbc_batches(cellcounts_global_t * cct_context){
	int x1;
	for(x1=0;x1<CELLBC_BATCH_NUMBER+2; x1++){
		char fname[MAX_FILE_NAME_LENGTH+200];
		SUBreadSprintf(fname, MAX_FILE_NAME_LENGTH+200,"%s/temp-cellcounts-%06d-%03d.tmpbin",cct_context->temp_file_dir,getpid(), x1);
		cct_context -> batch_files[x1] = fopen(fname,"wb");
		cellCounts_init_lock(cct_context -> batch_file_locks+x1, 0);
	}
	int umfpi;
	if(cct_context->input_mode == GENE_INPUT_BCL)for(umfpi=1; umfpi<=4; umfpi++){
		char fname [MAX_FILE_NAME_LENGTH+20];
		char * ftype = "R1";
		if(3==umfpi && !cct_context->is_dual_index) continue;
		if(2==umfpi) ftype = "I1";
		if(3==umfpi) ftype = "I2";
		if(4==umfpi) ftype = "R2";
		SUBreadSprintf(fname, MAX_FILE_NAME_LENGTH+20, "UnassignedReads%03d_%s.fastq.gz", cct_context ->current_dataset_no, ftype);
		parallel_gzip_writer_init(cct_context-> fastq_unassigned_writer+(umfpi-1), fname, cct_context->total_threads);
	}
	cellCounts_init_lock(&cct_context->fastq_unassigned_lock, 0);
	return 0;
}

// we will never use spin-lock.
void cellCounts_init_lock(cellCounts_lock_t * lock, int is_spin){
	pthread_mutex_init(lock, NULL);
}

void cellCounts_destroy_lock(cellCounts_lock_t * lock){
	//if(lock -> is_spin_lock) pthread_spin_destroy(&lock->spinlock);
	//	else
	 pthread_mutex_destroy(lock);
}

/*
int cellCounts_lock_occupy(cellCounts_lock_t * lock){
	//if(lock -> is_spin_lock)
	//	return pthread_spin_lock(&lock->spinlock);
	//else
	return pthread_mutex_lock(lock);
}

int cellCounts_lock_release(cellCounts_lock_t * lock){
	//if(lock -> is_spin_lock) return pthread_spin_unlock(&lock->spinlock);
	//else
	return pthread_mutex_unlock(lock);
}
*/
#define cellCounts_lock_occupy pthread_mutex_lock
#define cellCounts_lock_release pthread_mutex_unlock

int cellCounts_load_context(cellcounts_global_t * cct_context){
	int rv = 0;
	cellCounts_init_lock(&cct_context -> input_dataset_lock, 1 || (cct_context -> input_mode == GENE_INPUT_BCL));

	if(cct_context -> input_mode == GENE_INPUT_BCL){
		rv = rv || geinput_open_bcl(cct_context -> input_dataset_name , & cct_context -> input_dataset , cct_context -> reads_per_chunk, cct_context -> total_threads);
		if(!rv)cct_context -> is_dual_index = cct_context -> input_dataset.bcl_input.is_dual_index;
	} else if(cct_context -> input_mode == GENE_INPUT_SCRNA_FASTQ)
		rv = rv || geinput_open_scRNA_fqs(cct_context -> input_dataset_name , & cct_context -> input_dataset , cct_context -> reads_per_chunk, cct_context -> total_threads);
	else if(cct_context -> input_mode == GENE_INPUT_SCRNA_BAM)
		rv = rv || geinput_open_scRNA_BAM(cct_context -> input_dataset_name , & cct_context -> input_dataset , cct_context -> reads_per_chunk, cct_context -> total_threads);
	else 	rv = -1;

	rv = rv || load_offsets(& cct_context -> chromosome_table, cct_context -> index_prefix);
	rv = rv || determine_total_index_blocks(cct_context);

	int bitmap_size = (4096 / 8)*1024*1024 *2; // the last "*2" is for the extended exon regions
	rv = rv || ((cct_context -> exonic_region_bitmap = calloc(bitmap_size, 1))==NULL);
	rv = rv || cellCounts_load_base_value_indexes(cct_context);
	rv = rv || cellCounts_load_scRNA_tables(cct_context);
	rv = rv || cellCounts_load_annotations(cct_context);
	rv = rv || cellCounts_open_cellbc_batches(cct_context);

	return rv;
}

void * delete_file_thread(void * arg);
int cellCounts_destroy_context(cellcounts_global_t * cct_context){
	int x1;
	pthread_join(cct_context ->thread_delete_files,NULL);
	for(x1=0;x1<CELLBC_BATCH_NUMBER+2; x1++)
		cellCounts_destroy_lock(cct_context -> batch_file_locks+x1);
	cellCounts_destroy_lock(&cct_context -> input_dataset_lock);

	if(cct_context -> is_BAM_and_FQ_out_generated){
		HashTableDestroy(cct_context->sample_BAM_writers);
		cellCounts_destroy_lock(&cct_context->fastq_unassigned_lock);
		int umfpi;
		if(cct_context->input_mode == GENE_INPUT_BCL)
			for(umfpi = 0; umfpi < 4; umfpi++)
				if(umfpi!=2 || cct_context->is_dual_index)
					parallel_gzip_writer_close(cct_context->fastq_unassigned_writer+umfpi);
	}

	geinput_close(&cct_context -> input_dataset);
	destroy_offsets(&cct_context->chromosome_table);
	HashTableDestroy(cct_context->sample_sheet_table);
	HashTableDestroy(cct_context->lineno1B_to_sampleno1B_tab);
	ArrayListDestroy(cct_context->sample_id_to_name);
	ArrayListDestroy(cct_context->sample_barcode_list);
	ArrayListDestroy(cct_context->cell_barcodes_array);
	ArrayListDestroy(cct_context->all_features_array);
	HashTableDestroy(cct_context->gene_name_table);
	HashTableDestroy(cct_context->cell_barcode_head_tail_table);
	HashTableDestroy(cct_context->chromosome_exons_table);
	gvindex_destory(cct_context->value_index);
	free(cct_context -> cmd_rebuilt);
	free(cct_context -> value_index);
	free(cct_context -> exonic_region_bitmap);
	free(cct_context -> features_sorted_chr);
	free(cct_context -> features_sorted_geneid);
	free(cct_context -> features_sorted_start);
	free(cct_context -> features_sorted_stop);
	free(cct_context -> features_sorted_strand);
	free(cct_context -> block_end_index);
	free(cct_context -> block_min_start);
	free(cct_context -> block_max_end);
	free(cct_context -> gene_name_array);
	free(cct_context -> unistr_buffer_space);

	print_in_box(80,0,0,"");
	print_in_box(80,2,0,"");
	SUBREADputs("");
	return 0;
}

void cellCounts_go_chunk_nextchunk(cellcounts_global_t * cct_context){
	if(cct_context -> input_mode == GENE_INPUT_BCL)
		cacheBCL_go_chunk_end(&cct_context -> input_dataset.bcl_input);
	cct_context -> running_processed_reads_in_chunk=0;
}

void cellCounts_clean_context_after_chunk(cellcounts_global_t * cct_context) {
	cct_context -> running_processed_reads_in_chunk = 0;
	cct_context -> processed_reads_in_chunk = 0;
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

typedef struct{
	int thread_id;
	int block_start;
	int block_end;
	HashTable * result_tab;
	int * small_side_ordered_event_ids, * large_side_ordered_event_ids;
	chromosome_event_t * event_space;
	cellcounts_global_t * cct_context;
} AT_context_t;


#define MIN_EVENT_SUPPORT_NO 2

typedef struct{
	char current_cigar_decompress[CORE_MAX_CIGAR_STR_LEN + 1];
	char cigar [CORE_MAX_CIGAR_STR_LEN];

	unsigned short chimeric_sections;
	unsigned int out_poses[CIGAR_PERFECT_SECTIONS];
	short out_lens[CIGAR_PERFECT_SECTIONS];
	char out_cigars[CIGAR_PERFECT_SECTIONS][60];
	char out_strands[CIGAR_PERFECT_SECTIONS];

	char additional_information[CORE_ADDITIONAL_INFO_LENGTH + 1];
	voting_location_t * raw_result;

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
	voting_location_t ** out_raws;
} cellCounts_output_context_t;

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

#define SCORING_MAX_QUALITY_MAPPING  10000000ll

srInt_64 cellCounts_calculate_pos_weight_1sec(cellcounts_global_t * cct_context, unsigned int pos, int len){
	unsigned int pos_i;
	srInt_64 ret = 10ll;
	for(pos_i = pos+1; pos_i <= pos+len ; pos_i ++){ //because exons are 1-based. 
		int exonic_map_byte= pos_i /8;
		int exonic_map_bit = pos_i %8;
		if (cct_context ->exonic_region_bitmap [exonic_map_byte] & (1<<exonic_map_bit))
			return SCORING_MAX_QUALITY_MAPPING;
		exonic_map_byte +=  (4096 / 8)*1024*1024;
		if (cct_context ->exonic_region_bitmap [exonic_map_byte] & (1<<exonic_map_bit))
			ret = 13ll;
	}
	return ret;
}

srInt_64 cellCounts_calculate_pos_weight(cellcounts_global_t * cct_context, unsigned int pos, char * cigar){ 
	int tmpi=0, nch;
	srInt_64 max_weight = 10;
	while(0!=(nch = *(cigar++))){
		if(isdigit(nch)){
			tmpi = tmpi*10+nch-'0';
		}else{
			int toadd_chro = 0;
			if(nch=='M'){
				toadd_chro = tmpi;
				max_weight = max(max_weight, cellCounts_calculate_pos_weight_1sec(cct_context, pos, toadd_chro));
				if(max_weight >= SCORING_MAX_QUALITY_MAPPING)return max_weight;
			} else if(nch == 'D' || nch == 'N' || nch == 'S'){ // At this step, the head-soft-clipped bases were not offset to the mapping position. 
				toadd_chro = tmpi;
			}
			pos += toadd_chro;
			tmpi = 0;
		}
	}

	return max_weight;
}

#define MAX_HIT_NUMBER (1000*1000)

void cellCounts_find_hits_for_mapped_section(cellcounts_global_t * cct_context, int thread_no, char * chro_name, int section_begin_pos, int section_end_pos, int is_fragment_negative_strand, int * nhits){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	int start_reverse_table_index = section_begin_pos / REVERSE_TABLE_BUCKET_LENGTH;
	int end_reverse_table_index = (1+section_end_pos) / REVERSE_TABLE_BUCKET_LENGTH;

	//SUBREADprintf("CHRO:%p(%s)\n", chro_name,chro_name);
	fc_chromosome_index_info * this_chro_info = HashTableGet(cct_context -> chromosome_exons_table, chro_name);
	if(this_chro_info == NULL) {
		if(cct_context -> sam_chro_to_anno_chr_alias) {
			char * anno_chro_name = HashTableGet(cct_context -> sam_chro_to_anno_chr_alias, chro_name);
			if(anno_chro_name) this_chro_info = HashTableGet(cct_context -> chromosome_exons_table, anno_chro_name);
		}
		if(this_chro_info == NULL && memcmp(chro_name, "chr", 3)==0) {
			this_chro_info = HashTableGet(cct_context -> chromosome_exons_table, chro_name+3);
		}
		if(this_chro_info == NULL && strlen(chro_name)<=2) {
			char chro_name_buff[MAX_CHROMOSOME_NAME_LEN+5];
			strcpy(chro_name_buff, "chr");
			strcpy(chro_name_buff+3, chro_name);
			this_chro_info = HashTableGet(cct_context -> chromosome_exons_table, chro_name_buff);
		}
	}

	//SUBREADprintf("CHROINF %s=%p ; PossibleLen = %d\n", chro_name, this_chro_info, this_chro_info?this_chro_info-> chro_possible_length:-1);
	if(this_chro_info) {
		unsigned int search_start, search_end, search_block_id;
		assert(this_chro_info -> reverse_table_start_index);
		start_reverse_table_index = min(start_reverse_table_index, this_chro_info-> chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH);
		end_reverse_table_index = min(end_reverse_table_index, this_chro_info-> chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH+ 1);

		while(start_reverse_table_index<=end_reverse_table_index) {
			search_start = *(this_chro_info -> reverse_table_start_index +start_reverse_table_index);
			if(search_start<0xffffff00)break;
			start_reverse_table_index++;
		}
		//SUBREADprintf("REV_TAB RANGE: %d ~ %d for CHRO %s => SEARCH_START %d ; SEARCH_BLOCK_ID %d\n", start_reverse_table_index, end_reverse_table_index, chro_name, search_start, search_block_id);
		if(search_start>0xffffff00) return;
		search_end = this_chro_info -> chro_block_table_end;

		for(search_block_id=search_start; search_block_id<search_end; search_block_id++){
			//SUBREADprintf("TEST BLOCK %s : %d --  %d ~ %d ?? %d ~ %d\n", chro_name, search_block_id, cct_context -> block_min_start[search_block_id] , cct_context -> block_max_end[search_block_id], section_begin_pos, section_end_pos);
			if (cct_context -> block_min_start[search_block_id] > section_end_pos) break;
			if (cct_context -> block_max_end[search_block_id] < section_begin_pos) continue;

			int search_item_start = 0, search_item_end = cct_context -> block_end_index[search_block_id];
			if(search_block_id>0)search_item_start = cct_context -> block_end_index[search_block_id-1];
			//SUBREADprintf("REACHED SEARCH %s : %d : BLOCK=%d. ITEM START=%d ~ %d\n", chro_name, section_begin_pos, search_block_id, search_item_start, search_item_end);

			// search_item_id is the inner number of the exons.
			int search_item_id;
			for(search_item_id = search_item_start ; search_item_id < search_item_end; search_item_id++) {
				if (cct_context -> features_sorted_stop[search_item_id] >= section_begin_pos) {
					if (cct_context -> features_sorted_start[search_item_id] > section_end_pos) break;
					// there is an overlap >=1 between read and feature.
					// the overlap length is min(end_r, end_F) - max(start_r, start_F) + 1
					
					int is_strand_ok = (is_fragment_negative_strand == cct_context -> features_sorted_strand[search_item_id]);

					if(is_strand_ok){
						if((*nhits) >= thread_context -> hits_number_capacity - 1) {
							thread_context -> hits_number_capacity = max(10, thread_context -> hits_number_capacity*2);
							thread_context -> hits_start_pos = realloc(thread_context -> hits_start_pos , sizeof(int) * thread_context -> hits_number_capacity);
							thread_context -> hits_length = realloc(thread_context -> hits_length, sizeof(short) * thread_context -> hits_number_capacity);
							thread_context -> hits_chro = realloc(thread_context -> hits_chro, sizeof(char *) * thread_context -> hits_number_capacity);
							thread_context -> hits_indices = realloc(thread_context -> hits_indices, sizeof(srInt_64) * thread_context -> hits_number_capacity);
						}

						if((*nhits) <= MAX_HIT_NUMBER - 1) {
							thread_context -> hits_indices[(*nhits)] = search_item_id;
							(*nhits)++;
						} else {
							SUBREADprintf("ERROR: the read overlapped with more than %d features.\n", (*nhits));
							return ;
						}
					}
				} 
			}
		}
	}
}

#define SCRNA_READ_NAME_SPLIT_CHAR '|'

int cellCounts_scan_read_name_str(cellcounts_global_t * cct_context, char * rbin, char * read_name, char ** sample_seq, char ** sample_qual, char ** BC_seq, char ** BC_qual, char ** UMI_seq, char ** UMI_qual, char ** lane_str, char ** RG, int * rname_trimmed_len){
	char * testi;
	int field_i=0;
	if(NULL == read_name && rbin) read_name = rbin + 36;
	for(testi = read_name +1; * testi; testi ++){
		if((*testi)== SCRNA_READ_NAME_SPLIT_CHAR || ((*testi)== ':' && cct_context -> input_mode == GENE_INPUT_BCL )){
			field_i++;
			if(field_i == 1) {
				if(rname_trimmed_len) (*rname_trimmed_len)=testi-read_name;
				if(BC_seq)(*BC_seq) = testi+1;
				if(UMI_seq)(*UMI_seq) = testi+1+cct_context -> known_cell_barcode_length;
			}else if(field_i == 2){
				if(BC_qual)(*BC_qual) = testi+1;
				if(UMI_qual)(*UMI_qual) = testi+1+cct_context -> known_cell_barcode_length;
			}else if(field_i == 3){
				*sample_seq = testi + 1;
				if(RG)(*RG) = *sample_seq;
			}else if(field_i == 4){
				if(sample_qual)(*sample_qual) = testi + 1;
			}else if(field_i == 5){
				(*lane_str) = testi + 1;
				if(memcmp(*lane_str, "@RgLater@", 9)==0) (*lane_str) += 9;
				break;
			}
		}
	}
	if(cct_context -> UMI_length <1){ // no locking is needed because it can be safely done many times.
		int umi_end_pos=0,nch;
		for(umi_end_pos=0; 0!=(nch = (*UMI_seq) [umi_end_pos]); umi_end_pos++) if(!isalpha(nch))break;
		if(umi_end_pos > MAX_UMI_LEN){
			SUBREADprintf("ERROR: the UMI length is abnormaly long (%d bases). This can be caused by an incorrect cell barcode file.\n", umi_end_pos);
		  	umi_end_pos = MAX_UMI_LEN;
			cct_context -> has_error = 1;
		}
		cct_context -> UMI_length = umi_end_pos;
	}

	return field_i;
}


unsigned int cellCounts_explain_read(cellcounts_global_t * cct_context, int thread_no, realignment_result_t * realigns, int read_len, char * read_name , char *read_text, char *qual,  int voting_loc_no, int is_negative_strand);

int cellCounts_get_cellbarcode_no(cellcounts_global_t * cct_context, int thread_no, char * cbc){
	//return -1;
	char tmpc [MAX_READ_NAME_LEN];
	int xx1;
	ArrayList * ret=NULL;

	for(xx1=0;xx1<3;xx1++){
		int xx2;
		if(xx1==1) ret = ArrayListCreate(100);

		if(xx1>0){
			tmpc[0] = (xx1==2)?'S':'F';
			for(xx2=0; xx2<cct_context -> known_cell_barcode_length/2 ; xx2++)
				tmpc[1+xx2] = cbc[2*xx2+xx1-1];
			tmpc[1+cct_context -> known_cell_barcode_length/2]=0;
		}else{
			memcpy(tmpc, cbc, cct_context -> known_cell_barcode_length);
			tmpc[cct_context -> known_cell_barcode_length]=0;
		}

		void *xrawarr = HashTableGet(cct_context -> cell_barcode_head_tail_table, tmpc);

		if(xx1 == 0){
			//if(xrawarr) SUBREADprintf("CAFE ? %p\n", xrawarr);
			srInt_64 xint = xrawarr - NULL;
			if(( xint & 0xFFFFFFFFF0000000llu)== IMPOSSIBLE_MEMORY_SPACE){
				int only_cell_id = xint - IMPOSSIBLE_MEMORY_SPACE;
				// no memory was allocated.
				return only_cell_id;
			}
		}else{
			ArrayList * rawarr = xrawarr;
			if(rawarr){
				int xx3,xx2, found;
				for(xx2=0; xx2<rawarr->numOfElements; xx2++){
					int bcno = ArrayListGet(rawarr, xx2)-NULL;
					found=0;
					for(xx3=0;xx3<ret -> numOfElements;xx3++){
						if(ArrayListGet(ret, xx3)==NULL+bcno){
							found=1;
							break;
						}
					}

					if(!found)ArrayListPush(ret, NULL+bcno);
				}
			}
		}
	}


	int tb1=-1;
	for(xx1=0; xx1<ret -> numOfElements; xx1++){
		int tbcn = ArrayListGet(ret,xx1)-NULL;
		char * known_cbc = ArrayListGet(cct_context -> cell_barcodes_array, tbcn);
		int hc = hamming_dist_ATGC_max2( known_cbc, cbc );

	//	cbc[16]=0; if(hc <=3)SUBREADprintf("TEST_CBC %s ~ %s = %d\n", known_cbc, cbc, hc);
		if(hc==1){
			tb1 = tbcn;
			break;
		}
	}
	//SUBREADprintf("CANDIDATE CELL BARCODES=%ld ; hit = %d\n", ret->numOfElements, tb1);
	ArrayListDestroy(ret);

	return tb1;
}

void cellCounts_build_read_bin(cellcounts_global_t * cct_context, int thread_no, char * rbin, char * read_name, int read_name_len, int trimmed_rname_len, int read_len, char * read_text, char * qual_text, char * chro_name, int chro_pos, int reporting_index, int multi_mapping_number, int this_multi_mapping_i, int editing_dist){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	char * cigar = NULL;
	int mapping_quality = 255;

	int flags = 4;
	if(reporting_index >=0){
		flags = thread_context -> reporting_flags[reporting_index];
		cigar = thread_context -> reporting_cigars[reporting_index];
		mapping_quality = thread_context -> reporting_mapq[reporting_index];
		if(this_multi_mapping_i>1) flags += 256;
	}

	int cigar_opts[CIGAR_PERFECT_SECTIONS], xk1, cover_length = 0;
	int cigar_opt_len = 0;
	if(cigar) cigar_opt_len = SamBam_compress_cigar(cigar, cigar_opts, & cover_length, CIGAR_PERFECT_SECTIONS);

	int record_length = 4 + 4 + 4 + 4 + 4 +  /* l_seq: */ 4 + 4 + 4 + 4 + /* read_name:*/ read_name_len + cigar_opt_len * 4 + (read_len + 1) /2 + read_len;

	int bin = SamBam_reg2bin(chro_pos -1, chro_pos-1+cover_length);
	int refID = -1;
	if(chro_name) refID = HashTableGet(cct_context -> chromosome_table.read_name_to_index, chro_name) - NULL - 1;
	int bin_mq_nl = (bin<<16) | (mapping_quality << 8) | read_name_len ;
	int fag_nc = (flags<<16) | cigar_opt_len;
	int nextRefID = -1, temp_len = 0, next_chro_pos=-1;
	chro_pos--;
	memcpy(rbin + 4 , & refID , 4);
	memcpy(rbin + 8 , & chro_pos , 4);
	memcpy(rbin +12 , & bin_mq_nl , 4);
	memcpy(rbin +16 , & fag_nc , 4);
	memcpy(rbin +20 , & read_len , 4);
	memcpy(rbin +24 , & nextRefID , 4);
	memcpy(rbin +28 , & next_chro_pos , 4);
	memcpy(rbin +32 , & temp_len , 4);
	memcpy(rbin +36 , read_name, read_name_len);
	memcpy(rbin +36 + read_name_len, cigar_opts, cigar_opt_len*4);
	SamBam_read2bin(read_text  , rbin +36 +read_name_len +cigar_opt_len*4);

	int basev = 36 +read_name_len +cigar_opt_len*4+(read_len + 1) /2;
	for(xk1=0; xk1<read_len; xk1++) *(rbin +basev+xk1) = qual_text[xk1]-33;
	
	if(multi_mapping_number>=0){
		rbin[record_length ++]='H';
		rbin[record_length ++]='I';
		rbin[record_length ++]='C';
		rbin[record_length ++]=this_multi_mapping_i;

		rbin[record_length ++]='N';
		rbin[record_length ++]='H';
		rbin[record_length ++]='C';
		rbin[record_length ++]=multi_mapping_number;
	}

	if(editing_dist>=0){
		rbin[record_length ++]='N';
		rbin[record_length ++]='M';
		rbin[record_length ++]='C';
		rbin[record_length ++]=editing_dist;
	}

	record_length -=4;
	memcpy(rbin     , & record_length , 4);
}

void cellCounts_write_one_read_bin(cellcounts_global_t * cct_context, int thread_no, FILE * binfp, int sample_no, int cellbarcode_no, char * umi_barcode, char * readbin, int nhits, int notmapped){
	int x1;
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;

	fwrite(&sample_no,4,1,binfp);
	if(0==notmapped){
		fwrite(&cellbarcode_no,4,1,binfp);
		if(nhits <1){
			srInt_64 zero_genes = 1LLU<<63;
			fwrite(&zero_genes, 8,1,binfp);
		}else if(nhits<2){
			srInt_64 gene_no = thread_context -> hits_indices[0];
			fwrite(&gene_no, 8,1,binfp);
		}else{
			srInt_64 total_genes = nhits + (1LLU<<63);
			fwrite(&total_genes, 8,1,binfp);
			for(x1=0;x1<nhits;x1++){
				srInt_64 gene_no = thread_context -> hits_indices[x1];
				fwrite(&gene_no, 8,1,binfp);
			}
		}
		fwrite(umi_barcode, cct_context->UMI_length, 1, binfp);
	}
	memcpy(&x1, readbin, 4);
	x1+=4;
	fwrite(readbin, x1, 1, binfp);
}

int cellCounts_get_sample_id(cellcounts_global_t * cct_context, char * sbc, int read_laneno){
	int x1;
	for(x1=0; x1 < cct_context -> sample_barcode_list -> numOfElements ; x1++ ){
		char ** lane_and_barcode = ArrayListGet(cct_context -> sample_barcode_list, x1);
		int sheet_lane_no = lane_and_barcode[0]-(char*)NULL;
		if(sheet_lane_no == LANE_FOR_ALL_LANES || read_laneno == sheet_lane_no){
			int sample_no = lane_and_barcode[1]-(char*)NULL;
			char * knownbar = lane_and_barcode[2];
			if(lane_and_barcode[3]){
				int hd = hamming_dist_ATGC_max1_2p( sbc, knownbar );
				if(hd<=2) return sample_no;
			}else{
				int hd = hamming_dist_ATGC_max1( sbc, knownbar );
				if(hd<=1) return sample_no;
			}
		}
	}
	return -1;
}

int cellCounts_parallel_gzip_writer_add_read_fqs_scRNA(parallel_gzip_writer_t**outfps, char * bambin, int thread_no, char * read_text_raw, char * qual_text_raw){ // the text-raw variables are the reads in their input form (not potentially reversed form)
	int reclen=0;
	parallel_gzip_writer_t * outR1fp = outfps[0];
	parallel_gzip_writer_t * outI1fp = outfps[1];
	parallel_gzip_writer_t * outI2fp = outfps[2];
	parallel_gzip_writer_t * outR2fp = outfps[3];

	memcpy(&reclen, bambin,4);
	int flag = 0, l_seq = 0, l_read_name = 0, n_cigar_ops = 0;
	memcpy(&l_read_name, bambin+12,1);
	memcpy(&n_cigar_ops, bambin+16,2);
	memcpy(&flag, bambin+18,2);
	memcpy(&l_seq, bambin+20,4);

	parallel_gzip_writer_add_text(outR2fp,"@",1,thread_no);
	parallel_gzip_writer_add_text(outR1fp,"@",1,thread_no);
	parallel_gzip_writer_add_text(outI1fp,"@",1,thread_no);
	if(outI2fp) parallel_gzip_writer_add_text(outI2fp,"@",1,thread_no);
	char * readname = bambin+36;
	parallel_gzip_writer_add_text(outR1fp,readname, 12,thread_no);
	parallel_gzip_writer_add_text(outR2fp,readname, 12,thread_no);
	parallel_gzip_writer_add_text(outI1fp,readname, 12,thread_no);
	if(outI2fp) parallel_gzip_writer_add_text(outI2fp,readname, 12,thread_no);

	parallel_gzip_writer_add_text(outR1fp,"\n",1,thread_no);
	parallel_gzip_writer_add_text(outR2fp,"\n",1,thread_no);
	parallel_gzip_writer_add_text(outI1fp,"\n",1,thread_no);
	if(outI2fp) parallel_gzip_writer_add_text(outI2fp,"\n",1,thread_no);

	//SUBREADprintf("WRITEFQ RNAME '%s'\n", bambin+36);
	char * R1seq = bambin+36+13;
	int R1len = 0;
	for(R1len=0; R1seq[R1len] && R1seq[R1len]!='|' ;R1len++);
	char * R1qual = R1seq + R1len + 1;
	parallel_gzip_writer_add_text(outR1fp,R1seq, R1len,thread_no);
	parallel_gzip_writer_add_text(outR1fp,"\n+\n",3,thread_no);
	parallel_gzip_writer_add_text_qual(outR1fp,R1qual, R1len,thread_no);
	parallel_gzip_writer_add_text(outR1fp,"\n",1,thread_no);

	char * I1seq = R1qual + R1len + 1;
	int I1I2len = 0;
	for(I1I2len=0; I1seq[I1I2len] && I1seq[I1I2len]!='|' ;I1I2len++);
	int I1len = I1I2len;
	if(outI2fp) I1len /=2;

	char * I2seq = I1seq + I1len;
	char * I1qual = I1seq + I1I2len + 1;
	char * I2qual = I1seq + I1I2len + I1len + 1;

	parallel_gzip_writer_add_text(outI1fp,I1seq, I1len,thread_no);
	parallel_gzip_writer_add_text(outI1fp,"\n+\n",3,thread_no);
	parallel_gzip_writer_add_text_qual(outI1fp,I1qual, I1len,thread_no);
	parallel_gzip_writer_add_text(outI1fp,"\n",1,thread_no);

	if(outI2fp){
		parallel_gzip_writer_add_text(outI2fp,I2seq, I1len,thread_no);
		parallel_gzip_writer_add_text(outI2fp,"\n+\n",3,thread_no);
		parallel_gzip_writer_add_text_qual(outI2fp,I2qual, I1len,thread_no);
		parallel_gzip_writer_add_text(outI2fp,"\n",1,thread_no);
	}

	char * oseq = read_text_raw;
	parallel_gzip_writer_add_text(outR2fp, oseq, l_seq,thread_no);
	parallel_gzip_writer_add_text(outR2fp,"\n+\n",3,thread_no);

	oseq = qual_text_raw;
	parallel_gzip_writer_add_text(outR2fp, oseq, l_seq,thread_no);
	parallel_gzip_writer_add_text(outR2fp,"\n",1,thread_no);
	return 0;
}



void cellCounts_vote_and_add_count(cellcounts_global_t * cct_context, int thread_no, char * read_name, int rlen, char * read_text, char * qual_text, char * raw_text, char * raw_qual, char * chro_name, int chro_pos, int reporting_index, int nhits, int multi_mapping_number, int this_multi_mapping_i, int editing_dist){
	int batch_no;
	
	char * sample_seq=NULL, *sample_qual=NULL, *BC_qual=NULL, *BC_seq=NULL, *UMI_seq=NULL, *UMI_qual=NULL, *lane_str=NULL, *RG=NULL, *testi;
	int rname_trimmed_len=0;
	cellCounts_scan_read_name_str(cct_context, NULL, read_name, &sample_seq, &sample_qual, &BC_seq, &BC_qual, &UMI_seq, &UMI_qual, &lane_str, &RG, &rname_trimmed_len);

	int sample_no = -1;
	if(cct_context -> input_mode == GENE_INPUT_SCRNA_BAM){
		sample_no = 1;  // Only one sample in the BAM mode. A sample may have multiple BAM files but they have the same sample_no.
				// Multiple input samples are mapped/counted in multiple C_cellCounts calls. Each call only does one sample.
	}else if(lane_str){
		int laneno = 0;
		for(testi = lane_str+1; *testi; testi++){
			if(!isdigit(*testi))break;
			laneno = laneno*10 + (*testi)-'0';
		}
		sample_no = cellCounts_get_sample_id(cct_context, sample_seq, laneno); 
	}else if(memcmp("input#", sample_seq,6)==0){
		int lineno = (sample_seq[6]-'0')*1000+(sample_seq[7]-'0')*100+(sample_seq[8]-'0')*10+(sample_seq[9]-'0') +1;
		sample_no = HashTableGet(cct_context -> lineno1B_to_sampleno1B_tab, NULL+lineno)-NULL;
	}else SUBREADprintf("Wrong read name: %s\n", read_name);
	if(0 && FIXLENstrcmp("R00004925364",read_name)==0){
		SUBREADprintf("WRITE-INTO %s:: %s : %d\n", read_name, chro_name, chro_pos) ;
	}
	int cell_barcode_no = cellCounts_get_cellbarcode_no(cct_context, thread_no, BC_seq);

	if(nhits > 1 && !cct_context -> allow_multi_overlapping_reads) nhits = 0;
	if(reporting_index >=0){
		if(cell_barcode_no>=0 && sample_no>0) batch_no = cell_barcode_no % CELLBC_BATCH_NUMBER;
		else if(sample_no>0) batch_no = CELLBC_BATCH_NUMBER;
	}else batch_no = CELLBC_BATCH_NUMBER+1;
	
	char readbin[READ_BIN_BUF_SIZE];
	cellCounts_build_read_bin(cct_context, thread_no, readbin, read_name, strlen(read_name), rname_trimmed_len, rlen, read_text, qual_text, chro_name, chro_pos, reporting_index, multi_mapping_number, this_multi_mapping_i, editing_dist);

	//if(batch_no < CELLBC_BATCH_NUMBER) SUBREADprintf("WRITE_FP_BTCH %d \n", batch_no);
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	if(sample_no>0){
		cellCounts_lock_occupy(cct_context -> batch_file_locks + batch_no);
		FILE * binfp = cct_context -> batch_files [ batch_no ];
		cellCounts_write_one_read_bin(cct_context, thread_no, binfp, sample_no, cell_barcode_no, UMI_seq, readbin, nhits, batch_no == CELLBC_BATCH_NUMBER+1);
		cellCounts_lock_release(cct_context -> batch_file_locks + batch_no);

		if(reporting_index >=0 && this_multi_mapping_i==1){
			thread_context -> mapped_reads_per_sample[ sample_no -1 ]++;
			thread_context -> reads_per_sample[ sample_no -1 ]++;
			if(nhits >0) thread_context -> assigned_reads_per_sample[ sample_no -1 ]++;
		}else thread_context -> reads_per_sample[ sample_no -1 ]++;
	}else thread_context -> reads_per_sample[ cct_context-> sample_sheet_table -> numOfElements ]++; // unassigned reads in NO+1:

	void * pps[6];
	void ** sample_bam_2fps = (void**)pps;
	if(sample_no<=0){
		pps[0]=NULL;
		pps[1]=cct_context -> fastq_unassigned_writer;
		pps[2]=cct_context -> fastq_unassigned_writer+1;
		if(cct_context -> is_dual_index)pps[3]=cct_context -> fastq_unassigned_writer+2;
		else pps[3]=NULL;
		pps[4]=cct_context -> fastq_unassigned_writer+3;
		pps[5]=&cct_context -> fastq_unassigned_lock;
	}else sample_bam_2fps = HashTableGet(cct_context -> sample_BAM_writers, NULL+(sample_no-1) + 1); // sample_id-1: 0,1,2,...

	if(GENE_INPUT_SCRNA_FASTQ != cct_context -> input_mode){
		parallel_gzip_writer_t **gz3fps = (parallel_gzip_writer_t **)sample_bam_2fps+1;
		cellCounts_parallel_gzip_writer_add_read_fqs_scRNA(gz3fps, readbin, thread_no, raw_text, raw_qual);
		// I2 always has the same length as I1, hence no test is required.
		if( gz3fps[0]-> thread_objs[thread_no].in_buffer_used >= PARALLEL_GZIP_TXT_BUFFER_SIZE - PARALLEL_GZIP_TXT_BUFFER_MARGIN ||
		    gz3fps[1]-> thread_objs[thread_no].in_buffer_used >= PARALLEL_GZIP_TXT_BUFFER_SIZE - PARALLEL_GZIP_TXT_BUFFER_MARGIN ||
		    gz3fps[3]-> thread_objs[thread_no].in_buffer_used >= PARALLEL_GZIP_TXT_BUFFER_SIZE - PARALLEL_GZIP_TXT_BUFFER_MARGIN ){
			parallel_gzip_zip_texts(gz3fps[0], thread_no, 0);
			parallel_gzip_zip_texts(gz3fps[1], thread_no, 0);
			if(gz3fps[2]) parallel_gzip_zip_texts(gz3fps[2], thread_no, 0);
			parallel_gzip_zip_texts(gz3fps[3], thread_no, 0);
			cellCounts_lock_occupy(sample_bam_2fps[5]);
			parallel_gzip_writer_flush(gz3fps[0], thread_no);
			parallel_gzip_writer_flush(gz3fps[1], thread_no);
			if(gz3fps[2]) parallel_gzip_writer_flush(gz3fps[2], thread_no);
			parallel_gzip_writer_flush(gz3fps[3], thread_no);
			cellCounts_lock_release(sample_bam_2fps[5]);
		}
	}
}

void cellCounts_summarize_entrez_hits(cellcounts_global_t * cct_context, int thread_no, int * nhits){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	if((*nhits) == 0){
		return;
	}else if((*nhits) == 1){
		thread_context -> hits_indices [0] = cct_context -> features_sorted_geneid[thread_context -> hits_indices [0]];
		return;
	}else{
		int scaned_idx, wrt_ptr=0;
		for(scaned_idx = 0; scaned_idx < *nhits; scaned_idx++)
			thread_context -> hits_indices [scaned_idx] = cct_context -> features_sorted_geneid[thread_context -> hits_indices [scaned_idx]];

		for(scaned_idx = 0; scaned_idx < *nhits; scaned_idx++){
			srInt_64 geneid = thread_context -> hits_indices [scaned_idx];
			int k2, found=0;
			for(k2=0; k2<wrt_ptr; k2++) if(thread_context -> hits_indices [k2]==geneid) found=1;
			if(!found) thread_context -> hits_indices [wrt_ptr++] = geneid;
		}
		(*nhits) = wrt_ptr;
	}
}

void cellCounts_write_read_in_batch_bin(cellcounts_global_t * cct_context, int thread_no, int reporting_index, char * read_name, char * read_text, char * qual_text, char * raw_text, char * raw_qual, int rlen){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	char * chro_name = NULL;
	int chro_pos = 0;
	if(reporting_index>=0){
		locate_gene_position(thread_context -> reporting_positions[reporting_index]+1, &cct_context -> chromosome_table, &chro_name, &chro_pos );
		if(!chro_name) reporting_index=-1;
	}
	if(reporting_index>=0){
		int is_negative = (thread_context -> reporting_flags[reporting_index] & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0;
		chro_pos += get_soft_clipping_length(thread_context -> reporting_cigars[reporting_index]);

		int nch, nchi, add_curs, tmpi=0, chro_curs=chro_pos, nhits=0;
		for(nchi=0; 0!=(nch= thread_context -> reporting_cigars[reporting_index][nchi]); nchi++){
			if(isdigit(nch)){
				tmpi = 10*tmpi + (nch-'0');
			}else{
				add_curs = 0;
				if(nch=='M'||nch=='D'||nch=='N') add_curs =tmpi;
				if(nch=='M')cellCounts_find_hits_for_mapped_section(cct_context, thread_no, chro_name, chro_curs, chro_curs+add_curs, is_negative,&nhits);
				tmpi=0;
			}
		}
		cellCounts_summarize_entrez_hits(cct_context, thread_no, &nhits);
		cellCounts_vote_and_add_count(cct_context, thread_no, read_name, rlen, read_text, qual_text, raw_text, raw_qual, chro_name, chro_pos, reporting_index, nhits, thread_context -> reporting_multi_alignment_no, thread_context -> reporting_this_alignment_no +1, thread_context -> reporting_editing_distance[reporting_index]);
	}else //unmapped
		cellCounts_vote_and_add_count(cct_context, thread_no, read_name, rlen, read_text, qual_text, raw_text, raw_qual, NULL, 0, -1, 0, 0, 0, -1);
}

int cellCounts_add_repeated_buffer(cellcounts_global_t * cct_context, unsigned int * repeated_buffer_position, char ** repeated_buffer_cigar, int * repeated_count, realignment_result_t * res1);
int cellCounts_fetch_next_read_pair(cellcounts_global_t * cct_context, int thread_no,int *read_len, char * read_name, char * read_text, char * qual_text, subread_read_number_t * read_no_in_chunk) ;
int cellCounts_do_voting(cellcounts_global_t * cct_context, int thread_no) ;

void * cellCounts_run_in_thread(void * params){
	void ** parameters = (void **)params;
	cellcounts_global_t * cct_context = (cellcounts_global_t * ) parameters[0];
	int thread_no = parameters[1]-NULL;
	int task = parameters[2]-NULL;
	int * ret_value_pointer = (int *)parameters[3];
	free(parameters);

	switch(task) {
		case STEP_VOTING:
			*ret_value_pointer = cellCounts_do_voting(cct_context, thread_no);
		break;

	}
	return NULL;
}

#define MAX_EVENT_NUMBER_INIT 200000

int cellCounts_prepare_context_for_align(cellcounts_global_t * cct_context, int thread_no, int task) {
	return 0;
}

int cellCounts_sort_junction_entry_table(cellcounts_global_t * cct_context);

int cellCounts_run_maybe_threads(cellcounts_global_t * cct_context, int task){
	int ret_value =0;

	int current_thread_no ;
	cellcounts_align_thread_t * thread_contexts = calloc(sizeof(cellcounts_align_thread_t) , cct_context->total_threads);
	int ret_values[64];

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
		cct_context -> loconf_map  +=thread_contexts[current_thread_no].loconf_map;
		cct_context -> hiconf_map  +=thread_contexts[current_thread_no].hiconf_map;

		if(STEP_VOTING == task) cellCounts_free_topKbuff(cct_context, current_thread_no);
		ret_value += *(ret_values + current_thread_no);
		int smpno;
		for(smpno = 0; smpno < cct_context-> sample_sheet_table -> numOfElements; smpno ++){
			cct_context -> mapped_reads_per_sample[smpno] += thread_contexts[current_thread_no].mapped_reads_per_sample[smpno];
			cct_context -> assigned_reads_per_sample[smpno] += thread_contexts[current_thread_no].assigned_reads_per_sample[smpno];
			cct_context -> reads_per_sample[smpno] += thread_contexts[current_thread_no].reads_per_sample[smpno];
		}
		cct_context -> reads_per_sample[smpno] += thread_contexts[current_thread_no].reads_per_sample[smpno]; //  for non-assigned
		if(ret_value)break;
	}
	//SUBREADprintf("HICONF MAPPING (SIMPLE) = %lld, LOWCONF MAPPING (ALL SUBREADS, NOT SIMPLE) = %lld\n", cct_context -> hiconf_map , cct_context -> loconf_map );

	free(thread_contexts);
	return ret_value;
}

typedef struct {
	unsigned int thread_bodytable_number;
	short thread_no;
} concatinating_events_record_t;

int cellCounts_add_repeated_buffer(cellcounts_global_t * cct_context, unsigned int * repeated_buffer_position, char ** repeated_buffer_cigar, int * repeated_count, realignment_result_t * res1){
	int x1, is_repeated = 0;

	char * r1_cigar = res1 -> cigar_string;
	unsigned int r1_pos = res1 -> first_base_position;

	for(x1 = 0; x1 < (*repeated_count); x1 ++)
		if(repeated_buffer_position[x1] == r1_pos && strcmp(repeated_buffer_cigar[x1], r1_cigar) == 0) {
			is_repeated = 1;
			break;
		}

	if((!is_repeated) && (*repeated_count) < MAX_ALIGNMENT_PER_ANCHOR * cct_context -> max_reported_alignments_per_read){
		repeated_buffer_position[*repeated_count] = r1_pos;
		strcpy(repeated_buffer_cigar[*repeated_count], r1_cigar);
		(*repeated_count) +=1;
	}

	return is_repeated;
}

static int is_dual_index, idx_offset, base_offset, sread_len, total_bin_len, rname_tail_pos;
int cellCounts_copy_bin_to_textread(cellcounts_global_t * cct_context, int readlane, unsigned char* readbin, char * read_name, char * seq, char *qual, int * psread_lens, subread_read_number_t rno){
	int  bii;
	if(sread_len < 1){
		is_dual_index = psread_lens[3]>0;
		idx_offset  = psread_lens[0];
		total_bin_len = psread_lens[0]+psread_lens[1]+psread_lens[2]+psread_lens[3];
		base_offset = psread_lens[1] + ( is_dual_index?psread_lens[2]:0) + idx_offset;
		rname_tail_pos = 16 +2*base_offset;
		sread_len = psread_lens[2+is_dual_index];
	}

	//SUBREADprintf("RLENs=%d, idx=%d, base=%d\n", total_bin_len, idx_offset, base_offset);
	#ifdef __MINGW32__
	SUBreadSprintf(read_name, 15, "R%011" PRIu64 ":", rno);
	#else
	SUBreadSprintf(read_name, 15, "R%011llu:", rno);
	#endif

	read_name[13+idx_offset]='|';
	read_name[14+2*idx_offset]='|';
	read_name[15+base_offset+idx_offset]='|';
	SUBreadSprintf(read_name + rname_tail_pos, MAX_READ_NAME_LEN +1- rname_tail_pos , "|@RgLater@L%03d" , readlane);

	for(bii = 0; bii < total_bin_len; bii++){
		unsigned int nch = readbin[bii];
		char nbase, nqual;
		if(nch > 0){
			int nch1 = nch & 0x3;
			nbase = (1413956417) >> (nch1*8); // 1413956417 = 'TGCA'
			nqual=33+(nch >>2);
			if(nqual <='!'){
				nbase='N';
				nqual='#';
			}
		}else{
			nqual='#';
			nbase='N';
		}
		if(nqual >= '/' && bii < base_offset ) nqual++;
		if(bii < idx_offset){
			read_name[13+bii] = nbase;
			read_name[14+idx_offset+bii]= nqual;
		}else if(bii < base_offset ){
			read_name[15+idx_offset+bii] = nbase;
			read_name[16+base_offset+bii]= nqual;
		}else{
			seq[bii - base_offset] = nbase;
			qual[bii - base_offset] = nqual;
		}
	}
	return sread_len;
}

int cellCounts_fetch_next_read_pair(cellcounts_global_t * cct_context, int thread_no,int *read_len, char * read_name, char * read_text, char * qual_text, subread_read_number_t * read_no_in_chunk) {
	int rl1=0;
	subread_read_number_t this_number = -1;
	gene_input_t * ginp1 = &cct_context -> input_dataset;

	if( ginp1 -> file_type == GENE_INPUT_BCL ){
		int * read_lengths = ginp1 -> bcl_input.single_read_lengths;
		cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
		if(thread_context -> bcl_input_local_cached<1){
			cellCounts_lock_occupy(&cct_context -> input_dataset_lock); 
			int new_reads = cacheBCL_next_readbin(&ginp1 -> bcl_input, thread_context -> bcl_input_local_readlane , thread_context -> bcl_input_local_readbin, BCL_READBIN_ITEMS_LOCAL, &thread_context -> bcl_input_local_start_no);
			if(new_reads)
				thread_context -> bcl_input_local_cached = thread_context -> bcl_input_local_filled = new_reads;
			else if(!cct_context -> running_processed_reads_in_chunk)
				cct_context -> running_processed_reads_in_chunk = (ginp1 -> bcl_input.read_no_in_chunk);
			cellCounts_lock_release(&cct_context -> input_dataset_lock); 
		}
		if(thread_context -> bcl_input_local_cached>0){
			int posnumb = thread_context -> bcl_input_local_filled - thread_context -> bcl_input_local_cached;
			this_number = thread_context -> bcl_input_local_start_no + posnumb;
			thread_context -> bcl_input_local_cached --;
			rl1 = cellCounts_copy_bin_to_textread(cct_context, thread_context -> bcl_input_local_readlane [posnumb], (unsigned char*)thread_context -> bcl_input_local_readbin[posnumb],
				read_name, read_text , qual_text, read_lengths, this_number);
		}
		else rl1=0;

	}else{
		cellCounts_lock_occupy(&cct_context -> input_dataset_lock); 
		if(cct_context -> running_processed_reads_in_chunk < cct_context -> reads_per_chunk) {
			rl1 = geinput_next_read_with_lock(ginp1, read_name, read_text , qual_text);

			if(rl1 > 0){
				this_number = cct_context -> running_processed_reads_in_chunk;
				cct_context -> running_processed_reads_in_chunk ++;
			}
		}
		cellCounts_lock_release(&cct_context -> input_dataset_lock); 
	}

	if(rl1>0 && this_number>=0) {
		*read_no_in_chunk = this_number;
		*read_len = rl1;
		read_text[rl1] = qual_text[rl1] = 0;
		return 0;
	} else {
		*read_no_in_chunk = -1;
		*read_len = -1;
		if(rl1 == -2) cct_context -> has_error=1;
		return 1;
	}
}

void cellCounts_update_top_three(cellcounts_global_t * cct_context, int * top_buffer_3i, int new_value){
	if(new_value > top_buffer_3i[cct_context -> max_distinct_top_vote_numbers - 1]){
		int x1;
		for(x1 = 0;x1 < cct_context -> max_distinct_top_vote_numbers ; x1++){
			if(new_value > top_buffer_3i[x1]){
				int x2;
				for(x2 = cct_context -> max_distinct_top_vote_numbers - 1 ; x2 > x1 ; x2 --){
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
#define _test_record_size	if(current_record_number >= current_record_size - 2){\
		current_record_size *= 1.5;\
		records=realloc(records, sizeof(scanning_events_record_t)*current_record_size);\
		if(NULL == records) return -1;\
	}\

#define _add_record	records[current_record_number].scanning_positons = body -> event_small_side;\
	records[current_record_number].thread_bodytable_number = xx1;\
	current_record_number++;\
	records[current_record_number].scanning_positons = body -> event_large_side;\
	records[current_record_number].thread_bodytable_number = xx1;\
	current_record_number++;\

typedef struct {
	unsigned int scanning_positons;
	unsigned int thread_bodytable_number;
} scanning_events_record_t;


int cellCounts_find_new_indels(cellcounts_global_t * cct_context, int thread_no, char * read_name, char * read_text, char * qual_text, int read_len, int voting_loc_no) {
	return 0;
}

int cellCounts_indel_recorder_copy(gene_vote_number_t * alnrec, gene_vote_number_t * votrec, int tolimax, int applied_subreads_per_strand, int * read_pos_move, char * read_name, unsigned int absloc){
	if(0 && FIXLENstrcmp("R00001325252", read_name)==0){
		int toli;
		SUBREADprintf("RAW Indel Record of %s mapped to %u\n", read_name, absloc);
		for(toli = 0; toli < tolimax; toli +=3)
			SUBREADprintf("  %d %d %d\n", votrec[toli], votrec[toli+1], votrec[toli+2]);
	}
	if(tolimax<=3){
		alnrec[0]=votrec[0];
		alnrec[1]=votrec[1];
		alnrec[2]=0;
		alnrec[3]=0;
		return 0;
	}

	char offsets[applied_subreads_per_strand];
	int sri,toli;
	memset(offsets,0x77, sizeof(char)*applied_subreads_per_strand);
	for(toli=0; toli<tolimax; toli+=3)
		for(sri=votrec[toli]; sri<=votrec[toli+1]; sri++) offsets[sri -1] = votrec[toli+2];

	for(sri=0; sri < applied_subreads_per_strand ; sri++){
		if(offsets[sri]!=0x77){
			*read_pos_move=offsets[sri];
			break;
		}
	}

	int offset_cur = 0, subread0 =-1, high_conf_index = 0;
	toli = 0;
	for(sri=0; ; sri++){
		if(offsets[sri]!=0x77){
			offsets[sri] -= (*read_pos_move);
			if(subread0 < 0){
				subread0 = sri+1;
				offset_cur = offsets[sri];
			}
			high_conf_index = offsets[sri];
		}
		if(offset_cur != offsets[sri] && subread0>0){
			alnrec[toli]= subread0;
			alnrec[toli +1]= sri; // sri is the "next section", which is the subread_last + 1
			alnrec[toli +2]= offset_cur;

			offset_cur = offsets[sri];
			subread0 = sri+1;
			toli +=3;
			if(toli >= MAX_INDEL_TOLERANCE*3){
				//SUBREADprintf("TOLIbreak %d\n", toli);
				break;
			}else alnrec[toli+3]=0;
		}
		if(offsets[sri]==0x77) subread0 =-1;
		if(sri == applied_subreads_per_strand-1){
			if(subread0 >0){
				alnrec[toli]= subread0;
				alnrec[toli +1]= sri+1;
				alnrec[toli +2]= offset_cur;
				if( toli +3 < MAX_INDEL_TOLERANCE*3 )alnrec[toli+3]=0;
				toli+=3;
			}
			break;
		}
	}
	if(0 && FIXLENstrcmp("R00000000057", read_name)==0){
		SUBREADprintf("Copying Indel\n");
		for(toli=0; alnrec[toli]; toli+=3){
			SUBREADprintf("  %d ~ %d : indels=%d\n", alnrec[toli], alnrec[toli+1], alnrec[toli+2]);
			if(toli >= MAX_INDEL_TOLERANCE*3)break;
		}
	}
	return high_conf_index;
}

int cellCounts_matchBin_chro(char * read_bin, int base_offset, gene_value_index_t * index, unsigned int pos, int test_len){
	int ret = 0;

	unsigned int idx_offset_byte, idx_offset_bit;
	unsigned int rbin_offset_byte, rbin_offset_bit;
	gvindex_baseno2offset_m(pos, index , idx_offset_byte, idx_offset_bit);
	if(idx_offset_byte >= index-> values_bytes)return 0;
	char idx_intv = index->values [idx_offset_byte];

	rbin_offset_byte = base_offset/4;
	rbin_offset_bit = (base_offset*2)%8;
	char read_intv = read_bin[rbin_offset_byte];
	int read_i;
	for(read_i = 0; read_i < test_len ; read_i ++){
		char tt = (idx_intv >> idx_offset_bit) & 3;
		char tv = (read_intv >> rbin_offset_bit) & 3;
		if(tt == tv)ret++;
		idx_offset_bit+=2;
		if(idx_offset_bit==8){
			idx_offset_byte++;
			if(idx_offset_byte == index-> values_bytes)return 0;
			idx_intv = index->values [idx_offset_byte];
			idx_offset_bit = 0;
		}
		rbin_offset_bit+=2;
		if(rbin_offset_bit==8){
			rbin_offset_byte ++;
			read_intv = read_bin[rbin_offset_byte];
			rbin_offset_bit =0;
		}
	}
	return ret;
}

#define MINM_INVALID_INDEL (-9999999)

int cellCounts_meet_in_the_middle(cellcounts_global_t * cct_context, int thread_no, unsigned int first_half_abs_pos, char * read_bin, int read_bin_base, int gap_on_read_len, int expected_indel_len, char * read_name, int * gap_mismatch){
	int second_half_start_in_read; 
	unsigned short first_half_matched[MAX_SCRNA_READ_LENGTH], second_half_matched[MAX_SCRNA_READ_LENGTH];
	int max_matched_bases = -99999, best_second_start_in_read = MINM_INVALID_INDEL;
	int x1, summ1=0;
	gene_value_index_t * current_value_index = cct_context->value_index;

	for(x1=0; x1<gap_on_read_len; x1++){
		int idx_value = cellCounts_get_index_int(current_value_index, first_half_abs_pos+ x1);
		int read_value = cellCounts_get_read_int(read_bin, read_bin_base+x1);
		first_half_matched[x1] = summ1; // "if second half starts at x1 (included x1), then how many matched in the first half?""
		summ1 += idx_value == read_value;
	}

	summ1 = 0;
	int indel_offset_first = max(0,(-expected_indel_len));
	for(x1=gap_on_read_len-1; x1>=indel_offset_first; x1--){
		int idx_value = cellCounts_get_index_int(current_value_index, first_half_abs_pos+ x1 +expected_indel_len);
		int read_value = cellCounts_get_read_int(read_bin, read_bin_base+x1); // "if second half starts at x1 (included x1), then how many matched in the second half?"
		summ1 += idx_value == read_value;
		second_half_matched[x1] = summ1;
	}

	if(1) for(second_half_start_in_read = indel_offset_first; second_half_start_in_read < gap_on_read_len; second_half_start_in_read++){
		int sum_here = first_half_matched[ second_half_start_in_read - indel_offset_first ] + second_half_matched[second_half_start_in_read];
		if(sum_here > max_matched_bases){
			max_matched_bases = sum_here ;
			best_second_start_in_read = second_half_start_in_read ;
		}
	}else for(second_half_start_in_read = max(0,(-expected_indel_len)); second_half_start_in_read < gap_on_read_len; second_half_start_in_read++){
		int first_half_length = second_half_start_in_read - max(0,(-expected_indel_len));
		int second_half_length = gap_on_read_len - second_half_start_in_read;
		unsigned int second_half_abs_pos = first_half_abs_pos + second_half_start_in_read + expected_indel_len;
		int first_half_matched = cellCounts_matchBin_chro(read_bin, read_bin_base ,  cct_context -> value_index, first_half_abs_pos, first_half_length);
		int second_half_matched = cellCounts_matchBin_chro(read_bin, read_bin_base + second_half_start_in_read,  cct_context -> value_index, second_half_abs_pos, second_half_length);
		int both_matched = first_half_matched + second_half_matched;
		if(both_matched>max_matched_bases){
			max_matched_bases = both_matched;
			best_second_start_in_read = second_half_start_in_read;
		}
	}
	(* gap_mismatch) = gap_on_read_len - max_matched_bases + min(0, expected_indel_len);
	if(0)SUBREADprintf("FOUND %d indel MEET at %d in %d gap ; %d matched,\n", expected_indel_len, best_second_start_in_read, gap_on_read_len, max_matched_bases );
	return best_second_start_in_read + min(0, expected_indel_len);
}

srInt_64 cellCounts_test_score(cellcounts_global_t * cct_context, int thread_no, char * read_name, int read_len, unsigned int abs_pos, char * cigar, int head_soft_clipped, int tail_soft_clipped, int all_matched_bases, int all_mismatched_bases){
	if(all_mismatched_bases > cct_context -> max_mismatching_bases_in_reads) return 0;
	return all_matched_bases*1000000llu / (1llu+all_mismatched_bases);
}

#define REVERSED_READ_BIN_OFFSET ( MAX_SCRNA_READ_LENGTH /4+1 )

// do: 
//   1, indel detection (meet-in-the-middle or Smith-Waterman)
//   2, build CIGAR
//   3, calculate matched/mismatched
//   4, calculate and save scores in array
srInt_64 cellCounts_explain_one_read(cellcounts_global_t * cct_context, int thread_no,char * read_name, char * read_bin, char * read_text, int read_len,  gene_vote_number_t all_subreads, gene_sc_vote_t * votetab, int vote_i, int vote_j){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	gene_vote_number_t indel_offsets [MAX_INDEL_TOLERANCE*3];
	int read_pos_move = 0, toli, in_cigar_readlen = 0, all_mismatched_bases = 0, all_matched_bases = 0, all_mapped_bases = 0;
	char newcigar[30];

	int rbin_offset_for_reversed = (votetab -> masks[vote_i][vote_j] & IS_NEGATIVE_STRAND)? REVERSED_READ_BIN_OFFSET :0;

	unsigned int abs_pos = votetab -> pos[vote_i][vote_j];
	int tolimax = votetab -> toli[vote_i][vote_j], all_indel_length=0;
	cellCounts_indel_recorder_copy(indel_offsets, votetab -> indel_recorder[vote_i][vote_j], tolimax, all_subreads, &read_pos_move, read_name, abs_pos);
	abs_pos += read_pos_move;
	int last_indel = indel_offsets[2], last_section_subread_no = indel_offsets[1]-1, head_soft_clipped = -1, last_mapped_base_in_read = 0;
	if(0 && FIXLENstrcmp("R00001325252", read_name)==0)SUBREADprintf("\n%s\n%s\n000V THREE_TOLI %d %d %d\n", read_name, read_text, indel_offsets[0], indel_offsets[1], indel_offsets[2]);
	thread_context -> reporting_cigars[thread_context -> reporting_count][0]=0;
	for(toli = 3; toli < tolimax; toli+=3){
		if(indel_offsets[toli]==0)break;
		int hiconf_vote_first = indel_offsets[toli]-1;
		int hiconf_vote_last = indel_offsets[toli+1]-1;
		int indel_offset = indel_offsets[toli+2];
		int indel_diff = indel_offset - last_indel;
		if(abs(indel_diff)>= cct_context -> max_indel_length) continue;
		int last_correct_base = find_subread_end(read_len, all_subreads , last_section_subread_no) - 8;
		int first_correct_base = find_subread_end(read_len, all_subreads , hiconf_vote_first) - 16 + 8;
		if(0 && FIXLENstrcmp("R00001325252", read_name)==0)
			SUBREADprintf("TEST_2POSES LCB=%d [sr %d] FCB=%d [sr %d] INCIGAR=%d\n", last_correct_base, last_section_subread_no, first_correct_base, hiconf_vote_first, in_cigar_readlen);

		if(last_correct_base < in_cigar_readlen) last_correct_base= in_cigar_readlen;

		unsigned int meet_start = abs_pos + last_correct_base + last_indel;
		if(head_soft_clipped <0){
			int first_mapped_base_in_read = find_subread_end(read_len, all_subreads , indel_offsets[0]-1);
			head_soft_clipped = cellCounts_find_soft_clipping(cct_context, thread_no, read_bin+rbin_offset_for_reversed, 0, abs_pos /* this can only happen if no indel is in read */, first_mapped_base_in_read , 0, first_mapped_base_in_read);
			if(head_soft_clipped > 0)SUBreadSprintf(thread_context -> reporting_cigars[thread_context -> reporting_count], 11,"%dS", head_soft_clipped );
			if(meet_start < head_soft_clipped + abs_pos ) meet_start= head_soft_clipped + abs_pos;
			if(head_soft_clipped > last_correct_base) last_correct_base = head_soft_clipped;
			in_cigar_readlen = head_soft_clipped;
		}


		int gap_mismatched = 0;
		int indel_pos = cellCounts_meet_in_the_middle(cct_context, thread_no, meet_start, read_bin + rbin_offset_for_reversed, last_correct_base, first_correct_base - last_correct_base, indel_diff, read_name, &gap_mismatched);
		if(indel_pos == MINM_INVALID_INDEL){
			indel_pos = (first_correct_base - last_correct_base)/2;
			all_indel_length += abs(indel_diff);
		}

		if(0 && FIXLENstrcmp("R00001325252", read_name)==0) SUBREADprintf("THREE_TOLI %d %d %d\nTESTING PARAM: LAST_CORR=%d ; FIRST_CORR=%d, GAP=%d. FOUND_INDEL_POS=%d\n", indel_offsets[toli], indel_offsets[toli+1], indel_offsets[toli+2], last_correct_base, first_correct_base, first_correct_base - last_correct_base, indel_pos);

		int section_matched = cellCounts_matchBin_chro(read_bin +rbin_offset_for_reversed , in_cigar_readlen , cct_context -> value_index, abs_pos + in_cigar_readlen + last_indel , last_correct_base - in_cigar_readlen);

		all_mismatched_bases += gap_mismatched + (last_correct_base - in_cigar_readlen - section_matched);
		all_matched_bases += section_matched + first_correct_base - last_correct_base - gap_mismatched + min(0, indel_diff);

		SUBreadSprintf(newcigar, 45, "%dM%dM%d%c%dM", last_correct_base - in_cigar_readlen , indel_pos, abs(indel_diff), indel_diff>0?'D':'I', first_correct_base - last_correct_base - indel_pos + min(0, indel_diff) );
		all_mapped_bases += first_correct_base - in_cigar_readlen + min(0, indel_diff);
		strcat(thread_context -> reporting_cigars[thread_context -> reporting_count], newcigar);
		in_cigar_readlen = first_correct_base;
		last_mapped_base_in_read = find_subread_end(read_len, all_subreads , hiconf_vote_last) - 16 + 9;;
		last_indel = indel_offset;
		last_section_subread_no = hiconf_vote_last;
	}

	if(head_soft_clipped <0){
		int first_mapped_base_in_read = find_subread_end(read_len, all_subreads , indel_offsets[0]-1) - 8;
		head_soft_clipped = cellCounts_find_soft_clipping(cct_context, thread_no, read_bin+rbin_offset_for_reversed,0, abs_pos /* this can only happen if no indel is in read */, first_mapped_base_in_read , 0, first_mapped_base_in_read);

		if(head_soft_clipped > 0)SUBreadSprintf(thread_context -> reporting_cigars[thread_context -> reporting_count], 11,"%dS", head_soft_clipped );
		in_cigar_readlen = head_soft_clipped;
		last_mapped_base_in_read = find_subread_end(read_len, all_subreads, last_section_subread_no) - 16 + 8;
	}

	int tail_soft_clipped = cellCounts_find_soft_clipping(cct_context, thread_no, read_bin+rbin_offset_for_reversed, last_mapped_base_in_read, abs_pos + last_mapped_base_in_read  + last_indel,  read_len - last_mapped_base_in_read , 1, 1);
	int section_matched = cellCounts_matchBin_chro(read_bin +rbin_offset_for_reversed , in_cigar_readlen , cct_context -> value_index, abs_pos + in_cigar_readlen + last_indel, read_len - in_cigar_readlen - tail_soft_clipped);
	all_mismatched_bases += (read_len - in_cigar_readlen - tail_soft_clipped  - section_matched);
	all_matched_bases += section_matched;

	SUBreadSprintf(newcigar, 11, "%dM", read_len - in_cigar_readlen - tail_soft_clipped);
	all_mapped_bases += read_len - in_cigar_readlen - tail_soft_clipped;
	strcat(thread_context -> reporting_cigars[thread_context -> reporting_count], newcigar);
	if(tail_soft_clipped){
		SUBreadSprintf(newcigar, 11, "%dS", tail_soft_clipped);
		strcat(thread_context -> reporting_cigars[thread_context -> reporting_count], newcigar);
	}

	char tmp_new_cigar[MAX_SCRNA_READ_LENGTH+20];
	int rebuilt_rlen = cellCounts_reduce_Cigar(thread_context -> reporting_cigars[thread_context -> reporting_count], tmp_new_cigar );

	strcpy(thread_context -> reporting_cigars[thread_context -> reporting_count], tmp_new_cigar);

	srInt_64 weight = cellCounts_calculate_pos_weight(cct_context, abs_pos, tmp_new_cigar);
	srInt_64 score = 0;
	if(rebuilt_rlen==read_len && all_mapped_bases >= cct_context -> min_mapped_length_for_mapped_read )score=cellCounts_test_score(cct_context, thread_no, read_name, read_len, abs_pos, thread_context -> reporting_cigars[thread_context -> reporting_count], head_soft_clipped, tail_soft_clipped, all_matched_bases, all_mismatched_bases)*weight;

	if(0 && (FIXLENstrcmp("R00000002478", read_name) == 0|| rebuilt_rlen!=read_len)){
		char posstr[100], posstr2[100];
		cellCounts_absoffset_to_posstr(cct_context, abs_pos , posstr);
		cellCounts_absoffset_to_posstr(cct_context, abs_pos  + in_cigar_readlen + last_indel , posstr2);
		SUBREADprintf("\nREAD_EXP %s\n%s\nMAPPED=%s [%s] %u; TESTINGPOS=%s ; CIGAR=%s -> %s / %d-bases\n MAPPED=%d ;  MATCH=%d ; MM=%d\nGOT Score=%lld ; at weight=%lld\n", read_name, read_text, posstr, rbin_offset_for_reversed ?"NEG":"POS", abs_pos , posstr2, thread_context -> reporting_cigars[thread_context -> reporting_count], tmp_new_cigar, rebuilt_rlen, all_mapped_bases, all_matched_bases, all_mismatched_bases , score , weight );
	}

	thread_context -> reporting_scores[thread_context -> reporting_count] = score;
	thread_context -> reporting_positions[thread_context -> reporting_count] = abs_pos;
	thread_context -> reporting_flags[thread_context -> reporting_count] = rbin_offset_for_reversed?SAM_FLAG_REVERSE_STRAND_MATCHED:0;
	thread_context -> reporting_mapq[thread_context -> reporting_count] = 40 - all_mismatched_bases;
	thread_context -> reporting_editing_distance[thread_context -> reporting_count] = all_mismatched_bases + all_indel_length;
	return score;
}

int sort_readscore_compare_LargeFirst(void * vp , int i , int j){
	void ** pp = vp;
	cellcounts_align_thread_t * thread_context = pp[0];
	int * sorting_index = pp[1];
	int idxI = sorting_index [i];
	int idxJ = sorting_index [j];
	if(thread_context -> reporting_scores[idxI] > thread_context -> reporting_scores[idxJ]) return -1; // large number first
	if(thread_context -> reporting_scores[idxI] < thread_context -> reporting_scores[idxJ]) return 1;
	return 0;
}

void sort_readscore_exchange(void * vp , int i , int j){
	void ** pp = vp;
	int * sorting_index = pp[1];
	int idxI = sorting_index [i];
	sorting_index [i] = sorting_index [j];
	sorting_index [j] = idxI ;
}

int cellCounts_select_and_write_alignments(cellcounts_global_t * cct_context, int thread_no, gene_sc_vote_t * votetab, char * read_name, char * read_text, char * read_bin, char * read_qual, int read_len, gene_vote_number_t all_subreads) {
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	int i,j,reverse_text_offset, distinct_vote_number_i;

	thread_context -> reporting_multi_alignment_no = 0;
	thread_context -> reporting_count=0;
	if(votetab -> max_vote >= cct_context -> min_votes_per_mapped_read){
		int top_distinct_vote_numbers[cct_context -> max_distinct_top_vote_numbers];
		memset(top_distinct_vote_numbers, 0 , cct_context -> max_distinct_top_vote_numbers * sizeof(int));

		for (i=0; i<GENE_SCRNA_VOTE_TABLE_SIZE; i++){
			for (j=0; j< votetab->items[i]; j++){
				int vv = votetab -> votes[i][j];
				if(vv>=cct_context -> min_votes_per_mapped_read)cellCounts_update_top_three(cct_context, top_distinct_vote_numbers, vv);
			}
		}

		for(distinct_vote_number_i = 0 ; distinct_vote_number_i < cct_context -> max_distinct_top_vote_numbers; distinct_vote_number_i ++){
			if(thread_context -> reporting_count >= cct_context -> max_voting_simples)break;
			int this_vote_N = top_distinct_vote_numbers[distinct_vote_number_i];
			if(this_vote_N < 1 || (top_distinct_vote_numbers[0] - this_vote_N > cct_context -> max_differential_from_top_vote_number )) break;

			for (i=0; i<GENE_SCRNA_VOTE_TABLE_SIZE; i++){
				if(thread_context -> reporting_count >= cct_context -> max_voting_simples)break;
				for (j=0; j< votetab->items[i]; j++){
					if(thread_context -> reporting_count >= cct_context -> max_voting_simples)break;

					int vv = votetab->votes[i][j];
					if(vv == this_vote_N){
						srInt_64 this_score = cellCounts_explain_one_read(cct_context, thread_no, read_name, read_bin, read_text, read_len, all_subreads, votetab, i, j);
//						SUBREADprintf("READQV %s : VOTE=%d  SCORE=%lld\n", read_name, vv, this_score);
						thread_context -> reporting_count ++;
						if(this_score >0)thread_context -> reporting_multi_alignment_no ++;
					}
				}
			}
		}
		thread_context -> reporting_multi_alignment_no = min(thread_context -> reporting_multi_alignment_no, cct_context -> max_reported_alignments_per_read);
	}

//	SUBREADprintf("READQV %s : VOTE %d , ALN %d\n\n", read_name, votetab->max_vote, thread_context -> reporting_multi_alignment_no);


	if(thread_context -> reporting_multi_alignment_no) {
		int sorting_index [thread_context -> reporting_count];
		for(distinct_vote_number_i = 0 ; distinct_vote_number_i < thread_context -> reporting_count; distinct_vote_number_i ++) sorting_index [distinct_vote_number_i ] = distinct_vote_number_i ;
		void * sorting_ptr[2];
		sorting_ptr[0] = thread_context;
		sorting_ptr[1] = sorting_index;
		//sort : large number first
		quick_sort(sorting_ptr , thread_context -> reporting_count , sort_readscore_compare_LargeFirst, sort_readscore_exchange); // The last many records are 0-score records. Only "reporting_multi_alignment_no" records ahead are worth writting (score > 0).

		for(thread_context -> reporting_this_alignment_no = 0 ; thread_context -> reporting_this_alignment_no < thread_context -> reporting_multi_alignment_no; thread_context -> reporting_this_alignment_no ++){
			int myno = sorting_index[thread_context -> reporting_this_alignment_no ];
			if(thread_context -> reporting_scores[ myno ] < 1)continue;
			if(thread_context -> reporting_this_alignment_no >= cct_context -> max_reported_alignments_per_read) break;
			reverse_text_offset = (thread_context -> reporting_flags[myno] & SAM_FLAG_REVERSE_STRAND_MATCHED)?MAX_SCRNA_READ_LENGTH+1:0;
			if(reverse_text_offset >0 && 0==read_qual[reverse_text_offset]){
				strcpy(read_qual+reverse_text_offset, read_qual);
				reverse_quality(read_qual+reverse_text_offset, read_len);
			}
			if(0 && FIXLENstrcmp("R00000002478", read_name)==0){
				char posstr[100];
				cellCounts_absoffset_to_posstr(cct_context, thread_context -> reporting_positions[myno] +1, posstr);
				SUBREADprintf("\nFINAL WRITTEN %s BY score=%lld\n====================\n\n",posstr,thread_context -> reporting_scores[ myno ]);
			}
			cellCounts_write_read_in_batch_bin(cct_context, thread_no, myno, read_name, read_text + reverse_text_offset, read_qual+reverse_text_offset, read_text , read_qual , read_len);
		}
	} else cellCounts_write_read_in_batch_bin(cct_context, thread_no, -1, read_name, read_text, read_qual , read_text, read_qual, read_len);

	return 0;
}

int cellCounts_genekey2int(char *key) {
	int ret = 0;

	char * keyc = key +16;
	for (; key < keyc; key++) {
		char c1 = *key;
		ret = (ret << 2) | base2int(c1);
	}
	return ret;
}


int cellCounts_process_copy_ptrs_to_votes_compare(void * arrp, int i, int j){
	int * Xptrs = ((temp_votes_per_read_t**)arrp)[1] -> votes;
	int * trying_subread_no = ((int**)arrp)[0];
	int diffv = Xptrs[ trying_subread_no [i] ] - Xptrs[ trying_subread_no [j] ] ;
	return diffv;
}

void cellCounts_process_copy_ptrs_to_votes_exchange(void * arrp, int i, int j){
	int * trying_subread_no = ((int**)arrp)[0];
	int tmpi = trying_subread_no [i];
	trying_subread_no [i] = trying_subread_no [j];
	trying_subread_no [j] = tmpi;
}


#define INDEL_SEGMENT_SIZE 5
#define VOTING_PRIME_NUMBER 66889

//#define _index_vote(key) (((unsigned int)(key))%GENE_SCRNA_VOTE_TABLE_SIZE)
#define _index_vote_tol(key) (((unsigned int)(key)/INDEL_SEGMENT_SIZE)%GENE_SCRNA_VOTE_TABLE_SIZE)
#define cellCounts_voting_update_topK(vt, vn, ij)  { \
	int x3, x3done=0;\
	for(x3= 0; x3 < cct_context -> max_voting_simples +1; x3++){\
		if(ij == vt->topK_IJ[x3]){\
			if(x3 > 0 && vn > vt->topK_votes[x3-1]){\
				int x4;\
				for(x4 = x3; x4< cct_context -> max_voting_simples; x4++){\
					vt->topK_votes[x4] = vt->topK_votes[x4+1];\
					vt->topK_IJ[x4] = vt->topK_IJ[x4+1];\
				}\
				vt->topK_votes[cct_context -> max_voting_simples]=0;\
			}else{\
				vt->topK_votes[x3] = vn;\
				x3done=1;\
			}\
			break;\
		}\
	}\
	if(!x3done)for(x3= 0; x3 < cct_context -> max_voting_simples +1; x3++){\
		if(vn > vt->topK_votes[x3]){\
			int x4;\
			if(vt->topK_votes[x3]>0)for(x4 = cct_context -> max_voting_simples ; x4 > x3; x4--){\
				vt->topK_votes[x4] = vt->topK_votes[x4-1];\
				vt->topK_IJ[x4] = vt->topK_IJ[x4-1];\
			}\
			vt->topK_votes[x3]=vn;\
			vt->topK_IJ[x3]=ij;\
			break; \
		}\
	}\
}


void cellCounts_process_copy_ptrs_to_votes(cellcounts_global_t * cct_context, int thread_no, temp_votes_per_read_t * ptrs, gene_sc_vote_t * vote, int applied_subreads_per_strand, char * read_name){
	int subreads = applied_subreads_per_strand*2;
	int x1, x2, trying_subread_no[subreads];
	for(x1=0; x1<subreads; x1++) trying_subread_no[x1]=x1;
	void * sort_arr[2];
	sort_arr [0] = trying_subread_no;
	sort_arr [1] = ptrs;
	quick_sort(sort_arr , subreads , cellCounts_process_copy_ptrs_to_votes_compare, cellCounts_process_copy_ptrs_to_votes_exchange);

	init_gene_vote(vote);
	int cct_indel_len = cct_context -> max_indel_length, cct_indel_neg = -cct_context -> max_indel_length;
	for(x1=0; x1<subreads; x1++){
		int myno = trying_subread_no[x1];
		int has_votes = ptrs->votes[myno];
		if(has_votes<1) continue;

		int mynoP1PStr=(myno%applied_subreads_per_strand)+1;
		int offset = ptrs->offsets[myno];
		int of_p_16 = offset + 16;
		int is_reversed = (myno >= applied_subreads_per_strand)?IS_NEGATIVE_STRAND:0;
		unsigned int * index_ptr = ptrs->start_location_in_index [myno];
		int vote_prime_sum = ptrs->votes[trying_subread_no[subreads-1]];

		for(x2 = 0 ; x2 < has_votes; x2++){
			unsigned int kv = index_ptr[ vote_prime_sum % has_votes ] - offset;
			vote_prime_sum += VOTING_PRIME_NUMBER;
			int iix, offsetX2, offsetX, datalen, datalen2, found = 0;
			offsetX = offsetX2 = _index_vote_tol(kv);
			datalen = datalen2 = vote -> items[offsetX];
			unsigned int * dat2, *dat;
			dat = dat2 = vote -> pos[offsetX];

			for(iix = 0; iix<=INDEL_SEGMENT_SIZE; iix = iix>0?-iix:(-iix+INDEL_SEGMENT_SIZE)){
				if(iix) {
					offsetX = _index_vote_tol(kv+iix);
					datalen = vote -> items[offsetX];
					if(!datalen)continue;
					dat = vote -> pos[offsetX];
				} else if(!datalen) continue;

				int itemidx;
				for (itemidx=0;itemidx<datalen;itemidx++){
					int dist0 = kv-dat[itemidx];
					if( dist0 >= cct_indel_neg  && dist0 <= cct_indel_len  && is_reversed == vote->masks[offsetX][itemidx]){
						int toli, tolimax=vote -> toli[offsetX][itemidx], known_indel=0;

						gene_vote_number_t * indelrec = vote -> indel_recorder[offsetX][itemidx];
						for(toli = 0; toli < tolimax; toli +=3){
							if( indelrec [toli+2] == dist0 ){
								if(indelrec [toli] > mynoP1PStr) indelrec [toli]  = mynoP1PStr;
								if(indelrec [toli +1] < mynoP1PStr) indelrec [toli +1]  = mynoP1PStr;
								known_indel=1;
								break;
							}
						}
					
						if(tolimax < MAX_INDEL_TOLERANCE*3 && !known_indel){
							indelrec [tolimax] = mynoP1PStr;
							indelrec [tolimax+1] = mynoP1PStr;
							indelrec [tolimax+2] = dist0;
							vote -> toli[offsetX][itemidx]+=3;
						}

						gene_vote_number_t test_max = (vote->votes[offsetX][itemidx]);
						test_max ++;
						vote -> votes[offsetX][itemidx] = test_max;
						if(vote->max_vote < test_max) vote->max_vote = test_max;

						if (offset < vote->coverage_start [offsetX][itemidx])
							vote->coverage_start [offsetX][itemidx] = offset;
						if (of_p_16 > vote->coverage_end [offsetX][itemidx])
							vote->coverage_end [offsetX][itemidx] = of_p_16;

						found = 1;
						break;
					}
				}
				if(found)break;
			}
			if((!found) && datalen2 < GENE_SCRNA_VOTE_SPACE){
				vote -> items[offsetX2] = datalen2+1;
				dat2[datalen2] = kv;
				vote -> masks[offsetX2][datalen2] = is_reversed;
				vote -> votes[offsetX2][datalen2] = 1;
				vote -> indel_recorder[offsetX2][datalen2][0] = mynoP1PStr;
				vote -> indel_recorder[offsetX2][datalen2][1] = mynoP1PStr;
				vote -> indel_recorder[offsetX2][datalen2][2] = 0;
				vote -> toli[offsetX2][datalen2] = 3;
				vote->coverage_start [offsetX2][datalen2] = offset;
				vote->coverage_end [offsetX2][datalen2] = of_p_16;

				if (vote->max_vote==0) vote->max_vote = 1;
			}
		}
	}
}



void simpleMode_cellCounts_process_copy_ptrs_to_votes(cellcounts_global_t * cct_context, int thread_no, temp_votes_per_read_t * ptrs, gene_sc_vote_t * vote, int applied_subreads_per_strand, char * read_name, int simple_mode){
	int subreads = applied_subreads_per_strand*2, space_limit = GENE_SCRNA_VOTE_SPACE - (simple_mode?1:0);
	int x1, x2, trying_subread_no[subreads];
	for(x1=0; x1<subreads; x1++) trying_subread_no[x1]=x1;
	void * sort_arr[2];
	sort_arr [0] = trying_subread_no;
	sort_arr [1] = ptrs;
	quick_sort(sort_arr , subreads , cellCounts_process_copy_ptrs_to_votes_compare, cellCounts_process_copy_ptrs_to_votes_exchange);
	if(0 && FIXLENstrcmp("R00004925364",read_name)==0){
		for(x1=0;x1<subreads ;x1++){
			int myno = trying_subread_no[x1];
			SUBREADprintf("SORTED-VOTES: %s #%d has %d votes ; subread_no=%d\n", read_name, x1, ptrs->votes[myno], myno);
		}
	}

	if(simple_mode){
		init_gene_vote(vote);
		for(x1=0;x1<cct_context -> max_voting_simples +1; x1++){
			vote -> topK_votes[x1]=0;
			vote -> topK_IJ[x1]=0xffffffff;
		}
	}
	int cct_indel_len = cct_context -> max_indel_length, cct_indel_neg = -cct_context -> max_indel_length;
	for(x1=0; x1<subreads; x1++){
		int myno = trying_subread_no[x1];
		int has_votes = ptrs->votes[myno];
		if(has_votes<1) continue;

		int mynoP1PStr=(myno%applied_subreads_per_strand)+1;
		int offset = ptrs->offsets[myno];
		int of_p_16 = offset + 16;
		int is_reversed = (myno >= applied_subreads_per_strand)?IS_NEGATIVE_STRAND:0;
		unsigned int * index_ptr = ptrs->start_location_in_index [myno];
		int vote_prime_sum = ptrs->votes[trying_subread_no[subreads-1]];
		int ignore_creation_voteloc = 0;

		//#warning "======== IGNORED SOME NEW CANDIDATE LOCATIONS BY (1) =================="
		if(0)if(x1 >= subreads - 5 && has_votes >10) {
			if(vote->max_vote >= 4) ignore_creation_voteloc =1;
			else if(has_votes > 20 && vote->max_vote >= 3) ignore_creation_voteloc =1;
		}

		for(x2 = 0 ; x2 < has_votes; x2++){
			unsigned int kv = index_ptr[ vote_prime_sum % has_votes ] - offset;
			vote_prime_sum += VOTING_PRIME_NUMBER;
			int iix, offsetX2, offsetX, datalen, datalen2, found = 0;
			offsetX = offsetX2 = _index_vote_tol(kv);
			datalen = datalen2 = vote -> items[offsetX];
			unsigned int * dat2, *dat;
			dat = dat2 = vote -> pos[offsetX];
			for(iix = 0; iix<=INDEL_SEGMENT_SIZE; iix = iix>0?-iix:(-iix+INDEL_SEGMENT_SIZE)) {
				if(iix) {
					offsetX = _index_vote_tol(kv+iix);
					datalen = vote -> items[offsetX];
					if(!datalen)continue;
					dat = vote -> pos[offsetX];
				} else if(!datalen) continue;

				int itemidx;
				for (itemidx=0;itemidx<datalen;itemidx++){
					int dist0 = kv-dat[itemidx];
					if( dist0 >= cct_indel_neg  && dist0 <= cct_indel_len  && is_reversed == vote->masks[offsetX][itemidx]){
						int toli, tolimax=vote -> toli[offsetX][itemidx], known_indel=0;

						gene_vote_number_t * indelrec = vote -> indel_recorder[offsetX][itemidx];
						for(toli = 0; toli < tolimax; toli +=3){
							if( indelrec [toli+2] == dist0 ){
								if(indelrec [toli] > mynoP1PStr) indelrec [toli]  = mynoP1PStr;
								if(indelrec [toli +1] < mynoP1PStr) indelrec [toli +1]  = mynoP1PStr;
								known_indel=1;
								break;
							}
						}
						
						if(0 && FIXLENstrcmp("R00004925364",read_name)==0){ 
							SUBREADprintf("  ADDING VOTE #%d : dist0=%d , rev=%d , known=%d , toli=%d tolirec = %d %d %d\n", myno, dist0 , is_reversed , known_indel, tolimax, indelrec [tolimax-3], indelrec [tolimax-2], indelrec [tolimax-1]);
						}

						if(tolimax < MAX_INDEL_TOLERANCE*3 && !known_indel){
							indelrec [tolimax] = mynoP1PStr;
							indelrec [tolimax+1] = mynoP1PStr;
							indelrec [tolimax+2] = dist0;
							vote -> toli[offsetX][itemidx]+=3;
						}

						gene_vote_number_t test_max = (vote->votes[offsetX][itemidx]);
						test_max ++;
						vote -> votes[offsetX][itemidx] = test_max;

						if(simple_mode && test_max > vote->topK_votes[cct_context -> max_voting_simples]){
							int myIJ = (offsetX << 16) | itemidx;
							cellCounts_voting_update_topK(vote, test_max , myIJ);
						}

						if(vote->max_vote < test_max) vote->max_vote = test_max;

						if (offset < vote->coverage_start [offsetX][itemidx])
							vote->coverage_start [offsetX][itemidx] = offset;
						if (of_p_16 > vote->coverage_end [offsetX][itemidx])
							vote->coverage_end [offsetX][itemidx] = of_p_16;

						found = 1;
						break;
					}
				}
				if(found)break;
			}
			if ( datalen2< space_limit && (!ignore_creation_voteloc) && (!found)){
				vote -> items[offsetX2] = datalen2+1;
				dat2[datalen2] = kv;
				vote -> masks[offsetX2][datalen2] = is_reversed;
				vote -> votes[offsetX2][datalen2] = 1;
				vote -> indel_recorder[offsetX2][datalen2][0] = mynoP1PStr;
				vote -> indel_recorder[offsetX2][datalen2][1] = mynoP1PStr;
				vote -> indel_recorder[offsetX2][datalen2][2] = 0;
				vote -> toli[offsetX2][datalen2] = 3;
				vote->coverage_start [offsetX2][datalen2] = offset;
				vote->coverage_end [offsetX2][datalen2] = of_p_16;

				if (vote->max_vote==0) vote->max_vote = 1;

				if(simple_mode){
					int x3;
					for(x3=0; x3< cct_context -> max_voting_simples+1;x3++) if(0==vote -> topK_votes[x3]){
						vote -> topK_votes[x3]=1;
						vote -> topK_IJ[x3] = (offsetX2 << 16) | datalen2;
					}
				}
			}
		}
	}
}
int cellCounts_simple_mode_highconf(cellcounts_global_t * cct_context, int thread_no, int applied_subreads, gene_sc_vote_t * vote, char * read_name){
	int x1, vlast = vote -> max_vote;
	for(x1 = 1; x1 < cct_context -> max_voting_simples +1; x1++){
		int vdiff = vlast - vote->topK_votes[x1];
		if(vdiff >= 3) return 1;

		vlast = vote->topK_votes[x1];
	}
	return 0;

/*
	int max_I = vote -> max_vote_IJ >> 16;
	int max_J = vote -> max_vote_IJ & 0xffff;

	if(0){
		int tolimax = vote -> toli[max_I][max_J];
		if(tolimax >3)return 0;
	}

	int SRp1_start = vote -> indel_recorder[max_I][max_J][0];
	int SRp1_end = vote -> indel_recorder[max_I][max_J][1];
	return SRp1_start <= 4 && SRp1_end >= applied_subreads - 4 - ( applied_subreads>=12?1:0 );
*/
}

int cellCounts_build_simple_mode_subread_masks(cellcounts_global_t * cct_context, int thread_no, int applied_subreads){
	if(applied_subreads<9) return 0;
	int last_sr_0B = applied_subreads-1;
	int sr_step = (last_sr_0B -1) * 10000 / 4 + 1, sri, ret=0;
	for(sri = 0 ; sri < last_sr_0B*10000+100; sri+= sr_step) ret |= 1<<( sri / 10000 );
	return ret;
}


int cellCounts_do_voting(cellcounts_global_t * cct_context, int thread_no) {
	subread_read_number_t current_read_number=0;
	char * read_text, * qual_text;
	char read_name[MAX_READ_NAME_LEN+1];
	char read_bin[REVERSED_READ_BIN_OFFSET * 2];
	int read_len=0;

	read_text = malloc(MAX_SCRNA_READ_LENGTH * 2+2);
	qual_text = malloc(MAX_SCRNA_READ_LENGTH * 2+2);

	temp_votes_per_read_t prefill_ptrs;
	gene_sc_vote_t * vote_me = malloc(sizeof(gene_sc_vote_t));

	if(vote_me==NULL) {
		SUBREADprintf("Cannot allocate voting memory.\n");
		return -1;
	}

	int index_gap_width = cct_context -> current_index -> index_gap;

	while(!cct_context -> has_error) {
		int subread_no;
		int is_reversed, applied_subreads = 0;

		cellCounts_fetch_next_read_pair(cct_context, thread_no,  &read_len, read_name, read_text, qual_text, &current_read_number);
		if(current_read_number < 0) break;
		if(read_len< 16) continue;

		int CR15GLS = (read_len - 15 - index_gap_width)<<16;
		int subread_step =  CR15GLS /(cct_context -> total_subreads_per_read -1);
		if(subread_step<(index_gap_width<<16))subread_step = index_gap_width<<16;
		applied_subreads = 1 + CR15GLS / subread_step;

		int building_rbin_offset = 0, read_text_rev_offset =0;
		for(is_reversed = 0; is_reversed<2; is_reversed++) {
			gehash_key_t subread_integer = 0;
			int last_vote_rpos = -16;
			for(subread_no=0; subread_no < applied_subreads ; subread_no++) {
				int subread_offset = ((subread_step * subread_no) >> 16);
				#define SHIFT_SUBREAD_INT(ii, pp) { int nch = read_text [pp+read_text_rev_offset]; ii = (ii << 2) | base2int( nch );}
				#define BUILD_RBIN  {  int new2b = subread_integer & 3;\
					int rbin_byte = building_rbin_offset + (last_vote_rpos +16)/4;\
					int rbin_bit =(last_vote_rpos +16)%4 *2;\
					if(rbin_bit ==0) read_bin[rbin_byte]=0;\
					read_bin[rbin_byte] |= new2b<<rbin_bit;  }
				for(; last_vote_rpos  < subread_offset ; last_vote_rpos ++){
					SHIFT_SUBREAD_INT(subread_integer , last_vote_rpos  +16);
					BUILD_RBIN;
				}
				prefill_votes(cct_context->current_index, &prefill_ptrs, applied_subreads, subread_integer, subread_offset, subread_no, is_reversed);
			}
			if(last_vote_rpos > read_len - 16)SUBREADprintf("ERROR: exceeded offset %d > %d\n", last_vote_rpos , read_len - 16);
			for(; last_vote_rpos  < read_len - 16 ; last_vote_rpos ++){
				SHIFT_SUBREAD_INT(subread_integer , last_vote_rpos  +16);
				BUILD_RBIN;
			}

			if(is_reversed) {
				cellCounts_process_copy_ptrs_to_votes(cct_context, thread_no, &prefill_ptrs, vote_me, applied_subreads, read_name);
#ifdef __MINGW32__
				if(current_read_number % 1000000 == 0 && current_read_number>0) print_in_box(80,0,0,"  Mapped : % 13" PRId64 " reads; time elapsed : % 5.1f mins\n", cct_context -> all_processed_reads_before_chunk + current_read_number, ( - cct_context -> program_start_time + miltime() ) / 60.);
#else
				if(current_read_number % 1000000 == 0 && current_read_number>0) print_in_box(80,0,0,"  Mapped : % 13lld reads; time elapsed : % 5.1f mins\n", cct_context -> all_processed_reads_before_chunk + current_read_number, ( - cct_context -> program_start_time + miltime() ) / 60.);
#endif
				cellCounts_select_and_write_alignments(cct_context, thread_no, vote_me, read_name, read_text, read_bin, qual_text, read_len, applied_subreads);
			} else {
				building_rbin_offset = REVERSED_READ_BIN_OFFSET;
				read_text_rev_offset = MAX_SCRNA_READ_LENGTH+1;
				strcpy(read_text+read_text_rev_offset, read_text);
				reverse_read(read_text+read_text_rev_offset, read_len, GENE_SPACE_BASE);
				qual_text[read_text_rev_offset] = 0;
			}
		}
	}

	free(vote_me);
	free(read_text);
	free(qual_text);

	return cct_context -> has_error;
}






#define MAKE_SUBREAD_OFFSET	if(subread_no == applied_subreads -1) subread_offset= read_len-16; else subread_offset= ((subread_step * subread_no) >> 16);

#define MAKE_SUBREAD_INTVAL	subread_integer =0; for(xk1 = 0; xk1 < 16; xk1++){\
					int rbin_byte = (xk1+subread_offset)/4 + is_reversed * REVERSED_READ_BIN_OFFSET;\
					int rbin_bit = (xk1+subread_offset) %4 *2;\
					unsigned int vtmp = ( read_bin[ rbin_byte ] >> rbin_bit )&3;\
					subread_integer |= vtmp<<(2*(15-xk1)); }

int simpleMode_cellCounts_do_voting(cellcounts_global_t * cct_context, int thread_no) {
	int xk1;
	subread_read_number_t current_read_number=0;
	char * read_text, * qual_text;
	char read_name[MAX_READ_NAME_LEN+1];
	char read_bin[REVERSED_READ_BIN_OFFSET * 2];
	int read_len=0;
	int subread_in_simple_mode_masks[SCRNA_SUBREADS_HARD_LIMIT];

	memset(subread_in_simple_mode_masks, 0xff, sizeof(int)*SCRNA_SUBREADS_HARD_LIMIT);
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	read_text = malloc(MAX_SCRNA_READ_LENGTH * 2+2);
	qual_text = malloc(MAX_SCRNA_READ_LENGTH * 2+2);

	temp_votes_per_read_t prefill_ptrs;
	gene_sc_vote_t * vote_me = malloc(sizeof(gene_sc_vote_t));

	if(vote_me==NULL) {
		SUBREADprintf("Cannot allocate voting memory.\n");
		return -1;
	}

	int index_gap_width = cct_context -> current_index -> index_gap;

	while(1) {
		int subread_no;
		int is_reversed, applied_subreads = 0;

		cellCounts_fetch_next_read_pair(cct_context, thread_no,  &read_len, read_name, read_text, qual_text, &current_read_number);
		if(current_read_number < 0) break;
		if(read_len< 16) continue;

		int CR15GLS = (read_len - 15 - index_gap_width)<<16;
		int subread_step =  CR15GLS /(cct_context -> total_subreads_per_read -1);
		if(subread_step<(index_gap_width<<16))subread_step = index_gap_width<<16;
		applied_subreads = 1 + CR15GLS / subread_step;

		if(1){
			for(is_reversed = 0; is_reversed<2; is_reversed++) {
				for(xk1=0; xk1<read_len ; xk1++){
					int nch;
					if(is_reversed){
						nch = read_text[read_len - 1 - xk1];
						if(nch == 'T') nch = 'A';
						else if(nch == 'C') nch = 'G';
						else if(nch == 'G') nch = 'C';
						else nch = 'T';
					}else nch = read_text[xk1];
					int rbin_byte = xk1/4 + is_reversed * REVERSED_READ_BIN_OFFSET;
					int rbin_bit = xk1%4 *2;
					if(rbin_bit == 0) read_bin[ rbin_byte ] = 0;
					read_bin[ rbin_byte ] |= (base2int(nch)) << rbin_bit;
				}
			}

			if(subread_in_simple_mode_masks[applied_subreads] < 0)
				subread_in_simple_mode_masks[applied_subreads] = cellCounts_build_simple_mode_subread_masks(cct_context, thread_no, applied_subreads);

			int written = 0, simple_mode;
			for(simple_mode = 1; simple_mode >=0; simple_mode --){
				int read_text_rev_offset =0;
				for(is_reversed = 0; is_reversed<2; is_reversed++) {
					for(subread_no=0; subread_no < applied_subreads ; subread_no++) {
						int subread_offset;
						gehash_key_t subread_integer;

						int this_subread_used_in_simple;
						if(is_reversed) this_subread_used_in_simple = subread_in_simple_mode_masks[applied_subreads] & (1<<( applied_subreads -1 - subread_no));
						else this_subread_used_in_simple = subread_in_simple_mode_masks[applied_subreads] & (1<<subread_no);
						if(( simple_mode && !this_subread_used_in_simple ) || ( 0==simple_mode  && this_subread_used_in_simple )){
							int my_no = subread_no+ is_reversed*applied_subreads;
							prefill_ptrs.votes[my_no] = 0;
							continue;
						}

						MAKE_SUBREAD_OFFSET;
						MAKE_SUBREAD_INTVAL;

						prefill_votes(cct_context->current_index, &prefill_ptrs, applied_subreads, subread_integer, subread_offset, subread_no, is_reversed);
					}

					if(is_reversed) {
						simpleMode_cellCounts_process_copy_ptrs_to_votes(cct_context, thread_no, &prefill_ptrs, vote_me, applied_subreads, read_name, simple_mode);

						if(current_read_number % 1000000 == 0) SUBREADprintf("Mapping and counting: %lld  ; %.2f mins\n", cct_context -> all_processed_reads_before_chunk + current_read_number, ( - cct_context -> program_start_time + miltime() ) / 60.);
						if((simple_mode && cellCounts_simple_mode_highconf(cct_context, thread_no, applied_subreads, vote_me, read_name)) || (simple_mode == 0)){
							cellCounts_select_and_write_alignments(cct_context, thread_no, vote_me, read_name, read_text, read_bin, qual_text, read_len, applied_subreads);
							thread_context -> hiconf_map += simple_mode;
							thread_context -> loconf_map += (0==simple_mode);
							written =1;
						}
					} else {
						read_text_rev_offset = MAX_SCRNA_READ_LENGTH+1;
						if(simple_mode){
							strcpy(read_text+read_text_rev_offset, read_text);
							reverse_read(read_text+read_text_rev_offset, read_len, GENE_SPACE_BASE);
							qual_text[read_text_rev_offset] = 0;
						}
					}
				}
				if(written) break;
			}
		}
	}

	free(vote_me);
	free(read_text);
	free(qual_text);

	return 0;
}


void cellCounts_absoffset_to_posstr(cellcounts_global_t * cct_context, unsigned int pos, char * res){
	char * ch;
	int off;
	locate_gene_position(pos, &cct_context -> chromosome_table, &  ch, &off);
	SUBreadSprintf(res, 100 , "%s:%u", ch, off);
}

void known_pointer_strcat(char * targ, char * src, char ** buf){
	int srclen = strlen(src);
	if( (*buf) == NULL){
		(*buf) = targ;
	}
	memcpy((*buf), src, srclen);
	(*buf) += srclen;
	(**buf) = 0;
}

int cellCounts_write_gene_list(cellcounts_global_t * cct_context){
	int xk1;
	char ofname[MAX_FILE_NAME_LENGTH + 100];
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 100,"%s.Annot",cct_context->output_prefix);
	FILE * fp_out = fopen( ofname , "w" );
	fprintf(fp_out,"GeneID\tChr\tStart\tEnd\tStrand\tLength\n");

	unsigned int * gene_exons_number = calloc(sizeof(unsigned int) , cct_context -> gene_name_table -> numOfElements);
	unsigned int * gene_exons_pointer = calloc(sizeof(unsigned int) , cct_context -> gene_name_table -> numOfElements);
	unsigned int * gene_exons_start = malloc(sizeof(unsigned int) * cct_context -> all_features_array -> numOfElements);
	unsigned int * gene_exons_end = malloc(sizeof(unsigned int) * cct_context -> all_features_array -> numOfElements);
	char ** gene_exons_chr = malloc(sizeof(char *) * cct_context -> all_features_array -> numOfElements);
	char * gene_exons_strand = malloc(cct_context -> all_features_array -> numOfElements);

	for(xk1 = 0; xk1 < cct_context -> all_features_array -> numOfElements; xk1++) {
		int gene_id = cct_context -> features_sorted_geneid[xk1];
		gene_exons_number[gene_id]++;
	}

	unsigned int accumulative_no = 0;
	unsigned longest_gene_exons = 0;
	for(xk1 = 0 ; xk1 < cct_context -> gene_name_table -> numOfElements; xk1++) {
		unsigned int this_gene_exons = gene_exons_number[xk1];
		longest_gene_exons = max(longest_gene_exons, this_gene_exons);
		gene_exons_number[xk1] = accumulative_no;
		accumulative_no += this_gene_exons;
	}

	for(xk1 = 0; xk1 < cct_context -> all_features_array -> numOfElements; xk1++) {
		int gene_id = cct_context -> features_sorted_geneid[xk1];
		int gene_write_ptr = gene_exons_number[gene_id] + gene_exons_pointer[gene_id];

		gene_exons_chr[gene_write_ptr] = cct_context -> features_sorted_chr[xk1];
		gene_exons_start[gene_write_ptr] = cct_context -> features_sorted_start[xk1]; 
		gene_exons_end[gene_write_ptr] = cct_context -> features_sorted_stop[xk1]; 
		gene_exons_strand[gene_write_ptr] = cct_context -> features_sorted_strand[xk1]; 

		gene_exons_pointer[gene_id]++;
	}

	char *is_occupied = malloc(longest_gene_exons);
	unsigned int * input_start_stop_list = malloc(longest_gene_exons * sizeof(int) * 2);
	unsigned int * output_start_stop_list = malloc(longest_gene_exons * sizeof(int) * 2);
	int disk_is_full = 0;

	char * out_chr_list = malloc(longest_gene_exons * (1+cct_context -> longest_chro_name) + 1), * tmp_chr_list = NULL;
	char * out_start_list = malloc(11 * longest_gene_exons + 1), * tmp_start_list = NULL;
	char * out_end_list = malloc(11 * longest_gene_exons + 1), * tmp_end_list = NULL;
	char * out_strand_list = malloc(2 * longest_gene_exons + 1), * tmp_strand_list = NULL;

	for(xk1 = 0 ; xk1 < cct_context -> gene_name_table -> numOfElements; xk1++) {
		int xk2;
		
		memset(is_occupied,0,gene_exons_pointer[xk1]);
		tmp_chr_list = NULL;
		tmp_start_list = NULL;
		tmp_end_list = NULL;
		tmp_strand_list = NULL;
		out_chr_list[0]=0;
		out_start_list[0]=0;
		out_end_list[0]=0;
		out_strand_list[0]=0;
		int gene_nonoverlap_len =0;

		unsigned char * gene_symbol = cct_context -> gene_name_array [xk1];
		for(xk2=0; xk2<gene_exons_pointer[xk1]; xk2++) {
			if(!is_occupied[xk2]) {
				int xk3;
				char * matched_chr = gene_exons_chr[xk2 + gene_exons_number[xk1]];
				char matched_strand = gene_exons_strand[xk2 + gene_exons_number[xk1]];

				memset(input_start_stop_list, 0, gene_exons_pointer[xk1] * sizeof(int) * 2);
				int gap_merge_ptr = 1;
				input_start_stop_list[0] = gene_exons_start[xk2 + gene_exons_number[xk1]];
				input_start_stop_list[1] = gene_exons_end[xk2 + gene_exons_number[xk1]] + 1;

				for(xk3 = xk2; xk3 < gene_exons_pointer[xk1]; xk3++)
				{
					if(xk3==xk2)continue;

					if((!is_occupied[xk3]) && strcmp(matched_chr, gene_exons_chr[xk3+gene_exons_number[xk1]])==0 && matched_strand == gene_exons_strand[xk3 + gene_exons_number[xk1]])
					{
						is_occupied[xk3]=1;
						input_start_stop_list[gap_merge_ptr*2] = gene_exons_start[xk3+gene_exons_number[xk1]]; 
						input_start_stop_list[gap_merge_ptr*2+1] = gene_exons_end[xk3+gene_exons_number[xk1]]+1;

						gap_merge_ptr++;
					}
				}

				{
						int merged_gaps = mergeIntervals(input_start_stop_list, output_start_stop_list, gap_merge_ptr);

						for(xk3=0; xk3<gap_merge_ptr; xk3++)
						{
							char numbbuf[12];
							known_pointer_strcat(out_chr_list, matched_chr, &tmp_chr_list);
							known_pointer_strcat(out_chr_list, ";", &tmp_chr_list);

							SUBreadSprintf(numbbuf,12,"%u;", input_start_stop_list[xk3 * 2]);
							known_pointer_strcat(out_start_list, numbbuf, &tmp_start_list);
							SUBreadSprintf(numbbuf,12,"%u;", input_start_stop_list[xk3 * 2 + 1] - 1);
							known_pointer_strcat(out_end_list, numbbuf, &tmp_end_list);
							SUBreadSprintf(numbbuf,12,"%c;", (matched_strand==1)?'-':( ( matched_strand==0 )? '+':'.'));
							known_pointer_strcat(out_strand_list, numbbuf, &tmp_strand_list);

						}
						for(xk3=0; xk3<merged_gaps; xk3++)
							gene_nonoverlap_len += output_start_stop_list[xk3 * 2 + 1] - output_start_stop_list[xk3 * 2];
				}
			}
		}
		#define _cut_tail(x) (x)[strlen(x)-1]=0

		_cut_tail(out_chr_list);
		_cut_tail(out_start_list);
		_cut_tail(out_end_list);
		_cut_tail(out_strand_list);

		int wlen = fprintf(fp_out, "%s\t%s\t%s\t%s\t%s\t%d\n", gene_symbol, out_chr_list, out_start_list, out_end_list, out_strand_list, gene_nonoverlap_len);
		if(wlen < 6)disk_is_full = 1;
	}

	free(is_occupied);
	free(input_start_stop_list);
	free(output_start_stop_list);
	free(out_chr_list);
	free(out_strand_list);
	free(out_start_list);
	free(out_end_list);

	free(gene_exons_number);
	free(gene_exons_pointer);
	free(gene_exons_chr);
	free(gene_exons_start);
	free(gene_exons_end);
	free(gene_exons_strand);
	fclose(fp_out);

	if(disk_is_full){
		SUBREADprintf("ERROR: disk is full; the count file cannot be generated.\n");
		unlink(ofname);
		return -1;
	}
	return 0;
}



int cellCounts_run_mapping(cellcounts_global_t * cct_context){
	int chunk_no = 0;

	cct_context -> current_index = (gehash_t*) malloc(sizeof(gehash_t));
	cct_context -> running_processed_reads_in_chunk=0;
	cct_context -> processed_reads_in_chunk=0;
	cct_context -> all_processed_reads_before_chunk = 0;
	sread_len = 0;

	while(1) {
		int ret = 0;

		for(cct_context->current_index_block_number = 0; cct_context->current_index_block_number < cct_context->total_index_blocks; cct_context->current_index_block_number++) {	   
			char tmp_fname[MAX_FILE_NAME_LENGTH+30];

			if(cct_context->total_index_blocks > 1 || chunk_no == 0) {	   
				SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+30, "%s.%02d.b.tab", cct_context->index_prefix, cct_context->current_index_block_number);
				print_in_box(80,0,0, "Load the %d-%s index block...",1+ cct_context->current_index_block_number, cct_context->current_index_block_number==0?"st":(cct_context->current_index_block_number==1?"nd":"th"));
				if(gehash_load(cct_context -> current_index, tmp_fname)) return -1;
				print_in_box(80,0,0, "The index block has been loaded. Now map the reads...");
				print_in_box(80,0,0, "");
				SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+30, "%s.%02d.b.array", cct_context->index_prefix, cct_context->current_index_block_number);
			}
			
			if(cct_context->total_index_blocks == cct_context->current_index_block_number + 1)
				cct_context -> is_final_voting_run = 1;
			else	cct_context -> is_final_voting_run = 0;
			
			ret = cellCounts_run_maybe_threads(cct_context, STEP_VOTING);
			cct_context -> processed_reads_in_chunk = cct_context -> running_processed_reads_in_chunk;
			int is_last_chunk = cct_context -> processed_reads_in_chunk < cct_context-> reads_per_chunk;
			
			if(cct_context->total_index_blocks > 1 || is_last_chunk)
				gehash_destory_fast(cct_context -> current_index);
			
			if(ret) break;
			if(!cct_context -> processed_reads_in_chunk) break;
		}

		cellCounts_go_chunk_nextchunk(cct_context);
		cct_context -> all_processed_reads_before_chunk += cct_context -> processed_reads_in_chunk ;

		if(ret) return ret;

		if(cct_context -> processed_reads_in_chunk < cct_context -> reads_per_chunk ||
		  (cct_context -> output_binfiles_are_full))
			// base value indexes loaded in the last circle are not destroyed and are used in writting the indel VCF.
			break;

		cellCounts_clean_context_after_chunk(cct_context);
		chunk_no++;
	}

	free(cct_context -> current_index);
	return 0;
}

#define CELLCOUNTS_BAMBLOCK_SIZE 60000
#define CELLCOUNTS_BAMBLOCK_COMP_NUMBER 1
#define CELLRANGER_MERGER_WORKER_BINSIZE 62000

struct scRNA_merge_batches_worker_task{
	int sample_id;
	int inbin_len;
	int inbin_number;
	int inbin_batch_start_offsets [CELLCOUNTS_BAMBLOCK_COMP_NUMBER ];
	srInt_64 block_number;
	char inbin[(READ_BIN_BUF_SIZE+CELLCOUNTS_BAMBLOCK_SIZE )*CELLCOUNTS_BAMBLOCK_COMP_NUMBER ];
};

struct scRNA_merge_batches_worker_current{
	struct scRNA_merge_batches_worker_task * task;
	char outbin[CELLRANGER_MERGER_WORKER_BINSIZE * CELLCOUNTS_BAMBLOCK_COMP_NUMBER ];
	int outbin_len[CELLCOUNTS_BAMBLOCK_COMP_NUMBER ];
	unsigned int crc32[CELLCOUNTS_BAMBLOCK_COMP_NUMBER ];

	z_stream strm;
};

struct cell_gene_umi_supp{
	int cellbc;
	srInt_64 gene_no;
	char umi[MAX_UMI_LEN];
	int supp_reads;
};

int cellCounts_hamming_max2_fixlen(char * u1, char * u2, int ulen){
	int x, ret=0;
	for(x=0; x<ulen; x++){
		if(u1[x]!=u2[x]) ret++;
		if(ret>1)return ret;
	}
	return ret;
}

#define ADD_count_hash(bc,gn,no)   HashTablePut(cellBCp0_genep0_P1_to_UMIs, NULL +1+(((1LLU*(bc))<<32)| (gn) ),  HashTableGet(   cellBCp0_genep0_P1_to_UMIs, NULL +1+(((1LLU*(bc))<<32)| (gn))) +(no) )
void cellCounts_do_one_batch_UMI_merge_one_cell(ArrayList* structs, int sec_start, int sec_end, int is_UMI_step2, HashTable * filtered_CGU_table, srInt_64 * remove_count){
	int x1;
	void ** app1 = structs -> appendix1;
	cellcounts_global_t * cct_context = app1[0];
	HashTable * cellBCp0_genep0_P1_to_UMIs = app1[2];
	int sample_id = app1[3]-NULL;

	if(is_UMI_step2){
		// NB: when this function is called, sec_end - sec_start MUST be >=2.
		for(x1 = sec_start; x1<sec_end; x1++) {
			struct cell_gene_umi_supp * str1 = ArrayListGet(structs, x1);
			if(x1 == sec_start){
				struct cell_gene_umi_supp * str2 = ArrayListGet(structs, sec_start+1);
				if(str1 -> supp_reads > str2 -> supp_reads){
					ADD_count_hash(str1->cellbc, str1->gene_no,1);
					continue;
				} else if(remove_count)(*remove_count)++;
			}

			str1 -> cellbc = -1;

			char replaced_key[40+MAX_UMI_LEN];
#ifdef __MINGW32__
			int keyptr = SUBreadSprintf(replaced_key, 40+MAX_UMI_LEN,"%d-%" PRId64 "-", str1 -> cellbc, str1 -> gene_no);
#else
			int keyptr = SUBreadSprintf(replaced_key, 40+MAX_UMI_LEN,"%d-%lld-", str1 -> cellbc, str1 -> gene_no);
#endif
			memcpy(replaced_key+keyptr, str1 -> umi, cct_context -> UMI_length);
			replaced_key[keyptr+cct_context -> UMI_length]=0;
			HashTablePut(filtered_CGU_table, strdup(replaced_key), NULL-1);
		}
	}else{
		ArrayList * accepted_list =NULL;
		HashTable * looktable = NULL;
		if(sec_end - sec_start >30){
			looktable = StringTableCreate((sec_end - sec_start)/5);
			HashTableSetDeallocationFunctions(looktable, free, (void (*)(void *value))ArrayListDestroy);
		}else accepted_list = ArrayListCreate(sec_end - sec_start);

		for(x1=sec_start; x1<sec_end; x1++){
			struct cell_gene_umi_supp * try_str = ArrayListGet(structs , x1);
			int x2, found = 0;
			ArrayList * test_accs;
			int hx;

			if(looktable){
				for(hx = 0; hx<2; hx++){
					char test_ky[MAX_UMI_LEN];
					test_ky[0] = hx?'S':'F';
					memcpy(test_ky +1, try_str -> umi + hx * cct_context -> UMI_length/2 , cct_context -> UMI_length/2);
					test_ky[1+cct_context -> UMI_length/2]=0;

					test_accs = HashTableGet(looktable, test_ky);
					if(!test_accs)continue;

					for(x2=0; x2<test_accs->numOfElements; x2++){
						struct cell_gene_umi_supp * acc_str = ArrayListGet(test_accs, x2);
						if(cellCounts_hamming_max2_fixlen(acc_str -> umi, try_str -> umi, cct_context -> UMI_length)<2){
							found=1;
							acc_str -> supp_reads += try_str -> supp_reads;
							try_str -> cellbc = -1;

							char replaced_key[55+MAX_UMI_LEN];
#ifdef __MINGW32__
							int keyptr = SUBreadSprintf(replaced_key, 55+MAX_UMI_LEN,"%d-%d-%" PRId64 "-", sample_id, try_str -> cellbc, try_str -> gene_no);
#else
							int keyptr = SUBreadSprintf(replaced_key, 55+MAX_UMI_LEN,"%d-%d-%lld-", sample_id, try_str -> cellbc, try_str -> gene_no);
#endif

							memcpy(replaced_key+keyptr, try_str -> umi, cct_context -> UMI_length);
							replaced_key[keyptr+cct_context -> UMI_length]=0;
							HashTablePut(filtered_CGU_table, strdup(replaced_key), acc_str -> umi);
							break;
						}
					}
					if(found)break;
				}
			}else{
				test_accs = accepted_list;

				for(x2=0; x2<test_accs->numOfElements; x2++){
					struct cell_gene_umi_supp * acc_str = ArrayListGet(test_accs, x2);
					if(cellCounts_hamming_max2_fixlen(acc_str -> umi, try_str -> umi, cct_context -> UMI_length)<2){
						found=1;
						acc_str -> supp_reads += try_str -> supp_reads;
						try_str -> cellbc = -1;

						char replaced_key[55+MAX_UMI_LEN];
#ifdef __MINGW32__
						int keyptr = SUBreadSprintf(replaced_key, 55+MAX_UMI_LEN,"%d-%d-%" PRId64 "-", sample_id, try_str -> cellbc, try_str -> gene_no);
#else
						int keyptr = SUBreadSprintf(replaced_key, 55+MAX_UMI_LEN,"%d-%d-%lld-", sample_id, try_str -> cellbc, try_str -> gene_no);
#endif
						memcpy(replaced_key+keyptr, try_str -> umi, cct_context -> UMI_length);
						replaced_key[keyptr+cct_context -> UMI_length]=0;
						HashTablePut(filtered_CGU_table, strdup(replaced_key), acc_str -> umi);
						break;
					}
				}
			}
			if(!found){
				if(looktable){
					for(hx = 0; hx<2; hx++){
						char test_ky[MAX_UMI_LEN];
						test_ky[0] = hx?'S':'F';
						memcpy(test_ky +1, try_str -> umi + hx * cct_context -> UMI_length/2 , cct_context -> UMI_length/2);
						test_ky[1+cct_context -> UMI_length/2]=0;
						test_accs = HashTableGet(looktable, test_ky);
						if(!test_accs){
							test_accs = ArrayListCreate(10);
							HashTablePut(looktable, strdup(test_ky), test_accs);
						}
						ArrayListPush(test_accs, try_str);
					}
				}else ArrayListPush(accepted_list, try_str);
			}
		}

		if(looktable)HashTableDestroy(looktable);
		else ArrayListDestroy(accepted_list);
	}
}

void cellCounts_do_one_batch_UMI_merge_one_step(ArrayList* structs, int is_UMI_step2, HashTable * filtered_CGU_table, srInt_64 * remove_count){
	void ** app1 = structs -> appendix1;
	cellcounts_global_t * cct_context = app1[0];
	HashTable * cellBCp0_genep0_P1_to_UMIs = app1[2];
	srInt_64 x1, sec_start = 0;
	srInt_64 old_sec_key = -1;

	for(x1=0; x1<=structs -> numOfElements; x1++){
		srInt_64 sec_key = -1;
		int is_umi_changed = 0;

		struct cell_gene_umi_supp * str1 =NULL;
		if(x1<structs -> numOfElements){
			str1 = ArrayListGet(structs, x1);
			if(str1 -> cellbc <0) continue;
			sec_key = str1 -> cellbc;
			sec_key = sec_key << 32;
			if(is_UMI_step2 && sec_key == old_sec_key){
				struct cell_gene_umi_supp * strold = ArrayListGet(structs, sec_start);
				is_umi_changed = memcmp(strold -> umi, str1-> umi, cct_context-> UMI_length);
			}else if(!is_UMI_step2) sec_key = sec_key | str1 -> gene_no;
				// gene_no itself is 64-bit, but it is nearly impossible to have two neighbouring
				// structures that have the same last 32-bit of gene_no.
		}

		if( (x1>sec_start && sec_key!=old_sec_key) || is_umi_changed){ // when x1 == numOfElements, sec_key is -1. If old_sec_key is also -1, no item is included in the list. If old_sec_key is >=0, the last sec is processed.
			struct cell_gene_umi_supp * str0 = ArrayListGet(structs, sec_start);
			if(x1 - sec_start>1 && str0->cellbc>=0) cellCounts_do_one_batch_UMI_merge_one_cell(structs, sec_start, x1, is_UMI_step2, filtered_CGU_table, remove_count);
			else if(is_UMI_step2 && str0->cellbc>=0) ADD_count_hash(str0->cellbc,str0->gene_no,1);

			sec_start = x1;
		}
		old_sec_key = sec_key;
	}
}


int cellCounts_do_one_batch_sort_compare(void * ar, int l, int r){
	void ** arr = ar;
	void ** bin_ptrs = arr[0];
	cellcounts_global_t * cct_context = arr[1];

	char * Lptr = bin_ptrs[l];
	char * Rptr = bin_ptrs[r];
	srInt_64 Lgenes=0, Rgenes=0;
	memcpy(&Lgenes, Lptr+8, 8);
	memcpy(&Rgenes, Rptr+8, 8);
	if(Lgenes & (1LLU<<63))Lgenes=Lgenes & 0x7fffffffllu; else Lgenes=0;
	if(Rgenes & (1LLU<<63))Rgenes=Rgenes & 0x7fffffffllu; else Rgenes=0;
	srInt_64 Lpos= ((0LLU+*(int*)(Lptr+16+Lgenes*8+cct_context->UMI_length+4))<<32) | *(unsigned int*)(Lptr+16+Lgenes*8+cct_context->UMI_length+4+4);
	srInt_64 Rpos= ((0LLU+*(int*)(Rptr+16+Rgenes*8+cct_context->UMI_length+4))<<32) | *(unsigned int*)(Rptr+16+Rgenes*8+cct_context->UMI_length+4+4);
	if(Lpos>Rpos)return 1;
	if(Lpos<Rpos)return -1;
	return 0;
}

void cellCounts_do_one_batch_sort_exchange(void * ar, int l, int r){
	void ** arr = ar;
	void ** bin_ptrs = arr[0];
	void * tp = bin_ptrs[l];
	bin_ptrs[l]=bin_ptrs[r];
	bin_ptrs[r]=tp;
}

void cellCounts_do_one_batch_sort_merge(void * ar, int start, int items, int items2){
	void ** arr = ar;
	void ** bin_ptrs = arr[0];
	bin_ptrs +=start;

	void ** tmp = malloc(sizeof(void*)*(items2+items));
	int i1_cursor=0, i2_cursor=items, wptr=0;
	while(1){
		if(i1_cursor == items && i2_cursor == items + items2 )break;
		int select_items_1 = (i2_cursor == items + items2) || (i1_cursor < items && cellCounts_do_one_batch_sort_compare(ar, start+ i1_cursor,start + i2_cursor) <= 0);
		if(select_items_1) tmp[wptr++] = bin_ptrs[i1_cursor++];
		else tmp[wptr++] = bin_ptrs[i2_cursor++];
	}
	memcpy(bin_ptrs, tmp, sizeof(void*)*(items2+items));
	free(tmp);
}

int cellCounts_do_one_batch_tab_to_struct_list_compare(void * L_elem, void * R_elem, ArrayList * me){
	struct cell_gene_umi_supp *L = L_elem, *R = R_elem;
	void ** app1 = me -> appendix1;
	cellcounts_global_t * cct_context = app1[0];
	int sort_by_geneid_then_umi = app1[1] - NULL;

	if(L->cellbc > R->cellbc) return 1;
	if(L->cellbc < R->cellbc) return -1;

	if(sort_by_geneid_then_umi){
		if(L->gene_no>R->gene_no) return 1;
		if(L->gene_no<R->gene_no) return -1;
	}else{
		int umicmps = memcmp(L->umi, R->umi, cct_context -> UMI_length);
		if(umicmps) return umicmps;
	}

	if(L->supp_reads < R->supp_reads) return 1;
	if(L->supp_reads > R->supp_reads) return -1; // reversed by # supp reads

	if(sort_by_geneid_then_umi){
		int umicmps = memcmp(L->umi, R->umi, cct_context -> UMI_length);
		if(umicmps) return umicmps;
	}else{

		if(L->gene_no>R->gene_no) return 1;
		if(L->gene_no<R->gene_no) return -1;
	}
	return 0;
}

void cellCounts_do_one_batch_tab_to_struct_list(void *ky, void *val, HashTable * tab){
	int supp_reads = val-NULL;
	ArrayList ** cell_gene_umi_list = tab -> appendix1;
	int UMI_length = tab -> counter1;

	struct cell_gene_umi_supp * new_item = malloc(sizeof(struct cell_gene_umi_supp));
	char * kyptr = ky;
	int sample_id = atoi(kyptr); // one-based sample id
	for(; '-' != *kyptr; kyptr++);
	kyptr++;
	new_item -> cellbc = atoi(kyptr);
	for(; '-' != *kyptr; kyptr++);
	kyptr++;
	new_item -> gene_no = atoll(kyptr);
	for(; '-' != *kyptr; kyptr++);
	memcpy(new_item->umi, kyptr+1, UMI_length);
	new_item -> supp_reads = supp_reads;
	if(sample_id<1)SUBREADprintf("WRONG SAMPLE ID: %d from '%s'\n", sample_id, (char*)ky);
	ArrayListPush(cell_gene_umi_list[sample_id-1], new_item);
}

void cellCounts_do_one_batch_write_UMIs(void * vcell_gene, void * vumis, HashTable * me){
	FILE * fp = me->appendix1;
	vcell_gene --;
	fwrite(&vcell_gene,1,8,fp);
	fwrite(&vumis,1,8,fp);
}


void * cellCounts_merge_batches_worker(void * vp){
	void **vpp = vp;
	cellcounts_global_t * cct_context = vpp[0];
	worker_master_mutex_t * worker_mut  = vpp[1];
	int my_worker_id = vpp[2] - NULL;
	struct scRNA_merge_batches_worker_current * my_current_job = vpp[3];
	free(vp);

	int Z_DEFAULT_MEM_LEVEL = 8;
	worker_thread_start(worker_mut, my_worker_id);
	while(1){
		if(worker_wait_for_job(worker_mut, my_worker_id)) break;
		if(!cct_context -> is_BAM_and_FQ_out_generated) continue;


		struct scRNA_merge_batches_worker_task * current_input = my_current_job -> task;
		int current_blk;
		for(current_blk =0; current_blk < current_input -> inbin_number ; current_blk ++){
			char * inbin_blk = current_input -> inbin + current_input -> inbin_batch_start_offsets[current_blk ];
			int inblock_size = -1;
			if(current_blk == current_input -> inbin_number -1) inblock_size= current_input -> inbin_len - current_input -> inbin_batch_start_offsets[current_blk ];
			else if(CELLCOUNTS_BAMBLOCK_COMP_NUMBER>1) inblock_size= current_input  -> inbin_batch_start_offsets[current_blk +1] - current_input -> inbin_batch_start_offsets[current_blk ];

			deflateInit2(&my_current_job -> strm , SAMBAM_COMPRESS_LEVEL_NORMAL, Z_DEFLATED, SAMBAM_GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
			my_current_job -> strm.avail_in = inblock_size ;
			my_current_job -> strm.next_in = (unsigned char*)inbin_blk ;
			my_current_job -> strm.avail_out = CELLRANGER_MERGER_WORKER_BINSIZE;
			my_current_job -> strm.next_out = (unsigned char*)(my_current_job -> outbin + CELLRANGER_MERGER_WORKER_BINSIZE * current_blk);

			deflate(&my_current_job -> strm, Z_FINISH);
			my_current_job -> outbin_len [current_blk] = CELLRANGER_MERGER_WORKER_BINSIZE-my_current_job -> strm.avail_out;
			my_current_job -> crc32 [current_blk] = SamBam_CRC32(inbin_blk, inblock_size);
			deflateEnd(&my_current_job -> strm);
		}
	}
	return NULL;
}

int cellCounts_make_barcode_bam_bin(cellcounts_global_t * cct_context, char * rbin, char * new_rbin, int binlen, char * fixedbc_seq, char * fixedumi_seq, srInt_64 gene_no, srInt_64 * genes) {
	char * cellbc_seq=NULL,*umi_seq=NULL, * cellbc_qual=NULL,*umi_qual=NULL, *sample_seq=NULL, *sample_qual=NULL, *lane_str=NULL;
	int rname_trimmed_len=0;
	cellCounts_scan_read_name_str(cct_context, rbin, NULL, & sample_seq, & sample_qual, & cellbc_seq, & cellbc_qual, & umi_seq, & umi_qual, &lane_str, NULL, &rname_trimmed_len);
	int new_rbin_len = 0, n_cigar_op =0, l_read_name=0, l_seq=0;

	memcpy(new_rbin, rbin, 36);
	new_rbin_len += 36;

	memcpy(&n_cigar_op, rbin+16,2);
	memcpy(&l_seq, rbin+20,4);
	l_read_name=((unsigned char*)rbin)[12];
	new_rbin[12] = rname_trimmed_len+1;
	memcpy(new_rbin+new_rbin_len, rbin+36, rname_trimmed_len);
	new_rbin[36+rname_trimmed_len]=0;
	new_rbin_len+= rname_trimmed_len+1;
	memcpy(new_rbin+new_rbin_len, rbin +36 + l_read_name, 4*n_cigar_op + l_seq + (l_seq+1)/2);
	new_rbin_len += 4*n_cigar_op + l_seq + (l_seq+1)/2;
	char * ext_bin_ptr = rbin + 36 + l_read_name +4*n_cigar_op + l_seq + (l_seq+1)/2;

#ifndef DO_STARSOLO_THING
	int CR_found=0, CB_found=0, CY_found=0, UR_found=0, UY_found=0, UB_found=0;
	while(ext_bin_ptr < rbin+binlen+4){
		char * tagstr = NULL; 
		int taglen = 0;
		if(ext_bin_ptr[0]=='C' && ext_bin_ptr[1]=='R' && ext_bin_ptr[2]=='Z'){
			CR_found = 1;
			tagstr = cellbc_seq;
			taglen = cct_context -> known_cell_barcode_length;
		}else if(ext_bin_ptr[0]=='C' && ext_bin_ptr[1]=='B' && ext_bin_ptr[2]=='Z'){
			CB_found = 1;
			tagstr = fixedbc_seq;
			taglen = cct_context -> known_cell_barcode_length;
		}else if(ext_bin_ptr[0]=='C' && ext_bin_ptr[1]=='Y' && ext_bin_ptr[2]=='Z'){
			CY_found = 1;
			tagstr = cellbc_qual;
			taglen = cct_context -> known_cell_barcode_length;
		}else if(ext_bin_ptr[0]=='U' && ext_bin_ptr[1]=='R' && ext_bin_ptr[2]=='Z'){
			UR_found = 1;
			tagstr = umi_seq;
			taglen = cct_context -> UMI_length;
		}else if(ext_bin_ptr[0]=='U' && ext_bin_ptr[1]=='B' && ext_bin_ptr[2]=='Z'){
			UB_found = 1;
			tagstr = fixedumi_seq;
			taglen = cct_context -> UMI_length;
		}else if(ext_bin_ptr[0]=='U' && ext_bin_ptr[1]=='Y' && ext_bin_ptr[2]=='Z'){
			UY_found = 1;
			tagstr = umi_qual;
			taglen = cct_context -> UMI_length;
		}
	
		if(tagstr){
			new_rbin[new_rbin_len++]=*(ext_bin_ptr++);
			new_rbin[new_rbin_len++]=*(ext_bin_ptr++);
			new_rbin[new_rbin_len++]=*(ext_bin_ptr++);
			int taglenold = strlen(ext_bin_ptr);
			memcpy(new_rbin+new_rbin_len,tagstr, taglen);
			*(new_rbin+new_rbin_len+taglen)=0;
			ext_bin_ptr += taglenold+1;
			new_rbin_len += taglen+1;
		}else{
			int content_len = SAP_pairer_skip_tag_body_len(ext_bin_ptr);
			memcpy(new_rbin + new_rbin_len, ext_bin_ptr, content_len );
			new_rbin_len += content_len;
			ext_bin_ptr += content_len;
		}
	}
	if(!CR_found){
		new_rbin[new_rbin_len++]='C';new_rbin[new_rbin_len++]='R';new_rbin[new_rbin_len++]='Z';
		memcpy(new_rbin+new_rbin_len, cellbc_seq, cct_context -> known_cell_barcode_length);
		*(new_rbin+new_rbin_len+cct_context -> known_cell_barcode_length)=0;
		new_rbin_len += cct_context -> known_cell_barcode_length+1;
	}
	if(fixedbc_seq && !CB_found){
		new_rbin[new_rbin_len++]='C';new_rbin[new_rbin_len++]='B';new_rbin[new_rbin_len++]='Z';
		memcpy(new_rbin+new_rbin_len, fixedbc_seq, cct_context -> known_cell_barcode_length);
		*(new_rbin+new_rbin_len+cct_context -> known_cell_barcode_length)=0;
		new_rbin_len += cct_context -> known_cell_barcode_length+1;
	}
	if(!CY_found){
		new_rbin[new_rbin_len++]='C';new_rbin[new_rbin_len++]='Y';new_rbin[new_rbin_len++]='Z';
		memcpy(new_rbin+new_rbin_len, cellbc_qual, cct_context -> known_cell_barcode_length);
		*(new_rbin+new_rbin_len+cct_context -> known_cell_barcode_length)=0;
		new_rbin_len += cct_context -> known_cell_barcode_length+1;
	}

	if(!UR_found){
		new_rbin[new_rbin_len++]='U';new_rbin[new_rbin_len++]='R';new_rbin[new_rbin_len++]='Z';
		memcpy(new_rbin+new_rbin_len, umi_seq, cct_context -> UMI_length);
		*(new_rbin+new_rbin_len+cct_context -> UMI_length)=0;
		new_rbin_len += cct_context -> UMI_length+1;
	}
	if(fixedumi_seq && !UB_found){
		new_rbin[new_rbin_len++]='U';new_rbin[new_rbin_len++]='B';new_rbin[new_rbin_len++]='Z';
		memcpy(new_rbin+new_rbin_len, fixedumi_seq, cct_context -> UMI_length);
		*(new_rbin+new_rbin_len+cct_context -> UMI_length)=0;
		new_rbin_len += cct_context -> UMI_length+1;
	}
	if(!UY_found){
		new_rbin[new_rbin_len++]='U';new_rbin[new_rbin_len++]='Y';new_rbin[new_rbin_len++]='Z';
		memcpy(new_rbin+new_rbin_len, umi_qual, cct_context -> UMI_length);
		*(new_rbin+new_rbin_len+cct_context -> UMI_length)=0;
		new_rbin_len += cct_context -> UMI_length+1;
	}
#endif

	new_rbin_len-=4;
	memcpy(new_rbin, &new_rbin_len,4);
	return new_rbin_len;
}

void cellCounts_do_one_batch_write_extend_rbin(cellcounts_global_t * cct_context, char * rbin, int binlen, FILE * fp, char * fixedbc_seq, char * fixedumi_seq, srInt_64 gene_no, srInt_64 * genes){
	char new_rbin[ binlen + 150 ]; // removed barcodes/qual from read names, add them to extra fields if they weren't there. Gene names are not put here.
	int new_rbin_len = cellCounts_make_barcode_bam_bin( cct_context, rbin, new_rbin, binlen, fixedbc_seq, fixedumi_seq, gene_no, genes );
	fwrite(new_rbin, 1, new_rbin_len+4, fp);
}

#ifdef __MINGW32__
#define ADD_key_FMT1 "%d-%d-%" PRId64 "-%s"
#else
#define ADD_key_FMT1 "%d-%d-%lld-%s"
#endif
#define ADD_key_struct { char my_key [50+MAX_UMI_LEN]; \
	SUBreadSprintf(my_key, 50+MAX_UMI_LEN,ADD_key_FMT1, sample_id, cell_no, gene_no, UMI_str); \
	srInt_64 supp_reads = HashTableGet(supp_reads_SCGU, my_key)-NULL; \
	if(1>supp_reads) HashTablePut(supp_reads_SCGU, strdup(my_key), NULL+1); \
	else HashTablePutReplaceEx(supp_reads_SCGU, my_key, NULL+supp_reads+1, 0,0,0); }

void * cellCounts_do_one_batch(void * paramsp1){
	srInt_64 x1;
	void ** params = paramsp1;
	cellcounts_global_t * cct_context = params[0];
	ArrayList * file_size_list = params[2];
	char *temp_dir = cct_context -> temp_file_dir;
	free(paramsp1);
	int me_max_Rbin_len = 0;
	int me_max_genes = 0;
	char ** bin_ptrs = malloc(sizeof(char*) * 1500000), * batch_content=NULL;
	int bin_ptr_size = 1500000;
	srInt_64 removed_UMIs = 0;
	while(1){
		int this_batch_no = -1;
		cellCounts_lock_occupy(&cct_context -> input_dataset_lock);
		//cellCounts_lock_occupy(&cct_context -> input_dataset_lock);
		if(cct_context -> do_one_batch_runner_current < CELLBC_BATCH_NUMBER +1){
			int this_batch_sorted_idx = (cct_context -> do_one_batch_runner_current ++);
			srInt_64 this_batch_size_and_no = ArrayListGet(file_size_list, file_size_list->numOfElements-1 -this_batch_sorted_idx)-NULL;
			this_batch_no = (int)(this_batch_size_and_no&0xfffffllu);
		}
		if(me_max_genes > cct_context -> barcode_batched_max_genes) cct_context -> barcode_batched_max_genes = me_max_genes;
		if(me_max_Rbin_len > cct_context-> barcode_batched_max_Rbin_len) cct_context-> barcode_batched_max_Rbin_len = me_max_Rbin_len;
		cellCounts_lock_release(&cct_context -> input_dataset_lock);
		//cellCounts_lock_release(&cct_context -> input_dataset_lock);
		if(0>this_batch_no)break;
		char tmp_fname[MAX_FILE_NAME_LENGTH+80];
		SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+80, "%s/temp-cellcounts-%06d-%03d.tmpbin", temp_dir, getpid(), this_batch_no);
		FILE * fp = fopen(tmp_fname, "rb");
		fseeko(fp, 0, SEEK_END);
		srInt_64 batch_fsize = ftello(fp);
		fseeko(fp, 0, SEEK_SET);
		if(batch_content==NULL) batch_content = malloc(batch_fsize);
		srInt_64 batch_content_len = fread(batch_content, 1, batch_fsize, fp);
		fclose(fp);
		if(batch_content_len!=batch_fsize){
			SUBREADprintf("ERROR: Cannot load file at once: %d!\n", this_batch_no);
			return NULL;
		}

		HashTable * supp_reads_SCGU = StringTableCreate(500000);
		HashTableSetDeallocationFunctions(supp_reads_SCGU, free, NULL);
		srInt_64 scanptr = 0;
		int rbin_no = 0;
		char UMI_str[MAX_UMI_LEN+1];

		while(scanptr < batch_content_len-1){
			int cell_no=0, sample_id=0;
			srInt_64 gene_no=0;
			if(bin_ptr_size<=rbin_no){
				bin_ptr_size = bin_ptr_size*2;
				bin_ptrs = realloc(bin_ptrs, sizeof(char*)*bin_ptr_size);
			}
			bin_ptrs[rbin_no] = batch_content+scanptr;
			memcpy(&sample_id, batch_content+scanptr, 4);
			scanptr += 4; // sample_ID 
			memcpy(&cell_no, batch_content+scanptr, 4);
			scanptr += 4; // cellbarcode_NO
			memcpy(&gene_no, batch_content+scanptr, 8);
			scanptr += 8; // gene_id 
			if(gene_no & (1LLU<<63)){
				int genes = (int)(gene_no & 0x7fffffffllu);
				if(genes > me_max_genes)me_max_genes=genes;
				memcpy(UMI_str, batch_content+scanptr+8*genes, cct_context -> UMI_length);
				UMI_str[cct_context -> UMI_length]=0;

				for(x1=0; x1<genes; x1++){
					memcpy(&gene_no, batch_content+scanptr, 8);
					scanptr += 8;
					ADD_key_struct;
				}
			}else{
				UMI_str[cct_context -> UMI_length]=0;
				memcpy(UMI_str, batch_content+scanptr, cct_context -> UMI_length);
				ADD_key_struct;
			}
			scanptr += cct_context -> UMI_length ; // UMI str
			int rbinlen = 0;
			memcpy(&rbinlen, batch_content+scanptr, 4);
			if(me_max_Rbin_len < rbinlen) me_max_Rbin_len = rbinlen;
			scanptr += rbinlen +4; // read_bin
			rbin_no++;
		}
		ArrayList ** cell_gene_umi_list = malloc(sizeof(void*)*cct_context -> sample_sheet_table -> numOfElements);
		for(x1 =0; x1< cct_context -> sample_sheet_table -> numOfElements; x1++){
			cell_gene_umi_list[x1]=ArrayListCreate(2000000);
			ArrayListSetDeallocationFunction(cell_gene_umi_list[x1], free);
		}
		supp_reads_SCGU -> appendix1 = cell_gene_umi_list;
		supp_reads_SCGU -> appendix2 = cct_context;
		supp_reads_SCGU -> counter1 = cct_context -> UMI_length;
		HashTableIteration(supp_reads_SCGU, cellCounts_do_one_batch_tab_to_struct_list);
		HashTable * filtered_SCGU_table = StringTableCreate(max(10000,cell_gene_umi_list[0] -> numOfElements / 10));
		HashTableSetDeallocationFunctions(filtered_SCGU_table, free, NULL);

		fp = fopen(tmp_fname, "wb");
		for(x1 = 0; x1 < cct_context -> sample_sheet_table -> numOfElements; x1++){
			HashTable * cellbcP0_to_geneno0B_P1_to_UMIs = HashTableCreate(500000);

			void * app1[3];
			cell_gene_umi_list[x1] -> appendix1 = app1;
			app1[0] = cct_context;
			app1[1] = NULL+1;
						// 1 : sorted by cell_bc, then gene, then supported_reads, then UMIstr (this is for step1 UMI merging)
						// 0 : sorted by cell_bc, then UMIstr, then supported_reads, then gene (this is for step2 UMI merging)
						// supported_reads : large -> small; the other: small -> large
			ArrayListSort(cell_gene_umi_list[x1],  cellCounts_do_one_batch_tab_to_struct_list_compare);
			cellCounts_do_one_batch_UMI_merge_one_step(cell_gene_umi_list[x1], 0, filtered_SCGU_table, NULL);

			app1[1] = NULL+0;
			app1[2] = cellbcP0_to_geneno0B_P1_to_UMIs;
			ArrayListSort(cell_gene_umi_list[x1], cellCounts_do_one_batch_tab_to_struct_list_compare);

			cellCounts_do_one_batch_UMI_merge_one_step(cell_gene_umi_list[x1], 1, filtered_SCGU_table, &removed_UMIs);

			cellbcP0_to_geneno0B_P1_to_UMIs -> appendix1 = fp;

			fwrite(&cellbcP0_to_geneno0B_P1_to_UMIs -> numOfElements,1,8,fp);
			HashTableIteration(cellbcP0_to_geneno0B_P1_to_UMIs, cellCounts_do_one_batch_write_UMIs);
			HashTableDestroy(cellbcP0_to_geneno0B_P1_to_UMIs);
		}

		void * sort_base[2];
		sort_base[0] = bin_ptrs;
		sort_base[1] = cct_context;
		merge_sort(sort_base, rbin_no, cellCounts_do_one_batch_sort_compare, cellCounts_do_one_batch_sort_exchange, cellCounts_do_one_batch_sort_merge);

		for(x1 = 0; x1 < rbin_no; x1++){
			char * binptr = bin_ptrs[x1];
			int cellid =0, sampleid = 0;
			srInt_64 gene_no =0, genes = 0, geneno_0 = 0;
			char * umi, * glist_ptr =NULL;
			memcpy(&sampleid, binptr, 4);
			memcpy(&cellid, binptr+4, 4);
			memcpy(&gene_no, binptr+8, 8);
			if(gene_no & (1LLU<<63)){
				glist_ptr =binptr + 16;
				genes = (int)(gene_no & 0x7fffffff);
				memcpy(&geneno_0, binptr+16, 8);
			}
			umi = binptr + 16 + 8*genes;
			char SCGU_key [40+MAX_UMI_LEN];

#ifdef __MINGW32__
			int keyptr = SUBreadSprintf(SCGU_key, 40+MAX_UMI_LEN,"%d-%d-%" PRId64 "-", sampleid, cellid,  (gene_no & (1LLU<<63))? geneno_0: gene_no);
#else
			int keyptr = SUBreadSprintf(SCGU_key, 40+MAX_UMI_LEN,"%d-%d-%lld-", sampleid, cellid,  (gene_no & (1LLU<<63))? geneno_0: gene_no);
#endif
			memcpy(SCGU_key+keyptr, umi, cct_context -> UMI_length);
			SCGU_key[keyptr+ cct_context -> UMI_length] = 0;

			char * new_UMI = HashTableGet(filtered_SCGU_table, SCGU_key);
			if(new_UMI) umi = new_UMI;
			if(umi == NULL-1) umi="-----------------------------------------";
			fwrite(&sampleid, 1, 4, fp);
			fwrite(&cellid, 1, 4, fp);
			fwrite(&gene_no, 1, 8, fp);
			if(gene_no & (1LLU<<63)) fwrite( glist_ptr, 1, 8*genes, fp );
			fwrite(umi,1, cct_context -> UMI_length, fp);
			int binlen;

			memcpy(&binlen, binptr+16+8*genes+cct_context -> UMI_length,4 );
			char * new_cellbc = NULL;
			if(cellid>=0)new_cellbc = ArrayListGet(cct_context -> cell_barcodes_array, cellid);
			cellCounts_do_one_batch_write_extend_rbin(cct_context, binptr+16+8*genes+cct_context ->UMI_length, binlen, fp, new_cellbc, umi[0]=='-'?NULL:umi, gene_no, (srInt_64*)glist_ptr);
		}
		fclose(fp);
		HashTableDestroy(supp_reads_SCGU);
		HashTableDestroy(filtered_SCGU_table);
		for(x1 =0; x1< cct_context -> sample_sheet_table -> numOfElements; x1++)ArrayListDestroy(cell_gene_umi_list[x1]);
		free(cell_gene_umi_list);
	}
	free(batch_content);
	free(bin_ptrs);
	return NULL + removed_UMIs;
}

#ifdef DO_STARSOLO_THING
#define DO_CREATE_BAI_FOR_BAM 0
#else
#define DO_CREATE_BAI_FOR_BAM 1
#endif

void cellCounts_save_BAM_result(cellcounts_global_t * cct_context, struct scRNA_merge_batches_worker_current * finished_job){
	if(!finished_job -> task)return;
	if(cct_context -> is_BAM_and_FQ_out_generated){
		int sample_id = finished_job -> task -> sample_id;
		void ** fps = HashTableGet(cct_context -> sample_BAM_writers, NULL+sample_id);
		simple_bam_writer * wtr = (*fps);
		int inbin_pos = 0, outblocki = 0;
		int nextoffset = -1;
		if(CELLCOUNTS_BAMBLOCK_COMP_NUMBER>1)nextoffset = finished_job -> task -> inbin_batch_start_offsets[1];
		int block_number_this = finished_job -> task -> block_number - finished_job -> task -> inbin_number +1;
		if(DO_CREATE_BAI_FOR_BAM) while(inbin_pos < finished_job -> task -> inbin_len){
			int binlen = 0;
			if(outblocki < finished_job -> task -> inbin_number -1 &&  inbin_pos == nextoffset ){
				outblocki ++;
				if(outblocki < finished_job -> task -> inbin_number -1 && CELLCOUNTS_BAMBLOCK_COMP_NUMBER>1)nextoffset = finished_job -> task -> inbin_batch_start_offsets[outblocki +1];
				block_number_this = finished_job -> task -> block_number - (finished_job -> task -> inbin_number -1 - outblocki);
			}
			binlen=*(int*)(finished_job -> task -> inbin+inbin_pos);
			simple_bam_writer_update_index(wtr, finished_job -> task -> inbin+inbin_pos, binlen, block_number_this, inbin_pos);
			inbin_pos += 4+binlen;
		}

		for(outblocki=0; outblocki < finished_job -> task -> inbin_number ; outblocki ++){
			int block_number_this = finished_job -> task -> block_number - (finished_job -> task -> inbin_number -1 - outblocki);
			int inblock_size = -1;
			if(outblocki  == finished_job -> task -> inbin_number -1) inblock_size = finished_job -> task -> inbin_len - finished_job -> task -> inbin_batch_start_offsets[outblocki];
			else if(CELLCOUNTS_BAMBLOCK_COMP_NUMBER>1)inblock_size = finished_job -> task -> inbin_batch_start_offsets[outblocki +1] -  finished_job -> task -> inbin_batch_start_offsets[outblocki];
			simple_bam_write_compressed_block(wtr, finished_job -> outbin + CELLRANGER_MERGER_WORKER_BINSIZE*outblocki , finished_job -> outbin_len[outblocki ], inblock_size, finished_job -> crc32[outblocki], block_number_this);
		}
	}
	finished_job -> task = NULL;
}

void cellCounts_merged_write_sparse_unique_genes(void * ky, void * va, HashTable * tab){
	HashTable * unique_geneno1B_tab = tab -> appendix1;
	HashTable * used_cellnoP1_tab = tab -> appendix2;

	int cellbcP1 = ky-NULL;
	if(used_cellnoP1_tab && !HashTableGet(used_cellnoP1_tab, NULL+cellbcP1))return;
	HashTable * g2u = va;
	ArrayList * g2ul = HashTableKeys(g2u);
	int x1;
	for(x1=0; x1<g2ul->numOfElements; x1++){
		void *geneno1B_ptr = ArrayListGet(g2ul,x1);
		if(!HashTableGet(unique_geneno1B_tab, ArrayListGet(g2ul,x1))) HashTablePut(unique_geneno1B_tab, geneno1B_ptr, NULL+1);
		tab -> counter1 += HashTableGet(g2u, geneno1B_ptr)?1:0;
	}
	ArrayListDestroy(g2ul);
}



int cellCounts_merged_write_sparse_matrix(cellcounts_global_t * cct_context, HashTable * cellP1_to_geneP1_to_umis_tab, HashTable * cellnoP1_to_umis_tab, ArrayList * used_cell_barcodes, int sample_index, char * tabtype, ArrayList * loaded_features, HashTable * sorted_order_p1_to_i_p1_tab){
	int x1,x2;

	char ofname[MAX_FILE_NAME_LENGTH + 100];
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 100,"%s.scRNA.%03d.%s.summary",cct_context->output_prefix, sample_index+1,tabtype);
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 100,"%s.scRNA.%03d.%s.BCtab",cct_context->output_prefix, sample_index+1,tabtype);
	FILE * ofp_bcs = fopen( ofname , "w" );
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 100,"%s.scRNA.%03d.%s.GENEtab",cct_context->output_prefix, sample_index+1,tabtype);
	FILE * ofp_genes = fopen( ofname , "w" );
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 100,"%s.scRNA.%03d.%s.spmtx",cct_context->output_prefix, sample_index+1,tabtype);
	FILE * ofp_mtx = fopen( ofname , "w" );
	fprintf(ofp_mtx,"%%%%MatrixMarket matrix coordinate integer general\n");

	HashTable * used_cellnoP1_tab = ArrayListToLookupTable_Int(used_cell_barcodes);
	HashTable * unique_NZ_geneno1B_table = HashTableCreate(10000);
	cellP1_to_geneP1_to_umis_tab -> counter1 = 0;
	cellP1_to_geneP1_to_umis_tab -> appendix1 = unique_NZ_geneno1B_table;
	cellP1_to_geneP1_to_umis_tab -> appendix2 = used_cellnoP1_tab;
	HashTableIteration(cellP1_to_geneP1_to_umis_tab, cellCounts_merged_write_sparse_unique_genes);
	srInt_64 total_lines = cellP1_to_geneP1_to_umis_tab -> counter1;
	ArrayList * unique_NZ_genenosP1_list = HashTableKeys(unique_NZ_geneno1B_table);
	HashTableDestroy(unique_NZ_geneno1B_table);
	HashTableDestroy(used_cellnoP1_tab);
	ArrayListSort(unique_NZ_genenosP1_list, NULL);

	#ifdef __MINGW32__
	fprintf(ofp_mtx, "%" PRId64 " %" PRId64 " %" PRId64 "\n", unique_NZ_genenosP1_list -> numOfElements , used_cell_barcodes -> numOfElements,  total_lines );
	#else
	fprintf(ofp_mtx, "%lld %lld %lld\n", unique_NZ_genenosP1_list -> numOfElements , used_cell_barcodes -> numOfElements,  total_lines );
	#endif

	for(x2=0; x2 < unique_NZ_genenosP1_list -> numOfElements; x2++){
		int gene_index_0B = ArrayListGet(unique_NZ_genenosP1_list, x2) - NULL-1;
		char* gene_name = (char*)cct_context -> gene_name_array [gene_index_0B];
		fprintf(ofp_genes,"%s\n", gene_name);
	}

	for(x1 = 0; x1 < used_cell_barcodes -> numOfElements; x1++){
		srInt_64 cellno = ArrayListGet(used_cell_barcodes, x1)-NULL;
		char * cellbc_seq = ArrayListGet(cct_context -> cell_barcodes_array, cellno);
		fprintf(ofp_bcs,"%s\n", cellbc_seq);
	}

	for(x1 = 0; x1 < used_cell_barcodes -> numOfElements; x1++){
		srInt_64 cellno = ArrayListGet(used_cell_barcodes, x1)-NULL;
		HashTable * geneno1B_to_UMIs = HashTableGet(cellP1_to_geneP1_to_umis_tab, NULL+1+cellno);

		for(x2=0; x2 < unique_NZ_genenosP1_list -> numOfElements; x2++){
			int geneno1B = ArrayListGet(unique_NZ_genenosP1_list, x2)-NULL;
			int this_umis = HashTableGet(geneno1B_to_UMIs, NULL+geneno1B) -NULL;
			if(this_umis>0)fprintf(ofp_mtx,"%d %d %d\n", x2+1, x1+1, this_umis);
		}
	}
	ArrayListDestroy(unique_NZ_genenosP1_list);
	fclose(ofp_bcs);
	fclose(ofp_genes);
	fclose(ofp_mtx);

	return 0;
}

void cellCounts_merged_45K_to_90K_sum_SUM_Level2(void * GeneNo1B, void * vUMIs, HashTable * m2){
	HashTable * summed_gene_to_umis = m2 -> appendix1;
	HashTablePut(summed_gene_to_umis, GeneNo1B, vUMIs + (HashTableGet(summed_gene_to_umis, GeneNo1B)-NULL));
}

void cellCounts_merged_45K_to_90K_sum_SUM(void * keyCellNoP1, void * Vgno_umi_tab, HashTable * me){
	HashTable * summed_gene_to_umis  = me -> appendix1;
	HashTable * bcid_look_tab = me -> appendix2;
	HashTable * geneno1B_to_UMIs_tab = Vgno_umi_tab;
	if(!HashTableGet(bcid_look_tab, keyCellNoP1))return;
	geneno1B_to_UMIs_tab -> appendix1 = summed_gene_to_umis;
	HashTableIteration(geneno1B_to_UMIs_tab ,cellCounts_merged_45K_to_90K_sum_SUM_Level2 );
}

void cellCounts_merged_45K_to_90K_sum_WRT(void * kyGeneID, void * valUMIs, HashTable * me){
	cellcounts_global_t * cct_context = me -> appendix1;
	
	FILE * ofp = me -> appendix2;

	unsigned char * gene_name = cct_context -> gene_name_array[ kyGeneID - NULL-1 ];
	fprintf(ofp, "%s\t%u\n", gene_name, (unsigned int) (valUMIs-NULL));
}

void cellCounts_merged_45K_to_90K_sum(cellcounts_global_t * cct_context, HashTable * cellP1_geneP1_UMIs_tab, ArrayList * bcid_P0_arr, int sample_no, ArrayList * loaded_features, HashTable * sorted_index_p1_to_i_p1_tab){
	HashTable * summed_gene_to_umis = HashTableCreate( 3+cellP1_geneP1_UMIs_tab->numOfElements/6 );
	HashTable * bcid_look_tab = ArrayListToLookupTable_Int(bcid_P0_arr);
	cellP1_geneP1_UMIs_tab -> appendix1 = summed_gene_to_umis;
	cellP1_geneP1_UMIs_tab -> appendix2 = bcid_look_tab;
	cellP1_geneP1_UMIs_tab -> appendix3 = cct_context;
	HashTableIteration( cellP1_geneP1_UMIs_tab, cellCounts_merged_45K_to_90K_sum_SUM );

	char ofname[MAX_FILE_NAME_LENGTH + 100];
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 100,"%s.scRNA.%03d.AmbSum",cct_context->output_prefix, sample_no+1);
	FILE * write_fp = fopen(ofname,"w");
	fprintf(write_fp,"GeneID\tUMIs\n");
	summed_gene_to_umis -> appendix1 = cct_context;
	summed_gene_to_umis -> appendix2 = write_fp;
	void * vp2[2];
	vp2[0]=loaded_features;
	vp2[1]=sorted_index_p1_to_i_p1_tab;
	summed_gene_to_umis -> appendix3 = vp2;
	summed_gene_to_umis -> counter1 = sample_no;
	HashTableIteration( summed_gene_to_umis, cellCounts_merged_45K_to_90K_sum_WRT );
	HashTableDestroy(bcid_look_tab);
	HashTableDestroy(summed_gene_to_umis);
	fclose(write_fp);
}

void cellCounts_merged_write_nozero_geneids_WRT(void *k, void *v, HashTable* me){
	FILE * fp = me->appendix1;
	cellcounts_global_t * cct_context = me->appendix2;
	unsigned char* gene_symbol = cct_context -> gene_name_array [k-NULL-1];
	fprintf(fp, "%s\n", gene_symbol);
}

void cellCounts_merged_write_nozero_geneids(cellcounts_global_t * cct_context, HashTable * no0genes, int samplenno, ArrayList * loaded_features, HashTable * sorted_order_p1_to_i_p1_tab){
	char ofname[MAX_FILE_NAME_LENGTH + 100];
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 100,"%s.scRNA.%03d.no0Genes", cct_context->output_prefix, samplenno+1);
	FILE * fp = fopen( ofname , "w" );
	no0genes -> appendix1 =fp;
	void * tv2[2];
	no0genes -> appendix2 = cct_context;
	tv2[0]=loaded_features;
	tv2[1]=sorted_order_p1_to_i_p1_tab;
	no0genes -> appendix3 =tv2;
	HashTableIteration(no0genes, cellCounts_merged_write_nozero_geneids_WRT);
	fclose(fp);
}

void cellCounts_merged_to_tables_write_build_UMIcount_in(void * ky, void * val, HashTable * tab){
	tab -> counter1 += (val-NULL);
}

void cellCounts_merged_to_tables_write_build_UMIcounts(void * ky, void * val, HashTable * tab){
	HashTable * cellbcP1_to_umis_tab = tab -> appendix1;
	int cell_no = ky-NULL-1;
	HashTable * geneP1_to_counts_tab = val;

	geneP1_to_counts_tab -> counter1 = 0;
	HashTableIteration(geneP1_to_counts_tab, cellCounts_merged_to_tables_write_build_UMIcount_in);
	HashTablePut(cellbcP1_to_umis_tab, NULL+1+cell_no, NULL+geneP1_to_counts_tab -> counter1);
}

#define MIN_UMIS_FOR_CANDIDATE_RESCUE 500
#define SCRNA_AMBIENT_RESCURE_MEDIAN_FRACTION 0.01
#define MAX_CANDIDATE_CELLS 20000
void cellCounts_merged_ambient_rescure(cellcounts_global_t * cct_context, HashTable * cellP1_to_geneP1_to_umis_tab, HashTable * cellnoP1_to_umis_tab, ArrayList * this_sample_45k_90k_barcode_no_P0, ArrayList * this_sample_ambient_rescure_candi, ArrayList * highconf_cellbc_list){
	ArrayList * sorted_bcno_p1 = HashTableSortedIndexes( cellnoP1_to_umis_tab, 1);
	HashTable * highconf_cellbc_list_tab = ArrayListToLookupTable_Int(highconf_cellbc_list);
	srInt_64 x1, high_conf_cells = 0;
	for(x1=0; x1 < sorted_bcno_p1 -> numOfElements; x1++){
		void * this_bc_pnt = ArrayListGet(sorted_bcno_p1 ,  x1);
		if(HashTableGet(highconf_cellbc_list_tab, this_bc_pnt)) high_conf_cells = x1+1;
		else break; // assuming that all high-umi barcodes are high-confident, this makes x1 being the # of total high-confidence barcodes.	
	}
	if(high_conf_cells >0){
		srInt_64 median_umis = HashTableGet(cellnoP1_to_umis_tab, ArrayListGet(sorted_bcno_p1 ,  (high_conf_cells-1)/2))-NULL;
		srInt_64 median_umis_001_cut = (srInt_64)(median_umis *1. *SCRNA_AMBIENT_RESCURE_MEDIAN_FRACTION +0.50000001);
		for(x1=0; x1 < sorted_bcno_p1 -> numOfElements; x1++){
			void * this_bc_pnt_p1 = ArrayListGet(sorted_bcno_p1 ,  x1);
			if(HashTableGet(highconf_cellbc_list_tab, this_bc_pnt_p1)){
				continue; // it is in high-conf list
			}
			srInt_64 this_bc_umis = HashTableGet(cellnoP1_to_umis_tab, this_bc_pnt_p1) - NULL;
			if(this_bc_umis < median_umis_001_cut) break;
			if(this_bc_umis < MIN_UMIS_FOR_CANDIDATE_RESCUE) break;
			if(x1 >= 45000) break;

			// The MAX_CANDIDATE_CELLS was added on 31 JAN 2023 when I realised that Cell Ranger only tests 20,000 candidates.
			if(this_sample_ambient_rescure_candi -> numOfElements < MAX_CANDIDATE_CELLS)
				ArrayListPush(this_sample_ambient_rescure_candi, this_bc_pnt_p1-1);
		}
	}
	for(x1=45000; x1 < sorted_bcno_p1 -> numOfElements; x1++){
		if(x1 >= 90000) break;
		ArrayListPush(this_sample_45k_90k_barcode_no_P0, ArrayListGet(sorted_bcno_p1 ,  x1)-1 );
	}
	ArrayListDestroy(sorted_bcno_p1);
	HashTableDestroy(highconf_cellbc_list_tab);
}



#define SCRNA_BOOTSTRAP_HIGH_INDEX 30
#define SCRNA_BOOTSTRAP_SAMPLING_TIMES 100

int cellCounts_merged_bootstrap_a_sample(cellcounts_global_t * cct_context, HashTable * cellP1_to_geneP1_to_umis_tab, HashTable * cellnoP1_to_umis_tab, ArrayList * highconf_cellbc_list){
	ArrayList * sorted_idx = HashTableSortedIndexes( cellnoP1_to_umis_tab, 1);
	srInt_64 x2, x1;
	float cellCounts_umi_cutoff = cct_context -> umi_cutoff;

	#define SCRNA_IDX_PRIME_NUMBER_BIG 11218439llu;
	srInt_64 this_total = 0, seed_rand = sorted_idx -> numOfElements/2;

	int last_umi_no= -1;
	if(cellCounts_umi_cutoff >= 0.0){
		for(x1 = 0; x1 < sorted_idx -> numOfElements ; x1++){
			void * cellbc_p1_ptr = ArrayListGet(sorted_idx,x1);
			srInt_64 this_umis = HashTableGet(cellnoP1_to_umis_tab, cellbc_p1_ptr )-NULL;
			if(this_umis >= cellCounts_umi_cutoff-0.1){
				ArrayListPush(highconf_cellbc_list, ArrayListGet( sorted_idx, x1 ) - 1 );
				last_umi_no = this_umis;
			}else break;	// #UMI-sorted so no need to scan more
		}
	}else{
		for(x1 = 0; x1 < SCRNA_BOOTSTRAP_SAMPLING_TIMES; x1++){
			ArrayList * resampled_list_of_umis = ArrayListCreate( sorted_idx->numOfElements );
			for(x2 = 0; x2 < sorted_idx -> numOfElements ; x2++){
				seed_rand %= sorted_idx -> numOfElements;
				void * cellbc_p1_ptr = ArrayListGet(sorted_idx, seed_rand);
				seed_rand += SCRNA_IDX_PRIME_NUMBER_BIG;
				srInt_64 this_umis = HashTableGet( cellnoP1_to_umis_tab, cellbc_p1_ptr )-NULL;
				ArrayListPush(resampled_list_of_umis,NULL+this_umis);
			}
			ArrayListSort( resampled_list_of_umis, NULL );
			srInt_64 UMIs_30th_div10 = ArrayListGet(resampled_list_of_umis, resampled_list_of_umis -> numOfElements - SCRNA_BOOTSTRAP_HIGH_INDEX) -NULL;
			UMIs_30th_div10 = (srInt_64)(UMIs_30th_div10*1./10 + 0.500000001);

			for(x2 =0; x2< resampled_list_of_umis -> numOfElements; x2++){
				srInt_64 lli = resampled_list_of_umis -> numOfElements -1 -x2;
				srInt_64 this_umis = ArrayListGet(resampled_list_of_umis, lli)-NULL;
				if(this_umis >= UMIs_30th_div10) this_total ++;
				else break;
			}
			ArrayListDestroy(resampled_list_of_umis);
		}
		double total_f = this_total*1. / SCRNA_BOOTSTRAP_SAMPLING_TIMES;
		this_total = (int)(total_f + 0.500000001);

		void * last_ptr =NULL;
		for(x1 = 0; x1 < min(sorted_idx -> numOfElements, this_total) ; x1++){
			last_ptr = ArrayListGet( sorted_idx, x1 );
			ArrayListPush(highconf_cellbc_list, last_ptr - 1 );
		}
		last_umi_no = HashTableGet(cellnoP1_to_umis_tab ,last_ptr)-NULL;
	}
	ArrayListDestroy(sorted_idx);
	return last_umi_no;
}

void cellCounts_merged_to_tables_write(cellcounts_global_t * cct_context, HashTable ** cellP1_to_geneP1_to_umis, ArrayList * loaded_features, srInt_64 nexons){
	char ofname[MAX_FILE_NAME_LENGTH + 20];
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 20,"%s.scRNA.SampleTable",cct_context->output_prefix);
	FILE * sample_tab_fp = fopen( ofname , "w" );
	int x1;

	fprintf(sample_tab_fp,"SampleName\tUMICutoff\tTotalReads\tMappedReads\tAssignedReads\tIndex\n");
	for(x1 = 0; x1 < cct_context -> sample_sheet_table -> numOfElements ; x1++){
		ArrayList * high_confid_barcode_index_list = ArrayListCreate(20000);
		ArrayList * this_sample_ambient_rescure_candi = ArrayListCreate(10000);
		ArrayList * this_sample_45k_90k_barcode_no_P0 = ArrayListCreate(90000 - 45000 + 100);

		HashTable * cellbcP1_to_umis_tab = HashTableCreate(cellP1_to_geneP1_to_umis[x1] -> numOfElements);
		cellP1_to_geneP1_to_umis[x1] -> appendix1 = cellbcP1_to_umis_tab;
		HashTableIteration(cellP1_to_geneP1_to_umis[x1], cellCounts_merged_to_tables_write_build_UMIcounts);

		int applied_umi_cut = -1;
		if(cellbcP1_to_umis_tab -> numOfElements >0) applied_umi_cut = cellCounts_merged_bootstrap_a_sample(cct_context, cellP1_to_geneP1_to_umis[x1], cellbcP1_to_umis_tab, high_confid_barcode_index_list);
		cct_context -> applied_umi_cut[x1] = applied_umi_cut;
		cellCounts_merged_ambient_rescure(cct_context, cellP1_to_geneP1_to_umis[x1], cellbcP1_to_umis_tab, this_sample_45k_90k_barcode_no_P0, this_sample_ambient_rescure_candi, high_confid_barcode_index_list);

		int umi_cutoff = cct_context -> applied_umi_cut[x1];
		char * this_sample_name = ArrayListGet( cct_context -> sample_id_to_name, x1);
#ifdef __MINGW32__
		fprintf(sample_tab_fp,"%s\t%d\t%" PRId64 "\t%" PRId64 "\t%" PRId64 "\t%d\n", this_sample_name, umi_cutoff,  cct_context -> reads_per_sample[x1],  cct_context -> mapped_reads_per_sample[x1],  cct_context -> assigned_reads_per_sample[x1] ,x1+1);
#else
		fprintf(sample_tab_fp,"%s\t%d\t%lld\t%lld\t%lld\t%d\n", this_sample_name, umi_cutoff, cct_context -> reads_per_sample[x1],  cct_context -> mapped_reads_per_sample[x1],  cct_context -> assigned_reads_per_sample[x1], x1+1);
#endif
		srInt_64 xk1;
		HashTable * sorted_order_p1_to_i_p1_tab = HashTableCreate(nexons/4);
		for(xk1 = 0; xk1 < nexons ; xk1++){
			fc_feature_info_t * feature1 = ArrayListGet(loaded_features, xk1);
			HashTablePut(sorted_order_p1_to_i_p1_tab, NULL+ feature1->sorted_order+1 , NULL+xk1+1 );
		}

		if(cct_context -> report_excluded_barcodes){
			ArrayList * all_barcode_index_list = HashTableKeys(cellbcP1_to_umis_tab);
			for(xk1=0; xk1<all_barcode_index_list->numOfElements; xk1++) all_barcode_index_list->elementList[xk1]--;
			cellCounts_merged_write_sparse_matrix(cct_context, cellP1_to_geneP1_to_umis[x1], cellbcP1_to_umis_tab, all_barcode_index_list, x1, "RawOut", loaded_features, sorted_order_p1_to_i_p1_tab);
			ArrayListDestroy(all_barcode_index_list);
		}
		cellCounts_merged_write_sparse_matrix(cct_context, cellP1_to_geneP1_to_umis[x1], cellbcP1_to_umis_tab, high_confid_barcode_index_list, x1, "HighConf",  loaded_features, sorted_order_p1_to_i_p1_tab);
		cellCounts_merged_write_sparse_matrix(cct_context, cellP1_to_geneP1_to_umis[x1], cellbcP1_to_umis_tab, this_sample_ambient_rescure_candi, x1, "RescCand",  loaded_features, sorted_order_p1_to_i_p1_tab);
		cellCounts_merged_45K_to_90K_sum( cct_context, cellP1_to_geneP1_to_umis[x1], this_sample_45k_90k_barcode_no_P0, x1 , loaded_features, sorted_order_p1_to_i_p1_tab);
		HashTable * no0genes = HashTableCreate(10000);
		cellP1_to_geneP1_to_umis[x1] -> appendix1 = no0genes;
		cellP1_to_geneP1_to_umis[x1] -> appendix2 = NULL;
		HashTableIteration(cellP1_to_geneP1_to_umis[x1], cellCounts_merged_write_sparse_unique_genes);
		cellCounts_merged_write_nozero_geneids(cct_context, no0genes, x1, loaded_features, sorted_order_p1_to_i_p1_tab);

		HashTableDestroy(no0genes);
		ArrayListDestroy(this_sample_ambient_rescure_candi);
		ArrayListDestroy(this_sample_45k_90k_barcode_no_P0);
		ArrayListDestroy(high_confid_barcode_index_list);
		HashTableDestroy(cellbcP1_to_umis_tab);
		HashTableDestroy(sorted_order_p1_to_i_p1_tab);
	}

	fclose(sample_tab_fp);

}

int bin_file_exists(char *filename) {
	struct stat   buffer;   
	return (stat (filename, &buffer) == 0);
}

void * delete_file_thread(void * arg){
	void ** ptrs =arg;
	srInt_64 *current_sorting_key = NULL;
	cellcounts_global_t * cct_context = NULL;
	current_sorting_key = ptrs[0];
	cct_context = ptrs[1];
	
	while(1){
		int all_closed = 1;
		int x1;
		for(x1=0; x1<CELLBC_BATCH_NUMBER +2; x1++){
			if(NULL== current_sorting_key || current_sorting_key[x1] == 0x7fffffffffffffffLLU){
				char tmp_fname[MAX_FILE_NAME_LENGTH+50];
				SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+50, "%s/temp-cellcounts-%06d-%03d.tmpbin", cct_context -> temp_file_dir, getpid(), x1);
				if(bin_file_exists(tmp_fname)) unlink(tmp_fname);
			}else all_closed = 0;
		}
		if(all_closed ) break;
		sleep(2);
	}
	if(current_sorting_key)free(current_sorting_key);
	return NULL;
}

int cellCounts_do_cellbc_batches(cellcounts_global_t * cct_context){
	pthread_t *threads = malloc(sizeof(pthread_t)* cct_context-> total_threads);
	int sample_i,xk2,xk1,compress_workers = max(1, cct_context-> total_threads-1);
	HashTable * cellnoP1_to_genenoP1_to_UMIs[cct_context -> sample_sheet_table -> numOfElements];
	struct scRNA_merge_batches_worker_task * task_buffers = malloc(sizeof(struct scRNA_merge_batches_worker_task) * (1+compress_workers)* cct_context->sample_sheet_table -> numOfElements);
	int current_filling_worker_per_sample [cct_context-> sample_sheet_table -> numOfElements];
	struct scRNA_merge_batches_worker_current * worker_current_jobs = calloc(sizeof(struct scRNA_merge_batches_worker_current), compress_workers);

	ArrayList * file_size_list = ArrayListCreate(CELLBC_BATCH_NUMBER +1);
	for(xk1=0; xk1<CELLBC_BATCH_NUMBER +2; xk1++){
		if(xk1<CELLBC_BATCH_NUMBER +1){
			srInt_64 batchsize = ftello( cct_context -> batch_files[xk1]);
			ArrayListPush(file_size_list, NULL+( batchsize<<20 | xk1));
		}
		fclose(cct_context -> batch_files[xk1]);
	}
	ArrayListSort(file_size_list, NULL);

	srInt_64 block_numbers_current [cct_context-> sample_sheet_table -> numOfElements];
	for(xk1=0; xk1< cct_context-> sample_sheet_table -> numOfElements; xk1++){
		cellnoP1_to_genenoP1_to_UMIs[xk1] = HashTableCreate(10000);
		HashTableSetDeallocationFunctions(cellnoP1_to_genenoP1_to_UMIs[xk1], NULL,(void (*) (void*))HashTableDestroy);
		current_filling_worker_per_sample[xk1] = 0;
		task_buffers[xk1].inbin_len = 0;
		task_buffers[xk1].inbin_number = 0;
		task_buffers[xk1].inbin_batch_start_offsets[0]=0;
		block_numbers_current[xk1] = 0;
	}

	for(xk1=0; xk1<compress_workers+1; xk1++)for(xk2 = 0; xk2 < cct_context-> sample_sheet_table -> numOfElements; xk2++) task_buffers[xk1*cct_context->sample_sheet_table -> numOfElements + xk2].sample_id = xk2+1;

	for(xk1=0; xk1< cct_context-> total_threads; xk1++){
		void ** vpp = malloc(sizeof(void*)*3);
		vpp[0] = cct_context;
		vpp[1] = NULL + xk1;
		vpp[2] = file_size_list;
		pthread_create(threads + xk1, NULL, cellCounts_do_one_batch, vpp);
	}

	srInt_64 removed_umis = 0;
	for(xk1=0; xk1<cct_context-> total_threads; xk1++){
		void * pret = NULL;
		pthread_join(threads[xk1], &pret);
		removed_umis += (pret - NULL);
	}
	//SUBREADprintf("After processing batches, %lld UMIs were removed in step2 of UMI merging.\n", removed_umis);
	print_in_box(80,0,0,"");
	if(cct_context -> input_mode== GENE_INPUT_BCL){
		srInt_64 all_extracted_reads = 0;
		for(xk1 = 0; xk1 < cct_context-> sample_sheet_table -> numOfElements+1; xk1++) 
			all_extracted_reads += cct_context-> reads_per_sample[xk1];

		for(xk1 = 0; xk1 < cct_context-> sample_sheet_table -> numOfElements; xk1++) {
			srInt_64 extracted_reads = cct_context-> reads_per_sample[xk1];
			char * sample_name = ArrayListGet(cct_context-> sample_id_to_name, xk1);
#ifdef __MINGW32__
			print_in_box(81,0,0,"  % 13" PRId64 " (%4.1f%%%%) reads were assigned to %s.\n", extracted_reads, extracted_reads*100./all_extracted_reads, sample_name);
#else
			print_in_box(81,0,0,"  %'13lld (%4.1f%%%%) reads were assigned to %s.\n", extracted_reads, extracted_reads*100./all_extracted_reads, sample_name);
#endif
		}
		print_in_box(80,0,0,"");

#ifdef __MINGW32__
		if(cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements] < 0.005*all_extracted_reads) print_in_box(81,0,0,"  % 13" PRId64 "(%4.0f%%%%) reads were assigned to samples in total.", all_extracted_reads - cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements], 100.-cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements]*100./all_extracted_reads);
		else print_in_box(81,0,0,"  % 13" PRId64 " (%4.1f%%%%) reads were assigned to samples in total.", all_extracted_reads - cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements], 100.-cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements]*100./all_extracted_reads);
#else
		if(cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements] < 0.005*all_extracted_reads) print_in_box(81,0,0,"  %'13lld (%4.0f%%%%) reads were assigned to samples in total.", all_extracted_reads - cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements], 100.-cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements]*100./all_extracted_reads);
		else print_in_box(81,0,0,"  %'13lld (%4.1f%%%%) reads were assigned to samples in total.", all_extracted_reads - cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements], 100.-cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements]*100./all_extracted_reads);
#endif


	print_in_box(80,0,0,"");
	}
	print_in_box(80,0,0,"Generate UMI count tables...");
	ArrayListDestroy(file_size_list);

	worker_master_mutex_t worker_mut;
	worker_master_mutex_init(&worker_mut, max(1, cct_context-> total_threads-1));

	for(xk1=0; xk1<max(1, cct_context-> total_threads-1); xk1++){
		void ** vpp = malloc(sizeof(void*)*4);
		vpp[0] = cct_context;
		vpp[1] = &worker_mut;
		vpp[2] = NULL + xk1;
		vpp[3] = worker_current_jobs + xk1;
		pthread_create(threads + xk1, NULL, cellCounts_merge_batches_worker, vpp);
	}

	FILE * input_fps[CELLBC_BATCH_NUMBER+2];
	char * last_rbin_buffer[CELLBC_BATCH_NUMBER+1];
	srInt_64 * current_sorting_key = malloc(sizeof(srInt_64)*(CELLBC_BATCH_NUMBER+2));
	
	for(xk1=0; xk1< CELLBC_BATCH_NUMBER+2; xk1++){
		char tmp_fname[MAX_FILE_NAME_LENGTH+80];
		SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+80, "%s/temp-cellcounts-%06d-%03d.tmpbin", cct_context -> temp_file_dir, getpid(), xk1);
		input_fps[xk1] = fopen(tmp_fname,"rb");
		if(xk1 == CELLBC_BATCH_NUMBER+1)break;

		srInt_64 section1_items=0;
		for(sample_i = 0; sample_i < cct_context -> sample_sheet_table -> numOfElements; sample_i++){
			size_t frret = fread(&section1_items,1, 8, input_fps[xk1]);
			for(xk2 = 0; xk2 < section1_items; xk2++){
				srInt_64 cellbcP0_geneno0B=0, umis=0;
				frret += fread(&cellbcP0_geneno0B,1,8,input_fps[xk1]);
				frret += fread(&umis,1,8,input_fps[xk1]);

				int cellbc_no = cellbcP0_geneno0B>>32;
				int gene_no0B = (int)(cellbcP0_geneno0B&0xffffffffu);
				HashTable *gene_tab = HashTableGet(cellnoP1_to_genenoP1_to_UMIs[sample_i], NULL+cellbc_no+1);
				if(gene_tab==NULL){
					gene_tab = HashTableCreate(300);
					HashTablePut(cellnoP1_to_genenoP1_to_UMIs[sample_i], NULL+cellbc_no+1, gene_tab); 
				}
				HashTablePut(gene_tab, NULL+gene_no0B+1 , NULL+umis);
			}
		}
		last_rbin_buffer[xk1] = malloc( cct_context -> barcode_batched_max_genes *8 + cct_context -> barcode_batched_max_Rbin_len + 4 + MAX_UMI_LEN + 16 + 10000);
		int rlen = fread(last_rbin_buffer[xk1], 1, 16, input_fps[xk1]);
		if(rlen >0){
			int binlen = 0;
			srInt_64 genes = 0;
			memcpy(&genes, last_rbin_buffer[xk1]+8, 8);
			if(genes & (1LLU<<63))genes = genes & 0x7fffffff;
			else genes= 0;
			
			size_t frret = fread(last_rbin_buffer[xk1]+16, 1, 8*genes+ cct_context -> UMI_length + 4, input_fps[xk1]);
			memcpy(&binlen, last_rbin_buffer[xk1] +16 +8*genes+ cct_context -> UMI_length  , 4);
			frret += fread(last_rbin_buffer[xk1] + 16+ 8*genes+ cct_context -> UMI_length + 4, 1, binlen, input_fps[xk1]);

			srInt_64 sorting_key = *(int*)(last_rbin_buffer[xk1] + 16 +8*genes+cct_context -> UMI_length +4);
			sorting_key = sorting_key << 32;
			sorting_key |= *(int*)(last_rbin_buffer[xk1] + 16+ 8*genes+cct_context -> UMI_length +8);
			current_sorting_key[xk1] = sorting_key;
		}else{
			fclose(input_fps[xk1]);
			free(last_rbin_buffer[xk1]);
			current_sorting_key[xk1] = 0x7fffffffffffffffLLU;
		}
	}

	current_sorting_key[CELLBC_BATCH_NUMBER+1]=0;
	void * del_ptrs[2];
	del_ptrs[0] = current_sorting_key;
	del_ptrs[1] = cct_context ;
	pthread_create(&cct_context->thread_delete_files, NULL, delete_file_thread, del_ptrs);
	int current_worker = 0;
	while(1){
		int selected_fp_no = 0;
		srInt_64 selected_fp_key = current_sorting_key[0];
		for(xk1=1; xk1<CELLBC_BATCH_NUMBER+1; xk1++){
			if(current_sorting_key[xk1] < selected_fp_key){
				selected_fp_key = current_sorting_key[xk1] ;
				selected_fp_no = xk1;
			}
		}
		if(selected_fp_key == 0x7fffffffffffffffLLU) break;

		int sample_id = 0, binlen = 0;
		srInt_64 genes = 0;
		memcpy(&sample_id, last_rbin_buffer[selected_fp_no], 4);
		memcpy(&genes, last_rbin_buffer[selected_fp_no]+8, 8);
		if(genes & (1LLU<<63)) genes = genes & 0x7fffffff;
		else genes = 0;
		memcpy(&binlen,last_rbin_buffer[selected_fp_no]+16+8*genes+cct_context -> UMI_length,4);
		struct scRNA_merge_batches_worker_task * tofill = task_buffers+(current_filling_worker_per_sample[sample_id-1] * cct_context->sample_sheet_table -> numOfElements +sample_id-1);
		memcpy(tofill->inbin + tofill-> inbin_len, last_rbin_buffer[selected_fp_no]+16+8*genes+cct_context -> UMI_length, binlen + 4);
		tofill -> inbin_len += (binlen + 4);
		if(tofill -> inbin_number ==0) tofill -> inbin_number =1;
		//SUBREADprintf("ADDING BLOCKKK = %d  WKR = %d  IT THINK IT'S %d ; GENES=%d\n", tofill -> inbin_len, current_worker, tofill -> sample_id, genes);
		if(tofill-> inbin_len > CELLCOUNTS_BAMBLOCK_SIZE * CELLCOUNTS_BAMBLOCK_COMP_NUMBER){
			master_wait_for_job_done(&worker_mut, current_worker);
			struct scRNA_merge_batches_worker_current * my_finished_job = worker_current_jobs+current_worker;
			cellCounts_save_BAM_result(cct_context, my_finished_job);
			my_finished_job -> task = tofill;
			my_finished_job -> task -> block_number = (block_numbers_current[sample_id-1]++);
			master_notify_worker(&worker_mut, current_worker);

			current_filling_worker_per_sample[sample_id-1] ++;
			if(current_filling_worker_per_sample[sample_id-1] == compress_workers +1) current_filling_worker_per_sample[sample_id-1] = 0;
			tofill = task_buffers+(current_filling_worker_per_sample[sample_id-1] * cct_context->sample_sheet_table -> numOfElements +sample_id-1);
			tofill -> inbin_len = 0;
			tofill -> inbin_number = 0;
			tofill -> inbin_batch_start_offsets[0]=0;
			
			current_worker ++;
			if(current_worker == compress_workers) current_worker=0;
		}else if( tofill -> inbin_len - tofill -> inbin_batch_start_offsets[ tofill -> inbin_number - 1 ] > CELLCOUNTS_BAMBLOCK_SIZE ){
			tofill -> inbin_batch_start_offsets[tofill -> inbin_number]= tofill -> inbin_len;
			tofill -> inbin_number++;
			block_numbers_current[sample_id-1]++;
		}

		int rlen = fread(last_rbin_buffer[selected_fp_no], 1, 16, input_fps[selected_fp_no]);
		if(rlen >0){
			int binlen = 0;
			srInt_64 genes = 0;
			memcpy(&genes, last_rbin_buffer[selected_fp_no]+8, 8);
			if(genes & (1LLU<<63))genes = genes & 0x7fffffff;
			else genes= 0;
			size_t frret = fread(last_rbin_buffer[selected_fp_no]+16, 1, 8*genes+ cct_context -> UMI_length + 4, input_fps[selected_fp_no]);
			memcpy(&binlen, last_rbin_buffer[selected_fp_no] +16 +8*genes+ cct_context -> UMI_length  , 4);
			frret += fread(last_rbin_buffer[selected_fp_no] + 16+ 8*genes+ cct_context -> UMI_length + 4, 1, binlen, input_fps[selected_fp_no]);
			srInt_64 sorting_key = *(int*)(last_rbin_buffer[selected_fp_no] + 16+8*genes +cct_context -> UMI_length +4);
			sorting_key = sorting_key << 32;
			sorting_key |= *(int*)(last_rbin_buffer[selected_fp_no] + 16 +8*genes+cct_context -> UMI_length +8);
			current_sorting_key[selected_fp_no] = sorting_key;
		} else {
			fclose(input_fps[selected_fp_no]);
			free(last_rbin_buffer[selected_fp_no]);
			current_sorting_key[selected_fp_no] = 0x7fffffffffffffffLLU;
		}
	}

	for(xk1=0; xk1<cct_context -> sample_sheet_table -> numOfElements; xk1++){
		struct scRNA_merge_batches_worker_task * tofill = task_buffers+(current_filling_worker_per_sample[xk1] * cct_context-> sample_sheet_table -> numOfElements +xk1);
		if(tofill->inbin_len<1) continue;

		master_wait_for_job_done(&worker_mut, current_worker);
		struct scRNA_merge_batches_worker_current * my_finished_job = worker_current_jobs+current_worker;
		cellCounts_save_BAM_result(cct_context, my_finished_job);
		my_finished_job -> task = tofill;
		my_finished_job -> task -> block_number = (block_numbers_current[xk1]++);
		master_notify_worker(&worker_mut, current_worker);
		current_worker ++;
		if(current_worker == compress_workers) current_worker=0;
	}
	for(xk1=0; xk1<compress_workers; xk1++){
		struct scRNA_merge_batches_worker_current * my_finished_job = worker_current_jobs+current_worker;
		if(my_finished_job -> task)master_wait_for_job_done(&worker_mut, current_worker);
		cellCounts_save_BAM_result(cct_context, my_finished_job);

		current_worker ++;
		if(current_worker == compress_workers) current_worker=0;
	}

	for(xk1 = 0; xk1 < 1+compress_workers; xk1++) for(xk2 = 0; xk2 < cct_context->sample_sheet_table -> numOfElements;xk2++) {
		task_buffers[ xk1 * cct_context->sample_sheet_table -> numOfElements  + xk2 ].inbin_len = 0;
		task_buffers[ xk1 * cct_context->sample_sheet_table -> numOfElements  + xk2 ].inbin_number = 0;
	}
	current_worker = 0;
	FILE * notmapped_fp = input_fps[CELLBC_BATCH_NUMBER+1];
	while(1){
		int sample_id = 0, binlen = 0;
		int rlen = fread(&sample_id, 1, 4, notmapped_fp);
		if(rlen < 4) break;
		struct scRNA_merge_batches_worker_task * tofill = task_buffers+(current_filling_worker_per_sample[sample_id -1] * cct_context->sample_sheet_table -> numOfElements +sample_id-1);
		size_t frret = fread(&binlen, 1, 4, notmapped_fp);
		char old_bin[binlen+4];
		memcpy(old_bin, &binlen,4);
		frret += fread(old_bin+4, 1, binlen, notmapped_fp);
		int new_binlen = cellCounts_make_barcode_bam_bin(cct_context, old_bin, tofill -> inbin + tofill -> inbin_len, binlen, NULL, NULL, -1, NULL);
		tofill -> inbin_len += 4+ new_binlen; // block size: Total length of the alignment record, excluding this field. Then, the alignment record.

		if(tofill -> inbin_number ==0) tofill -> inbin_number =1;
		if(tofill-> inbin_len > CELLCOUNTS_BAMBLOCK_SIZE * CELLCOUNTS_BAMBLOCK_COMP_NUMBER){
			struct scRNA_merge_batches_worker_current * my_finished_job = worker_current_jobs+current_worker;
			if(my_finished_job -> task)master_wait_for_job_done(&worker_mut, current_worker);
			cellCounts_save_BAM_result(cct_context, my_finished_job);
			my_finished_job -> task = tofill;
			my_finished_job -> task -> block_number = (block_numbers_current[sample_id-1]++);
			master_notify_worker(&worker_mut, current_worker);

			current_filling_worker_per_sample[sample_id-1] ++;
			if(current_filling_worker_per_sample[sample_id-1] == compress_workers +1) current_filling_worker_per_sample[sample_id-1] = 0;
			tofill = task_buffers+(current_filling_worker_per_sample[sample_id-1] * cct_context-> sample_sheet_table -> numOfElements +sample_id-1);
			tofill -> inbin_len = 0;
			tofill -> inbin_number = 0;
			tofill -> inbin_batch_start_offsets[0]=0;
			
			current_worker ++;
			if(current_worker == compress_workers) current_worker=0;
		}else if( tofill -> inbin_len - tofill -> inbin_batch_start_offsets[ tofill -> inbin_number - 1 ] > CELLCOUNTS_BAMBLOCK_SIZE ){
			tofill -> inbin_batch_start_offsets[tofill -> inbin_number]= tofill -> inbin_len;
			block_numbers_current[sample_id-1]++;
			tofill -> inbin_number++;
		}
	}
	fclose(notmapped_fp);
	current_sorting_key[CELLBC_BATCH_NUMBER+1] = 0x7fffffffffffffffLLU;

	for(xk1=0; xk1<cct_context -> sample_sheet_table -> numOfElements; xk1++){
		struct scRNA_merge_batches_worker_task * tofill = task_buffers+(current_filling_worker_per_sample[xk1] * cct_context-> sample_sheet_table -> numOfElements +xk1);
		if(tofill->inbin_len<1) continue;

		master_wait_for_job_done(&worker_mut, current_worker);
		struct scRNA_merge_batches_worker_current * my_finished_job = worker_current_jobs+current_worker;
		cellCounts_save_BAM_result(cct_context, my_finished_job);
		my_finished_job -> task = tofill;
		my_finished_job -> task -> block_number = (block_numbers_current[xk1]++);
		master_notify_worker(&worker_mut, current_worker);
		current_worker ++;
		if(current_worker == compress_workers) current_worker=0;
	}

	for(xk1=0; xk1<compress_workers; xk1++){
		master_wait_for_job_done(&worker_mut, current_worker);
		struct scRNA_merge_batches_worker_current * my_finished_job = worker_current_jobs+current_worker;
		cellCounts_save_BAM_result(cct_context, my_finished_job);

		current_worker ++;
		if(current_worker == compress_workers) current_worker=0;
	}

	terminate_workers(&worker_mut);
	free(task_buffers);
	free(worker_current_jobs);

	for(xk1=0; xk1< compress_workers; xk1++) pthread_join(threads[xk1],NULL);

	worker_master_mutex_destroy(&worker_mut);

	if(!strstr(cct_context->index_prefix,".InternalTest.Rsr.step5@rand"))cellCounts_merged_to_tables_write(cct_context , cellnoP1_to_genenoP1_to_UMIs , cct_context -> all_features_array, cct_context -> all_features_array->numOfElements);

	for(xk1=0; xk1<cct_context -> sample_sheet_table -> numOfElements; xk1++)
		HashTableDestroy(cellnoP1_to_genenoP1_to_UMIs[xk1]);
	return 0;
}

int cellCounts_run_counting(cellcounts_global_t * cct_context){
	int ret= 0;
	ret = ret || cellCounts_do_cellbc_batches(cct_context);
	ret = ret || cellCounts_write_gene_list(cct_context);
	return ret;
}

void cellCounts_finalise_error_run(cellcounts_global_t * cct_context){
	void * argsv[2];
	argsv[0]=NULL;
	argsv[1]=cct_context;
	delete_file_thread(argsv);
}

#ifdef MAKE_STANDALONE
	#define cellCounts_main main
#endif
int cellCounts_main(int argc, char** argv){
	setlocale(LC_ALL,"");
	cellcounts_global_t * cct_context = calloc(sizeof(cellcounts_global_t),1);
	cct_context -> program_start_time = miltime();

	int ret = 0;
	ret = ret || cellCounts_args_context(cct_context, argc, argv);
	ret = ret || cellCounts_load_context(cct_context);
	ret = ret || cellCounts_run_mapping(cct_context);
	ret = ret || cellCounts_run_counting(cct_context);
	ret = ret || cellCounts_destroy_context(cct_context);
	if(cct_context -> has_error) cellCounts_finalise_error_run(cct_context);

	free(cct_context);
	if(ret) SUBREADprintf("cellCounts terminates with errors.\n");
	return ret;
}
