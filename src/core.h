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
  
  
//Core.c is the totally refactored body of the current subread, subjunc, subfusion and subindel programs.
//This program runs in this way:
// 1, it loads reads from FASTQ, FASTA, PLAIN or SAM files
// 2, it change the read and the quality strings according to FF/FR/RR/RF orders, such that every read is in the "F" manner.
//    Namely, the second read in a pair is reversed.
// 3, it builds voting table for each read.
// 4, it performs the first iteration.
// 5, it performs the second iteration.
// 6, it reports results regarding the reads. 
// 7, if the input files contain too many reads (e.g., greater than 14M in total), step 1 ~ 6 are repeated until all reads are finalised.
// 8, it reports global results.

//Core.c maintains the following context:
// 1, the global context.
// 2, the overall context for each module.
// 3, the local context for voting each read for each module.
// The structure of each context is named as
//     Context_ModuleName_Scope

//Thread context is maintained for each thread.

// There should be a function for each module to store/load its global context.

#ifndef _SUBREAD_CORE_H_
#define _SUBREAD_CORE_H_
#include <pthread.h> 
#include "subread.h"
#include "sambam-file.h"

//#define is_target(s) (memcmp((s),"chr901_242237_242737_0:0:0_0:0:0_3761", 22)==0)

#define is_target(s) 0

//#define _global_retrieve_voting_context(global_context, pair_number) (global_context->chunk_vote_records[pair_number])

//#define _global_retrieve_alignment_ptr(global_context, pair_number, is_second_read, best_read_id) ((global_context->input_reads.is_paired_end_reads?(global_context -> chunk_alignment_records + (2*pair_number+is_second_read)* global_context->config.multi_best_reads + best_read_id):(global_context -> chunk_alignment_records+pair_number*global_context->config.multi_best_reads + best_read_id)))

//#define _global_retrieve_subjunc_ptr(global_context, pair_number, is_second_read, best_read_id) ((global_context->input_reads.is_paired_end_reads?(global_context -> chunk_subjunc_records + (2*pair_number+is_second_read)* global_context->config.multi_best_reads + best_read_id):(global_context -> chunk_subjunc_records+pair_number*global_context->config.multi_best_reads + best_read_id)))

//#define _global_retrieve_big_margin_ptr(global_context,pair_number, is_second_read)  (is_second_read?(global_context->big_margin_record + (2*pair_number+is_second_read)* global_context->config.big_margin_record_size):(global_context->big_margin_record + pair_number * global_context->config.big_margin_record_size))

#define mark_gapped_read(res) (res)-> result_flags|= CORE_IS_GAPPED_READ;

#define CORE_MAX_CIGAR_STR_LEN 110
#define CORE_ADDITIONAL_INFO_LENGTH 400
#define READPAIRS_FOR_CALC_EXPT_TLEN 1000

typedef struct{
	subread_read_number_t fragments;
	unsigned int maximum_interval_length;
	unsigned int expected_items;
	HashTable * entry_table;
} bucketed_table_t;

typedef struct{
	int capacity;
	int items;
	unsigned int keyed_bucket;
	unsigned int maximum_interval_length;

	unsigned int * positions;
	void ** details;
} bucketed_table_bucket_t;


typedef struct{
	unsigned long long capacity;
	subread_read_number_t fragments;
	subread_read_number_t * fragment_numbers;
} fragment_list_t;




typedef struct{
	int is_paired_end_reads;
	int is_internal_error;
	gene_input_t first_read_file;
	gene_input_t second_read_file;
	unsigned long long first_read_file_size;
	subread_lock_t input_lock;
	double avg_read_length;
} read_input_t;

typedef struct{
	char read_name[MAX_READ_NAME_LEN];
	unsigned short flags;
	char chro_name[MAX_CHROMOSOME_NAME_LEN];
	unsigned int location;
	unsigned short map_quality;
	char cigar[CORE_MAX_CIGAR_STR_LEN];
	char other_chro_name[MAX_CHROMOSOME_NAME_LEN];
	unsigned int other_location;
	long long int tlen;
	int rlen;
	char read_text[MAX_READ_LENGTH];
	char qual_text[MAX_READ_LENGTH];
	char additional_columns[CORE_ADDITIONAL_INFO_LENGTH];
} output_read_buffer_t;

typedef struct{
	int fragment_number_in_chunk;
	int multi_mapping_locations;
	int this_mapping_location;
	output_read_buffer_t r1;
	output_read_buffer_t r2;
} output_fragment_buffer_t;	

typedef struct{
	// running_scheme
	int all_threads;
	int is_first_iteration_running;
	int is_second_iteration_running;
	int is_third_iteration_running;
	int use_memory_buffer;
	float memory_use_multiplex;
	char temp_file_prefix[MAX_FILE_NAME_LENGTH];
	subread_read_number_t reads_per_chunk;
	int fast_run;

	// input_scheme
	char first_read_file[MAX_FILE_NAME_LENGTH * MAX_SCRNA_FASTQ_FILES * 3];
	char second_read_file[MAX_FILE_NAME_LENGTH];
	char exon_annotation_file[MAX_FILE_NAME_LENGTH];
	char exon_annotation_file_screen_out[MAX_FILE_NAME_LENGTH];
	char exon_annotation_alias_file[MAX_FILE_NAME_LENGTH];
	int exon_annotation_file_type;
	char exon_annotation_gene_id_column[MAX_READ_NAME_LEN];
	char exon_annotation_feature_name_column[MAX_READ_NAME_LEN];


	int do_remove_neighbour_for_scRNA;
	int scRNA_input_mode;
	short read_trim_5;
	short read_trim_3;
	int is_first_read_reversed;
	int is_second_read_reversed;
	int space_type;
	int is_methylation_reads;
	int phred_score_format;
	int is_SAM_file_input;
	int is_gzip_fastq;


	// reporting scheme
	char read_group_id[MAX_FILE_NAME_LENGTH];
	char read_group_txt[MAX_FILE_NAME_LENGTH];
	char output_prefix[MAX_FILE_NAME_LENGTH];
	int report_sam_file;
	int report_no_unpaired_reads;
	int min_mapped_fraction;
	int max_mismatch_exonic_reads;
	int max_mismatch_junction_reads;
	int ignore_unmapped_reads;
	int report_unmapped_using_mate_pos;
	int report_multi_mapping_reads;
	int downscale_mapping_quality;
	int is_BAM_input;
	int is_BAM_output;
	int is_input_read_order_required;
	int sort_reads_by_coordinates;
	int convert_color_to_base;
	int SAM_extra_columns;
	int report_multiple_best_in_pairs;
	unsigned int multi_best_reads;
	unsigned int reported_multi_best_reads;

	// basic voting
	char index_prefix[MAX_FILE_NAME_LENGTH];
	int top_scores;
	int max_vote_combinations;
	int max_vote_simples;
	int max_vote_number_cutoff;
	int total_subreads;
	int minimum_subread_for_first_read;
	int minimum_subread_for_second_read;
	float minimum_exonic_subread_fraction;
	int minimum_pair_distance;
	int maximum_pair_distance;
	int restrected_read_order;
	int show_soft_cliping;
	int max_indel_length;
	int expected_pair_distance;
	int ambiguous_mapping_tolerance;
	int use_hamming_distance_break_ties;
	int use_quality_score_break_ties;
	int big_margin_record_size;
	int PE_predominant_weight;
	int no_TLEN_preference;

	// subjunc
	int entry_program_name;
	int experiment_type;
	int do_breakpoint_detection;
	int do_big_margin_filtering_for_junctions;
	int do_big_margin_filtering_for_reads;
	int limited_tree_scan;
	int use_hamming_distance_in_exon;
	unsigned int maximum_intron_length;
	int high_quality_base_threshold;
	int max_insertion_at_junctions;
	int check_donor_at_junctions;

	// subfusion
	int do_fusion_detection;
	int do_long_del_detection;
	int do_structural_variance_detection;
	int prefer_donor_receptor_junctions;
	int more_accurate_fusions;
	int maximum_translocation_length;
	int maximum_colocating_distance;

	// indel
	char do_superlong_indel_detection;
	char extending_search_indels;
	int k_mer_length;
	int realignment_minimum_variant_distance;
	int reassembly_start_read_number;
	int reassembly_key_length;
	int reassembly_subread_length;
	int reassembly_window_multiplex;
	int reassembly_tolerable_voting;
	int reassembly_window_alleles;
	int init_max_event_number;
	int use_dynamic_programming_indel;
	int use_bitmap_event_table;
	int maximise_sensitivity_indel;
	int flanking_subread_indel_mismatch;
	int DP_penalty_create_gap;
	int DP_penalty_extend_gap;
	int DP_match_score;
	int DP_mismatch_penalty;

} configuration_t;

#define CORE_IS_GT_AG_DONORS 1
#define CORE_NOTFOUND_DONORS 2
#define CORE_IS_STRAND_JUMPED 4
#define CORE_IS_NEGATIVE_STRAND 8
#define CORE_IS_FULLY_EXPLAINED 16
#define CORE_IS_BREAKEVEN 32
#define CORE_IS_GAPPED_READ 64
#define CORE_IS_PAIRED_END 128
#define CORE_TOO_MANY_MISMATCHES 256 

#define CORE_CIGAR_OPT_M 0
#define CORE_CIGAR_OPT_S 1
#define CORE_CIGAR_OPT_D 2
#define CORE_CIGAR_OPT_I 3
#define CORE_CIGAR_OPT_B 4
#define CORE_CIGAR_OPT_N 5
#define CORE_CIGAR_OPT_BB 6
#define CORE_CIGAR_OPT_NN 7

#define CORE_PROGRAM_SUBREAD 100
#define CORE_PROGRAM_SUBJUNC 200
#define CORE_PROGRAM_SUBINDEL 1000

#define CORE_EXPERIMENT_DNASEQ 1000
#define CORE_EXPERIMENT_RNASEQ 2000

#define PRINT_BOX_WRAPPED 4
#define PRINT_BOX_NOCOLOR_FOR_COLON 2
#define PRINT_BOX_CENTER 1

typedef struct{
	unsigned int event_small_side;	// the last base before the event
	unsigned int event_large_side;	// the first base after the event
					// for exon-exon junctions, these two numbers are the last and first base in the exons. (ZERO-BASED INDEX)
					// for indels, the smaller side is the base before in inserted/deleted bases; the larger side is the base after them.
					// for a fusion, the two sides are the base before the fusion point (the fusion point is not a base, but a vertical bar between two bases -- two sides.)
					//
					//
					//					The 'V' marked bases are the two sides.
					//                                      V              V
					//    Exon-exon junction:    100 .......s|gt........ag|s....... 300	chrX
					//					 ^            ^   These two bars are the junction point.
					//
					//
					//
					//    DNA-fusions:
					//
					//    small_side ===>  .|.  <=== large_side		normal arrangement (N/Y)
					//
					//    small_side ===>  .|
					//    large_side ===>  .|				strand jumpped (N/N)
					//
					//    		        |.  <=== small_side		strand jumpped (Y/Y)  
					//    		        |.  <=== large_side
					//
					//                      ^ This is fusion point. It is not a base, but a vertical bar between two sides. 
					//
					//
					//
					// note: the fusion point is larger than event_small_side if small_side_increasing_coordinate == N; it is smaller than event_small_side if small_side_increasing_coordinate == Y;
					//       the fusion point is smaller than event_large_side if large_side_increasing_coordinate == Y; it is larger than event_large_side if large_side_increasing_coordinate == N. 

	short indel_length;
	short junction_flanking_left;
	short junction_flanking_right;

	char indel_at_junction;
	char is_negative_strand;	// this only works to junction detection, according to 'GT/AG' or 'CT/AC' donors. This only applys to junctions.
	char is_strand_jumped;		// "strand jumped" means that the left and right sides are on different strands. This only applys to fusions.
	char is_donor_found_or_annotation;		// only for junctions: GT/AG is found at the location. 1: found, 0:not found: 64: from annotation (thus unknown)
						// Also, if "is_strand_jumped" is true, all coordinates (e.g., splicing points, cover_start, cover_end, etc) are on "reversed read" view.

	char small_side_increasing_coordinate;
	char large_side_increasing_coordinate;	// normal exon-exon junctions must have N and Y on small/large increasing coordinates.

	//char is_ambiguous;
	char connected_next_event_distance;	// the distance (negative or positive) to the next event in the table. For example, if the cigar string is 10M3I1M1I10M, event "3I" will have 1 here .
	char connected_previous_event_distance;	// the distance (negative or positive) to the next event in the table. For example, if the cigar string is 10M3I1M1I10M, event "1I" will have 1 here.

	//char inserted_bases[(1+MAX_INSERTION_LENGTH) / 4 + 1];
	char * inserted_bases;
	unsigned short supporting_reads;
	unsigned short anti_supporting_reads;
	unsigned short final_counted_reads;
	unsigned short final_reads_mismatches;
	unsigned char event_type;
	unsigned int global_event_id;
	float event_quality;

	unsigned long long critical_read_id;
	int critical_supporting_reads;
} chromosome_event_t;



typedef struct
{

	unsigned int selected_position;
	short result_flags;
	short read_length;
	// 4 bytes
	gene_vote_number_t selected_votes;
	gene_vote_number_t used_subreads_in_vote;
	unsigned char noninformative_subreads_in_vote;
	// this coverage is the range on reads, in point of view of "main piece" strand (i.e., "is_negative_strand")
	char indels_in_confident_coverage; 
	char is_fully_covered;

	gene_vote_number_t selected_indel_record [MAX_INDEL_SECTIONS*3 + 1];
	unsigned short confident_coverage_start; 
	unsigned short confident_coverage_end;

	short subread_quality; 

} mapping_result_t;

typedef struct{
	mapping_result_t * mapping_result;
	unsigned int first_base_position;
	char cigar_string[CORE_MAX_CIGAR_STR_LEN];
	chromosome_event_t * supporting_chromosome_events[MAX_EVENTS_IN_READ];
	short flanking_size_left[MAX_EVENTS_IN_READ];
	short flanking_size_right[MAX_EVENTS_IN_READ];
	char crirical_support[MAX_EVENTS_IN_READ];
	
	char first_base_is_jumpped;
	short final_mismatched_bases;
	short final_matched_bases;
	short best_second_diff_bases;
	short realign_flags;
	short final_quality;
	short chromosomal_length;
	short MAPQ_adjustment;
	int known_junction_supp;
	int final_penalty;
} realignment_result_t;

#define BUCKETED_TABLE_INIT_ITEMS 3
#define FRAGMENT_LIST_INIT_ITEMS 3


typedef struct
{
	short split_point;
	gene_vote_number_t minor_votes;
	char  double_indel_offset;
	char indel_at_junction;
	char small_side_increasing_coordinate;
	char large_side_increasing_coordinate;
	unsigned int minor_position;
	// this coverage is the range on reads, in point of view of "main piece" strand
	unsigned short minor_coverage_start; 
	unsigned short minor_coverage_end; 

} subjunc_result_t;

// THREE_TOP_UPDATE dictates the number of top-votes.

typedef struct {
	int is_vote_t_item;
	int item_index_i;
	int item_index_j;
	unsigned int mapping_position;
	int major_half_votes;
	unsigned short read_start_base;
}simple_mapping_t;

typedef struct{
	simple_mapping_t * r1_loc;
	simple_mapping_t * r2_loc;

	unsigned long long score_adj;
} vote_combination_t;

typedef struct{
	unsigned long long maxinum_read_number;
	int result_chunk_fd;
	int junction_chunk_fd;
	int big_margin_chunk_fd;

	mapping_result_t * result_chunk_addr;
	subjunc_result_t * junction_chunk_addr;
	unsigned short * big_margin_chunk_addr;

	subread_lock_t resize_lock;
} bigtable_t;


// this cached item is for an entire read (not read pair, not one best result)
typedef struct {
	int status;
	subread_read_number_t pair_number;
	int is_second_read;

	unsigned short big_margin_data[9*3];
	mapping_result_t * alignment_res;
	subjunc_result_t * subjunc_res;

} bigtable_cached_result_t;


typedef struct{
	void * junction_tmp_r1;
	void * junction_tmp_r2;

	void * alignment_tmp_r1;
	void * alignment_tmp_r2;

	void * vote_simple_1_buffer;
	void * vote_simple_2_buffer;

	void * comb_buffer;
} topK_buffer_t;

typedef struct{
	unsigned long long all_correct_PE_reads;
	int thread_id;
	pthread_t thread;

	// module_thread_context is different to module_context, though some contents in module_thread_context could be also in module context (e.g., junction tables)
	// modules functions may deside to use objects in which context.
	void * module_thread_contexts[5];
	gene_value_index_t * current_value_index;

	subread_lock_t output_lock;
	unsigned int all_mapped_reads;
	unsigned int not_properly_pairs_wrong_arrangement;
	unsigned int not_properly_pairs_different_chro;
	unsigned int not_properly_different_strands;
	unsigned int not_properly_pairs_TLEN_wrong;
	unsigned int all_unmapped_reads;
	unsigned int not_properly_pairs_only_one_end_mapped;
	unsigned int all_multimapping_reads;
	unsigned int all_uniquely_mapped_reads;

	topK_buffer_t topKbuff;
} thread_context_t;


typedef struct{
	// basic running configuration
	configuration_t config;

	// for the index
	gehash_t * current_index;
	gene_value_index_t * current_value_index;
	gene_value_index_t all_value_indexes[100];
	int index_block_number;
	int current_index_block_number;
	int will_remove_input_file;
	int is_phred_warning;
	int is_final_voting_run;

	// global locks
	subread_lock_t thread_initial_lock;

	// global FILE*
	SamBam_Writer * output_bam_writer;
	FILE * output_sam_fp;
	FILE * long_insertion_FASTA_fp;
	char * output_sam_inner_buffer;
	int output_sam_is_full;

	// running contexts
	void * module_contexts[5]; 
	thread_context_t * all_thread_contexts;
	subread_read_number_t last_written_fragment_number;
	int need_merge_buffer_now;
	read_input_t input_reads;
	bigtable_t bigtable;
	int rebuilt_command_line_size;
	char * rebuilt_command_line;

	subread_lock_t bigtable_lock;
	subread_lock_t output_lock;
	int bigtable_cache_size;
	FILE * bigtable_cache_file_fp;
	long long bigtable_cache_file_loaded_fragments_begin;
	long long bigtable_cache_file_fragments;
	bigtable_cached_result_t * bigtable_cache;
	char *bigtable_cache_malloc_ptr, *bigtable_cache_malloc_ptr_junc;
	unsigned int bigtable_chunked_fragments;
	int bigtable_dirty_data;
	
	gene_offset_t chromosome_table;
	double start_time;
	double align_start_time;
	double timecost_load_index;
	double timecost_voting;
	double timecost_before_realign;
	double timecost_for_realign;

	unsigned long long expected_TLEN_sum;
	unsigned int expected_TLEN_read_numbers;
	unsigned long long all_processed_reads;
	unsigned long long all_correct_PE_reads;
	unsigned int all_junctions;
	unsigned int all_fusions;
	unsigned int all_indels;

	unsigned int all_mapped_reads;
	unsigned int not_properly_pairs_wrong_arrangement;
	unsigned int not_properly_pairs_different_chro;
	unsigned int not_properly_different_strands;
	unsigned int not_properly_pairs_TLEN_wrong;
	unsigned int all_unmapped_reads;
	unsigned int not_properly_pairs_only_one_end_mapped;
	unsigned int all_multimapping_reads;
	unsigned int all_uniquely_mapped_reads;

	unsigned long long current_circle_start_abs_offset_file1;
	gene_inputfile_position_t current_circle_start_position_file1;
	gene_inputfile_position_t current_circle_start_position_file2;
	gene_inputfile_position_t current_circle_end_position_file1;
	gene_inputfile_position_t current_circle_end_position_file2;
	subread_read_number_t processed_reads_in_chunk;
	subread_read_number_t running_processed_reads_in_chunk;

	// sunfusion structural variance
	bucketed_table_t funky_table_BC;
	bucketed_table_t funky_table_DE;
	fragment_list_t funky_list_A;
	fragment_list_t funky_list_DE;

	bucketed_table_t breakpoint_table_P;
	bucketed_table_t breakpoint_table_QR;
	bucketed_table_t breakpoint_table_YZ;

	bucketed_table_t translocation_result_table;
	bucketed_table_t inversion_result_table;
	topK_buffer_t topKbuff;

	// per chunk parameters
	subread_read_number_t read_block_start;
	char * exonic_region_bitmap;
	HashTable * sam_chro_to_anno_chr_alias;
	HashTable * annotation_chro_table;
} global_context_t;


#define MODULE_INDEL_ID 0
#define MODULE_JUNCTION_ID 1

#define STEP_VOTING 10
#define STEP_ITERATION_ONE 20
#define STEP_ITERATION_TWO 30
#define STEP_ITERATION_THREE 40
#define STEP_WRITE_CHUNK_RESULTS 90

#define MEMORY_OPTIMISATION_LAPTOP 10
#define MEMORY_OPTIMISATION_DESKTOP 20
#define MEMORY_OPTIMISATION_SUPERCOMPUTER 100

// This function initialise the module contexts that are needed according to the running configuration.
int init_modules(global_context_t * context);

// This function undo any memory allocation and close opened files in init_modules.
int destroy_modules(global_context_t * context);

// This function shows welcome messages and shows the current configuration to stdout.
// This function also checks the legality of the parameters in context->config
int print_configuration(global_context_t * context);

// This function opens all related files except the index. It however counts the piece number of the index.
int load_global_context(global_context_t * context);

// This function deallocates memory and closes files that were allowcated and opened in load_global_context.
int destroy_global_context(global_context_t * context);


// This function processes the read chunks (24M reads each), including three steps: 1, voting; 2, iteration one; 3, iteration two.
int read_chunk_circles(global_context_t * context);


// write the final results in each module. Note that the SAM file is writen in read_chunk_circles, after each read chunk is processed.
int write_final_results(global_context_t * context);

// this function is called by an interface to start the main process with a customized opt parser.
int core_main(int argc , char ** argv, int (parse_opts (int , char **, global_context_t * )));

// decompress a cigar string to cigar
// cigar_len is the maximum length of decompressed cigar
// bincigar is the maximum length of compressed cigar
// it returns the length of decompressed cigar, or -1 if buffer is too short.
int bincigar2cigar(char * cigar, int cigar_len, char * bincigar, int bincigar_max_len, int read_len);

// compress a cigar string to bincigar
// bincigar_len is the maximum length of compressed cigar
// it returns the length of compressed cigar, or -1 if buffer is too short.
int cigar2bincigar(char *cigar, char *bincigar, int bincigar_len);

// print the logo of our program
void print_subread_logo();

// print a line in the box
void print_in_box(int line_width, int is_boundary, int options, char * pattern,...);

// find the value index covering this read
// it returns NULL if no index is found.
gene_value_index_t * find_current_value_index(global_context_t * global_context, unsigned int pos, int len);

// generate the time string 
void char_strftime(char * tbuf);

int term_strncpy(char * dst, char * src, int max_dst_memory);

int is_result_in_PE(mapping_result_t * aln);

void core_version_number(char * program);

mapping_result_t * _global_retrieve_alignment_ptr(global_context_t * global_context, subread_read_number_t pair_number, int is_second_read, int best_read_id);
subjunc_result_t * _global_retrieve_subjunc_ptr(global_context_t * global_context, subread_read_number_t pair_number, int is_second_read, int best_read_id);
unsigned short * _global_retrieve_big_margin_ptr(global_context_t * global_context, subread_read_number_t pair_number, subread_read_number_t is_second_read);

// This assumes the first part of Cigar has differet strandness to the main part of the cigar.
// Pos is the LAST WANTED BASE location before the first strand jump (split by 'b' or 'n').
// The first base in the read actually has a larger coordinate than Pos. 
unsigned int reverse_cigar(unsigned int pos, char * cigar, char * new_cigar);

int chimeric_cigar_parts(global_context_t * global_context , unsigned int sel_pos, int is_first_section_negative_strand, int is_first_section_reversed, char * in_cigar, unsigned int * out_poses, char ** out_cigars, char * out_strands, int read_len, short * out_read_lens, char * read_name);

void warning_file_limit();
void quick_sort(void * arr,int arr_size, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r));
void basic_sort(void * arr,int arr_size, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r));

// L_Minus_R should return -1, 0 or 1 when L<R, L==R or L>R.
// The result is from Small to Large.
void merge_sort(void * arr, int arr_size, int L_Minus_R (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge_SmallFirst(void * arr, int start, int items, int items2));

void absoffset_to_posstr(global_context_t * global_context, unsigned int pos, char * res);

void test_PE_and_same_chro(global_context_t * global_context , unsigned int pos1, unsigned int pos2, int * is_PE_distance, int * is_same_chromosome , int rlen1, int rlen2);

int FIXLENstrcmp(char * fixed_len, char * rname);

int is_valid_digit(char * optarg, char * optname);
int is_valid_digit_range(char * optarg, char * optname, int min, int max_inc);
int is_valid_float(char * optarg, char * optname);
int exec_cmd(char * cmd, char * outstr, int out_limit);
int is_pos_in_annotated_exon_regions(global_context_t * global_context, unsigned int pos);
char * get_sam_chro_name_from_alias(HashTable * tab, char * anno_chro);
void subread_rebuild_cmd(int argc, char ** argv, global_context_t * global_context);
#endif
