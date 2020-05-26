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
  
  
#ifndef SUBREAD_CORE_INDEL_H_
#define SUBREAD_CORE_INDEL_H_

#include "subread.h"
#include "hashtable.h"
#include "core.h"

// chromosome events can be indels, junctions or fusions.
// if it is an insertion event, event_large_site = event_small_site+1.

//#define MAX_EVENT_ENTRIES_PER_SITE 5
//#define MAX_EVENT_ENTRIES_PER_SITE 12 
//
#define EVENT_ENTRIES_INIT_SIZE (9) 
#define MAX_EVENT_ENTRIES_PER_SITE 9
#define CHRO_EVENT_TYPE_REMOVED 0
#define CHRO_EVENT_TYPE_INDEL 8
#define CHRO_EVENT_TYPE_LONG_INDEL 16 
#define CHRO_EVENT_TYPE_POTENTIAL_INDEL 32 
#define CHRO_EVENT_TYPE_JUNCTION 64 
#define CHRO_EVENT_TYPE_FUSION 128
#define CHRO_EVENT_TYPE_SNP 256 

#define EVENT_SEARCH_BY_SMALL_SIDE 10
#define EVENT_SEARCH_BY_LARGE_SIDE 20
#define EVENT_SEARCH_BY_BOTH_SIDES 30


#define REASSEMBLY_WINDOW_LENGTH 350

//#define is_target_window_X(x) ((x + 1) * REASSEMBLY_WINDOW_LENGTH / 2 >= (10734463 % BASE_BLOCK_LENGTH) && (x- 1) * REASSEMBLY_WINDOW_LENGTH /2-1 <= (10734463%BASE_BLOCK_LENGTH) )
#define is_target_window_X(x) 0
//#define MAXIMUM_EVENT_NUMBER 300000


typedef struct{
	int is_precisely_called;
	unsigned int source_left_side;	// the base BEFORE the translocated sequence.
	unsigned int target_left_side;  // tge base BEFORE the inserted translocated sequence.
	unsigned int length;

	unsigned int event_P_number;
	unsigned int event_Q_number;
	unsigned int event_R_number;

	int is_inv;
	unsigned int all_sup_P;
	unsigned int max_sup_QR;
} translocation_result_t;

typedef struct{
	int is_precisely_called;

	unsigned int event_Y_rough_small_abs;
	unsigned int event_Z_rough_large_abs;

	unsigned int small_side;	// the base BEFORE the reversed sequence
	unsigned int length;

	unsigned int event_Y_number;	// event_no in the event space.
	unsigned int event_Z_number;

	unsigned int all_sup_D;
	unsigned int max_sup_E;
} inversion_result_t;

struct reassmebly_window_allele
{
	char rebuilt_window[8000];
	float allele_quality;
	int rebuilt_size;
};

typedef struct{
	gehash_t * voting_indexes;
	char * chro_name;
	unsigned long long int * start_keys;
	short * start_offsets;

	unsigned int * read_no_counter;
	unsigned int block_start_linear_pos;
	HashTable * read_sequence_table;
	HashTable * read_position_table;
	HashTable * read_quality_table;
	gene_vote_t * vote_list;
	gene_vote_t * vote_list_rectify;
	short * read_rectify_space;

	char rebuilt_window[2500];
	int rebuilt_window_size;


	struct reassmebly_window_allele * final_alleles;

	unsigned int used_read_ids[2000];
	int used_read_number;


	int search_cost;
	int total_matched_bases;
	int max_matched_bases;
	unsigned int window_start_pos;
} reassembly_by_voting_block_context_t;



typedef struct{
	HashTable ** de_bruijn_graphs;
	char * chro_name;
	unsigned long long int * start_keys;
	short * start_offsets;

	unsigned int block_start_linear_pos;
} reassembly_block_context_t;

#define EVENT_BODY_LOCK_BUCKETS 14929


typedef struct{
	HashTable * event_entry_table;
	unsigned int total_events;
	unsigned int current_max_event_number;
	chromosome_event_t * event_space_dynamic;
	HashTable * local_reassembly_pileup_files;
	subread_lock_t event_body_locks[EVENT_BODY_LOCK_BUCKETS];

	short ** dynamic_align_table;
	char ** dynamic_align_table_mask;
} indel_context_t;

typedef struct{
	HashTable * event_entry_table;
	unsigned int total_events;
	unsigned int current_max_event_number;
	chromosome_event_t * event_space_dynamic;
	unsigned short * final_counted_reads_array;
	unsigned short * final_reads_mismatches_array;

	short ** dynamic_align_table;
	char ** dynamic_align_table_mask;
} indel_thread_context_t;

int init_indel_tables(global_context_t * context);
int destroy_indel_module(global_context_t * context);
int init_indel_thread_contexts(global_context_t * global_context, thread_context_t * thread_context, int task);
int sort_global_event_table(global_context_t * global_context);
int load_known_junctions(global_context_t * global_context);
int finalise_indel_and_junction_thread(global_context_t * global_context, thread_context_t * thread_contexts, int task);
int find_new_indels(global_context_t * global_context, thread_context_t * thread_context, int pair_number, char * read_name, char * read_text, char * qual_text, int read_len, int is_second_read, int best_read_id);
int write_indel_final_results(global_context_t * context);
int search_event(global_context_t * global_context,HashTable * event_table, chromosome_event_t * event_space, unsigned int pos, int search_type, unsigned char event_type, chromosome_event_t ** return_buffer);

void set_alignment_result(global_context_t * global_context, int pair_number, int is_second_read, int best_read_id, unsigned int position, int votes, gene_vote_number_t * indel_record, short best_cover_start, short best_cover_end, int is_negative_strand, int is_PE, unsigned int minor_position, unsigned int minor_votes, unsigned int minor_coverage_start, unsigned int minor_coverage_end, unsigned int split_point, int inserted_bases, int is_strand_jumped, int is_GT_AG_donors, int used_subreads_in_vote, int noninformative_subreads_in_vote, int major_indel_offset, int minor_indel_offset, int main_hamming, int minor_hamming, int main_quality, int minor_quality);

void put_new_event(HashTable * event_table, chromosome_event_t * new_event , int event_no);
void remove_neighbour(global_context_t * global_context);
int build_local_reassembly(global_context_t *global_context , thread_context_t *thread_context , int pair_number, char * read_name_1 , char * read_text_1 ,char * qual_text_1 , int read_len_1, int read_len_2, int is_second_read, int best_read_id, int is_paired_unmapped, mapping_result_t * current_res, mapping_result_t * mate_res);
int finalise_long_insertions(global_context_t * global_context);

// This function sets the global context with default values.
void init_global_context(global_context_t * context);

int write_local_reassembly(global_context_t *global_context, HashTable *pileup_fp_table, unsigned int anchor_pos, char * read_name , char * read_text ,char * qual_text , int read_len, int is_anchor_certain);

int finalise_long_insertions_by_hashtable(global_context_t * global_context);

void destroy_pileup_table(HashTable* local_reassembly_pileup_files);

chromosome_event_t * reallocate_event_space(global_context_t* global_context,thread_context_t* thread_context,int event_no);

int there_are_events_in_range(char * bitmap, unsigned int pos, int sec_len);

int anti_supporting_read_scan(global_context_t * global_context);

int core_dynamic_align(global_context_t * global_context, thread_context_t * thread_context, char * read, int read_len, unsigned int begin_position, char * movement_buffer, int expected_offset, char * read_name);

void init_core_temp_path(global_context_t * context);

chromosome_event_t * local_add_indel_event(global_context_t * global_context, thread_context_t * thread_context, HashTable * event_table, char * read_text, unsigned int left_edge, int indels, int score_supporting_read_added, int is_ambiguous, int mismatched_bases,int * old_event_id);

void print_indel_table(global_context_t * global_context);
int sort_junction_entry_table(global_context_t * global_context);
void mark_event_bitmap(unsigned char * bitmap, unsigned int pos);
int check_event_bitmap(unsigned char * bitmap, unsigned int pos);
#endif
