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
  
  
#ifndef SUBREAD_CORE_JUNCTION_H_
#define SUBREAD_CORE_JUNCTION_H_
#include "subread.h"
#include "hashtable.h"
#include "core.h"

//#warning "======== REMOVE *2000 ============="
#define REALIGN_TOTAL_TRIES (50)

#define FUNKY_FRAGMENT_A	1	// same strand and gapped (0<gap<tra_len)
#define FUNKY_FRAGMENT_BC	2	// very far far away (>=tra_len) or chimeric.
#define FUNKY_FRAGMENT_DE	4	// tlen < tra_len and strand jumpped
#define NOT_FUNKY		0	// normal fragment 
#define FUNKY_COLOCATION_TOLERANCE 500
#define BREAK_POINT_MAXIMUM_TOLERANCE 80 
#define S12_LIST_CAPACITY 100


// as the python sub-string rule: start is the first wanted base and end is the first unwanted base.
typedef struct{
	short read_pos_start;
	short read_pos_end;
	unsigned int abs_offset_for_start;
	// jumped from the "main piece" view
	char is_strand_jumped;
	char is_connected_to_large_side;

	chromosome_event_t * event_after_section; 
} perfect_section_in_read_t;

typedef struct{
	// result context

	//unsigned char  back_search_confirmed_sections;
	//unsigned char front_search_confirmed_sections;
	// NOTE THAT EVERYTHING IN back_search_junctions IS BACKWARD.
	// 1, ORDER OF EXONS ARE BACKWARD
	// 2, "ABS_OFFSET_FOR_START" ARE ACTUALLY AT END OF SECTIONS
	//perfect_section_in_read_t  back_search_junctions[MAX_EVENTS_IN_READ];
	//perfect_section_in_read_t front_search_junctions[MAX_EVENTS_IN_READ];

	// middle result context
	unsigned char tmp_search_sections;
	perfect_section_in_read_t  tmp_search_junctions [MAX_EVENTS_IN_READ];
	char current_is_strand_jumped;

	perfect_section_in_read_t result_back_junctions [MAX_ALIGNMENT_PER_ANCHOR][MAX_EVENTS_IN_READ];
	perfect_section_in_read_t result_front_junctions [MAX_ALIGNMENT_PER_ANCHOR][MAX_EVENTS_IN_READ];
	int result_back_junction_numbers[MAX_ALIGNMENT_PER_ANCHOR];
	int result_front_junction_numbers[MAX_ALIGNMENT_PER_ANCHOR];
	int all_back_alignments;
	int all_front_alignments;
	int known_junctions;
	unsigned int total_tries;

//	unsigned int tmp_jump_length;
//	unsigned int best_jump_length;
	
	// for the BEST record
	// they are not restored
	int best_matching_bases;
	int best_second_match_diff;
	int second_best_matching_bases;
	int best_indel_penalty;
	int tmp_total_matched_bases;
	int tmp_indel_penalty;
	int is_currently_tie;
	int best_is_complex;
	int best_support_as_simple;
	int best_min_unsupport_as_simple;
	int best_min_support_as_complex;
	int best_is_pure_donor_found_explain;

	// for the "current" stack (they are restored after poping from the stack)
	int tmp_support_as_simple;
	int tmp_min_unsupport;
	int tmp_min_support_as_complex;
	int tmp_is_pure_donor_found_explain;

	// input context
	int full_read_len;
	int is_fully_covered;
	char * full_read_text;
	char * full_qual_text;
	char * read_name;
	int is_confirmed_section_negative_strand;
	subread_read_number_t pair_number;
	int is_second_read;
	int best_read_id;
}explain_context_t;

typedef struct{
} junction_context_t;

typedef struct
{
	int read_len_1;
	int read_len_2;
	char * read_text_1;
	char * read_text_2;
	char is_negative_strand;
} new_junction_context_t;

void new_explain_try_replace(global_context_t* global_context, thread_context_t * thread_context, explain_context_t * explain_context, int remainder_len, int search_to_back);

int init_junction_tables(global_context_t * context);
int destroy_junction_tables(global_context_t * context);
int process_voting_junction(global_context_t * global_context, thread_context_t * thread_context, subread_read_number_t pair_number, gene_vote_t * vote_1, gene_vote_t * vote_2, char * read_name_1, char * read_name_2, char * read_text_1, char * read_text_2,  int read_len_1, int read_len_2, int is_negative_strand,  gene_vote_number_t v1_all_subreads, gene_vote_number_t v2_all_subreads);
int init_junction_thread_contexts(global_context_t * global_context, thread_context_t * thread_context, int task);
int finalise_junction_thread(global_context_t * global_context, thread_context_t * thread_context, int task);
unsigned int explain_read(global_context_t * global_context, thread_context_t * thread_context, realignment_result_t * realigns, subread_read_number_t pair_number,int read_len, char * read_name , char *read_text, char *qual, int is_second_read, int best_read_id, int is_negative_strand);
int write_junction_final_results(global_context_t * global_context);

// back_search_read_tail IS THE EXACT VERY SURE POSITION IN THE READ
// back_search_tail_position IS THE EXACT VERY SURE POSITION ON CHROMOSOME CORRESPONDING TO back_search_read_tail
int do_explain_back_search(global_context_t * global_context, thread_context_t * thread_context, explain_context_t * explain_context, char * read_text, char * read_qual, int back_search_read_tail, unsigned int back_search_tail_position);
int do_explain_front_search(global_context_t * global_context, thread_context_t * thread_context, explain_context_t * explain_context, char * read_text, char * read_qual, int front_search_read_head, unsigned int front_search_head_position);

unsigned int finalise_explain_CIGAR(global_context_t * global_context, thread_context_t * thread_context, explain_context_t * explain_context, realignment_result_t * realigns);

void search_events_to_back(global_context_t * global_context, thread_context_t * thread_context, explain_context_t * explain_context, char * read_text , char * qual_text, unsigned int read_tail_abs_offset, short read_tail_pos, short sofar_matched, int suggested_movements, int do_not_jump);

void search_events_to_back(global_context_t * global_context, thread_context_t * thread_context, explain_context_t * explain_context, char * read_text , char * qual_text, unsigned int read_tail_abs_offset, short read_tail_pos, short sofar_matched, int suggested_movements, int do_not_jump);

void find_new_junctions(global_context_t * global_context, thread_context_t * thread_context, subread_read_number_t pair_number, char * read_name, char * read_text, char * qual_text, int read_len, int is_second_read, int best_read_id);

int donor_score(global_context_t * global_context, thread_context_t * thread_context, unsigned int left_virtualHead_abs_offset, unsigned int right_virtualHead_abs_offset, int left_indels, int right_indels, int normally_arranged, int guess_start, int guess_end,  char * read_text, int read_len, int * final_split_point, int * is_GT_AG_strand, int * is_donor_found, int * inserted_bases, int * small_side_inc_coor, int * large_side_inc_coor, char *read_name);

int donor_jumped_score(global_context_t * global_context, thread_context_t * thread_context, unsigned int left_virtualHead_abs_offset, unsigned int right_virtualHead_abs_offset, int guess_start, int guess_end,  char * read_text, int read_len, int is_left_half_negative, int is_right_half_negative, int is_left_part_on_left_as_reversed, int * final_split_point, int * is_GT_AG_strand, int * is_donor_found, int * small_side_inc_coor, int * large_side_inc_coor);

int write_fusion_final_results(global_context_t * global_context);


int is_ambiguous_voting(global_context_t * global_context, subread_read_number_t pair_number, int is_second_read, int max_vote, int max_start,int max_end, int read_len, int is_negative);
void core_search_short_exons(global_context_t * global_context, thread_context_t * thread_context, char * read_text, char * qualityb0, int rl, unsigned int P1_Pos, unsigned int P2_Pos, short read_coverage_start, short read_coverage_end);


void core_fragile_junction_voting(global_context_t * global_context, thread_context_t * thread_context, char * rname, char * read, char * qual, unsigned int full_rl, int negative_strand, int color_space, unsigned int low_border, unsigned int high_border, gene_vote_t *vote_p1);

int is_funky_fragment(global_context_t * global_context, char * rname1, char * chr1, unsigned int pos1, int rlen1, int is_1_negative, char * cigar1, char * seq1, char * rname2, char * chr2, unsigned int pos2, int rlen2, int is_2_negative, char * cigar2, char * seq2, int tlen_removed_intron);

void finalise_structural_variances(global_context_t * global_context);

void debug_show_event(global_context_t* global_context, chromosome_event_t * event);
void get_event_two_coordinates(global_context_t * global_context, unsigned int event_no, char ** small_chro, int * small_pos, unsigned int * small_abs, char ** large_chro,  int * large_pos, unsigned int * large_abs);
int get_offset_maximum_chro_pos(global_context_t * global_context, thread_context_t * thread_context, unsigned int linear);
#endif
