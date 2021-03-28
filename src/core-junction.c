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
  
  
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "subread.h"
#include "sublog.h"
#include "gene-value-index.h"
#include "gene-algorithms.h"
#include "input-files.h"
#include "core.h"
#include "core-indel.h"
#include "core-junction.h"
#include "core-bigtable.h"

#define TTTSNAME "V0112_0155:7:1308:1308:136442"

#define CLUSTER_ALIGNMENT_DONOR_R1_MAPPED 2
#define CLUSTER_ALIGNMENT_DONOR_R2_MAPPED 4
#define CLUSTER_ALIGNMENT_DONOR_NEGATIVE_STRAND 1


unsigned int abs32uint(unsigned int x){
	if(x > 0x7fffffff) x = (0xffffffff - x) + 1;
	return x;
}

typedef struct{
	unsigned int piece_main_abs_offset;
	unsigned int piece_minor_abs_offset;
	int piece_main_masks;
	short piece_main_coverage_start;
	short piece_main_coverage_end;

	short piece_main_hamming_match;
	short piece_main_read_quality;
	short piece_minor_hamming_match;
	short piece_minor_read_quality;
	int piece_minor_score;
	short intron_length;

	gene_vote_number_t *piece_main_indel_record;
	unsigned short piece_main_indels;
	unsigned short piece_minor_indel_offset;
	gene_vote_number_t piece_main_votes;
	gene_vote_number_t piece_minor_votes;

	short piece_minor_coverage_start;
	short piece_minor_coverage_end;
	short split_point;
	char inserted_bases;
	char is_GT_AG_donors;
	char is_donor_found_or_annotation;
	char is_strand_jumped;
	char is_break_even;

	//unsigned long long int Score_H;
	//unsigned long long int Score_L;
} select_junction_record_t;


void debug_show_event(global_context_t* global_context, chromosome_event_t * event){
	char outpos1[100], outpos2[100];
	absoffset_to_posstr(global_context, event -> event_small_side, outpos1);
	absoffset_to_posstr(global_context, event -> event_large_side, outpos2);
	SUBREADprintf("Event between %s and %s\n", outpos1, outpos2);
}

int get_offset_maximum_chro_pos(global_context_t * global_context, thread_context_t * thread_context, unsigned int linear){
	gene_offset_t * chros =& global_context -> chromosome_table;
	int n = 0;
	int total_offsets = chros -> total_offsets;

	int LL = 0, RR = total_offsets-1;

	while(1){
		if(LL >= RR-1) break;
		int MM = (LL+RR)/2;
		if( linear > chros->read_offsets[MM]) LL = MM;
		else if(linear < chros->read_offsets[MM]) RR = MM;
		else break;
	}

	n = max(0, LL - 2);

	for (; n < chros -> total_offsets; n++) {
		if (chros->read_offsets[n] > linear) {
			int ret;
			unsigned int last_linear = 0;
			if(n==0)
				ret = chros->read_offsets[0] - chros -> padding  *2 +16;
			else{
				ret = ( chros->read_offsets[n] - chros->read_offsets[n-1]  ) - chros -> padding  *2 +16;
				last_linear =  chros->read_offsets[n-1];
			}
			linear -= last_linear;
			if(linear < chros -> padding || linear >= chros -> padding + ret) return -1;
			return ret;
		}
	}
	return -2;
}


// read_head_abs_pos is the offset of the FIRST WANTED base.
void search_events_to_front(global_context_t * global_context, thread_context_t * thread_context, explain_context_t * explain_context, char * read_text , char * qual_text, unsigned int read_head_abs_offset, short remainder_len, short sofar_matched, int suggested_movement, int do_not_jump)
{
	short tested_read_pos;
	//				#warning "SUBREAD_151 REMOVE THIS ASSERTION! "
	//				if(remainder_len >= 102)SUBREADprintf("FATAL:%d\n", remainder_len );
	//				assert(remainder_len < 102);

	HashTable * event_table = NULL;
	chromosome_event_t * event_space = NULL;

	if(thread_context) {
		event_table = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	} else {
		event_table = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}

	gene_value_index_t * value_index = thread_context?thread_context->current_value_index:global_context->current_value_index ;

	if((global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection)|| there_are_events_in_range(event_table->appendix1, read_head_abs_offset, remainder_len)) {
	
		int event_search_method;
		if((global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection))
			event_search_method = EVENT_SEARCH_BY_BOTH_SIDES;
		else
			event_search_method = EVENT_SEARCH_BY_SMALL_SIDE;
	
		// tested_read_pos is the index of the first base unwanted!
		
	
		int move_start = do_not_jump?0:global_context -> config.realignment_minimum_variant_distance;
		if(suggested_movement) move_start = suggested_movement-1;
		int is_junction_scanned = 0;
	
		if(MAX_EVENTS_IN_READ - 1 > explain_context -> tmp_search_sections)
		for(tested_read_pos = move_start ; tested_read_pos <= remainder_len; tested_read_pos++)
		{
			int xk1, matched_bases_to_site;
			chromosome_event_t *site_events[MAX_EVENT_ENTRIES_PER_SITE+1];

			int jump_penalty = 0;

			unsigned potential_event_pos;
			if(explain_context -> current_is_strand_jumped)
				potential_event_pos = read_head_abs_offset - tested_read_pos +1;
			else
				potential_event_pos = read_head_abs_offset + tested_read_pos -1;
			if(!check_event_bitmap(  event_table->appendix1, potential_event_pos )) continue;

			int search_types =  CHRO_EVENT_TYPE_INDEL | CHRO_EVENT_TYPE_JUNCTION | CHRO_EVENT_TYPE_FUSION;
			int site_events_no = search_event(global_context, event_table , event_space , potential_event_pos, event_search_method , search_types , site_events);


			if(0 && FIXLENstrcmp("R010442852", explain_context -> read_name) == 0)
			{
				SUBREADprintf("FOUND THE EVENT FRONT:%d at %u\n", site_events_no, potential_event_pos);
				if(site_events_no)
					SUBREADprintf("EVENT0_type = %d\n", site_events[0]->event_type);
			}

			if(!site_events_no)continue;

			unsigned int tested_chro_begin;
			if(explain_context -> current_is_strand_jumped)
				tested_chro_begin = read_head_abs_offset - tested_read_pos + 1;
			else
				tested_chro_begin = read_head_abs_offset;

			matched_bases_to_site = match_chro(read_text, value_index, tested_chro_begin, tested_read_pos, explain_context -> current_is_strand_jumped, global_context -> config.space_type);

			/*
			#warning "========= COMMENT TWO LINES ===================="
			SUBREADprintf("MBASETOSITE=%d, tested_read_pos=%d\n", matched_bases_to_site, tested_read_pos);
			SUBREADprintf("TXT=%s, tested_read_pos=%d\n", read_text, tested_chro_begin);
			*/

			int this_round_junction_scanned = 0;

			if(0 && FIXLENstrcmp("R010442852", explain_context -> read_name) == 0)
				SUBREADprintf("F_JUMP?  match=%d / tested=%d\n", matched_bases_to_site , tested_read_pos);

			//#warning "========= remove - 2000 from next line ============="
			if(explain_context -> total_tries < REALIGN_TOTAL_TRIES && tested_read_pos >0 && ( matched_bases_to_site*10000/tested_read_pos > 9000 - 2000 || global_context->config.maximise_sensitivity_indel) )
				for(xk1 = 0; xk1 < site_events_no ; xk1++)
				{
					chromosome_event_t * tested_event = site_events[xk1];

					if(explain_context -> is_fully_covered && tested_event -> event_type == CHRO_EVENT_TYPE_FUSION && tested_event -> event_large_side - tested_event -> event_small_side > MAX_DELETION_LENGTH){
						continue;
					}
					//if(explain_context -> pair_number == 23)
					if(0 && FIXLENstrcmp("R010442852", explain_context -> read_name) == 0){
						SUBREADprintf("F_JUMP?%d > %d    %s (%u) ; SEARCH_TAG=%u\n", (1+matched_bases_to_site)*10000 / tested_read_pos , 9000, read_text, tested_chro_begin, potential_event_pos);
						debug_show_event(global_context, tested_event);

					}

					// note that these two values are the index of the first wanted base.
					unsigned int new_read_head_abs_offset;

					if((global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection) && tested_event -> event_type == CHRO_EVENT_TYPE_INDEL)
					{
						if(explain_context ->current_is_strand_jumped){
							if(potential_event_pos == tested_event-> event_small_side) continue; 
						}else{
							if(potential_event_pos == tested_event-> event_large_side) continue; 
						}
					}
					if( tested_event -> event_type != CHRO_EVENT_TYPE_INDEL){
						this_round_junction_scanned = 1;
					}

					if((global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection))// && tested_event->event_type == CHRO_EVENT_TYPE_FUSION)
						new_read_head_abs_offset = (potential_event_pos == tested_event -> event_large_side)?tested_event -> event_small_side:tested_event -> event_large_side;
					else
						new_read_head_abs_offset = tested_event -> event_large_side;


					short new_remainder_len = remainder_len - tested_read_pos + min(0, tested_event->indel_length) - tested_event -> indel_at_junction;

	//				#warning "SUBREAD_151 REMOVE THIS ASSERTION! "
				//	assert(new_remainder_len < 102);

					if(new_remainder_len>0)
					{
						//if(explain_context -> pair_number==2074) printf("JUMPPED IN.\n");

						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].read_pos_end = explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].read_pos_start + tested_read_pos;
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].event_after_section = tested_event;
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].is_connected_to_large_side = (potential_event_pos == tested_event -> event_large_side);
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].read_pos_start = tested_read_pos - min(0, tested_event -> indel_length) + tested_event -> indel_at_junction;
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].abs_offset_for_start = new_read_head_abs_offset;
					

						if(tested_event->event_type == CHRO_EVENT_TYPE_FUSION) jump_penalty = 2;

						int current_is_jumped = explain_context -> current_is_strand_jumped;
						int current_sup_as_complex = explain_context -> tmp_min_support_as_complex;
						int current_sup_as_simple = explain_context -> tmp_support_as_simple;
						//int current_unsup_as_simple = explain_context -> tmp_min_unsupport;
						int current_pure_donor_found = explain_context -> tmp_is_pure_donor_found_explain;

						explain_context -> tmp_support_as_simple = tested_event -> supporting_reads;
						explain_context -> tmp_min_support_as_complex = min((tested_event -> is_donor_found_or_annotation & 64)?0x7fffffff:tested_event -> supporting_reads,explain_context -> tmp_min_support_as_complex);
						explain_context -> tmp_min_unsupport = min(tested_event -> anti_supporting_reads,explain_context -> tmp_min_unsupport);
						explain_context -> tmp_is_pure_donor_found_explain = explain_context -> tmp_is_pure_donor_found_explain && tested_event -> is_donor_found_or_annotation;
						explain_context -> tmp_indel_penalty += ( tested_event -> event_type == CHRO_EVENT_TYPE_INDEL );

						if(tested_event -> event_type == CHRO_EVENT_TYPE_FUSION && tested_event -> is_strand_jumped)
							explain_context -> current_is_strand_jumped = !explain_context -> current_is_strand_jumped;

						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].is_strand_jumped = explain_context -> current_is_strand_jumped;
						explain_context -> tmp_search_sections ++;

						explain_context -> total_tries ++;
						search_events_to_front(global_context, thread_context, explain_context, read_text + tested_event -> indel_at_junction + tested_read_pos -  min(0, tested_event->indel_length), qual_text + tested_read_pos -  min(0, tested_event->indel_length), new_read_head_abs_offset, new_remainder_len, sofar_matched + matched_bases_to_site - jump_penalty, tested_event -> connected_next_event_distance, 0);
						explain_context -> tmp_search_sections --;

						explain_context -> current_is_strand_jumped = current_is_jumped;
						explain_context -> tmp_indel_penalty -= ( tested_event -> event_type == CHRO_EVENT_TYPE_INDEL );
						explain_context -> tmp_min_support_as_complex = current_sup_as_complex;
						explain_context -> tmp_support_as_simple = current_sup_as_simple;
						//explain_context -> tmp_min_unsupport = current_unsup_as_simple;
						explain_context -> tmp_is_pure_donor_found_explain = current_pure_donor_found;
					}
					//if(global_context ->config.limited_tree_scan) break;
				}
			if( (global_context ->config.limited_tree_scan) && explain_context -> full_read_len <= EXON_LONG_READ_LENGTH) break;
			is_junction_scanned = max(is_junction_scanned, this_round_junction_scanned);
		}
	}
	//#warning "SUBREAD_151 REMOVE THE ASSERT! "
	//assert( remainder_len< 102 );
	int whole_section_matched = match_chro(read_text , value_index, explain_context -> current_is_strand_jumped?read_head_abs_offset - remainder_len +1:read_head_abs_offset, remainder_len , explain_context -> current_is_strand_jumped, global_context -> config.space_type);
 
	explain_context -> tmp_total_matched_bases = whole_section_matched + sofar_matched ;	

	new_explain_try_replace(global_context, thread_context, explain_context, remainder_len, 0);
}

void new_explain_try_replace(global_context_t* global_context, thread_context_t * thread_context, explain_context_t * explain_context, int remainder_len, int search_to_back)
{
	int is_better_result = 0, is_same_best = 0;

	if(0 && FIXLENstrcmp("simulated.11420793", explain_context->read_name)==0){
		SUBREADprintf("TRY_REPLACE : %s has best=%d, b_evn=%d, tscore=%d, t_evn=%d\n", explain_context->read_name, explain_context -> best_matching_bases ,  explain_context -> best_is_complex , explain_context-> tmp_total_matched_bases,  explain_context -> tmp_search_sections );
	}

	if(explain_context -> best_matching_bases - explain_context -> best_indel_penalty < explain_context-> tmp_total_matched_bases - explain_context -> tmp_indel_penalty)
	{
		is_better_result = 1;
		explain_context -> best_is_complex = explain_context -> tmp_search_sections ;
		explain_context -> is_currently_tie = 0;
		explain_context -> best_support_as_simple = explain_context -> tmp_support_as_simple;
		explain_context -> best_min_unsupport_as_simple = explain_context -> tmp_min_unsupport;
		explain_context -> best_min_support_as_complex = explain_context -> tmp_min_support_as_complex;
		explain_context -> best_is_pure_donor_found_explain = explain_context -> tmp_is_pure_donor_found_explain;
		explain_context -> second_best_matching_bases = max(explain_context -> second_best_matching_bases, explain_context -> best_matching_bases); 
		explain_context -> best_matching_bases = explain_context-> tmp_total_matched_bases ;
		explain_context -> best_indel_penalty = explain_context -> tmp_indel_penalty;
	}
	else if(explain_context -> best_matching_bases - explain_context -> best_indel_penalty == explain_context-> tmp_total_matched_bases - explain_context -> tmp_indel_penalty)
	{
		// only gapped explainations are complex counted.
		explain_context -> best_is_complex +=  explain_context -> tmp_search_sections;
		explain_context -> second_best_matching_bases = explain_context -> best_matching_bases;
		explain_context -> best_indel_penalty = explain_context -> tmp_indel_penalty;

		if(0 && FIXLENstrcmp("R010442852", explain_context -> read_name) == 0){
			SUBREADprintf("complexity: curr=%d, new=%d   ;   sections=%d\n", explain_context->best_min_support_as_complex, explain_context -> tmp_min_support_as_complex, explain_context -> tmp_search_sections );
		}
		if(explain_context -> best_is_complex > 1)
		{
			// is complex now!
			if(explain_context -> tmp_search_sections == 0)
			{
				if(explain_context -> tmp_min_unsupport >explain_context->best_min_support_as_complex){
					is_better_result = 1;
					explain_context->best_min_support_as_complex =explain_context -> tmp_min_unsupport;
					explain_context -> best_is_pure_donor_found_explain = explain_context -> tmp_is_pure_donor_found_explain;
					explain_context -> is_currently_tie = 0;
				}
				else if(explain_context -> tmp_min_unsupport == explain_context->best_min_support_as_complex)
				{
					explain_context -> is_currently_tie = 1;
					is_same_best = 1;
				}
			}
			else{
				if(explain_context -> tmp_min_support_as_complex  >explain_context->best_min_support_as_complex){
					is_better_result = 1;
					explain_context -> best_min_support_as_complex =explain_context -> tmp_min_support_as_complex;
					explain_context -> best_is_pure_donor_found_explain = explain_context -> tmp_is_pure_donor_found_explain;
					explain_context -> is_currently_tie = 0;
				}
				else if(explain_context -> tmp_min_support_as_complex  == explain_context->best_min_support_as_complex){
					explain_context -> is_currently_tie = 1;
					is_same_best = 1;
				}
			}

		}
		else
		{
			// this branch is reached ONLY if the last best is ONE-gapped (50M3D50M) and the current best is ungapped (100M)!
			if(explain_context -> best_is_pure_donor_found_explain)
			{
				if(explain_context -> best_min_unsupport_as_simple >= explain_context -> best_support_as_simple+2)
				{
					is_better_result = 1;
					explain_context -> best_min_support_as_complex = explain_context -> best_min_unsupport_as_simple;
					explain_context -> best_is_pure_donor_found_explain = explain_context -> tmp_is_pure_donor_found_explain;
					explain_context -> is_currently_tie = 0;
				}
			}
	//#warning "======= MAKE if(0) IS CORRECT BEFORE RELEASE ======"
			else if(0)
				if(explain_context -> best_min_unsupport_as_simple >= explain_context -> best_support_as_simple)
				{
					is_better_result = 1;
					explain_context -> best_min_support_as_complex = explain_context -> best_min_unsupport_as_simple;
					explain_context -> best_is_pure_donor_found_explain = explain_context -> tmp_is_pure_donor_found_explain;
					explain_context -> is_currently_tie = 0;
				}
		}
	}
	else return;

	if(is_better_result || is_same_best){
		if(search_to_back){
			explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].read_pos_start =  0;
		}else{
			explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].read_pos_end = explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].read_pos_start + remainder_len;
			explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].event_after_section = NULL;
		}
	}

	if(0 && FIXLENstrcmp("R010442852", explain_context -> read_name) == 0){
		SUBREADprintf("TRY_REPLACE_DESICION TO %s: BETTER=%d, SAME=%d ; CURRENT : %d secs ; NEWBEST : %d secs\n", search_to_back?"BACK":"FRONT", is_better_result, is_same_best, search_to_back? explain_context -> result_back_junction_numbers[0]:explain_context -> result_front_junction_numbers[0] ,explain_context -> tmp_search_sections);
		int xx1;
		for(xx1 = 0; xx1 < explain_context -> tmp_search_sections;xx1++){
			SUBREADprintf("  Event : %d ~ %d in read\n", explain_context -> tmp_search_junctions[xx1].read_pos_start, explain_context -> tmp_search_junctions[xx1].read_pos_end);
			if(explain_context -> tmp_search_junctions[xx1].event_after_section){
				SUBREADprintf("    ");
				debug_show_event(global_context, explain_context -> tmp_search_junctions[xx1].event_after_section);
			}
		}
	}

	if(is_better_result)
	{
		if(search_to_back){
			explain_context -> all_back_alignments = 1;
			explain_context -> result_back_junction_numbers[0] = explain_context -> tmp_search_sections +1;
			// checked: memory boundary
			memcpy(explain_context -> result_back_junctions[0], explain_context -> tmp_search_junctions , sizeof(perfect_section_in_read_t) * (explain_context -> tmp_search_sections +1)); 
	
		}else{
			explain_context -> all_front_alignments = 1;
			explain_context -> result_front_junction_numbers[0] = explain_context -> tmp_search_sections +1;
			// checked: memory boundary
			memcpy(explain_context -> result_front_junctions[0], explain_context -> tmp_search_junctions , sizeof(perfect_section_in_read_t) * (explain_context -> tmp_search_sections +1)); 
		}

	}else if(is_same_best){
		if(search_to_back && explain_context -> all_back_alignments < MAX_ALIGNMENT_PER_ANCHOR){
			explain_context -> result_back_junction_numbers[explain_context -> all_back_alignments] = explain_context -> tmp_search_sections +1;

			// checked: memory boundary
			memcpy(explain_context -> result_back_junctions[explain_context -> all_back_alignments], explain_context -> tmp_search_junctions , sizeof(perfect_section_in_read_t) * (explain_context -> tmp_search_sections +1)); 
			explain_context -> all_back_alignments ++;
		}else if((!search_to_back) && explain_context -> all_front_alignments < MAX_ALIGNMENT_PER_ANCHOR){
			explain_context -> result_front_junction_numbers[explain_context -> all_front_alignments] = explain_context -> tmp_search_sections +1;

			// checked: memory boundary
			memcpy(explain_context -> result_front_junctions[explain_context -> all_front_alignments], explain_context -> tmp_search_junctions , sizeof(perfect_section_in_read_t) * (explain_context -> tmp_search_sections +1)); 
			explain_context -> all_front_alignments ++;
		}
	}
}



void new_explain_try_replace_xe(global_context_t* global_context, thread_context_t * thread_context, explain_context_t * explain_context, int remainder_len, int search_to_back)
{
	int is_better_result = 0, is_same_best = 0;

	//SUBREADprintf("TRYING SET %s %s : Matched_bases : %d -> %d ; SECS : %d -> %d\n", explain_context -> read_name, search_to_back?"BACK":"FRONT", explain_context -> best_matching_bases, explain_context-> tmp_total_matched_bases, search_to_back? explain_context -> result_back_junction_numbers[0]:explain_context -> result_front_junction_numbers[0] ,explain_context -> tmp_search_sections );

	if(explain_context -> best_matching_bases < explain_context-> tmp_total_matched_bases)
	{
		is_better_result = 1;
		explain_context -> best_is_complex = explain_context -> tmp_search_sections ;
		explain_context -> is_currently_tie = 0;
		explain_context -> best_support_as_simple = explain_context -> tmp_support_as_simple;
		explain_context -> best_min_unsupport_as_simple = explain_context -> tmp_min_unsupport;
		explain_context -> best_min_support_as_complex = explain_context -> tmp_min_support_as_complex;
		explain_context -> best_is_pure_donor_found_explain = explain_context -> tmp_is_pure_donor_found_explain;
		explain_context -> second_best_matching_bases = max(explain_context -> second_best_matching_bases, explain_context -> best_matching_bases); 
		explain_context -> best_matching_bases = explain_context-> tmp_total_matched_bases ;

	} else if(explain_context -> tmp_search_sections < search_to_back? explain_context -> result_back_junction_numbers[0]:explain_context -> result_front_junction_numbers[0] - 1 && explain_context -> best_matching_bases == explain_context-> tmp_total_matched_bases) {
		is_better_result = 1;
		explain_context -> best_is_complex = explain_context -> tmp_search_sections ;
		explain_context -> is_currently_tie = 0;
		explain_context -> best_support_as_simple = explain_context -> tmp_support_as_simple;
		explain_context -> best_min_unsupport_as_simple = explain_context -> tmp_min_unsupport;
		explain_context -> best_min_support_as_complex = explain_context -> tmp_min_support_as_complex;
		explain_context -> best_is_pure_donor_found_explain = explain_context -> tmp_is_pure_donor_found_explain;
		explain_context -> second_best_matching_bases = max(explain_context -> second_best_matching_bases, explain_context -> best_matching_bases); 
		explain_context -> best_matching_bases = explain_context-> tmp_total_matched_bases ;
	} else if( 0 && explain_context -> best_matching_bases == explain_context-> tmp_total_matched_bases) {
		// only gapped explainations are complex counted.
		explain_context -> best_is_complex +=  explain_context -> tmp_search_sections;
		explain_context -> second_best_matching_bases = explain_context -> best_matching_bases;

		if(0 && FIXLENstrcmp("R010442852", explain_context -> read_name) == 0){
			SUBREADprintf("complexity: curr=%d, new=%d   ;   sections=%d\n", explain_context->best_min_support_as_complex, explain_context -> tmp_min_support_as_complex, explain_context -> tmp_search_sections );
		}
		if(explain_context -> best_is_complex > 1)
		{
			// is complex now!
			if(explain_context -> tmp_search_sections == 0)
			{
				if(explain_context -> tmp_min_unsupport >explain_context->best_min_support_as_complex){
					is_better_result = 1;
					explain_context->best_min_support_as_complex =explain_context -> tmp_min_unsupport;
					explain_context -> best_is_pure_donor_found_explain = explain_context -> tmp_is_pure_donor_found_explain;
					explain_context -> is_currently_tie = 0;
				}
				else if(explain_context -> tmp_min_unsupport == explain_context->best_min_support_as_complex)
				{
					explain_context -> is_currently_tie = 1;
					is_same_best = 1;
				}
			} else {
				if(explain_context -> tmp_min_support_as_complex  >explain_context->best_min_support_as_complex){
					is_better_result = 1;
					explain_context -> best_min_support_as_complex =explain_context -> tmp_min_support_as_complex;
					explain_context -> best_is_pure_donor_found_explain = explain_context -> tmp_is_pure_donor_found_explain;
					explain_context -> is_currently_tie = 0;
				}
				else if(explain_context -> tmp_min_support_as_complex  == explain_context->best_min_support_as_complex){
					explain_context -> is_currently_tie = 1;
					is_same_best = 1;
				}
			}

		} else {
			// this branch is reached ONLY if the last best is ONE-gapped (50M3D50M) and the current best is ungapped (100M)!
			if(explain_context -> best_is_pure_donor_found_explain)
			{
				if(explain_context -> best_min_unsupport_as_simple >= explain_context -> best_support_as_simple+2)
				{
					is_better_result = 1;
					explain_context -> best_min_support_as_complex = explain_context -> best_min_unsupport_as_simple;
					explain_context -> best_is_pure_donor_found_explain = explain_context -> tmp_is_pure_donor_found_explain;
					explain_context -> is_currently_tie = 0;
				}
			}
		}
	}
	else return;

	if(is_better_result || is_same_best){
		if(search_to_back){
			explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].read_pos_start =  0;
		}else{
			explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].read_pos_end = explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].read_pos_start + remainder_len;
			explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].event_after_section = NULL;
		}
	}

	if(0 && FIXLENstrcmp("R010442852", explain_context -> read_name) == 0){
		SUBREADprintf("RNAME=%s TRY_REPLACE_DESICION TO %s: BETTER=%d, SAME=%d ; CURRENT : %d secs ; NEWBEST : %d secs MM_bases : %d -> %d\n", explain_context -> read_name, search_to_back?"BACK":"FRONT", is_better_result, is_same_best, search_to_back? explain_context -> result_back_junction_numbers[0]:explain_context -> result_front_junction_numbers[0] ,explain_context -> tmp_search_sections, explain_context -> best_matching_bases , explain_context-> tmp_total_matched_bases);
		int xx1;
		for(xx1 = 0; xx1 < explain_context -> tmp_search_sections;xx1++){
			SUBREADprintf("  Event : %d ~ %d in read\n", explain_context -> tmp_search_junctions[xx1].read_pos_start, explain_context -> tmp_search_junctions[xx1].read_pos_end);
			if(explain_context -> tmp_search_junctions[xx1].event_after_section){
				SUBREADprintf("    ");
				debug_show_event(global_context, explain_context -> tmp_search_junctions[xx1].event_after_section);
			}
		}
	}

	if(is_better_result)
	{
		if(search_to_back){
			explain_context -> all_back_alignments = 1;
			explain_context -> result_back_junction_numbers[0] = explain_context -> tmp_search_sections +1;
			// checked: memory boundary
			memcpy(explain_context -> result_back_junctions[0], explain_context -> tmp_search_junctions , sizeof(perfect_section_in_read_t) * (explain_context -> tmp_search_sections +1)); 
	
		}else{
			explain_context -> all_front_alignments = 1;
			explain_context -> result_front_junction_numbers[0] = explain_context -> tmp_search_sections +1;
			// checked: memory boundary
			memcpy(explain_context -> result_front_junctions[0], explain_context -> tmp_search_junctions , sizeof(perfect_section_in_read_t) * (explain_context -> tmp_search_sections +1)); 
		}

	}else if(is_same_best && 0){
		if(search_to_back && explain_context -> all_back_alignments < MAX_ALIGNMENT_PER_ANCHOR){
			explain_context -> result_back_junction_numbers[explain_context -> all_back_alignments] = explain_context -> tmp_search_sections +1;

			// checked: memory boundary
			memcpy(explain_context -> result_back_junctions[explain_context -> all_back_alignments], explain_context -> tmp_search_junctions , sizeof(perfect_section_in_read_t) * (explain_context -> tmp_search_sections +1)); 
			explain_context -> all_back_alignments ++;
		}else if((!search_to_back) && explain_context -> all_front_alignments < MAX_ALIGNMENT_PER_ANCHOR){
			explain_context -> result_front_junction_numbers[explain_context -> all_front_alignments] = explain_context -> tmp_search_sections +1;

			// checked: memory boundary
			memcpy(explain_context -> result_front_junctions[explain_context -> all_front_alignments], explain_context -> tmp_search_junctions , sizeof(perfect_section_in_read_t) * (explain_context -> tmp_search_sections +1)); 
			explain_context -> all_front_alignments ++;
		}
	}
}


// read_tail_abs_offset is actually the offset of the base next to the last base in read tail.
// read_tail_pos is the FIRST UNWANTED BASE, after the read.
void search_events_to_back(global_context_t * global_context, thread_context_t * thread_context, explain_context_t * explain_context, char * read_text , char * qual_text, unsigned int read_tail_abs_offset, short read_tail_pos, short sofar_matched, int suggested_movement, int do_not_jump)
{
	short tested_read_pos;

	HashTable * event_table = NULL;
	chromosome_event_t * event_space = NULL;
	
	if(thread_context) {
		event_table = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	} else {
		event_table = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}

	gene_value_index_t * value_index = thread_context?thread_context->current_value_index:global_context->current_value_index ;
	if( there_are_events_in_range(event_table -> appendix2, read_tail_abs_offset - read_tail_pos, read_tail_pos) || global_context ->  config.do_fusion_detection ||global_context ->  config.do_long_del_detection ){
		int event_search_method;
		if((global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection))
			event_search_method = EVENT_SEARCH_BY_BOTH_SIDES;
		else
			event_search_method = EVENT_SEARCH_BY_LARGE_SIDE;
	
	
		int is_junction_scanned = 0;
		// minimum perfect section length is 1
		// tested_read_pos is the first WANTED BASE in section.
		int move_start = read_tail_pos - (do_not_jump?0:global_context -> config.realignment_minimum_variant_distance);
		if(suggested_movement) move_start = read_tail_pos - suggested_movement + 1;
	
	//#warning ">>>>>>>>>>>>>> COMMENT THIS <<<<<<<<<<<<<<<<<<<<<"
	//printf("OCT27-STEP-BKIN : %s , STT=%d, %u, %d\n", explain_context -> read_name, move_start, read_tail_abs_offset, read_tail_pos);

		if(MAX_EVENTS_IN_READ - 1> explain_context -> tmp_search_sections)
		for(tested_read_pos =  move_start; tested_read_pos >=0;tested_read_pos --)
		{
			int xk1, matched_bases_to_site;
			int jump_penalty = 0;
			chromosome_event_t *site_events[MAX_EVENT_ENTRIES_PER_SITE];

			int potential_event_pos;

			if(explain_context -> current_is_strand_jumped)
				potential_event_pos = read_tail_abs_offset + ( read_tail_pos - tested_read_pos);
			else
				potential_event_pos = read_tail_abs_offset - ( read_tail_pos - tested_read_pos);
	

			if(!check_event_bitmap(  event_table->appendix2, potential_event_pos )) continue;
			int search_types = CHRO_EVENT_TYPE_INDEL | CHRO_EVENT_TYPE_JUNCTION | CHRO_EVENT_TYPE_FUSION;
			int site_events_no = search_event(global_context, event_table , event_space , potential_event_pos, event_search_method , search_types, site_events);
			//#warning ">>>>>>>>>>>>>> COMMENT THIS <<<<<<<<<<<<<<<<<<<<<"
			//printf("OCT27-STEP-BKIN-SR: %s at %u, FOUND=%d\n" , explain_context -> read_name,potential_event_pos,site_events_no);

			if(!site_events_no)continue;

			unsigned int tested_chro_begin;
			if(explain_context -> current_is_strand_jumped)
				tested_chro_begin = read_tail_abs_offset + 1;
			else
				tested_chro_begin = read_tail_abs_offset - (read_tail_pos - tested_read_pos);

			matched_bases_to_site = match_chro(read_text + tested_read_pos, value_index, tested_chro_begin , read_tail_pos - tested_read_pos, explain_context -> current_is_strand_jumped, global_context -> config.space_type);

			int this_round_junction_scanned = 0;

			//#warning ">>>>>>>>>>>>>>>> REMOVE IT <<<<<<<<<<<<<<<<<<<<<<"
			//printf("OCT27-STEPSB-JB-%s: test %u = %d events; TEST=%d > 7000 : MA=%d; %s ; %u = %u - (%d - %d) ; LEV=%d\n", explain_context -> read_name, potential_event_pos, site_events_no, (read_tail_pos<=tested_read_pos)?(-1234):( matched_bases_to_site*10000/(read_tail_pos - tested_read_pos)) , matched_bases_to_site, read_text + tested_read_pos, potential_event_pos, read_tail_abs_offset, read_tail_pos, tested_read_pos, explain_context -> tmp_search_sections);
			//#warning "========= remove - 2000 from next line ============="
			if(explain_context -> total_tries < REALIGN_TOTAL_TRIES && (read_tail_pos>tested_read_pos) && ( matched_bases_to_site*10000/(read_tail_pos - tested_read_pos) > 9000 - 2000 || global_context->config.maximise_sensitivity_indel) )
				for(xk1 = 0; xk1 < site_events_no ; xk1++)
				{
					chromosome_event_t * tested_event = site_events[xk1];

					if(explain_context -> is_fully_covered && tested_event -> event_type == CHRO_EVENT_TYPE_FUSION && tested_event -> event_large_side - tested_event -> event_small_side > MAX_DELETION_LENGTH){
						continue;
					}

					if((global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection) && tested_event -> event_type == CHRO_EVENT_TYPE_INDEL)
					{
						if(explain_context->current_is_strand_jumped){
							if(potential_event_pos == tested_event-> event_large_side) continue; 
						}else{
							if(potential_event_pos == tested_event-> event_small_side) continue; 
						}
					}
					if( tested_event -> event_type != CHRO_EVENT_TYPE_INDEL){
						this_round_junction_scanned = 1;
					}


					if(0 && strcmp("S_chr901_565784_72M8D28M", explain_context -> read_name) == 0)
						SUBREADprintf("B_JUMP?%d > %d TLEN=%d \n", (1+matched_bases_to_site)*10000 / (read_tail_pos - tested_read_pos) , 9000, read_tail_pos - tested_read_pos);
					// note that read_tail_pos is the first unwanted base.
					int new_read_tail_pos = tested_read_pos;
					if(tested_event->event_type == CHRO_EVENT_TYPE_INDEL) new_read_tail_pos +=  min(0, tested_event -> indel_length);
					// note that read_tail_abs_offset is the first unwanted base.
					unsigned int new_read_tail_abs_offset;

					if((global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection))// && tested_event->event_type == CHRO_EVENT_TYPE_FUSION)
					{
						new_read_tail_abs_offset = (potential_event_pos == tested_event -> event_small_side)? tested_event -> event_large_side : tested_event -> event_small_side;
						if(tested_event->is_strand_jumped + explain_context -> current_is_strand_jumped == 1)
							new_read_tail_abs_offset--;
						else
							new_read_tail_abs_offset++;
					}
					else
						new_read_tail_abs_offset = tested_event -> event_small_side + 1;

					new_read_tail_pos -= tested_event -> indel_at_junction;

					if(new_read_tail_pos>0)
					{
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].read_pos_start = tested_read_pos;
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].event_after_section = tested_event;
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].is_connected_to_large_side = (potential_event_pos == tested_event -> event_small_side);
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].read_pos_end = tested_read_pos + min(0, tested_event->indel_length) - tested_event -> indel_at_junction;
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].abs_offset_for_start = new_read_tail_abs_offset; 

						if(tested_event->event_type == CHRO_EVENT_TYPE_FUSION) jump_penalty = 2;

						int current_is_jumped = explain_context -> current_is_strand_jumped ;
						int current_sup_as_complex = explain_context -> tmp_min_support_as_complex;
						int current_sup_as_simple = explain_context -> tmp_support_as_simple;
						int current_pure_donor_found = explain_context -> tmp_is_pure_donor_found_explain;

						explain_context -> tmp_support_as_simple = tested_event -> supporting_reads;
						explain_context -> tmp_min_support_as_complex = min((tested_event -> is_donor_found_or_annotation & 64)?0x7fffffff:tested_event -> supporting_reads,explain_context -> tmp_min_support_as_complex);
						explain_context -> tmp_min_unsupport = min(tested_event -> anti_supporting_reads,explain_context -> tmp_min_unsupport);
						explain_context -> tmp_is_pure_donor_found_explain = explain_context -> tmp_is_pure_donor_found_explain && tested_event -> is_donor_found_or_annotation;
						explain_context -> tmp_indel_penalty += ( tested_event -> event_type == CHRO_EVENT_TYPE_INDEL );

						if(tested_event -> event_type == CHRO_EVENT_TYPE_FUSION && tested_event -> is_strand_jumped)
							explain_context -> current_is_strand_jumped = !explain_context -> current_is_strand_jumped;
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].is_strand_jumped = explain_context -> current_is_strand_jumped;
						explain_context -> tmp_search_sections ++;
						explain_context -> total_tries ++;

						search_events_to_back(global_context, thread_context, explain_context, read_text , qual_text, new_read_tail_abs_offset , new_read_tail_pos, sofar_matched + matched_bases_to_site - jump_penalty, tested_event -> connected_previous_event_distance, 0);

						explain_context -> tmp_search_sections --;
						explain_context -> tmp_indel_penalty -= ( tested_event -> event_type == CHRO_EVENT_TYPE_INDEL );
						explain_context -> current_is_strand_jumped = current_is_jumped;
						explain_context -> tmp_min_support_as_complex = current_sup_as_complex;
						explain_context -> tmp_support_as_simple = current_sup_as_simple;
						explain_context -> tmp_is_pure_donor_found_explain = current_pure_donor_found;
					}
				}
			if(( global_context ->config.limited_tree_scan) && explain_context -> full_read_len <= EXON_LONG_READ_LENGTH) break;
			this_round_junction_scanned = max(this_round_junction_scanned, is_junction_scanned);
		} 
	}
	int whole_section_matched = match_chro(read_text , value_index, read_tail_abs_offset - (explain_context -> current_is_strand_jumped?-1:read_tail_pos), read_tail_pos , explain_context -> current_is_strand_jumped, global_context -> config.space_type);
 
	explain_context -> tmp_total_matched_bases = whole_section_matched + sofar_matched ;	

	new_explain_try_replace(global_context, thread_context, explain_context, 0, 1);
}

int init_junction_tables(global_context_t * context)
{
	fraglist_init(&context -> funky_list_A);
	fraglist_init(&context -> funky_list_DE);

	bktable_init(&context -> funky_table_BC, FUNKY_COLOCATION_TOLERANCE * 2, 10000000);
	bktable_init(&context -> funky_table_DE, FUNKY_COLOCATION_TOLERANCE * 2, 10000000);

	bktable_init(&context -> breakpoint_table_P, 2 * context -> config.maximum_pair_distance, 1000000);
	bktable_init(&context -> breakpoint_table_QR, 2 * BREAK_POINT_MAXIMUM_TOLERANCE, 1000000);
	bktable_init(&context -> breakpoint_table_YZ, 2 * context -> config.maximum_pair_distance, 1000000);

	bktable_init(&context -> translocation_result_table, 2*BREAK_POINT_MAXIMUM_TOLERANCE, 1000000);
	bktable_init(&context -> inversion_result_table, 2*BREAK_POINT_MAXIMUM_TOLERANCE, 1000000);
	return 0;
}

int destroy_junction_tables(global_context_t * context)
{
	fraglist_destroy(&context -> funky_list_A);
	fraglist_destroy(&context -> funky_list_DE);

	bktable_destroy(&context -> funky_table_BC);
	bktable_destroy(&context -> funky_table_DE);
	bktable_destroy(&context -> breakpoint_table_P);
	bktable_destroy(&context -> breakpoint_table_QR);
	bktable_destroy(&context -> breakpoint_table_YZ);

	HashTableIteration(context -> inversion_result_table.entry_table , bktable_free_ptrs);
	bktable_destroy(&context -> inversion_result_table);

	HashTableIteration(context -> translocation_result_table.entry_table , bktable_free_ptrs);
	bktable_destroy(&context -> translocation_result_table);

	return 0;
}
int init_junction_thread_contexts(global_context_t * global_context, thread_context_t * thread_context, int task)
{
	return 0;
}

void insert_big_margin_record(global_context_t * global_context , unsigned short * big_margin_record, unsigned char votes, short read_pos_start, short read_pos_end, int read_len, int is_negative)
{

	if( global_context->config.big_margin_record_size<3) return;

	unsigned short read_pos_start_2 = (is_negative?read_len -read_pos_end:read_pos_start) ;
	unsigned short read_pos_end_2 = (is_negative?read_len -read_pos_start:read_pos_end);

	int xk1;
	for(xk1=0; xk1< global_context->config.big_margin_record_size / 3; xk1++)
	{
		if( votes >= big_margin_record[xk1*3])
			break;
	}
	if(xk1< global_context->config.big_margin_record_size / 3)
	{
		int xk2;
		for(xk2 = global_context->config.big_margin_record_size-4; xk2 >= xk1*3; xk2--)
			big_margin_record[xk2 + 3] = big_margin_record[xk2];
		big_margin_record[xk1*3+0] = votes;
		big_margin_record[xk1*3+1] = read_pos_start_2;
		big_margin_record[xk1*3+2] = read_pos_end_2;
	}
}

int is_PE_distance(global_context_t * global_context, unsigned int pos1,  unsigned int pos2, int rlen1, int rlen2, int is_negative_R1, int is_negative_R2)
{
	long long int dist = pos2;
	dist -= pos1;

	is_negative_R1 = (is_negative_R1>0)?1:0;
	is_negative_R2 = (is_negative_R2>0)?1:0;

	if(pos1 > pos2) dist -= rlen1;
	else if(pos1 < pos2) dist += rlen2;
	else dist += max(rlen2, rlen1);

	if(abs(dist) > global_context->config.maximum_pair_distance || abs(dist)<global_context->config.minimum_pair_distance) return 0;

	if(is_negative_R1 != is_negative_R2) return 0;
	if(pos1 > pos2 && !is_negative_R1) return 0;
	if(pos1 < pos2 && is_negative_R1) return 0;
	return 1;
}


#define MAX_VOTE_TOLERANCE 1
//returns 1 if the vote number is not significantly higher than the vote numbers in the vote list. 
int test_small_minor_votes(global_context_t * global_context, int minor_i, int minor_j, int major_i, int major_j , gene_vote_t * votes, int read_len)
{
	int is_small_margin_minor = 0;
	long long dist = votes -> pos[minor_i][minor_j];
	dist -= votes -> pos[major_i][major_j];

	if(abs(dist)> global_context->config.maximum_intron_length)
	{
		int iii, jjj;
		for(iii=0; iii<GENE_VOTE_TABLE_SIZE; iii++)
		{
			for(jjj = 0; jjj < votes->items[iii]; jjj++)
			{
				if(iii == minor_i && jjj == minor_j) continue;
				// "2" is the tolerance.
				if(votes -> votes[minor_i][minor_j] - votes -> votes[iii][jjj] >=1) continue;

				int minor_coverage_start = votes -> coverage_start[minor_i][minor_j] ;
				int minor_coverage_end = votes -> coverage_end[minor_i][minor_j] ;

				int other_coverage_start = votes -> coverage_start[iii][jjj];
				int other_coverage_end = votes -> coverage_end[iii][jjj];

				int minor_negative = votes -> masks[minor_i][minor_j] & IS_NEGATIVE_STRAND;
				int other_negative = votes -> masks[iii][jjj] & IS_NEGATIVE_STRAND;

				if(minor_negative) {
					int ttt = read_len - minor_coverage_end;
					minor_coverage_end = read_len - minor_coverage_start;
					minor_coverage_start = ttt;
				}

				if(other_negative){
					int ttt = read_len - other_coverage_end;
					other_coverage_end = read_len - other_coverage_start;
					other_coverage_start = ttt;
				}

				if(abs(minor_coverage_end - other_coverage_end) < 7 && abs(minor_coverage_start - other_coverage_start)<7)
					is_small_margin_minor = 1;

				if(is_small_margin_minor) break;
			}
			if(is_small_margin_minor) break;
		}
	}	
	return is_small_margin_minor;
}


// function test_junction_minor returns 1 if the current anchor and current_vote[i][j] are not good mates in terms of junction reads:
// for example, if the distance is too far, if the coverered region overlapped or if the two mapped parts in the read are reversely arranged (expect in fusion detection) 
int test_junction_minor(global_context_t * global_context, thread_context_t * thread_context, gene_vote_t * votes, int vote_i, int vote_j, int i, int j, long long int dist)
{
	if(abs(dist)> global_context->config.maximum_intron_length) return 1; 
	if(votes -> coverage_start[vote_i][vote_j] == votes -> coverage_start[i][j])return 2;
	if(votes -> coverage_end[vote_i][vote_j]   == votes -> coverage_end[i][j])return 3;

	//SUBREADprintf( " COV_IN_READ: %d ~ %d         CHRO_POS: %u ~ %u \n", votes -> coverage_start[vote_i][vote_j] , votes -> coverage_start[i][j] , votes -> pos[vote_i][vote_j] ,  votes -> pos[i][j] );
	if(votes -> coverage_start[vote_i][vote_j] > votes -> coverage_start[i][j])
	{
		if(votes -> pos[vote_i][vote_j] < votes -> pos[i][j])return 4;
	}
	else
	{
		if(votes -> pos[vote_i][vote_j] > votes -> pos[i][j])return 5;
	}

	return 0;
}

void update_top_three(global_context_t * global_context, int * top_buffer_3i, int new_value){
	if(new_value > top_buffer_3i[global_context -> config.top_scores - 1]){
		int x1;
		for(x1 = 0;x1 < global_context -> config.top_scores ; x1++){
			if(new_value > top_buffer_3i[x1]){
				int x2;
				for(x2 = global_context -> config.top_scores - 1 ; x2 > x1 ; x2 --){
					top_buffer_3i[x2] = top_buffer_3i[x2-1];
				}
				top_buffer_3i[x1] = new_value;
				break;
			}else if(new_value == top_buffer_3i[x1]) break;
		}
	}
}



int comb_sort_compare(void * Vcomb_buffer, int i, int j){
	vote_combination_t * comb_buffer = (vote_combination_t *)Vcomb_buffer;
	return comb_buffer[i].score_adj - comb_buffer[j].score_adj;
}

void comb_sort_exchange(void * Vcomb_buffer, int i, int j){
	vote_combination_t * comb_buffer = (vote_combination_t *)Vcomb_buffer;
	vote_combination_t tmpv;
	memcpy(&tmpv, comb_buffer + i, sizeof(vote_combination_t));
	memcpy(comb_buffer + i, comb_buffer + j, sizeof(vote_combination_t));
	memcpy(comb_buffer + j, &tmpv, sizeof(vote_combination_t));
}

void comb_sort_merge(void * Vcomb_buffer, int start, int items, int items2){
	vote_combination_t * comb_buffer = (vote_combination_t *)Vcomb_buffer;
	vote_combination_t * merge_target = malloc(sizeof(vote_combination_t) * (items + items2));

	int items1_cursor = start, items2_cursor = start + items, x1;

	for(x1=0; x1 < items+items2; x1++){
		int select_items_1 = (items1_cursor < items + start && comb_sort_compare(comb_buffer, items1_cursor, items2_cursor) <=0) || (items2_cursor == start + items + items2);
		if(select_items_1){
			memcpy(merge_target+x1, comb_buffer+items1_cursor, sizeof(vote_combination_t));
			items1_cursor++;
		}else{
			memcpy(merge_target+x1, comb_buffer+items2_cursor, sizeof(vote_combination_t));
			items2_cursor++;
		}

	}

	memcpy(comb_buffer + start, merge_target, (items+items2) * sizeof(vote_combination_t));
	free(merge_target);

}

int is_better_inner(global_context_t * global_context, thread_context_t * thread_context, subjunc_result_t * junc_res, int old_intron_length,  gene_vote_number_t Vote_minor, int coverage_minor_length, int intron)
{
	if( Vote_minor > junc_res -> minor_votes ||
	  (Vote_minor ==junc_res -> minor_votes && coverage_minor_length > junc_res -> minor_coverage_end - junc_res -> minor_coverage_start) ||
	  (Vote_minor ==junc_res -> minor_votes && coverage_minor_length ==junc_res -> minor_coverage_end - junc_res -> minor_coverage_start && intron < old_intron_length))
		return 1;
	else    return 0;
}


#define COVERAGE_STAB_NUMBER 100
int test_fully_covered(global_context_t * global_context, gene_vote_t *  vote, int read_length){
	int i,j,xk1,xk2;
	char local_strands[COVERAGE_STAB_NUMBER];
	unsigned int local_locations[COVERAGE_STAB_NUMBER];
	unsigned long long local_coverage[COVERAGE_STAB_NUMBER];
	int used_stabs = 0;

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
	{
		for (j=0; j< vote->items[i]; j++)
		{
			if(vote -> votes[i][j]>2 && used_stabs < COVERAGE_STAB_NUMBER)
			{
				int is_fresh = 1;
				int is_negative = (vote -> masks[i][j] & IS_NEGATIVE_STRAND)?1:0; 
				for(xk1=0; xk1<used_stabs; xk1++){
					if(local_strands[xk1] == is_negative){
						long long dist = vote -> pos[i][j];
						dist -= local_locations[xk1];
						if(abs(dist) < MAX_DELETION_LENGTH)
						{
							is_fresh=0;
							break;
						}
					}
				}

				if(is_fresh){
					local_strands[used_stabs]=is_negative;
					local_locations[used_stabs]= vote -> pos[i][j];
					local_coverage[used_stabs] = 0;
					used_stabs++;
				}
			}
		}
	}
	if(!used_stabs) return 0;

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
	{
		for (j=0; j< vote->items[i]; j++)
		{
			if(vote -> votes[i][j]>=1)
			{
				int is_negative = (vote -> masks[i][j] & IS_NEGATIVE_STRAND)?1:0; 
				for(xk1=0; xk1<used_stabs; xk1++){
					if(local_strands[xk1] == is_negative){
						long long dist = vote -> pos[i][j];
						dist -= local_locations[xk1];
						if(abs(dist) < MAX_DELETION_LENGTH)
						{
							for(xk2 = vote -> coverage_start[i][j] * 64 / read_length; xk2 <= 
							    vote -> coverage_end[i][j] * 64 / read_length; xk2++){
								local_coverage[xk1] |= 1llu<<xk2;
							}  
						}
					}
				}
			}
		}
	}

	for(xk1=0; xk1<used_stabs; xk1++){
		int covered = 0;
		for(xk2 = 0; xk2<64; xk2++){
			covered += ( local_coverage[xk1] & (1llu<<xk2) )?1:0;
		}
		//SUBREADprintf("COVERAGE LEVEL=%d\n", covered);

		if(covered > 54){
			return 1;
		}
	}

	return 0;
}

int is_long_del_high_quality(global_context_t * global_context, thread_context_t * thread_context, int p1_start, int p1_end, int p2_start, int p2_end, int read_len, int p1_votes, int p2_votes){
	if(p1_votes < 3 || p2_votes < 3) return 0;
	if( min( p1_start, p2_start ) > 10 ) return 0;
	if( read_len - max(p1_end, p2_end)  > 10 ) return 0;
	return 1;
}
#define SE_READ_IN_KNOWN_EXON_REWARD 1

void copy_vote_to_alignment_res(global_context_t * global_context, thread_context_t * thread_context, mapping_result_t * align_res, subjunc_result_t * junc_res, gene_vote_t * current_vote, int vote_i, int vote_j, int curr_read_len, char * read_name, char * curr_read_text, int used_subreads_in_vote, int noninformative_subreads_in_vote, subread_read_number_t pair_number, int is_second_read, int * is_fully_covered)
{
	int vv = current_vote -> votes[vote_i][vote_j];
	if(global_context->config.scRNA_input_mode && !global_context -> input_reads.is_paired_end_reads) vv += SE_READ_IN_KNOWN_EXON_REWARD *is_pos_in_annotated_exon_regions(global_context, current_vote -> pos[vote_i][vote_j]);
	align_res -> selected_position = current_vote -> pos[vote_i][vote_j];
	align_res -> selected_votes = vv;
	align_res -> indels_in_confident_coverage = indel_recorder_copy(align_res -> selected_indel_record, current_vote -> indel_recorder[vote_i][vote_j]);
	align_res -> confident_coverage_end = current_vote -> coverage_end[vote_i][vote_j];
	align_res -> confident_coverage_start = current_vote -> coverage_start[vote_i][vote_j];
	align_res -> result_flags = (current_vote -> masks[vote_i][vote_j] & IS_NEGATIVE_STRAND)?(CORE_IS_NEGATIVE_STRAND):0;
	align_res -> used_subreads_in_vote = used_subreads_in_vote;
	align_res -> noninformative_subreads_in_vote = noninformative_subreads_in_vote;
	align_res -> is_fully_covered = *is_fully_covered ;

	if(global_context -> config.do_breakpoint_detection)
	{
		int i,j;

		// iterate all the anchors we have found in step 1:
		for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		{
			for (j=0; j< current_vote->items[i]; j++)
			{
				if(i == vote_i && j == vote_j) continue;
				if(align_res -> selected_votes < current_vote -> votes[i][j]) continue;	// major half must be the anchor

				long long int dist = current_vote -> pos[vote_i][vote_j]; 
				dist -= current_vote -> pos[i][j];

				int is_strand_jumpped = (current_vote -> masks[vote_i][vote_j] & IS_NEGATIVE_STRAND)!=(current_vote -> masks[i][j] & IS_NEGATIVE_STRAND);
				if((global_context-> config.do_fusion_detection || global_context-> config.do_long_del_detection) && (*is_fully_covered) && (dist > MAX_DELETION_LENGTH || is_strand_jumpped)) continue; 

				if((global_context-> config.do_fusion_detection || global_context-> config.do_long_del_detection)){
					// function test_small_minor_votes returns 1 if the vote number is not significantly
					// higher than the vote numbers in the vote list. 
					//#warning "SUBREAD_151 =========== THE TWO LINES SHOULD BE UNCOMMENTED IN RELEASED VERSION ==== WE COMMENT IT FOR A BETTER FUSION SENSITIVITY BUT ONLY FOR TEST ==================="
					if(1){
						int small_minor_bigmargin = test_small_minor_votes(global_context , i, j, vote_i, vote_j, current_vote, curr_read_len);
						if(small_minor_bigmargin) continue;
					}
				}else{
					// function test_junction_minor returns 1 if the current anchor and current_vote[i][j] 
					// are not good mates in terms of junction reads:
					//
					// for example, if the distance is too far, if the coverered region overlapped or 
					// if the covered region has a wrong arrangement to their relative positions.
					int test_minor_res = test_junction_minor(global_context, thread_context, current_vote, vote_i, vote_j, i, j, dist);
//#warning " ============ DEBUG 1 ==================== "
					if(0 && FIXLENstrcmp("R002403247", read_name) == 0) {
						char posout2[100];
						char posout1[100];
						absoffset_to_posstr(global_context, current_vote -> pos[vote_i][vote_j], posout1);
						absoffset_to_posstr(global_context, current_vote -> pos[i][j], posout2);
						SUBREADprintf("SMALL_MARGIN=%d at %s ~ %s\n", test_minor_res, posout1, posout2);
					}
					//	SUBREADprintf("TMR=%d (V=%d)\n", test_minor_res, current_vote -> votes[i][j]);
					if(test_minor_res)continue;
				}

				int is_better = is_better_inner(global_context, thread_context,
							junc_res, abs32uint(current_vote -> pos[vote_i][vote_j] - junc_res -> minor_position), current_vote -> votes[i][j], current_vote -> coverage_end[i][j] - current_vote -> coverage_start[i][j],
							abs32uint(current_vote -> pos[vote_i][vote_j] - current_vote -> pos[i][j]));

				int replace_minor = 0, minor_indel_offset = 0, inserted_bases = 0, is_GT_AG_donors = 0, is_donor_found_or_annotation = 0, final_split_point = 0, major_indels = 0, small_side_increasing_coordinate = 0, large_side_increasing_coordinate = 0;

				if(0 && FIXLENstrcmp("R002403247", read_name) == 0) 
				{
					char posout[100];
					absoffset_to_posstr(global_context, current_vote -> pos[i][j], posout);
					SUBREADprintf("IBT=%d (V=%d , OV=%d) at %s\n", is_better, current_vote -> votes[i][j], junc_res -> minor_votes, posout);
					SUBREADprintf("IBT OLD_INTRON=%d, INTRON=%d\n", abs32uint(current_vote -> pos[vote_i][vote_j] - junc_res -> minor_position),
							abs32uint(current_vote -> pos[vote_i][vote_j] - current_vote -> pos[i][j])
						);
				}

				if(is_better){
					// Determine the splicing point of the fusion or the junction
					// If the splicing point is determined, then set replace_minor = 1
					if(is_strand_jumpped){
						if(!global_context -> config.do_fusion_detection) continue;

						int minor_cover_end_as_reversed = (current_vote -> masks[i][j] & IS_NEGATIVE_STRAND)? current_vote -> coverage_end[i][j]:(curr_read_len - current_vote -> coverage_start[i][j]);
						int minor_cover_start_as_reversed = (current_vote -> masks[i][j] & IS_NEGATIVE_STRAND)? current_vote -> coverage_start[i][j]:(curr_read_len - current_vote -> coverage_end[i][j]);
						int main_cover_end_as_reversed = (current_vote -> masks[vote_i][vote_j] & IS_NEGATIVE_STRAND)?current_vote -> coverage_end[vote_i][vote_j]:(curr_read_len - current_vote -> coverage_start[vote_i][vote_j]);
						int main_cover_start_as_reversed = (current_vote -> masks[vote_i][vote_j] & IS_NEGATIVE_STRAND)?current_vote -> coverage_start[vote_i][vote_j]:(curr_read_len - current_vote -> coverage_end[vote_i][vote_j]);


						int overlapped ;
						if(main_cover_start_as_reversed > minor_cover_start_as_reversed)
							overlapped = minor_cover_end_as_reversed - main_cover_start_as_reversed;
						else
							overlapped = main_cover_end_as_reversed - minor_cover_start_as_reversed;

						if(overlapped > 14) continue;


						int guess_start_as_reversed = (main_cover_start_as_reversed > minor_cover_start_as_reversed)?
									 (minor_cover_end_as_reversed - 15): (main_cover_end_as_reversed - 15);

						int guess_end_as_reversed = (main_cover_start_as_reversed > minor_cover_start_as_reversed)?
									 (main_cover_start_as_reversed + 15): (minor_cover_start_as_reversed + 15);

						int is_small_half_negative = 0 != ((current_vote -> pos[vote_i][vote_j]>current_vote -> pos[i][j]?current_vote -> masks[i][j]:current_vote -> masks[vote_i][vote_j])&IS_NEGATIVE_STRAND); 
						int is_large_half_negative = !is_small_half_negative;

						int is_small_half_on_left_as_reversed = (main_cover_start_as_reversed > minor_cover_start_as_reversed) + (current_vote -> pos[vote_i][vote_j]> current_vote -> pos[i][j]) !=1;
						// small half on left(as reversed)  ===  small half on right (as 'forward' form of the read, i.e., the raw FASTQ form for read_A and reversed FASTQ form for read_B)

						unsigned int small_half_abs_offset = min(current_vote -> pos[i][j], current_vote -> pos[vote_i][vote_j]);
						unsigned int large_half_abs_offset = max(current_vote -> pos[i][j], current_vote -> pos[vote_i][vote_j]);

						// curr_read_text is the 'reversed' form of the read. I.e., the reversed FASTQ form for read_A and the raw FASTQ form for read_B.
						replace_minor = donor_jumped_score(global_context, thread_context, small_half_abs_offset, large_half_abs_offset,
									max(0, guess_start_as_reversed) , min( guess_end_as_reversed, curr_read_len),  curr_read_text,
									curr_read_len, is_small_half_negative, is_large_half_negative, is_small_half_on_left_as_reversed,
									& final_split_point, & is_GT_AG_donors, & is_donor_found_or_annotation, &small_side_increasing_coordinate, &large_side_increasing_coordinate);

						if( 0 && 1018082 == pair_number)
						{
							print_votes(current_vote, global_context -> config.index_prefix);
							SUBREADprintf("JUMP_001018082  NORMAL=%d  SMALL_NEG=%d  LARGE_NEG=%d,  SMALL_ABS=%u  LARGE_ABS=%u,  REPLACE=%d,   INCS=%d %d\n" ,  is_small_half_on_left_as_reversed, is_small_half_negative, is_large_half_negative, small_half_abs_offset, large_half_abs_offset, replace_minor, small_side_increasing_coordinate, large_side_increasing_coordinate);
						}


						// Now "final_split_point" is the read offset on the 'reversed' form of the read. It needs to be changed to (read_len - final_split_point) if the major half is on negative strand. 

						if(replace_minor>0) replace_minor += current_vote -> votes[i][j] * 100000;

					}
					else
					{

						int overlapped ;
						if(current_vote -> coverage_start[vote_i][vote_j] > current_vote -> coverage_start[i][j])
							overlapped = current_vote -> coverage_end[i][j] - current_vote -> coverage_start[vote_i][vote_j];
						else
							overlapped = current_vote -> coverage_end[vote_i][vote_j] - current_vote -> coverage_start[i][j];

						if(0 && FIXLENstrcmp("R000404427", read_name) == 0) 
						{
							SUBREADprintf("OVL=%d, DIST=%u\n", overlapped, (unsigned int)abs(dist));
						}

						if(overlapped > 14) continue;
						if(abs(dist)<6) continue;

						int guess_start = (current_vote -> coverage_start[vote_i][vote_j] > current_vote -> coverage_start[i][j])?
									 (current_vote -> coverage_end[i][j] - 8): (current_vote -> coverage_end[vote_i][vote_j] - 8);

						int guess_end = (current_vote -> coverage_start[vote_i][vote_j] < current_vote -> coverage_start[i][j])?
									 (current_vote -> coverage_start[i][j] + 8): (current_vote -> coverage_start[vote_i][vote_j] + 8);

						if((global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection) && !(current_vote -> masks[vote_i][vote_j] & IS_NEGATIVE_STRAND))
							// if for fusion, the current read must have been reversed.
							// hence, it is now changed to "main half" view.
							reverse_read(curr_read_text, curr_read_len, global_context -> config.space_type);

						int left_indel_offset=0,  right_indel_offset=0;
						int kx2;

						int normally_arranged = 1!=(current_vote -> coverage_start[vote_i][vote_j] > current_vote -> coverage_start[i][j]) + (current_vote -> pos[vote_i][vote_j] > current_vote -> pos[i][j]);

						if(curr_read_len > EXON_LONG_READ_LENGTH){

							int kx1;
							gene_vote_number_t * indel_recorder = current_vote -> indel_recorder[vote_i][vote_j];
							for(kx1=0; kx1<MAX_INDEL_SECTIONS; kx1++)
							{
								if(!indel_recorder[kx1*3]) break;
								major_indels += indel_recorder[kx1*3+2];
							}

							for(kx2=0; kx2<MAX_INDEL_SECTIONS; kx2++)
							{
								if(!current_vote -> indel_recorder[i][j][kx2*3]) break;
								minor_indel_offset += (current_vote -> indel_recorder[i][j][kx2*3+2]);
							}

							if(current_vote -> pos[vote_i][vote_j]  <  current_vote -> pos[i][j])
							{
								left_indel_offset=major_indels;
								right_indel_offset=minor_indel_offset;
							}
							else
							{
								right_indel_offset=major_indels;
								left_indel_offset=minor_indel_offset;

							}

							// the section having a smaller coordinate will have indel_offset !=0
							// the section having a larger coordiname MUST HAVE indel_offset == 0
							right_indel_offset=0;
						}

						if(is_long_del_high_quality( global_context, thread_context, current_vote -> coverage_start[i][j], current_vote -> coverage_end[i][j], current_vote -> coverage_start[vote_i][vote_j], current_vote -> coverage_end[vote_i][vote_j], curr_read_len, current_vote -> votes[i][j], current_vote -> votes[vote_i][vote_j])|| ! global_context -> config.do_long_del_detection)
						replace_minor = donor_score(global_context, thread_context, min(current_vote -> pos[vote_i][vote_j],
									current_vote -> pos[i][j]),max(current_vote -> pos[vote_i][vote_j] ,
									current_vote -> pos[i][j]), left_indel_offset, right_indel_offset, normally_arranged,
									max(0, guess_start), min( guess_end, curr_read_len), curr_read_text, curr_read_len,
									& final_split_point, & is_GT_AG_donors, & is_donor_found_or_annotation, & inserted_bases, &small_side_increasing_coordinate, &large_side_increasing_coordinate, read_name);
						else replace_minor = 0;

						// Now "final_split_point" is the read offset on the 'reversed' form of the read (I.e., the reversed FASTQ form for read_A and the raw FASTQ form for read_B.) if do_fusion_detection AND if the main half is on negative strand.
						// However, because the final_split_point is ALWAYS on the form where the major half can be mapped, final_split_point will never be changed.

						if(replace_minor>0) replace_minor += current_vote -> votes[i][j] * 100000;

						if(0 && ( FIXLENstrcmp("R006232475", read_name) == 0 ) )
							SUBREADprintf("NOJUMP_DONORs=%d   LOC=%u\n", replace_minor , current_vote -> pos[i][j]);
						if((global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection) && !(current_vote -> masks[vote_i][vote_j] & IS_NEGATIVE_STRAND))
							// changed back.
							reverse_read(curr_read_text, curr_read_len, global_context -> config.space_type);
					}
				}

				if(0 && ( FIXLENstrcmp("R006232475", read_name) == 0 ) )
				{
					char posout[100];
					absoffset_to_posstr(global_context, current_vote -> pos[i][j], posout);
					SUBREADprintf("TEST MINOR: POS=%s, REPLACE=%d\n", posout, replace_minor);
				}

				if(replace_minor){// && (replace_minor > current_piece_minor_score)){
					junc_res -> minor_position = current_vote -> pos[i][j];
					junc_res -> minor_votes = current_vote -> votes[i][j];

					junc_res -> minor_coverage_start = current_vote -> coverage_start[i][j];
					junc_res -> minor_coverage_end   = current_vote -> coverage_end  [i][j];

					junc_res -> double_indel_offset = (minor_indel_offset & 0xf)|((major_indels & 0xf)<<4);
					junc_res -> split_point = final_split_point;

					
					if(0 && 1018082 == pair_number)
					{
						SUBREADprintf("REPLACED: LOC %u, INCS=%d %d\n", junc_res -> minor_position, small_side_increasing_coordinate, large_side_increasing_coordinate);
					}

					junc_res -> small_side_increasing_coordinate = small_side_increasing_coordinate;
					junc_res -> large_side_increasing_coordinate = large_side_increasing_coordinate;
					junc_res -> indel_at_junction = inserted_bases;

					align_res -> result_flags &=~0x3;
					if( (!is_donor_found_or_annotation) || is_GT_AG_donors > 2) align_res -> result_flags |= 3;
					else	align_res -> result_flags = is_GT_AG_donors? (align_res -> result_flags|CORE_IS_GT_AG_DONORS):(align_res  -> result_flags &~CORE_IS_GT_AG_DONORS);

					align_res -> result_flags = is_strand_jumpped? (align_res -> result_flags|CORE_IS_STRAND_JUMPED):(align_res -> result_flags &~CORE_IS_STRAND_JUMPED);
				}
			}
		}

		if(0 && memcmp("V0112_0155:7:1101:1173:2204", read_name, 26) == 0) 
		{
			char leftpos[100], rightpos[100];
			absoffset_to_posstr(global_context, current_vote -> pos[vote_i][vote_j]  , leftpos);
			absoffset_to_posstr(global_context, junc_res -> minor_position, rightpos);
			SUBREADprintf("READ=%s, MAJOR=%s, MINOR=%s\n", read_name, leftpos, rightpos);
		}


		// This block runs after the minor half of this anchor is fully determined.
		// If the minor half is a fusion and there is a strand jump, move the minor half coverage to the major half strand.
		if(align_res -> result_flags & CORE_IS_STRAND_JUMPED)
		{
			// If "is_strand_jumped" is true, all coordinates so far are on the best voted strands (must be differnet strands, namely they're very likely to be overlapped). 
			int tmpv = junc_res -> minor_coverage_start;
			junc_res -> minor_coverage_start = curr_read_len - junc_res -> minor_coverage_end;
			junc_res -> minor_coverage_end = curr_read_len - tmpv;

			// Split_point is now the "negative strand read" view. It has to be changed to "main piece" view
			junc_res -> split_point = (align_res -> result_flags & CORE_IS_NEGATIVE_STRAND)?
							junc_res -> split_point :
							(curr_read_len - junc_res -> split_point);
		}
	}
}


void simple_PE_and_same_chro(global_context_t * global_context , simple_mapping_t * r1, simple_mapping_t * r2 , int * is_PE_distance, int * is_same_chromosome , int rlen1, int rlen2){
	test_PE_and_same_chro(global_context, r1 -> mapping_position, r2 -> mapping_position, is_PE_distance, is_same_chromosome, rlen1, rlen2);
}


#define MAX_CLUSTER_ELEMENTS 7

struct cluster_element{
	unsigned int initial_position;
	char cluster_members;
	char from_second_read[MAX_CLUSTER_ELEMENTS];
	int i_list[MAX_CLUSTER_ELEMENTS], j_list[MAX_CLUSTER_ELEMENTS];
};

int add_cluster_member(struct cluster_element * cl , int i, int j, int is_second_read){
	if(cl->cluster_members < MAX_CLUSTER_ELEMENTS){
		cl->i_list[(int)cl->cluster_members] = i;
		cl->j_list[(int)cl->cluster_members] = j;
		cl->from_second_read[(int)cl->cluster_members] = is_second_read;
		cl->cluster_members++;
	}
	return cl->cluster_members;
}

int is_same_cluster( global_context_t * global_context, struct cluster_element * cl , unsigned int pos){
	long long int test_pos = pos;
	test_pos -= cl -> initial_position;
	if(abs(test_pos) < global_context -> config.maximum_intron_length)
		return 1;
	return 0;
}

int process_voting_junction_PE_topK(global_context_t * global_context, thread_context_t * thread_context, subread_read_number_t pair_number, gene_vote_t * vote_1, gene_vote_t * vote_2, char * read_name_1, char * read_name_2, char * read_text_1, char * read_text_2, int read_len_1, int read_len_2, int is_negative_strand, gene_vote_number_t v1_all_subreads, gene_vote_number_t v2_all_subreads);
int align_cluster(global_context_t * global_context, thread_context_t * thread_context, struct cluster_element * this_cluster, char * read_name_1, char * read_name_2, char * read_text_1, char * read_text_2, int read_len_1, int read_len_2, int is_negative_strand,  gene_vote_t * vote_1, gene_vote_t * vote_2, int * this_score, int * ii_path, int * jj_path, int * masks, int * path_len, int * R1R2_mapped);

void simple_copy_vote_to_result( mapping_result_t * align_res, gene_vote_t * current_vote, int vote_i, int vote_j, int used_subreads_in_vote, int noninformative_subreads_in_vote, int score){
	align_res -> selected_position = current_vote -> pos[vote_i][vote_j];
	align_res -> selected_votes = score;
	align_res -> indels_in_confident_coverage = indel_recorder_copy(align_res -> selected_indel_record, current_vote -> indel_recorder[vote_i][vote_j]);
	align_res -> confident_coverage_end = current_vote -> coverage_end[vote_i][vote_j];
	align_res -> confident_coverage_start = current_vote -> coverage_start[vote_i][vote_j];
	align_res -> result_flags = (current_vote -> masks[vote_i][vote_j] & IS_NEGATIVE_STRAND)?(CORE_IS_NEGATIVE_STRAND):0;
	align_res -> used_subreads_in_vote = used_subreads_in_vote;
	align_res -> noninformative_subreads_in_vote = noninformative_subreads_in_vote;
}

int process_voting_junction_PE_juncs( global_context_t * global_context, thread_context_t * thread_context, subread_read_number_t pair_number, gene_vote_t * vote_1, gene_vote_t * vote_2, char * read_name_1, char * read_name_2, char * read_text_1, char * read_text_2, int read_len_1, int read_len_2, int is_negative_strand, gene_vote_number_t v1_all_subreads, gene_vote_number_t v2_all_subreads ){
	int current_cluster_number = 0,max_clusters = global_context -> config.max_vote_simples * 2;
	int i,j, is_second_read, tested_votes, x1;

	struct cluster_element * cluster_buffer = malloc(max_clusters * sizeof(struct cluster_element));
	int max_cluster_size_r1 = 0, max_cluster_size_r2 = 0;

	for( tested_votes = max(vote_1 -> max_vote, vote_2 -> max_vote); tested_votes > 0; tested_votes--) {
		for(is_second_read = 0 ; is_second_read < 1 + global_context -> input_reads.is_paired_end_reads; is_second_read ++) {
			gene_vote_t * current_vote = is_second_read?vote_2:vote_1;
			int * max_cluster_size = is_second_read?(&max_cluster_size_r2):(&max_cluster_size_r1);
			for (i=0; i<GENE_VOTE_TABLE_SIZE; i++) {
				for (j=0; j< current_vote->items[i]; j++) {
					if(current_vote->votes[i][j]!=tested_votes) continue;
					int is_added = 0;

					for(x1 = 0; x1 < current_cluster_number ; x1++){
						if(is_same_cluster(global_context, cluster_buffer+x1, current_vote->pos[i][j])){
							int new_size =add_cluster_member(cluster_buffer+x1, i, j, is_second_read);
							(*max_cluster_size) = max(*max_cluster_size, new_size);
							is_added = 1;
						}
						if(is_added)break;
					}
					if(current_cluster_number < max_clusters && !is_added){
						cluster_buffer[current_cluster_number].initial_position = current_vote->pos[i][j];
						cluster_buffer[current_cluster_number].i_list[0] = i;
						cluster_buffer[current_cluster_number].j_list[0] = j;
						cluster_buffer[current_cluster_number].from_second_read[0] = is_second_read;
						cluster_buffer[current_cluster_number].cluster_members = 1;
						current_cluster_number++;
					}
				}
			}
		}
	}
	
	if(1 || max_cluster_size_r1 == 3 || max_cluster_size_r2 == 3 ) // if there are 3-section clusters then parse it, else go to the regular procedure. There is a upper-limit to the sections to avoid fragile mapping.
	{
		for(x1 = 0 ; x1 < current_cluster_number ; x1++){
			int this_score = -1, path_len = -1, R1R2_mapped = 0;
			int this_ii_path[ MAX_CLUSTER_ELEMENTS ], this_jj_path[ MAX_CLUSTER_ELEMENTS ], this_masks [ MAX_CLUSTER_ELEMENTS ];
			align_cluster(global_context, thread_context, cluster_buffer + x1, read_name_1, read_name_2, read_text_1, read_text_2, read_len_1, read_len_2, is_negative_strand, vote_1, vote_2, &this_score, this_ii_path, this_jj_path, this_masks, &path_len, &R1R2_mapped);

			if(0 && FIXLENstrcmp("R00000003493", read_name_1)==0)
				SUBREADprintf("REAE_TEST : R12MAP=%d, PATHLEN=%d, SCORE=%d\n", R1R2_mapped, path_len, this_score);

			if(this_score > 0){
				if(( R1R2_mapped & CLUSTER_ALIGNMENT_DONOR_R1_MAPPED) && ( R1R2_mapped & CLUSTER_ALIGNMENT_DONOR_R2_MAPPED)){
					for(i = 0; i < global_context -> config.multi_best_reads; i++) {
						mapping_result_t * old_result_R1 = _global_retrieve_alignment_ptr(global_context, pair_number, 0, i);
						mapping_result_t * old_result_R2 = _global_retrieve_alignment_ptr(global_context, pair_number, 1, i);
						short old_score_R1 = old_result_R1 -> selected_votes;
						short old_score_R2 = old_result_R2 -> selected_votes;

						if( old_score_R1 < this_score || old_score_R2 < this_score ){

							for(j = global_context -> config.multi_best_reads - 2; j>=i; j--){
								mapping_result_t * shifted_result_R1 = _global_retrieve_alignment_ptr(global_context, pair_number, 0, j);
								mapping_result_t * shifted_result_R2 = _global_retrieve_alignment_ptr(global_context, pair_number, 1, j);
								mapping_result_t * target_result_R1 = _global_retrieve_alignment_ptr(global_context, pair_number, 0, j+1);
								mapping_result_t * target_result_R2 = _global_retrieve_alignment_ptr(global_context, pair_number, 1, j+1);
								memcpy( target_result_R1, shifted_result_R1 , sizeof(mapping_result_t));
								memcpy( target_result_R2, shifted_result_R2 , sizeof(mapping_result_t) );

							}

							int  best_R1_i = -1, best_R1_j = - 1 , highest_vote_R1 = -1, highest_vote_R2 = -1,  best_R2_i = -2, best_R2_j = -2;

							for(j = 0; j < path_len ; j++){
								if( this_masks[j] & CLUSTER_ALIGNMENT_DONOR_R1_MAPPED ){
									if( highest_vote_R1 < vote_1 -> votes [  this_ii_path [j] ] [  this_jj_path [j] ] ){
										best_R1_i = this_ii_path [j] ;
										best_R1_j = this_jj_path [j] ;
										highest_vote_R1 =  vote_1 -> votes [  this_ii_path [j] ] [  this_jj_path [j] ] ;
									}
								}else{
									if( highest_vote_R2 < vote_2 -> votes [  this_ii_path [j] ] [  this_jj_path [j] ] ){
										best_R2_i = this_ii_path [j] ;
										best_R2_j = this_jj_path [j] ;
										highest_vote_R2 =  vote_2 -> votes [  this_ii_path [j] ] [  this_jj_path [j] ] ;
									}
								}
								//SUBREADprintf("MASK=%d\n",  this_masks[j]);
							}

							//SUBREADprintf("IJ: R1=%d,%d  R2=%d,%d  MASK=%d\n", best_R1_i,best_R1_j,best_R2_i,best_R2_j);

							simple_copy_vote_to_result( old_result_R1, vote_1, best_R1_i, best_R1_j , v1_all_subreads, vote_1 -> noninformative_subreads, this_score);
							simple_copy_vote_to_result( old_result_R2, vote_2, best_R2_i, best_R2_j , v2_all_subreads, vote_2 -> noninformative_subreads, this_score);
							break;
						}
					}
				} else if(  R1R2_mapped & ( CLUSTER_ALIGNMENT_DONOR_R1_MAPPED | CLUSTER_ALIGNMENT_DONOR_R2_MAPPED ) ) {
					int is_R2_mapped = ( R1R2_mapped & CLUSTER_ALIGNMENT_DONOR_R2_MAPPED)?1:0;
					for(i = 0; i < global_context -> config.multi_best_reads; i++) {
						mapping_result_t * old_result_R = _global_retrieve_alignment_ptr(global_context, pair_number, is_R2_mapped, i);
						short old_score_R = old_result_R -> selected_votes;



						if( old_score_R < this_score ){

							for(j = global_context -> config.multi_best_reads - 2; j>=i; j--){
								mapping_result_t * shifted_result_R = _global_retrieve_alignment_ptr(global_context, pair_number, is_R2_mapped, j);
								mapping_result_t * target_result_R = _global_retrieve_alignment_ptr(global_context, pair_number, is_R2_mapped, j+1);
								memcpy( target_result_R, shifted_result_R , sizeof(mapping_result_t));

							}

							int  best_R_i = -1, best_R_j = - 1 , highest_vote_R = -1;
							gene_vote_t * this_vote = is_R2_mapped?vote_2:vote_1;

							for(j = 0; j < path_len ; j++){
								if( highest_vote_R < this_vote -> votes [  this_ii_path [j] ] [  this_jj_path [j] ] ){
									best_R_i = this_ii_path [j] ;
									best_R_j = this_jj_path [j] ;
									highest_vote_R =  this_vote -> votes [  this_ii_path [j] ] [  this_jj_path [j] ] ;
								}
								//SUBREADprintf("MASK=%d\n",  this_masks[j]);
							}

							//SUBREADprintf("IJ: R1=%d,%d  R2=%d,%d  MASK=%d\n", best_R1_i,best_R1_j,best_R2_i,best_R2_j);

							simple_copy_vote_to_result( old_result_R, this_vote, best_R_i, best_R_j , v1_all_subreads, this_vote -> noninformative_subreads, this_score);
							break;
						}


					}
				}
			}
			
			/*if(highest_score > 0){
				if(this_score >0){
					if(this_score > highest_score){
						highest_score = this_score;
						highest_occurance = 1;
						memcpy(best_ii_path, this_ii_path, sizeof(int)*path_len);
						memcpy(best_jj_path, this_jj_path, sizeof(int)*path_len);
						best_path_len = path_len;
					}else if(this_score == highest_score)
						highest_occurance ++;
				}
			}*/
		}


		// call new junctions from the path
		// then put the alignment into the best list.

	}else{
		return process_voting_junction_PE_topK(global_context, thread_context, pair_number, vote_1, vote_2, read_name_1, read_name_2, read_text_1, read_text_2, read_len_1, read_len_2, is_negative_strand, v1_all_subreads, v2_all_subreads);
	}

	free(cluster_buffer);
	return 0;
}


int compare_cluster_elements (void * arr, int l, int r){
	int * ii_array = ((void **)arr)[0];
	int * jj_array = ((void **)arr)[1];
	int * second_vote = ((void **)arr)[2];

	if(second_vote[l] != second_vote[r])
		return second_vote[l] - second_vote[r];

	gene_vote_t * vote_1 = ((void **)arr)[3];
	gene_vote_t * vote_2 = ((void **)arr)[4];
	

	gene_vote_t * this_vote_L = second_vote[l]?vote_2:vote_1;
	gene_vote_t * this_vote_R = second_vote[r]?vote_2:vote_1;

	return this_vote_L->coverage_start[ii_array[l]][jj_array[l]] - this_vote_R -> coverage_start[ii_array[r]][jj_array[r]];
}

void exchange_cluster_elements (void * arr, int l, int r){
	int * ii_array = ((void **)arr)[0];
	int * jj_array = ((void **)arr)[1];
	int * second_vote = ((void **)arr)[2];

	int ti;
	ti = ii_array[l];
	ii_array[l] = ii_array[r];
	ii_array[r]=ti;

	ti = jj_array[l];
	jj_array[l] = jj_array[r];
	jj_array[r]=ti;

	ti = second_vote[l];
	second_vote[l] = second_vote[r];
	second_vote[r]=ti;
}

#define NEW_EXTEND_SCAN_INTRON_LONGEST 6000
#define NEW_EXTEND_SCAN_EXON_SHORTEST 14

int find_path(global_context_t * global_context, thread_context_t * thread_context, int start_element_i, int target_element_i, int * ii_array, int * jj_array, int * is_second_vote_array,  gene_vote_t * vote_1, gene_vote_t * vote_2, char * read_name_1, char * read_name_2, char * read_text_1, char * read_text_2, int read_len_1, int read_len_2, int is_negative_strand, int * this_mask , int * exon_last_base);
int find_donor_receptor(global_context_t * global_context, thread_context_t * thread_context, char * rname, char * rtext, int rlen, int start_coverage, int end_coverage, unsigned int start_pos, unsigned int end_pos, int indels_in_start, int v1, int v2, int * misma_bases, int * matched_bases, int * is_negative_donor);
int extend_uncovered_region_juncs(global_context_t * global_context, thread_context_t * thread_context, char * rname, char * rtext, int rlen, int scan_to_tail, unsigned int scan_chro_start, int scan_read_start, unsigned short expect_donor, int * mismatched_bases_after_start, int * first_exon_last_base, unsigned int * first_exon_first_base, int * ret_mismatched_bases, int * is_negative_donor){

	if(  scan_to_tail  ) assert( scan_read_start < rlen - NEW_EXTEND_SCAN_EXON_SHORTEST );
	else	assert( scan_read_start >= NEW_EXTEND_SCAN_EXON_SHORTEST);

	gene_value_index_t * value_index = thread_context?thread_context->current_value_index:global_context->current_value_index;
	int x1, best_best_score = -1, best_best_occurance = 0;

	unsigned long long matching_target = 0, rolling_bases = 0;

	for(x1 = 0 ; x1 < 8 ; x1++){
		int nch = scan_to_tail? rtext[ rlen - 2 - x1 ] :  rtext[ 10 - x1 ] ;
		matching_target = ( matching_target << 8 ) | nch;
	}
		if(0 && FIXLENstrcmp("V0112_0155:7:1101:13762:2349#ACTTGA", rname ) == 0 )
			SUBREADprintf("TAG=%016llX\n",matching_target);

	for(x1 = 0; x1 < NEW_EXTEND_SCAN_INTRON_LONGEST ; x1++){
		int best_last_exon_base = -1, matched_in_the_uncovered_gap = -1, mismatched_bases = -1, extended_should_mismatch = -1;
		unsigned int scan_cursor = scan_chro_start ;
		if(scan_to_tail) scan_cursor+=x1;else scan_cursor-=x1;
		unsigned long long nch = gvindex_get( value_index, scan_cursor );
		if(scan_to_tail)
			rolling_bases = (rolling_bases >>  8) | nch << 56;
		else
			rolling_bases = nch | ( rolling_bases << 8 );

		//SUBREADprintf("MATCH:%016llX,%016llX\n", rolling_bases, matching_target);
		if(rolling_bases == matching_target){
			//SUBREADprintf("PNTT-M\n");
			if(scan_to_tail) {
				best_last_exon_base = find_donor_receptor(global_context, thread_context, rname, rtext, rlen, scan_read_start, rlen - 2 - 7,  scan_chro_start, scan_cursor - (rlen - 2) , 0, 0,0, &mismatched_bases , &matched_in_the_uncovered_gap, is_negative_donor);
				if(best_last_exon_base>0)
					extended_should_mismatch = match_chro( rtext + best_last_exon_base , value_index, scan_chro_start + best_last_exon_base, rlen - best_last_exon_base, 0, global_context->config.space_type);
			} else {
				best_last_exon_base = find_donor_receptor(global_context, thread_context, rname, rtext, rlen, 10, scan_read_start, scan_cursor - 3 , scan_chro_start, 0, 0,0, &mismatched_bases , &matched_in_the_uncovered_gap,is_negative_donor);
				if(best_last_exon_base>0)
					extended_should_mismatch = match_chro( rtext, value_index, scan_chro_start, best_last_exon_base , 0, global_context->config.space_type);
			}

		}

		if(best_last_exon_base >0 && extended_should_mismatch < ( scan_to_tail?( rlen - best_last_exon_base - 4  ):(  best_last_exon_base - 4  )) && mismatched_bases < 2 ){
			int this_score;
			if(scan_to_tail) this_score = rlen - scan_read_start - mismatched_bases;
			else	this_score = scan_read_start - mismatched_bases;
			if(best_best_score < this_score){
				best_best_score = this_score;
				(*mismatched_bases_after_start) = mismatched_bases;
				(*first_exon_last_base) = best_last_exon_base;
				(*first_exon_first_base) = scan_to_tail?( scan_cursor - (rlen - 2) ) : ( scan_cursor - 3 );
				(*ret_mismatched_bases) = mismatched_bases;
				best_best_occurance = 1;

			}else if( best_best_score == this_score ) best_best_occurance ++;
		}


		if(0 && (!scan_to_tail) && best_last_exon_base >0 && extended_should_mismatch < best_last_exon_base - 4 && mismatched_bases < 2){
			char out1pos[100], out2pos[100];
			absoffset_to_posstr(global_context, scan_chro_start + best_last_exon_base, out1pos);
			absoffset_to_posstr(global_context, scan_cursor - (rlen - 2) + best_last_exon_base, out2pos);
			SUBREADprintf("LIMMISMA: %d < %d - 4\t\tfor %s\n" , extended_should_mismatch,   best_last_exon_base ,rname);

			SUBREADprintf("HEAD MATCH: %s - %s : MM=%d ; SPLIT=%d\t%s\n",out1pos, out2pos, mismatched_bases, best_last_exon_base, rname);

			SUBREADprintf("R =%s\nS1=", rtext);
			int x2;
			for(x2 = 0; x2 <  rlen; x2++){
				if(x2 > best_last_exon_base + 16) SUBREADprintf(" ");
				else{
					int nch =  gvindex_get( value_index, scan_cursor - 3 + x2 );
					SUBREADprintf("%c", nch);
				}
			}
			SUBREADprintf("\nS2=");

			for(x2 = 0; x2 <  rlen; x2++){
				if(x2 > best_last_exon_base + 16 ) SUBREADprintf(" ");
				else{
					int nch =  gvindex_get( value_index, scan_chro_start +x2);
					SUBREADprintf("%c", nch);
				}
			}
			SUBREADprintf("\n   ");

			for(x2 = 0; x2 <  rlen; x2++){
				if(x2 < best_last_exon_base ) SUBREADprintf(" ");
				else if( x2 > best_last_exon_base + 1 ) SUBREADprintf(" ");
				else SUBREADprintf("|");
			}
			SUBREADprintf("\n\n");
		}
		if(0 && scan_to_tail && best_last_exon_base >0 && extended_should_mismatch < rlen - best_last_exon_base - 4 && mismatched_bases < 2){

			SUBREADprintf("LIMMISMA: %d < %d - 4\t\tfor %s\n" , extended_should_mismatch,  (rlen - best_last_exon_base ),rname);
			char out1pos[100], out2pos[100];
			absoffset_to_posstr(global_context, scan_chro_start + best_last_exon_base, out1pos);
			absoffset_to_posstr(global_context, scan_cursor - (rlen - 2) + best_last_exon_base, out2pos);
			SUBREADprintf("TAIL MATCH: %s - %s : MM=%d ; SPLIT=%d\t%s\n",out1pos, out2pos, mismatched_bases, best_last_exon_base, rname);

			SUBREADprintf("R =%s\nS1=", rtext);
			int x2;
			for(x2 = 0; x2 <  rlen; x2++){
				if(x2 < scan_read_start - 16) SUBREADprintf(" ");
				else{
					int nch =  gvindex_get( value_index, x2 + scan_chro_start);
					SUBREADprintf("%c", nch);
				}
			}
			SUBREADprintf("\nS2=");

			for(x2 = 0; x2 <  rlen; x2++){
				if(x2 < best_last_exon_base - 16 ) SUBREADprintf(" ");
				else{
					int nch =  gvindex_get( value_index, scan_cursor - (rlen - 2)+x2);
					SUBREADprintf("%c", nch);
				}
			}
			SUBREADprintf("\n   ");

			for(x2 = 0; x2 <  rlen; x2++){
				if(x2 < best_last_exon_base ) SUBREADprintf(" ");
				else if( x2 > best_last_exon_base + 1 ) SUBREADprintf(" ");
				else SUBREADprintf("|");
			}
			SUBREADprintf("\n\n");
		}
	}
	if(0&&best_best_occurance>0 && best_best_score>0)
		SUBREADprintf("OCCR=%d : SCR=%d\n", best_best_occurance, best_best_score);
	if(best_best_occurance == 1) return best_best_score;
	return -1;
}

void simple_add_junction( global_context_t * global_context, thread_context_t * thread_context, unsigned int left_edge_wanted, unsigned int right_edge_wanted, int indel_at_junction, int is_negative_donors ){
	char * chro_name_left, *chro_name_right;
	int chro_pos_left,chro_pos_right;

	locate_gene_position( left_edge_wanted , &global_context -> chromosome_table, &chro_name_left, &chro_pos_left);
	locate_gene_position( right_edge_wanted , &global_context -> chromosome_table, &chro_name_right, &chro_pos_right);
	if((!( global_context-> config.do_fusion_detection ||  global_context-> config.do_long_del_detection) ) && chro_name_right!=chro_name_left) return;

	//insert event
	HashTable * event_table = NULL;
	chromosome_event_t * event_space = NULL;
	if(thread_context)
	{
			event_table = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_entry_table;
			event_space = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}
	else
	{
			event_table = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_entry_table;
			event_space = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}
	chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];
	chromosome_event_t * found = NULL;
	int found_events = search_event(global_context, event_table, event_space, left_edge_wanted , EVENT_SEARCH_BY_SMALL_SIDE,  CHRO_EVENT_TYPE_INDEL | CHRO_EVENT_TYPE_JUNCTION | CHRO_EVENT_TYPE_FUSION, search_return);

	if(found_events)
	{
			int kx1;
			for(kx1 = 0; kx1 < found_events ; kx1++)
			{
					if(search_return[kx1] -> event_large_side == right_edge_wanted)
					{
							found = search_return[kx1];
							break;
					}
			}
	}
	
	if(found) found -> supporting_reads ++;
	else
	{
			int event_no;


			if(thread_context)
				event_no = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> total_events ++;
			else
				event_no = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) ->  total_events ++;


			chromosome_event_t * new_event = event_space+event_no;
			memset(new_event,0,sizeof(chromosome_event_t));
			new_event -> event_small_side = left_edge_wanted;
			new_event -> event_large_side = right_edge_wanted;
			new_event -> is_negative_strand= is_negative_donors;
			new_event -> event_type = CHRO_EVENT_TYPE_JUNCTION;
			new_event -> supporting_reads = 1;
			new_event -> indel_length = 0;
			new_event -> indel_at_junction = indel_at_junction;
			new_event -> is_donor_found_or_annotation = 1;
			new_event -> small_side_increasing_coordinate = 0;
			new_event -> large_side_increasing_coordinate = 1;
			put_new_event(event_table, new_event , event_no);
	}
}

int align_cluster(global_context_t * global_context, thread_context_t * thread_context, struct cluster_element * this_cluster, char * read_name_1, char * read_name_2, char * read_text_1, char * read_text_2, int read_len_1, int read_len_2, int is_negative_strand,  gene_vote_t * vote_1, gene_vote_t * vote_2, int * this_score, int * best_ii_path, int * best_jj_path, int * best_masks, int * best_path_length, int * R1R2_mapped){
	//int cluster_x1;

	//SUBREADprintf("\n === Cluster %s    %s  === \n", is_negative_strand?"NEG":"POS", read_name_1);
	//unsigned int min_frag_start = 0xffffffff;

	int ii_array[MAX_CLUSTER_ELEMENTS], jj_array[MAX_CLUSTER_ELEMENTS], is_second_vote_array[MAX_CLUSTER_ELEMENTS], dynamic_highest_mask[MAX_CLUSTER_ELEMENTS], x1;
	void * sort_pointers [5];

	for(x1 = 0 ; x1 < this_cluster->cluster_members; x1++){
		ii_array[x1] = this_cluster -> i_list[x1];
		jj_array[x1] = this_cluster -> j_list[x1];
		is_second_vote_array[x1] = this_cluster -> from_second_read[x1];
		
	}

	sort_pointers[0] = ii_array;
	sort_pointers[1] = jj_array;
	sort_pointers[2] = is_second_vote_array;
	sort_pointers[3] = vote_1;
	sort_pointers[4] = vote_2;

	basic_sort(sort_pointers, this_cluster->cluster_members, compare_cluster_elements, exchange_cluster_elements);

	if(0)
	for(x1 = 0 ; x1 < this_cluster->cluster_members; x1++){
		gene_vote_t * this_vote = is_second_vote_array[x1]?vote_2:vote_1;
		int ii = ii_array[x1];
		int jj = jj_array[x1];

		SUBREADprintf("   R%d %d - %d   POS=%u  VOTES=%d\n", 1+is_second_vote_array[x1], this_vote->coverage_start[ii][jj], this_vote->coverage_end[ii][jj], this_vote -> pos[ii][jj], this_vote->votes[ii][jj]);
	}

	int dynamic_highest_scores[MAX_CLUSTER_ELEMENTS], dynamic_last_exon[MAX_CLUSTER_ELEMENTS];
	char dynamic_highest_path[MAX_CLUSTER_ELEMENTS];
	memset(dynamic_highest_scores,0,sizeof(int)*MAX_CLUSTER_ELEMENTS);

	int target_element_i;

	for(target_element_i = 0; target_element_i < this_cluster->cluster_members; target_element_i++){
		gene_vote_t * this_vote = is_second_vote_array[target_element_i]?vote_2:vote_1;
		int ii = ii_array[target_element_i];
		int jj = jj_array[target_element_i];
		dynamic_highest_scores[target_element_i] =  this_vote->coverage_end[ii][jj] - this_vote->coverage_start[ii][jj];
		dynamic_highest_path[ target_element_i ] = -1;
	}

	int highest_score = -1;
	int highest_target_end = -1;
	for(target_element_i = 0; target_element_i < this_cluster->cluster_members; target_element_i++){
		int start_element_i;
		for(start_element_i = 0; start_element_i < this_cluster->cluster_members; start_element_i++){
			if(target_element_i <= start_element_i) continue;
			int this_mask = -1, breakpount_last_exon = -1;
			int increasing_score = find_path(global_context, thread_context, start_element_i, target_element_i, ii_array, jj_array, is_second_vote_array, vote_1, vote_2, read_name_1, read_name_2, read_text_1, read_text_2, read_len_1, read_len_2, is_negative_strand, &this_mask, &breakpount_last_exon);
			if(increasing_score >= 0 && increasing_score + dynamic_highest_scores[start_element_i] > dynamic_highest_scores[target_element_i]){
				dynamic_highest_path[ target_element_i ] = start_element_i;
				dynamic_highest_scores[target_element_i] = increasing_score + dynamic_highest_scores[start_element_i];
				dynamic_highest_mask[ target_element_i ] = this_mask;
				dynamic_last_exon[ target_element_i ] = breakpount_last_exon;
				if(  dynamic_highest_scores[target_element_i] > highest_score ){
					highest_score =  dynamic_highest_scores[target_element_i] ;
					highest_target_end = target_element_i;
				}
			}
		}
	}


	if(highest_target_end >=0 && highest_score > 160 - 159){
		int is_on_path [MAX_CLUSTER_ELEMENTS];
		memset(is_on_path,0,sizeof(int)*MAX_CLUSTER_ELEMENTS);

		gene_vote_t * last_vote = is_second_vote_array[ highest_target_end ]?vote_2:vote_1;
		int last_section_read_end = last_vote -> coverage_end[ ii_array [highest_target_end]  ] [ jj_array [highest_target_end]   ];
		int this_rlen = is_second_vote_array[ highest_target_end ]?read_len_2 : read_len_1;
		int this_votes = last_vote -> votes [ ii_array [highest_target_end]  ] [ jj_array [highest_target_end]   ];
		int tail_first_exon_last_base_in_read=-1, tail_mismatched_bases=-1;
		unsigned int tail_first_exon_first_base_on_chro, tail_mapped_section_pos;
		int front_first_exon_last_base_in_read=-1, front_mismatched_bases=-1;
		unsigned int front_first_exon_first_base_on_chro, front_mapped_section_pos ;
		int front_score = 0, tail_score = 0, front_negative_donor = 0, tail_negative_donor = 0;

		if(0 && last_section_read_end < this_rlen - NEW_EXTEND_SCAN_EXON_SHORTEST && this_votes > 4)
		{
			char * this_rname = is_second_vote_array[ highest_target_end ]?read_name_2:read_name_1;
			char * this_rtext = is_second_vote_array[ highest_target_end ]?read_text_2:read_text_1;
			int scan_to_tail = 1, mismatched_bases_after_start;
			tail_mapped_section_pos  =	last_vote -> pos[ ii_array [highest_target_end]  ] [ jj_array [highest_target_end]   ] + 
						last_vote -> current_indel_cursor[ ii_array [highest_target_end]  ] [ jj_array [highest_target_end]   ] ;
			if(0){
				char out1pos[100];
				absoffset_to_posstr(global_context, tail_mapped_section_pos, out1pos);
				SUBREADprintf("RN=%s\nSTART=%u, READ_START=%d, READ_FIRTS_BASE_POS=%s\n", this_rname, tail_mapped_section_pos , last_section_read_end, out1pos);
			}

			tail_score = extend_uncovered_region_juncs(global_context, thread_context, this_rname, this_rtext , this_rlen, scan_to_tail, tail_mapped_section_pos, last_section_read_end , -1, & mismatched_bases_after_start, &tail_first_exon_last_base_in_read, &tail_first_exon_first_base_on_chro, &tail_mismatched_bases, &tail_negative_donor);
		}
		(*best_path_length) = 0;
		if(highest_score>0) while(1){
			best_ii_path[(*best_path_length)] = ii_array[highest_target_end];
			best_jj_path[(*best_path_length)] = jj_array[highest_target_end];
			best_masks[ (*best_path_length) ] = dynamic_highest_mask[highest_target_end];

			if( dynamic_last_exon [ highest_target_end ] > 0 ) best_masks[ (*best_path_length) ] |= ( is_second_vote_array[ highest_target_end ]?CLUSTER_ALIGNMENT_DONOR_R2_MAPPED:CLUSTER_ALIGNMENT_DONOR_R1_MAPPED);

			if(  is_second_vote_array[ highest_target_end ] ) (*R1R2_mapped) |= CLUSTER_ALIGNMENT_DONOR_R2_MAPPED;
			else  (*R1R2_mapped) |= CLUSTER_ALIGNMENT_DONOR_R1_MAPPED;
			
			(*best_path_length)++;

			is_on_path[highest_target_end] = 1;
			if(  dynamic_highest_path[highest_target_end] == -1 ) break;
			highest_target_end = dynamic_highest_path[highest_target_end];		
		}

		gene_vote_t * first_vote = is_second_vote_array[ highest_target_end ]?vote_2:vote_1;
		int first_section_read_start = first_vote -> coverage_start [ ii_array [highest_target_end]  ] [ jj_array [highest_target_end] ] ;
		this_votes = first_vote  -> votes [ ii_array [highest_target_end]  ] [ jj_array [highest_target_end]   ];

		if(0 && first_section_read_start > NEW_EXTEND_SCAN_EXON_SHORTEST && this_votes > 4){
			char * this_rname = is_second_vote_array[ highest_target_end ]?read_name_2:read_name_1;
			char * this_rtext = is_second_vote_array[ highest_target_end ]?read_text_2:read_text_1;
			int scan_to_tail = 0, mismatched_bases_after_start;
			front_mapped_section_pos =  first_vote -> pos[ ii_array [highest_target_end]  ] [ jj_array [highest_target_end]   ];

			front_score = extend_uncovered_region_juncs(global_context, thread_context, this_rname, this_rtext , this_rlen, scan_to_tail, front_mapped_section_pos, first_section_read_start , -1, & mismatched_bases_after_start, &front_first_exon_last_base_in_read, &front_first_exon_first_base_on_chro, &front_mismatched_bases, &front_negative_donor);

		}
		
		if(0 && front_score>0 && tail_score>0){

			SUBREADprintf("\n>>> %s\n", read_name_1);

			for(target_element_i = 0; target_element_i < this_cluster->cluster_members; target_element_i++){
				SUBREADprintf("R%d\t", is_second_vote_array[target_element_i]+1);
			}
			SUBREADprintf("\n");

			for(target_element_i = 0; target_element_i < this_cluster->cluster_members; target_element_i++){
				gene_vote_t * this_vote = is_second_vote_array[target_element_i]?vote_2:vote_1;
				int ii = ii_array[target_element_i];
				int jj = jj_array[target_element_i];

				SUBREADprintf("%d%c%d\t", this_vote->coverage_start[ii][jj], is_on_path[target_element_i]?'=':'-', this_vote->coverage_end[ii][jj]);
			}
			SUBREADprintf("\n");

			for(target_element_i = 0; target_element_i < this_cluster->cluster_members; target_element_i++){
				SUBREADprintf("%d\t", dynamic_highest_scores[target_element_i]);
			}
			SUBREADprintf("\n");
			SUBREADprintf("Extra_scores = %d, %d\n", front_score, tail_score);
		}

		(*this_score) = highest_score + max(0, front_score) + max(0, tail_score);
		int applied_score_cut=0;
		if(((*R1R2_mapped) & CLUSTER_ALIGNMENT_DONOR_R1_MAPPED)&&( (*R1R2_mapped) & CLUSTER_ALIGNMENT_DONOR_R2_MAPPED ) )
			applied_score_cut = read_len_2 + read_len_1 - 70;
		else if((*R1R2_mapped) & CLUSTER_ALIGNMENT_DONOR_R1_MAPPED)
			applied_score_cut = read_len_1 - 30;
		else if((*R1R2_mapped) & CLUSTER_ALIGNMENT_DONOR_R1_MAPPED)
			applied_score_cut = read_len_2 - 30;

		if( (*this_score) >= applied_score_cut){
			for( x1 = 0;  x1 < MAX_CLUSTER_ELEMENTS; x1++){
				if(!is_on_path[x1]) continue;

				int x2, second_end = -1;
				for(x2 = x1 + 1; x2 < MAX_CLUSTER_ELEMENTS; x2++){
					if(is_on_path[x2]){
						second_end = x2;
						break;
					}
				}


				if(second_end > 0){
					if( dynamic_last_exon[second_end] >0 ){
						gene_vote_t * this_vote = is_second_vote_array[ x1 ]?vote_2:vote_1;
						unsigned int junction_small_side = this_vote -> pos[ ii_array[ x1 ]][ jj_array[ x1 ]] + 
											this_vote ->  current_indel_cursor[ ii_array[ x1 ]][ jj_array[ x1 ]] + dynamic_last_exon[second_end];

						unsigned int junction_large_side = this_vote -> pos[ ii_array[ second_end ]][ jj_array[ second_end ]] + dynamic_last_exon[second_end] + 1;

						if(0){
							char out1pos[100], out2pos[100];
							absoffset_to_posstr(global_context, junction_small_side, out1pos);
							absoffset_to_posstr(global_context, junction_large_side, out2pos);
							SUBREADprintf("CLUSTER_JUNCTION %s %s\n%s\n%s\n\n", out1pos, out2pos, read_text_1, read_text_2);
						}

						//#warning "SUBREAD_151 ============= MAKE SURE:  CHANGE '0' TO INSERTED BASES ================="
						simple_add_junction(global_context, thread_context, junction_small_side, junction_large_side, 0, (dynamic_highest_mask[ second_end ] & CLUSTER_ALIGNMENT_DONOR_NEGATIVE_STRAND)?1:0);
					}
				}
			}



			if(0 && front_mismatched_bases <1 && front_score >14){
				unsigned int junction_small_side = front_first_exon_first_base_on_chro + front_first_exon_last_base_in_read;
				unsigned int junction_large_side = front_mapped_section_pos + front_first_exon_last_base_in_read + 1;

				char out1pos[100], out2pos[100];
				absoffset_to_posstr(global_context, junction_small_side+1, out1pos);
				absoffset_to_posstr(global_context, junction_large_side+1, out2pos);
				//SUBREADprintf("FMB=%d\tFS=%d\nPOS=%s , %s\n\n", front_mismatched_bases, front_score, out1pos, out2pos);

				simple_add_junction(global_context, thread_context, junction_small_side, junction_large_side, 0, front_negative_donor);
			}

			if(0 && tail_mismatched_bases <1 && tail_score >14){
				unsigned int junction_small_side = tail_mapped_section_pos + tail_first_exon_last_base_in_read;
				unsigned int junction_large_side = tail_first_exon_first_base_on_chro + tail_first_exon_last_base_in_read; 

				char out1pos[100], out2pos[100];
				absoffset_to_posstr(global_context, junction_small_side+1, out1pos);
				absoffset_to_posstr(global_context, junction_large_side+1, out2pos);
				//SUBREADprintf("BMB=%d\tBS=%d\nPOS=%s , %s\n\n", tail_mismatched_bases, tail_score, out1pos, out2pos);



				simple_add_junction(global_context, thread_context, junction_small_side, junction_large_side, 0, tail_negative_donor);
			}
		}
	}
	return 0;
}

#define paired_donor_receptor_m2(s, c1, c2 ) ( s[0] == c1 && s[1] == c2 )

int is_paired_donor_receptor( char * small_bases, char * large_bases ){

	//SUBREADprintf("SITE1 = %c%c , SITE2 = %c%c\n", small_bases[0], small_bases[1], large_bases[0], large_bases[1]);
	//
	if ( paired_donor_receptor_m2( small_bases, 'G', 'T' ) && 
		 paired_donor_receptor_m2( large_bases, 'A', 'G' ) )
		return 1;

	if ( paired_donor_receptor_m2( small_bases, 'C', 'T' ) && 
		 paired_donor_receptor_m2( large_bases, 'A', 'C' ) )
		return 2;

	// http://www.ncbi.nlm.nih.gov/pmc/articles/PMC113136/
	//  the 99.24% of splice site pairs should be GT-AG,
	//  0.69% GC-AG,
	//  0.05% AT-AC
	//  and finally only 0.02% could consist of other types of non-canonical splice sites.

	// non-canonical : GC-AG (+) or CT-GC (-)
	if ( paired_donor_receptor_m2( small_bases, 'G', 'C' ) && 
		 paired_donor_receptor_m2( large_bases, 'A', 'G' ) )
		return 3;

	if ( paired_donor_receptor_m2( small_bases, 'C', 'T' ) && 
		 paired_donor_receptor_m2( large_bases, 'G', 'C' ) )
		return 4;
	

	// non-canonical : AT-AC (+) or GT-AT (-)
	if ( paired_donor_receptor_m2( small_bases, 'A', 'T' ) && 
		 paired_donor_receptor_m2( large_bases, 'A', 'C' ) )
		return 5;

	if ( paired_donor_receptor_m2( small_bases, 'G', 'T' ) && 
		 paired_donor_receptor_m2( large_bases, 'A', 'T' ) )
		return 6;
	

	return 0;
}

int find_donor_receptor(global_context_t * global_context, thread_context_t * thread_context, char * rname, char * rtext, int rlen, int start_coverage, int end_coverage, unsigned int start_pos, unsigned int end_pos, int indels_in_start, int v1, int v2, int * misma_bases, int * matched_bases, int * is_negative_donor){

	gene_value_index_t * value_index = thread_context?thread_context->current_value_index:global_context->current_value_index;
	int search_in_read_start = start_coverage - 8, search_in_read_end = end_coverage + 8;
	search_in_read_start = max(0, search_in_read_start);
	search_in_read_end = min( rlen, search_in_read_end );
	unsigned int search_in_chro_start = start_pos + indels_in_start + search_in_read_start;

	char chro_bases_startside[ search_in_read_end - search_in_read_start ], chro_bases_endside[search_in_read_end - search_in_read_start];

	int x1;

	for(x1 = 0; x1 < search_in_read_end - search_in_read_start; x1++){
		chro_bases_startside[x1] = gvindex_get( value_index, search_in_chro_start + x1 );
		chro_bases_endside[x1] = gvindex_get( value_index , end_pos + search_in_read_start + x1);
	}

	int insertion_in_between_i, best_testing_score = 500 * 1000;
	int best_insertion_in_between = -1, best_last_exon_base_in_start = -1;
	int applied_insertion_limit = global_context->config.max_insertion_at_junctions;
	for(insertion_in_between_i = 0; insertion_in_between_i <= applied_insertion_limit; insertion_in_between_i ++){
		int start_site_match [ search_in_read_end - search_in_read_start ], end_site_match[ search_in_read_end - search_in_read_start  ];
		int start_last_exon_base,  end_site_mismatches = 0, start_site_mismatches = 0;
		for(start_last_exon_base = 0 ; start_last_exon_base < search_in_read_end - search_in_read_start ; start_last_exon_base++){
			start_site_match[start_last_exon_base] = ( rtext[ search_in_read_start + start_last_exon_base ] == chro_bases_startside[start_last_exon_base] );
			int end_site_x = ( rtext[ search_in_read_start + start_last_exon_base] == chro_bases_endside[start_last_exon_base] );
			end_site_match[start_last_exon_base] = end_site_x;

			if(start_last_exon_base >=insertion_in_between_i )
				end_site_mismatches += !end_site_x;
		}

		for(start_last_exon_base = 0 ; start_last_exon_base < search_in_read_end - search_in_read_start - insertion_in_between_i ; start_last_exon_base++){
			end_site_mismatches -= (! end_site_match[start_last_exon_base + insertion_in_between_i] );
			start_site_mismatches += (! start_site_match[start_last_exon_base] );

			if(start_last_exon_base >= 2 && start_last_exon_base < search_in_read_end - search_in_read_start -insertion_in_between_i -2){


			if(0&& FIXLENstrcmp("V0112_0155:7:1101:12618:2466#ACTTGA", rname) == 0)
				SUBREADprintf("split=%d, ins=%d, MM=%d+%d \n", start_last_exon_base, insertion_in_between_i, start_site_mismatches, end_site_mismatches);


				if( (end_site_mismatches + start_site_mismatches) * 500 + insertion_in_between_i < best_testing_score ){
					int donor_paired_ret=is_paired_donor_receptor( chro_bases_startside + start_last_exon_base + 1, chro_bases_endside + insertion_in_between_i + start_last_exon_base - 1 );

					if( donor_paired_ret ) {
						best_insertion_in_between = insertion_in_between_i;
						best_last_exon_base_in_start = start_last_exon_base;
						best_testing_score = (end_site_mismatches + start_site_mismatches) * 500 + insertion_in_between_i;
						(*misma_bases) = (end_site_mismatches + start_site_mismatches);
						if(is_negative_donor) (*is_negative_donor) =(donor_paired_ret -1)%2;
						(*matched_bases) = end_coverage - start_coverage - insertion_in_between_i - (end_site_mismatches + start_site_mismatches);
					}
				}

			}
		}
	}


	if(0 && FIXLENstrcmp("V0112_0155:7:1101:12618:2466", rname)==0)
	{
		chro_bases_startside[x1] = 0;
		chro_bases_endside[x1] = 0;
		char sp1s[200];
		for(x1 =0; x1<200; x1++) sp1s[x1]=' ';
		sp1s[search_in_read_start] =0;

		char spE[200];
		for(x1 =0; x1<200; x1++) spE[x1]=' ';
		spE[search_in_read_start + best_last_exon_base_in_start] =0;

		char spBB[200];
		for(x1 =0; x1<200; x1++) spBB[x1]=' ';
		spBB[ best_insertion_in_between] =0;

		char out1pos[100];
		absoffset_to_posstr(global_context, search_in_chro_start, out1pos);

		if(0 && FIXLENstrcmp("chr14:105",out1pos)==0){
			SUBREADprintf("POS=%s\t\tINS=%d\t\t%s\n", out1pos, best_insertion_in_between, rname);
			SUBREADprintf("R= %s\nS1=%s%s\nS2=%s%s\n   %s|%s|\n\n", rtext, sp1s, chro_bases_startside, sp1s, chro_bases_endside, spE, spBB);
		}
	}

	if(best_last_exon_base_in_start>=0)
		return best_last_exon_base_in_start + search_in_read_start ;
	else return -1;
}

int find_path(global_context_t * global_context, thread_context_t * thread_context, int start_element_i, int target_element_i, int * ii_array, int * jj_array, int * is_second_vote_array,  gene_vote_t * vote_1, gene_vote_t * vote_2, char * read_name_1, char * read_name_2, char * read_text_1, char * read_text_2, int read_len_1, int read_len_2, int is_negative_strand, int * this_mask , int * exon_last_base){
	gene_vote_t * start_vote = is_second_vote_array[start_element_i]?vote_2:vote_1;
	gene_vote_t * end_vote = is_second_vote_array[target_element_i]?vote_2:vote_1;

	int start_coverage = start_vote->coverage_end[ ii_array[start_element_i] ][ jj_array[start_element_i] ];
	int end_coverage = end_vote->coverage_start[ ii_array[target_element_i] ][ jj_array[target_element_i] ];
	unsigned int start_pos =  start_vote->pos[ ii_array[start_element_i] ][ jj_array[start_element_i] ];
	unsigned int end_pos   =  end_vote->pos[ ii_array[target_element_i] ][ jj_array[target_element_i] ];
	int ret = -1;

	long long dist = start_pos;
	dist -= end_pos;
	(*this_mask)=0;
	if( abs(dist)<50000 ) {
		if(start_vote == end_vote){
			if(start_coverage < end_coverage + 9){
				char * this_read_name = is_second_vote_array[start_element_i]?read_name_2:read_name_1;
				char * this_read_text = is_second_vote_array[start_element_i]?read_text_2:read_text_1;
				int this_read_len = is_second_vote_array[start_element_i]?read_len_2:read_len_1, mismatched_bases = 0, matched_in_the_uncovered_gap = 0;
				if(start_pos < end_pos){
					int indels_in_start =  start_vote -> current_indel_cursor [ ii_array[start_element_i]] [ jj_array[start_element_i] ] , donor_receptor_neg_strand = -1;
					int best_last_base_in_start_exon = find_donor_receptor(global_context, thread_context, this_read_name, this_read_text, this_read_len, start_coverage, end_coverage, start_pos, end_pos, indels_in_start,  start_vote -> votes[  ii_array[start_element_i]] [ jj_array[start_element_i] ],   start_vote -> votes[  ii_array[target_element_i]] [ jj_array[target_element_i] ], &mismatched_bases , &matched_in_the_uncovered_gap, &donor_receptor_neg_strand);

					if(best_last_base_in_start_exon > 0 && mismatched_bases<1){
						ret = matched_in_the_uncovered_gap  +  end_vote->coverage_end[ ii_array[target_element_i] ][ jj_array[target_element_i] ] - end_coverage;
						(*this_mask) = donor_receptor_neg_strand? CLUSTER_ALIGNMENT_DONOR_NEGATIVE_STRAND : 0 ;

						if(0)SUBREADprintf("FROM %d-%d to %d-%d : INC=%d,  UNCOV=%d/%d\n",
									start_vote->coverage_start[ ii_array[start_element_i] ][ jj_array[start_element_i] ],
									start_vote->coverage_end[ ii_array[start_element_i] ][ jj_array[start_element_i] ],
									end_vote -> coverage_start[ ii_array[target_element_i] ][ jj_array[target_element_i] ],
									end_vote -> coverage_end[ ii_array[target_element_i] ][ jj_array[target_element_i] ], ret,
									matched_in_the_uncovered_gap , end_coverage - start_coverage);
					
						// # of matched bases, from the end of the "start" section to the end of the end section.
						*exon_last_base = best_last_base_in_start_exon;
					}
				}
			}
		}else{
			ret = end_vote->coverage_end[ ii_array[target_element_i] ][ jj_array[target_element_i] ] - end_vote->coverage_start[ ii_array[target_element_i] ][ jj_array[target_element_i] ] ;
			// if the two sections are on two reads, check the first base of the second read is after the first base of the first read.
		}
	}
	return ret;
}


int process_voting_junction_PE_topK(global_context_t * global_context, thread_context_t * thread_context, subread_read_number_t pair_number, gene_vote_t * vote_1, gene_vote_t * vote_2, char * read_name_1, char * read_name_2, char * read_text_1, char * read_text_2, int read_len_1, int read_len_2, int is_negative_strand, gene_vote_number_t v1_all_subreads, gene_vote_number_t v2_all_subreads)
{
	topK_buffer_t * topbuf = thread_context?&thread_context->topKbuff:&global_context ->topKbuff ;
	
	vote_combination_t * comb_buffer = (vote_combination_t *) topbuf -> comb_buffer;
	simple_mapping_t * vote_simple_1_buffer, * vote_simple_2_buffer;
	vote_simple_1_buffer =(simple_mapping_t *) topbuf -> vote_simple_1_buffer;
	vote_simple_2_buffer =(simple_mapping_t *) topbuf -> vote_simple_2_buffer;
	//memset(comb_buffer, 0 , sizeof(vote_combination_t) * global_context -> config.max_vote_combinations);

	int is_second_read,i,j;
	int third_highest_votes[2][9];
	int is_fully_covered_1 = 0;
	int is_fully_covered_2 = 0;

	for(is_second_read = 0 ; is_second_read < 1 + global_context -> input_reads.is_paired_end_reads; is_second_read ++)
	{
		gene_vote_t * current_vote = is_second_read?vote_2:vote_1;
		int *top_three_buff = third_highest_votes[is_second_read], i , j;
		int * is_fully_covered = is_second_read?&is_fully_covered_2:&is_fully_covered_1;
		int current_read_len = is_second_read?read_len_2:read_len_1;

		memset(top_three_buff, 0 , global_context -> config.top_scores * sizeof(int));

		if((global_context-> config.do_fusion_detection || global_context-> config.do_long_del_detection)){
			*is_fully_covered = test_fully_covered(global_context , current_vote, current_read_len);
		}



		for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		{
			for (j=0; j< current_vote->items[i]; j++){
				int vv = current_vote -> votes[i][j];
				if(global_context->config.scRNA_input_mode && !global_context -> input_reads.is_paired_end_reads)vv += SE_READ_IN_KNOWN_EXON_REWARD*is_pos_in_annotated_exon_regions(global_context, current_vote -> pos[i][j]);
				update_top_three(global_context, top_three_buff, vv);
			}
		}

		if(0 && FIXLENstrcmp("R00000003493",read_name_1)==0)SUBREADprintf("3N [R %d] =%d,%d,%d\n", 1+is_second_read, top_three_buff[0], top_three_buff[1], top_three_buff[2]);

		for(i = 0; i < global_context -> config.multi_best_reads; i++)
		{
			mapping_result_t * old_result = _global_retrieve_alignment_ptr(global_context, pair_number, is_second_read, i);
			if(old_result -> selected_votes>0)
			{
				update_top_three(global_context, top_three_buff, old_result -> selected_votes);
			}
		}
		if(0 && FIXLENstrcmp("R00000003493",read_name_1)==0)SUBREADprintf("3Q [R %d] =%d,%d,%d\n", 1+is_second_read, top_three_buff[0], top_three_buff[1], top_three_buff[2]);
	}
	

	int simple_record_numbers[2], third_k;
	
	for(is_second_read = 0 ; is_second_read < 1 + global_context -> input_reads.is_paired_end_reads; is_second_read ++)
	{ 
		int current_simple_number = 0;
		int current_read_len = is_second_read?read_len_2:read_len_1;
		// populate the two simple read lists
		for(third_k = 0 ; third_k < global_context -> config.top_scores; third_k ++)
		{
			if(current_simple_number >= global_context -> config.max_vote_simples)break;
			int this_vote_N = third_highest_votes [is_second_read][third_k];
			// only consider max_votes and max_votes - 1
			if(this_vote_N<1 || (third_highest_votes[is_second_read][0] - this_vote_N > global_context -> config.max_vote_number_cutoff )) break;

			simple_mapping_t * current_simple = is_second_read ? vote_simple_2_buffer: vote_simple_1_buffer;
			gene_vote_t * current_vote = is_second_read?vote_2:vote_1;
			for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
			{
				if(current_simple_number >= global_context -> config.max_vote_simples)break;
				for (j=0; j< current_vote->items[i]; j++)
				{
					if(current_simple_number >= global_context -> config.max_vote_simples)break;
					if(global_context->config.do_big_margin_filtering_for_junctions && third_k == 0 && current_vote->votes[i][j] >= third_highest_votes [is_second_read][global_context -> config.top_scores - 1])
						insert_big_margin_record(global_context , _global_retrieve_big_margin_ptr(global_context,pair_number, is_second_read), current_vote -> votes[i][j], current_vote -> coverage_start[i][j], current_vote -> coverage_end[i][j] , current_read_len, (current_vote -> masks[i][j] & IS_NEGATIVE_STRAND)?1:0);
					
					int vv = current_vote->votes[i][j];
					if(global_context->config.scRNA_input_mode && !global_context -> input_reads.is_paired_end_reads)vv += SE_READ_IN_KNOWN_EXON_REWARD*is_pos_in_annotated_exon_regions(global_context,  current_vote -> pos[i][j]);
					if(vv == this_vote_N && current_vote->votes[i][j] >= global_context->config.minimum_subread_for_second_read)
					{
						current_simple[current_simple_number].is_vote_t_item = 1;
						current_simple[current_simple_number].item_index_i = i;
						current_simple[current_simple_number].item_index_j = j;
						current_simple[current_simple_number].read_start_base = current_vote -> coverage_start[i][j];
						current_simple[current_simple_number].mapping_position = current_vote -> pos[i][j];
						current_simple[current_simple_number].major_half_votes = vv;

						current_simple_number ++;
					
					}
				}
			}

			for(i = 0; i < global_context -> config.multi_best_reads; i++)
			{
				mapping_result_t * old_result = _global_retrieve_alignment_ptr(global_context, pair_number, is_second_read, i);
				if(current_simple_number >= global_context -> config.max_vote_simples)break;
				if(old_result -> selected_votes == this_vote_N)
				{
					current_simple[current_simple_number].is_vote_t_item = 0;
					current_simple[current_simple_number].item_index_i = i;
					current_simple[current_simple_number].mapping_position = old_result -> selected_position;
					current_simple[current_simple_number].major_half_votes = old_result -> selected_votes;
					current_simple[current_simple_number].read_start_base = old_result -> confident_coverage_start;

					current_simple_number ++;
				}
			}

		}
		simple_record_numbers[is_second_read] = current_simple_number;
	}

	int used_comb_buffer = 0;
	//calculate all combinations	

	if(global_context -> input_reads.is_paired_end_reads){
		for(i = 0; i < simple_record_numbers[0]; i++){
			for(j = 0; j < simple_record_numbers[1]; j++){
				int target_index;
				int is_PE_distance = 0, is_same_chromosome = 0, is_both_exonic_regions = 0;

				if(max(vote_simple_1_buffer[i].major_half_votes, vote_simple_2_buffer[j].major_half_votes) < global_context->config.minimum_subread_for_first_read)continue;
				
				simple_PE_and_same_chro(global_context , vote_simple_1_buffer+i, vote_simple_2_buffer+j , &is_PE_distance, &is_same_chromosome , read_len_1, read_len_2);

				if((!is_PE_distance) && min(vote_simple_1_buffer[i].major_half_votes, vote_simple_2_buffer[j].major_half_votes) < global_context->config.minimum_subread_for_first_read)continue;
				if( global_context ->  exonic_region_bitmap && is_same_chromosome)is_both_exonic_regions =  is_pos_in_annotated_exon_regions(global_context, vote_simple_1_buffer[i].mapping_position + vote_simple_1_buffer[i].read_start_base ) && is_pos_in_annotated_exon_regions(global_context, vote_simple_2_buffer[j].mapping_position + vote_simple_2_buffer[j].read_start_base  ) ;

				int adjusted_weight;

				if(1){
					if (is_both_exonic_regions && is_PE_distance) adjusted_weight = 1800;
					else if(is_both_exonic_regions) adjusted_weight = 1300;
					else if(is_PE_distance) adjusted_weight = 1300;
					else if(is_same_chromosome) adjusted_weight = 1000;
					else adjusted_weight = 800;
				}else{
					if (is_both_exonic_regions) adjusted_weight = 1300;
					else if(is_PE_distance) adjusted_weight = 1300;
					else if(is_same_chromosome) adjusted_weight = 1000;
					else adjusted_weight = 800;
				}
					//int adjusted_weight = is_PE_distance?1600:(is_same_chromosome?1000:500);
				int adjusted_votes = (vote_simple_1_buffer[i].major_half_votes + vote_simple_2_buffer[j].major_half_votes) * adjusted_weight;

				for(target_index=0; target_index<used_comb_buffer; target_index++){
					if(comb_buffer[target_index].score_adj < adjusted_votes) break;
				}


				if(target_index < global_context -> config.max_vote_combinations){
					int move_i;

					for(move_i = min(used_comb_buffer, global_context -> config.max_vote_combinations - 1) ; move_i > target_index ; move_i --)
						//checked: memory boundary
						memcpy(comb_buffer + move_i, comb_buffer + move_i - 1 , sizeof(vote_combination_t) );

					comb_buffer[target_index].r1_loc = vote_simple_1_buffer+i;
					comb_buffer[target_index].r2_loc = vote_simple_2_buffer+j;
					comb_buffer[target_index].score_adj = adjusted_votes;

					if(used_comb_buffer < global_context -> config.max_vote_combinations)
						used_comb_buffer ++;
				}

			}
		}
	}

	mapping_result_t * alignment_tmp_r1, * alignment_tmp_r2;
	alignment_tmp_r1 = (mapping_result_t *) topbuf -> alignment_tmp_r1;
	alignment_tmp_r2 = (mapping_result_t *) topbuf -> alignment_tmp_r2;

	subjunc_result_t * junction_tmp_r2 , * junction_tmp_r1;
	junction_tmp_r1 = (subjunc_result_t *) topbuf -> junction_tmp_r1;
	junction_tmp_r2 = (subjunc_result_t *) topbuf -> junction_tmp_r2;

	memset(junction_tmp_r1, 0, sizeof(subjunc_result_t) * global_context->config.multi_best_reads);
	memset(junction_tmp_r2, 0, sizeof(subjunc_result_t) * global_context->config.multi_best_reads);

	memset(alignment_tmp_r1, 0, sizeof(mapping_result_t) * global_context->config.multi_best_reads);
	memset(alignment_tmp_r2, 0, sizeof(mapping_result_t) * global_context->config.multi_best_reads);

	int alignment_res_r1_cursor = 0, alignment_res_r2_cursor = 0;

	if(used_comb_buffer > 0){
		merge_sort(comb_buffer, used_comb_buffer, comb_sort_compare, comb_sort_exchange, comb_sort_merge);
		for(is_second_read = 0; is_second_read < 1 + global_context -> input_reads.is_paired_end_reads; is_second_read++){
			int current_read_len = is_second_read ? read_len_2:read_len_1;
			char * current_read_text = is_second_read ? read_text_2:read_text_1;
			int current_all_subreads = is_second_read ? v2_all_subreads:v1_all_subreads;
			mapping_result_t * current_alignment_tmp = is_second_read?alignment_tmp_r2:alignment_tmp_r1;
			int * current_r_cursor = is_second_read ? &alignment_res_r2_cursor:&alignment_res_r1_cursor;
			int * is_fully_covered = is_second_read?&is_fully_covered_2:&is_fully_covered_1;
			gene_vote_t * current_vote = is_second_read?vote_2:vote_1;

			subjunc_result_t * current_junction_tmp = NULL;
			if(global_context -> config.do_breakpoint_detection) current_junction_tmp = is_second_read?junction_tmp_r2:junction_tmp_r1;

			for(i = used_comb_buffer - 1; i >=0; i--){
				if((* current_r_cursor) >= global_context->config.multi_best_reads)break;

				// add the combination of comb_buffer[i] into the two mapping_result_t arrays
				simple_mapping_t * current_loc = is_second_read?comb_buffer[i].r2_loc:comb_buffer[i].r1_loc;
				assert(current_loc);
				unsigned int current_pos = current_loc->mapping_position;

				int is_exist = 0;
				for(j = 0; j < *current_r_cursor; j++)
				{
					if(current_alignment_tmp[j].selected_position == current_pos){
						is_exist = 1;
						break;
					}
				}
				//SUBREADprintf("CLLL BUF %d R_%d : %u ; EXIST %d. Written into the %d-th best location\n", i, 1+is_second_read,  current_loc->mapping_position, is_exist, *current_r_cursor);

				if(!is_exist){
					if(current_loc -> is_vote_t_item)
						copy_vote_to_alignment_res(global_context, thread_context, current_alignment_tmp + (*current_r_cursor), current_junction_tmp ? current_junction_tmp + (*current_r_cursor) : NULL, current_vote, current_loc -> item_index_i, current_loc -> item_index_j, current_read_len, read_name_1, current_read_text, current_all_subreads , current_vote -> noninformative_subreads, pair_number, is_second_read, is_fully_covered);
					else{
						//checked: memory boundary
						memcpy(current_alignment_tmp + (*current_r_cursor), _global_retrieve_alignment_ptr(global_context, pair_number, is_second_read, current_loc -> item_index_i), sizeof(mapping_result_t));
						if(current_junction_tmp)
							//checked: memory boundary
							memcpy(current_junction_tmp + (*current_r_cursor), _global_retrieve_subjunc_ptr(global_context, pair_number, is_second_read, current_loc -> item_index_i), sizeof(subjunc_result_t));
					}
					(*current_r_cursor)++;
				}
			}
		}
	}else{// if the one end is not mapped at all

		if(0 ==  simple_record_numbers[0])
			_global_retrieve_alignment_ptr(global_context, pair_number, 0, 0) -> noninformative_subreads_in_vote = vote_1 -> noninformative_subreads;
		if(global_context -> input_reads.is_paired_end_reads && 0 ==  simple_record_numbers[1])
			_global_retrieve_alignment_ptr(global_context, pair_number, 1, 0) -> noninformative_subreads_in_vote = vote_2 -> noninformative_subreads;

		if(simple_record_numbers[0]>0 || simple_record_numbers[1]>0)
		{
			// copy all the simple into the mapping_result_t

			for(is_second_read = 0; is_second_read < 1 + global_context -> input_reads.is_paired_end_reads; is_second_read++)
			{
				int * current_r_cursor = is_second_read ? &alignment_res_r2_cursor:&alignment_res_r1_cursor;

				int current_read_len = is_second_read ? read_len_2:read_len_1;
				char * current_read_text = is_second_read ? read_text_2:read_text_1;
				int current_all_subreads = is_second_read ? v2_all_subreads:v1_all_subreads;
				mapping_result_t * current_alignment_tmp = is_second_read?alignment_tmp_r2:alignment_tmp_r1;
				gene_vote_t * current_vote = is_second_read?vote_2:vote_1;
				int * is_fully_covered = is_second_read?&is_fully_covered_2:&is_fully_covered_1;

				subjunc_result_t * current_junction_tmp = NULL;
				if(global_context -> config.do_breakpoint_detection) current_junction_tmp = is_second_read?junction_tmp_r2:junction_tmp_r1;

				for(i = 0; i < simple_record_numbers[is_second_read]; i++){

					if((*current_r_cursor) >= global_context->config.multi_best_reads)break;

					simple_mapping_t * current_loc = is_second_read?vote_simple_2_buffer+i:vote_simple_1_buffer+i;

					if(current_loc -> major_half_votes < global_context->config.minimum_subread_for_first_read) continue;
					unsigned int current_pos = current_loc->mapping_position;

					int is_exist = 0;
					for(j = 0; j < *current_r_cursor; j++)
					{
						if(current_alignment_tmp[j].selected_position == current_pos){
							is_exist = 1;
							break;
						}
					}
					if(!is_exist){
						if(current_loc -> is_vote_t_item)
							copy_vote_to_alignment_res(global_context, thread_context, current_alignment_tmp + (*current_r_cursor), current_junction_tmp ? current_junction_tmp + (*current_r_cursor): NULL, current_vote, current_loc -> item_index_i, current_loc -> item_index_j, current_read_len, read_name_1, current_read_text, current_all_subreads , current_vote -> noninformative_subreads, pair_number, is_second_read, is_fully_covered);
						else{
							//checked:boundary
							memcpy(current_alignment_tmp + (*current_r_cursor), _global_retrieve_alignment_ptr(global_context, pair_number, is_second_read, current_loc -> item_index_i), sizeof(mapping_result_t));
							if(current_junction_tmp)
								//checked:boundary
								memcpy(current_junction_tmp + (*current_r_cursor), _global_retrieve_subjunc_ptr(global_context, pair_number, is_second_read, current_loc -> item_index_i), sizeof(subjunc_result_t));
						}

						(*current_r_cursor)++;
					}
				}
			}
		}
	}

	for(is_second_read = 0; is_second_read < 1 +  global_context -> input_reads.is_paired_end_reads; is_second_read++){
		int * current_r_cursor = is_second_read ? &alignment_res_r2_cursor:&alignment_res_r1_cursor;
		if((*current_r_cursor) > global_context->config.multi_best_reads){
			SUBREADprintf("ERROR: multi_best_locations excessed the boundary: %d > %d\n", (*current_r_cursor), global_context->config.multi_best_reads);
			return -1;
		}
	}

	for(is_second_read = 0; is_second_read < 1 +  global_context -> input_reads.is_paired_end_reads; is_second_read++)
	{
		int * current_r_cursor = is_second_read ? &alignment_res_r2_cursor:&alignment_res_r1_cursor;
		mapping_result_t * current_alignment_tmp = is_second_read?alignment_tmp_r2:alignment_tmp_r1;
		subjunc_result_t * current_junction_tmp = NULL;

		if(global_context -> config.do_breakpoint_detection) current_junction_tmp = is_second_read?junction_tmp_r2:junction_tmp_r1;

		for(i = 0; i < global_context->config.multi_best_reads ; i++){
			mapping_result_t * cur_res = _global_retrieve_alignment_ptr(global_context, pair_number, is_second_read, i);
			if( i < (*current_r_cursor))
				memcpy(cur_res, current_alignment_tmp + i, sizeof(mapping_result_t));
			else	cur_res -> selected_votes = 0;

			if(global_context -> config.do_breakpoint_detection) {
				subjunc_result_t * cur_junc =  _global_retrieve_subjunc_ptr(global_context, pair_number, is_second_read, i);
				if(i  < (*current_r_cursor))
					memcpy(cur_junc, current_junction_tmp + i , sizeof(subjunc_result_t));
				else	cur_junc -> minor_votes = 0;

			}
		}
	}

	return 0;
}


// seq1 and seq2 must be on the same strand!
// (seq2 is reversed)
// The second half of seq1 MUST BE the same as the first half of seq2 if the two reads have an overlapping part.
int is_gapped_as_funky(global_context_t * global_context, char * rname1, char * chr1, unsigned int pos1, int rlen1, int is_1_negative, char * cigar1, char * seq1, char * rname2, char * chr2, unsigned int pos2, int rlen2, int is_2_negative, char * cigar2, char * seq2, int tlen_removed_intron)
{
/*
	if(tlen_removed_intron >= rlen1 + rlen2) return 1;	// may be gapped.
	int try_overlapping;

	int best_matched_bases = 0;
	int best_overlapping_len = -1;

	int assumed_overlapping = rlen1+rlen2-tlen_removed_intron;
	for(try_overlapping = 0; try_overlapping < min(rlen1, rlen2); try_overlapping++)
	{
		int r1_start = rlen1 - try_overlapping;
		int r2_end = try_overlapping;
		int xk1;
		int all_matched = 0, all_mismatched = 0;
		for(xk1 = 0; xk1 < r2_end; xk1++){
			char r1ch = seq1[r1_start + xk1];
			char r2ch = seq2[xk1];
			if(r1ch==r2ch) all_matched++;
			else all_mismatched++;
		}

		if(all_mismatched <= 1 && try_overlapping == assumed_overlapping){
			// the assumed overlapping length is good enough.
			return 0;
		}
		if(all_mismatched <= 1 && all_matched > best_matched_bases){
			best_overlapping_len = try_overlapping;
			best_matched_bases = all_matched;
		}
	}

	if(best_overlapping_len <= 0)return 0;
	return assumed_overlapping 
*/
	return tlen_removed_intron > 600;
}

// the positions are not offset by adding the first soft clipping length. I.e., pos1 and pos2 may be smaller than those in the SAM files.
// seq1 and seq2 must be on the same strand!
// (seq2 is reversed)
int is_funky_fragment(global_context_t * global_context, char * rname1, char * chr1, unsigned int pos1, int rlen1, int is_1_negative, char * cigar1, char * seq1, char * rname2, char * chr2, unsigned int pos2, int rlen2, int is_2_negative, char * cigar2, char * seq2, int tlen_removed_intron)
{
	long long llraw_tlen = pos1;
	llraw_tlen -= pos2;
	if(llraw_tlen <0)
		llraw_tlen = -llraw_tlen;
	unsigned int raw_tlen = llraw_tlen;
	raw_tlen += max(rlen2, rlen1);

	//SUBREADprintf("CHRS=%p,%p,  POS=%u,%u,  RTLEN=%u\n", chr1, chr2, pos1, pos2, raw_tlen);

	if(chr1 != chr2) raw_tlen = 0;
	
	// note: the two pointers can be compared because they should be derived from the offset table.
	// Each chromosome name should have one and only one distinct char * pointer.
	if(chr1 == chr2 && raw_tlen <= global_context -> config.maximum_translocation_length && is_2_negative == is_1_negative)
	{
		if(is_gapped_as_funky(global_context, rname1, chr1, pos1, rlen1, is_1_negative, cigar1, seq1, rname2, chr2, pos2, rlen2, is_2_negative, cigar2, seq2, tlen_removed_intron))
			return FUNKY_FRAGMENT_A;
		else	return NOT_FUNKY;
	}
	else if( chr1 == chr2 && raw_tlen <= global_context -> config.maximum_translocation_length && is_2_negative != is_1_negative )
		return FUNKY_FRAGMENT_DE;
	else if( chr1 != chr2 || raw_tlen > global_context -> config.maximum_translocation_length) 
		return FUNKY_FRAGMENT_BC;

	return NOT_FUNKY;
}

int process_voting_junction(global_context_t * global_context, thread_context_t * thread_context, subread_read_number_t pair_number, gene_vote_t * vote_1, gene_vote_t * vote_2, char * read_name_1, char * read_name_2, char * read_text_1, char * read_text_2, int read_len_1, int read_len_2, int is_negative_strand, gene_vote_number_t v1_all_subreads, gene_vote_number_t v2_all_subreads){


	//#warning "FOR TESTING CLUSTER_BASED JUNCTION DETECTION ONLY!!."
	//return process_voting_junction_PE_juncs(global_context, thread_context, pair_number, vote_1, vote_2, read_name_1, read_name_2, read_text_1, read_text_2, read_len_1, read_len_2, is_negative_strand, v1_all_subreads, v2_all_subreads);
		return process_voting_junction_PE_topK(global_context, thread_context, pair_number, vote_1, vote_2, read_name_1, read_name_2, read_text_1, read_text_2, read_len_1, read_len_2, is_negative_strand, v1_all_subreads, v2_all_subreads);
		
}


unsigned int explain_read(global_context_t * global_context, thread_context_t * thread_context, realignment_result_t * final_realignments, subread_read_number_t pair_number, int read_len, char * read_name , char *read_text, char *qual_text, int is_second_read, int best_read_id, int is_negative_strand)
{
	explain_context_t explain_context;
	mapping_result_t *current_result = _global_retrieve_alignment_ptr(global_context, pair_number, is_second_read, best_read_id); 

	if(global_context -> config.do_big_margin_filtering_for_reads)
	{
		int current_repeated_times = is_ambiguous_voting(global_context, pair_number, is_second_read, current_result->selected_votes, current_result->confident_coverage_start, current_result->confident_coverage_end, read_len, (current_result->result_flags & CORE_IS_NEGATIVE_STRAND)?1:0);
		if( global_context -> config.do_big_margin_filtering_for_reads && current_repeated_times>1) return 0;
	}

	memset(&explain_context,0, sizeof(explain_context_t));
	explain_context.full_read_len = read_len;
	explain_context.is_fully_covered = current_result -> is_fully_covered ;
	explain_context.full_read_text = read_text;
	explain_context.full_qual_text = qual_text;
	explain_context.read_name = read_name;
	explain_context.is_confirmed_section_negative_strand = is_negative_strand ;
	explain_context.pair_number = pair_number;
	explain_context.is_second_read = is_second_read ;
	explain_context.best_read_id = best_read_id;
	explain_context.total_tries = 0;

	if(0 && FIXLENstrcmp("simulated.24700032", explain_context.read_name)==0)SUBREADprintf("BBFINAL %s SEL_POS=%u COV=%d - %d\n",  explain_context.read_name, current_result -> selected_position, current_result -> confident_coverage_start, current_result -> confident_coverage_end);

	unsigned int back_search_tail_position,front_search_start_position; 
	unsigned short back_search_read_tail, front_search_read_start;


	back_search_read_tail = min(explain_context.full_read_len , current_result -> confident_coverage_end );//- 5;
	back_search_tail_position = current_result -> selected_position + back_search_read_tail +  current_result -> indels_in_confident_coverage;

	//if( back_search_read_tail > 102)
	//SUBREADprintf("MAX back_search_read_tail : MIN %d , %d\n", explain_context.full_read_len , current_result -> confident_coverage_end);

	explain_context.tmp_search_junctions[0].read_pos_end = back_search_read_tail;
	explain_context.tmp_search_junctions[0].abs_offset_for_start = back_search_tail_position;

	explain_context.all_back_alignments = 0;
	explain_context.tmp_search_sections = 0;
	explain_context.best_indel_penalty =0;
	explain_context.best_matching_bases = -9999;
	explain_context.second_best_matching_bases = -9999;
	explain_context.tmp_indel_penalty = 0;
	explain_context.tmp_total_matched_bases = 0;
	explain_context.is_currently_tie = 0;
	explain_context.best_is_complex = 0;
	explain_context.best_support_as_simple = 0;
	explain_context.best_min_unsupport_as_simple = 0;
	explain_context.tmp_support_as_simple = 0;
	explain_context.tmp_min_support_as_complex = 999999;
	explain_context.tmp_min_unsupport = 999999;
	explain_context.tmp_is_pure_donor_found_explain = 1;
	explain_context.best_is_pure_donor_found_explain = 0;

	if(1) {
		front_search_read_start = back_search_read_tail > 8? back_search_read_tail - 8:0;
		front_search_start_position = back_search_tail_position>8?back_search_tail_position - 8:0;
	} else {
		//front_search_read_start = current_result -> confident_coverage_start + 5;
		front_search_read_start = min(explain_context.full_read_len , current_result -> confident_coverage_end);
		if(front_search_read_start > 2*global_context -> config.realignment_minimum_variant_distance) front_search_read_start -= 2*global_context -> config.realignment_minimum_variant_distance;
		else front_search_read_start = 0;
		front_search_start_position = current_result -> selected_position + front_search_read_start;
	}

	search_events_to_back(global_context, thread_context, &explain_context, read_text , qual_text, back_search_tail_position , back_search_read_tail, 0, 0, 1);
	int back_penalty = explain_context.best_indel_penalty;

	//int is_backsearch_tie = explain_context.is_currently_tie;
	int back_search_matches_diff = -9999;

	/*
 

	if(explain_context.back_search_confirmed_sections>0)
	{
		
		short last_section_length = explain_context.back_search_junctions[0].read_pos_end - explain_context.back_search_junctions[0].read_pos_start;
		
		front_search_read_start = explain_context.back_search_junctions[0].read_pos_start; 
		front_search_start_position = explain_context.back_search_junctions[0].abs_offset_for_start - last_section_length;

		int last_sec = explain_context.back_search_confirmed_sections-1;

		current_result -> selected_position = explain_context.back_search_junctions[last_sec].abs_offset_for_start - explain_context.back_search_junctions[last_sec].read_pos_end + explain_context.back_search_junctions[last_sec].read_pos_start;
		back_search_matches_diff = explain_context.best_matching_bases - explain_context.second_best_matching_bases;
 
		if(0 && memcmp(explain_context.read_name,  TTTSNAME, 26)==0)
		{
			int xk1;
			for(xk1 = 0; xk1 < explain_context.back_search_confirmed_sections; xk1++)
			{
				short pr_section_length = explain_context.back_search_junctions[xk1].read_pos_end - explain_context.back_search_junctions[xk1].read_pos_start;
				if(explain_context.back_search_junctions[xk1].event_after_section)
					SUBREADprintf("BACK_SECTIONS [%d], START IS %u; RPSS=%d ; RPED=%d ; LEN=%d ; EVENT is %u %u INDEL=%d\n", xk1, explain_context.back_search_junctions[xk1].abs_offset_for_start, explain_context.back_search_junctions[xk1].read_pos_start, explain_context.back_search_junctions[last_sec].read_pos_end, pr_section_length, explain_context.back_search_junctions[xk1].event_after_section->event_small_side, explain_context.back_search_junctions[xk1].event_after_section->event_large_side, explain_context.back_search_junctions[xk1].event_after_section->indel_length);
				else	SUBREADprintf("BACK_SECTIONS [%d], START IS %u; RPSS=%d ; RPED=%d ; LEN=%d\n", xk1, explain_context.back_search_junctions[xk1].abs_offset_for_start, explain_context.back_search_junctions[xk1].read_pos_start, explain_context.back_search_junctions[last_sec].read_pos_end, pr_section_length);
			}
		}

		//SUBREADprintf("DBI:%d - %d;\n", explain_context.best_matching_bases , explain_context.second_best_matching_bases);
	}
	else
	*/
	explain_context.all_front_alignments = 0;
	explain_context.tmp_search_sections = 0;
	explain_context.best_indel_penalty = 0;
	explain_context.best_matching_bases = -9999;
	explain_context.second_best_matching_bases = -9999;
	explain_context.tmp_total_matched_bases = 0;
	explain_context.tmp_indel_penalty = 0;

	explain_context.is_currently_tie = 0;
	explain_context.best_is_complex = 0;
	explain_context.best_support_as_simple = 0;
	explain_context.best_min_unsupport_as_simple = 0;
	explain_context.tmp_support_as_simple = 0;
	explain_context.tmp_min_support_as_complex = 999999;
	explain_context.tmp_min_unsupport = 999999;
	explain_context.tmp_is_pure_donor_found_explain = 1;
	explain_context.best_is_pure_donor_found_explain = 0;

	memset(explain_context.tmp_search_junctions, 0, sizeof(perfect_section_in_read_t ) * MAX_EVENTS_IN_READ);

	explain_context.tmp_search_junctions[0].read_pos_start = front_search_read_start;
	explain_context.tmp_search_junctions[0].abs_offset_for_start = front_search_start_position;

	if(0 && FIXLENstrcmp("R000002689",explain_context.read_name ) == 0)
		SUBREADprintf("Enter F_SEARCH: start=%u  read_pos=%d  REMAIN=%d\n", front_search_start_position, front_search_read_start,  read_len - front_search_read_start );


	short search_remain =  read_len - front_search_read_start;
	//#warning "SUBREAD_151 REMOVE THE ASSERT! "
	//if(search_remain >= 102)SUBREADprintf("FATAL: RLEN=%d, SEARCH=%d\n", read_len, front_search_read_start);
	//assert( search_remain < 102 );

	search_events_to_front(global_context, thread_context, &explain_context, read_text + front_search_read_start, qual_text + front_search_read_start, front_search_start_position, search_remain , 0, 0, 1);
	if(0 && FIXLENstrcmp("R_chr901_932716_91M1D9M",explain_context.read_name ) == 0)
		 SUBREADprintf("F_SEARCH has found %d result sets\n", explain_context.all_front_alignments);

	explain_context.best_indel_penalty += back_penalty;
	//int is_frontsearch_tie = explain_context.is_currently_tie;

	//SUBREADprintf("DFI:%d - %d;\n", explain_context.best_matching_bases , explain_context.second_best_matching_bases);
	int front_search_matches_diff = explain_context.best_matching_bases - explain_context.second_best_matching_bases;
	explain_context.best_second_match_diff = front_search_matches_diff + back_search_matches_diff;

	int realignment_number = finalise_explain_CIGAR(global_context, thread_context, &explain_context, final_realignments);

	if(0 && FIXLENstrcmp("SRR3439488.572382", explain_context.read_name)==0) 
		SUBREADprintf("TRYING_REALIGN:%s:%u\n", explain_context.read_name, explain_context.total_tries);

	return realignment_number;
}


void debug_clipping(global_context_t * global_context,  thread_context_t * thread_context, gene_value_index_t * current_value_index, char * read_text, unsigned int mapped_pos, int test_len,  int search_to_tail, int search_center, int number_of_clipped, char * read_name){

	//if(test_len>100)return;

	int xk1;

	SUBREADprintf("\n %s CENTER=%d, CLIPPED=%d, TLEN=%d    %s\n", read_name, search_center, number_of_clipped, test_len, search_to_tail?">>>>":"<<<<");

	for(xk1 = 0 ; xk1 < test_len ; xk1++)
	{
		char reference_base = gvindex_get(current_value_index, xk1 + mapped_pos);
		SUBREADprintf("%c", reference_base == read_text[xk1] ? '-':'#');
	}

	SUBREADprintf("\n");
	for(xk1 = 0 ; xk1 < test_len ; xk1++)
	{
		if(xk1 == search_center)
			SUBREADprintf("%c", search_to_tail?'>':'<');
		else SUBREADprintf(" ");
	}

	SUBREADprintf("\n");
	for(xk1 = 0 ; xk1 < test_len ; xk1++)
	{
		if( search_to_tail && xk1 >= test_len - number_of_clipped)
			SUBREADprintf("R");
		else if( (!search_to_tail) && xk1 <= number_of_clipped - 1)
			SUBREADprintf("L");
		else SUBREADprintf(" ");
	}

	SUBREADprintf("\n");

}

#define SOFT_CLIPPING_WINDOW_SIZE 5
#define SOFT_CLIPPING_MAX_ERROR   1

// it returns the number of bases to be clipped off.
int find_soft_clipping(global_context_t * global_context,  thread_context_t * thread_context, gene_value_index_t * current_value_index, char * read_text, unsigned int mapped_pos, int test_len,  int search_to_tail, int search_center)
{
	int base_in_window = 0;
	int added_base_index = 0, removed_base_index = 0;
	int search_start = 0;
	int matched_in_window = SOFT_CLIPPING_WINDOW_SIZE;
	int last_matched_base_index = -1, delta;

	if(search_to_tail)
	{
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

	for(added_base_index = search_start; added_base_index >= 0 && added_base_index < test_len; added_base_index += delta)
	{
		// add the new base
		char reference_base = gvindex_get(current_value_index, added_base_index + mapped_pos);

		if(0){
			char outpos1[100];
			absoffset_to_posstr(global_context, added_base_index + mapped_pos, outpos1);
			SUBREADprintf("CHMAT [%s] %s (%u) ref:read = %c:%c\n", search_to_tail?"T":"H" ,outpos1,  added_base_index + mapped_pos, reference_base,  read_text[added_base_index]);
		}
		int added_is_matched = (reference_base == read_text[added_base_index]);
		matched_in_window += added_is_matched;
		if(added_is_matched)
			last_matched_base_index = added_base_index;

		base_in_window ++;

		if(base_in_window > SOFT_CLIPPING_WINDOW_SIZE){
			removed_base_index = added_base_index - delta * SOFT_CLIPPING_WINDOW_SIZE;
			char removing_ref_base = gvindex_get(current_value_index, removed_base_index + mapped_pos);
			matched_in_window -= (removing_ref_base == read_text[removed_base_index]);
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

// read_head_abs_offset is the first WANTED base in read.
// If the first section in read is reversed, read_head_abs_offset is the LAST WANTED bases in this section. (the abs offset of the first base in the section is actually larger than read_head_abs_offset)
int final_CIGAR_quality(global_context_t * global_context, thread_context_t * thread_context, char * read_text, char * qual_text, int read_len, char * cigar_string, unsigned long read_head_abs_offset, int is_read_head_reversed, int * mismatched_bases, int covered_start, int covered_end, char * read_name, int * non_clipped_length, int *total_indel_length, int * matched_bases, int * chromosomal_length, int * full_section_clipped)
{
	int cigar_cursor = 0;
	int read_cursor = 0;
	unsigned int current_perfect_section_abs = read_head_abs_offset;
	int rebuilt_read_len = 0, total_insertion_length = 0;
	float all_matched_bases = 0;
	gene_value_index_t * current_value_index = thread_context?thread_context->current_value_index:global_context->current_value_index; 
	int current_reversed = is_read_head_reversed;
	int all_mismatched = 0;
	int is_First_M = 1, is_wrong_cigar = 0;
	int head_soft_clipped = -1, tail_soft_clipped = -1;
	unsigned int tmp_int = 0;

	//SUBREADprintf("Coverage : %d ~ %d\n", covered_start, covered_end);

	if(0){
		char posout1[100];
		int chro_max = get_offset_maximum_chro_pos(global_context,thread_context,read_head_abs_offset);
		absoffset_to_posstr(global_context, read_head_abs_offset, posout1);
		SUBREADprintf("READ %s : mapped to %s ; max_pos=%d\n", read_name,  posout1, chro_max);
	}

	while(1)
	{
		char nch = cigar_string[cigar_cursor++];
		if(!nch)break;
		if(isdigit(nch))
			tmp_int = tmp_int*10+(nch-'0');
		else{
			if(tmp_int == 0)is_wrong_cigar = 1;
			if(is_wrong_cigar) break;
			if(nch == 'M' || nch == 'S')
			{
				char *qual_text_cur;
				if(qual_text[0])qual_text_cur = qual_text+read_cursor;
				else	qual_text_cur = NULL;

				float section_qual;

				int is_Last_M = (cigar_string[cigar_cursor]==0);
				int has_clipping_this_section_head = 0, has_clipping_this_section_tail = 0;
				char * reversed_first_section_text = NULL;

				if(0){
					int is_head_in_chro = get_offset_maximum_chro_pos(global_context,thread_context, current_perfect_section_abs );
					int is_end_in_chro = get_offset_maximum_chro_pos(global_context,thread_context, current_perfect_section_abs + tmp_int );
					char posout1[100];
					char posout2[100];
					int chro_max = get_offset_maximum_chro_pos(global_context,thread_context, current_perfect_section_abs );
					absoffset_to_posstr(global_context, current_perfect_section_abs, posout1);
					absoffset_to_posstr(global_context, current_perfect_section_abs + tmp_int, posout2);
					SUBREADprintf("  %dM SECTION : mapped to %s ~ %s ; max_pos=%d ; Hin=%d, Ein=%d\n", tmp_int,  posout1, posout2, chro_max, is_head_in_chro, is_end_in_chro);
					SUBREADprintf("  %dM SECTION : Hin=%d, Ein=%d\n", tmp_int, is_head_in_chro, is_end_in_chro);
				}

				// find "J" sections if it is the first M
				if(is_First_M && global_context -> config.show_soft_cliping)
				{
					int adj_coverage_start = covered_start - read_cursor;

					if(current_reversed)
					{
						reversed_first_section_text = malloc(MAX_READ_LENGTH);
						memcpy(reversed_first_section_text, read_text, tmp_int);
						reverse_read(reversed_first_section_text, tmp_int,  global_context->config.space_type);

						head_soft_clipped = find_soft_clipping(global_context, thread_context, current_value_index, reversed_first_section_text, current_perfect_section_abs, tmp_int, 1, 0);
					}
					else
						head_soft_clipped = find_soft_clipping(global_context, thread_context, current_value_index, read_text, current_perfect_section_abs, tmp_int, 0, adj_coverage_start);
					//SUBREADprintf("SSHEAD:%d\n", head_soft_clipped);

					if(head_soft_clipped == tmp_int){
						(*full_section_clipped) = 1;
						head_soft_clipped = 0;
					}
					else has_clipping_this_section_head = 1;

					if(has_clipping_this_section_head){
						if( tmp_int - head_soft_clipped < 3 && head_soft_clipped > 1 ) (*full_section_clipped) = 1;
					}

					if(reversed_first_section_text)
						free(reversed_first_section_text);
					reversed_first_section_text = NULL;
				}
				if(is_Last_M && global_context -> config.show_soft_cliping)
				{
					int adj_coverage_end = covered_end - read_cursor;

					if(current_reversed)
					{
						reversed_first_section_text = malloc(MAX_READ_LENGTH);
						// checked: boundary
						memcpy(reversed_first_section_text, read_text + read_cursor, tmp_int);
						reverse_read(reversed_first_section_text, tmp_int,  global_context->config.space_type);
						tail_soft_clipped = find_soft_clipping(global_context, thread_context, current_value_index, reversed_first_section_text, current_perfect_section_abs, tmp_int, 0, tmp_int);
					}
					else
						tail_soft_clipped = find_soft_clipping(global_context, thread_context, current_value_index, read_text + read_cursor, current_perfect_section_abs, tmp_int, 1, adj_coverage_end);

					if(1 && FIXLENstrcmp("NS500643:556:HGTMTBGXB:4:13403:18179:8012", read_name)==0) 
						SUBREADprintf("SSTAIL:%d\n", tail_soft_clipped);

					if(1 && tail_soft_clipped == tmp_int){
						tail_soft_clipped = 0;
						if(full_section_clipped)(*full_section_clipped) = 1;
					} else has_clipping_this_section_tail = 1;

					if( has_clipping_this_section_tail ){
						if(tmp_int - tail_soft_clipped < 3 && tail_soft_clipped > 1) (*full_section_clipped) = 1;
					}

					if(reversed_first_section_text)
						free(reversed_first_section_text);
				}

				if(is_Last_M && is_First_M && tail_soft_clipped+head_soft_clipped >= tmp_int-1)
				{
					head_soft_clipped=0;
					tail_soft_clipped=0;
				}

				int mismatch_calculation_start = has_clipping_this_section_head?head_soft_clipped:0;
				int mismatch_calculation_end = has_clipping_this_section_tail?tail_soft_clipped:0;

				if(global_context -> config.space_type == GENE_SPACE_COLOR)
					section_qual =  match_base_quality_cs(current_value_index, read_text+read_cursor, current_perfect_section_abs, qual_text_cur, tmp_int, global_context->config.phred_score_format , mismatched_bases, &all_mismatched, global_context -> config.high_quality_base_threshold, mismatch_calculation_start, mismatch_calculation_end);
				else
					section_qual =  match_base_quality(current_value_index, read_text+read_cursor, current_perfect_section_abs, qual_text_cur, tmp_int, current_reversed, global_context->config.phred_score_format , mismatched_bases, &all_mismatched, global_context -> config.high_quality_base_threshold, mismatch_calculation_start, mismatch_calculation_end);
				all_matched_bases += section_qual;
				rebuilt_read_len += tmp_int;
				is_First_M=0;

				read_cursor += tmp_int;

				//move to the NEXT UNWANTED ABS OFFSET. 
				if(current_reversed)
					current_perfect_section_abs --;
				else
					current_perfect_section_abs += tmp_int;


			}
			else if(nch == 'I')
			{
				rebuilt_read_len += tmp_int;
				read_cursor += tmp_int;

				all_matched_bases += tmp_int;
				total_indel_length += tmp_int;
				total_insertion_length += tmp_int;
			}
			else if(nch == 'D')
			{
				total_indel_length ++;
				if(!current_reversed)
					current_perfect_section_abs += tmp_int;
			}
			else if(tolower(nch) == 'n')
			{
				total_indel_length ++;
				current_perfect_section_abs += tmp_int;
				if(nch == 'n') current_reversed = !current_reversed;
			}
			else if(tolower(nch) == 'b')
			{
				total_indel_length ++;
				current_perfect_section_abs -= tmp_int;
				if(nch == 'b') current_reversed = !current_reversed;
			}

			if(read_cursor>MAX_READ_LENGTH){
				SUBREADprintf("ERROR: Cigar section longer than read length: %d >= %d, '%s'\n", tmp_int , MAX_READ_LENGTH, cigar_string);
				is_wrong_cigar = 1;
			}

			tmp_int = 0;
		}
	}

	int my_non_clipped_length = read_len;
	my_non_clipped_length -= max(0,tail_soft_clipped);
	my_non_clipped_length -= max(0,head_soft_clipped);

	//#warning " ========== COMMENT THIS LINE !! ========="
	//printf("QCR ALL MM=%d, RBLEN=%d, MAPPED_LEN=%d ; CIGAR=%s\n", all_mismatched, rebuilt_read_len , my_non_clipped_length, cigar_string);
	
	if(is_wrong_cigar || rebuilt_read_len != read_len || my_non_clipped_length < global_context->config.min_mapped_fraction){
		(*mismatched_bases)=99999;
		all_matched_bases = 0;
		sprintf(cigar_string, "%dM", read_len);
	}
	else if((head_soft_clipped>0 || tail_soft_clipped>0))
	{
		char new_cigar_tmp[120];
		is_First_M=1;
		new_cigar_tmp[0]=0;
		cigar_cursor = 0;
		while(1)
		{
			char nch = cigar_string[cigar_cursor++];

			if(!nch)break;
			if(isdigit(nch))
				tmp_int = tmp_int*10+(nch-'0');
			else{
				char cigar_piece [30];
				cigar_piece[0]=0;

				if(nch == 'M')
				{
					char cigar_tiny [12];
					int is_Last_M = (cigar_string[cigar_cursor]==0);
					if(is_First_M && head_soft_clipped>0)
					{
						tmp_int -= head_soft_clipped;
						sprintf(cigar_tiny,"%dS",head_soft_clipped);
						strcat(cigar_piece, cigar_tiny);
					}
					if(is_Last_M && tail_soft_clipped>0)
					{
						tmp_int -= tail_soft_clipped;
					}
					sprintf(cigar_tiny,"%dM",tmp_int);
					strcat(cigar_piece, cigar_tiny);
					if(is_Last_M && tail_soft_clipped>0)
					{
						sprintf(cigar_tiny,"%dS",tail_soft_clipped);
						strcat(cigar_piece, cigar_tiny);
					}
					is_First_M = 0;
				}
				else
				{
					sprintf(cigar_piece, "%u%c", tmp_int, nch);
				}

				strcat(new_cigar_tmp, cigar_piece);
				tmp_int = 0;
			}
		}

		if(1 && FIXLENstrcmp("NS500643:556:HGTMTBGXB:4:13403:18179:8012", read_name)==0) 
			SUBREADprintf("NEW_CIGAR_2 : %s\n", new_cigar_tmp);
		strcpy(cigar_string, new_cigar_tmp);
	}

	if((*mismatched_bases) != 99999)
		(*mismatched_bases) = all_mismatched;

	(*non_clipped_length) = my_non_clipped_length;
	(*matched_bases) = my_non_clipped_length - all_mismatched - total_insertion_length;
	(*chromosomal_length) = current_perfect_section_abs - read_head_abs_offset + total_insertion_length;

	return max(0, (int)(all_matched_bases*60/my_non_clipped_length));
}

// this function also adds final_counting_reads in chromosome_events.
unsigned int finalise_explain_CIGAR(global_context_t * global_context, thread_context_t * thread_context, explain_context_t * explain_context, realignment_result_t * final_realignments)
{
	int xk1, front_i, back_i;
	char tmp_cigar[120];
	chromosome_event_t * to_be_supported [20];
	short flanking_size_left[20], flanking_size_right[20];
	int to_be_supported_count = 0;
	int is_junction_read = 0;
	int total_perfect_matched_sections = 0;

	mapping_result_t * result = _global_retrieve_alignment_ptr(global_context, explain_context->pair_number, explain_context->is_second_read, explain_context-> best_read_id); 
	result -> result_flags &= ~CORE_IS_FULLY_EXPLAINED;
	result -> result_flags &= ~CORE_IS_PAIRED_END;

	//SUBREADprintf("FINAL_CIGAR R1 %d[%d] = %p, FLAGS=%d\n", explain_context -> pair_number , explain_context-> best_read_id , result , result -> result_flags);
	// reverse the back_search result for every equally best alignment
	//

	for(back_i = 0; back_i < explain_context -> all_back_alignments; back_i++){
		if( explain_context -> result_back_junction_numbers[back_i] > MAX_EVENTS_IN_READ ){
			SUBREADprintf("ERROR: Too many cigar sections: %d > %d\n", explain_context -> result_back_junction_numbers[back_i] , MAX_EVENTS_IN_READ);
			return 0;
		}
		for(xk1=0; xk1<explain_context -> result_back_junction_numbers[back_i]/2; xk1++)
		{
			perfect_section_in_read_t tmp_exp;
			// checked: boundary
			memcpy(&tmp_exp, &explain_context -> result_back_junctions[back_i][xk1], sizeof(perfect_section_in_read_t));
			memcpy(&explain_context -> result_back_junctions[back_i][xk1],  &explain_context -> result_back_junctions[back_i][explain_context -> result_back_junction_numbers[back_i] - xk1 - 1] , sizeof(perfect_section_in_read_t));
			memcpy(&explain_context -> result_back_junctions[back_i][explain_context -> result_back_junction_numbers[back_i] - xk1 - 1] , &tmp_exp , sizeof(perfect_section_in_read_t));
		} 
	}

	// adding indel lengths in read lengths and relocate sections
	// note that the last section in back results has the same strand of the main piece.

	int is_cigar_overflow = 0, fusions_in_read = 0, final_alignment_number = 0;
	for(back_i = 0; back_i < explain_context -> all_back_alignments; back_i++){
		if(final_alignment_number >= MAX_ALIGNMENT_PER_ANCHOR)break;

		int is_first_section_negative = (result ->result_flags & CORE_IS_NEGATIVE_STRAND)?1:0; 
		for(xk1=0; xk1<explain_context -> result_back_junction_numbers[back_i]; xk1++)
		{
			int section_length = explain_context -> result_back_junctions[back_i][xk1].read_pos_end - explain_context -> result_back_junctions[back_i][xk1].read_pos_start; 
			if(0 && FIXLENstrcmp("simulated.11420793", explain_context->read_name)==0)SUBREADprintf("FINAL_EXPLAIN %s BACK_%d SEC_%d OLD_START=%d SEC_LENG=%d\n", explain_context->read_name, back_i, xk1,  explain_context -> result_back_junctions[back_i][xk1].abs_offset_for_start, section_length);
			unsigned int new_start_pos;

			if(explain_context -> result_back_junctions[back_i][xk1].is_strand_jumped)
				// the "strand_jumped" section do not need to move
				// however, the "abs_offset_for_start" is actually for the last base in this section.
				// this does not metter if we compare the reversed read to the chromosome.
				// "abs_offset_for_start" is the first UNWANTED base (smaller than the first WANTED base)
				new_start_pos = explain_context -> result_back_junctions[back_i][xk1].abs_offset_for_start +1;
			else
				// "abs_offset_for_start" is the first UNWANTED base. By subtracting the length, it becomes the first WANTED base.
				new_start_pos = explain_context -> result_back_junctions[back_i][xk1].abs_offset_for_start - section_length;

			explain_context -> result_back_junctions[back_i][xk1].abs_offset_for_start = new_start_pos;
			if(explain_context -> result_back_junctions[back_i][xk1].event_after_section
				&& explain_context -> result_back_junctions[back_i][xk1].event_after_section->is_strand_jumped) is_first_section_negative=!is_first_section_negative;
		}

		// build CIGAR
		for(front_i = 0; front_i < explain_context -> all_front_alignments; front_i++){
			if(final_alignment_number >= MAX_ALIGNMENT_PER_ANCHOR)break;

			to_be_supported_count = 0;
			tmp_cigar[0]=0;
			int known_junction_supp = 0;

			for(xk1 = 0; xk1 < explain_context -> result_back_junction_numbers[back_i] + explain_context -> result_front_junction_numbers[front_i] -1; xk1++)
			{
				char piece_cigar[25];
				int read_pos_start, read_pos_end;
				perfect_section_in_read_t * current_section, *next_section = NULL;

				int is_front_search = 0;
				if(xk1 >= explain_context -> result_back_junction_numbers[back_i] - 1) {
					current_section = &explain_context -> result_front_junctions[front_i][xk1 - explain_context -> result_back_junction_numbers[back_i] +1];
					if(xk1 - explain_context -> result_back_junction_numbers[back_i] +2 < explain_context -> result_front_junction_numbers[front_i])
						next_section = &explain_context -> result_front_junctions[front_i][xk1 - explain_context -> result_back_junction_numbers[back_i] +2];
					is_front_search = 1;
				} else {
					current_section = &explain_context -> result_back_junctions[back_i][xk1];
					if(xk1+1 <  explain_context ->  result_back_junction_numbers[back_i])
						next_section = &explain_context -> result_back_junctions[back_i][xk1+1];
				}


				if(xk1 == explain_context -> result_back_junction_numbers[back_i] - 1)
				     read_pos_start = explain_context -> result_back_junctions[back_i][xk1].read_pos_start;
				else read_pos_start = current_section -> read_pos_start;

				read_pos_end = current_section -> read_pos_end;
				chromosome_event_t *event_after = current_section -> event_after_section;

				sprintf(piece_cigar, "%dM", (read_pos_end - read_pos_start));
				total_perfect_matched_sections += (read_pos_end - read_pos_start);
				flanking_size_left[xk1] = (read_pos_end - read_pos_start);

				if(xk1<explain_context ->  result_back_junction_numbers[back_i] + explain_context ->  result_front_junction_numbers[front_i]  -2)
					assert(event_after);

				if(xk1>0)
					flanking_size_right[xk1-1] = (read_pos_end - read_pos_start);

				if(event_after)
				{
					if(event_after -> event_type == CHRO_EVENT_TYPE_INDEL)
						sprintf(piece_cigar+strlen(piece_cigar), "%d%c", abs(event_after->indel_length), event_after->indel_length>0?'D':'I');
					else if(event_after -> event_type == CHRO_EVENT_TYPE_JUNCTION||event_after -> event_type == CHRO_EVENT_TYPE_FUSION) {
						// the distance in CIGAR is the NEXT UNWANTED BASE of piece#1 to the FIRST WANTED BASE in piece#2
						int delta_one ;
						if(current_section -> is_strand_jumped + current_section -> is_connected_to_large_side == 1) delta_one = 1;
						else delta_one = -1;

						// if it is from front_search, the event side points to the first WANTED base of the next section; it should be moved to the last WANTED base the next section if the next section is jumped.
						if(next_section && (event_after -> is_strand_jumped + current_section -> is_strand_jumped==1))
						{
							if(is_front_search)
							{
								if(current_section -> is_connected_to_large_side)
									delta_one += (next_section->read_pos_end - next_section-> read_pos_start - 1);
								else
									delta_one -= (next_section->read_pos_end - next_section-> read_pos_start - 1);
							}
							else
							{
								if(current_section -> is_connected_to_large_side)
									delta_one += (next_section->read_pos_end - next_section-> read_pos_start - 1);
								else
									delta_one -= (next_section->read_pos_end - next_section-> read_pos_start - 1);
							}
						}
						
						char jump_mode = current_section -> is_connected_to_large_side?'B':'N';
						long long int movement = event_after -> event_large_side;
						movement -= event_after -> event_small_side - delta_one;
						if(1){
							if(jump_mode == 'B' && movement < 0){
								movement = - movement;
								jump_mode = 'N';
							}else if(jump_mode == 'N' && movement < 0){
								movement = - movement;
								jump_mode = 'B';
							}
						}
						
						if(event_after -> is_strand_jumped) jump_mode = tolower(jump_mode);
						fusions_in_read += (event_after -> event_type == CHRO_EVENT_TYPE_FUSION);
						sprintf(piece_cigar+strlen(piece_cigar), "%u%c", (int)movement, jump_mode);
						
						if(event_after -> indel_at_junction) sprintf(piece_cigar+strlen(piece_cigar), "%dI", event_after -> indel_at_junction);
						is_junction_read ++;
						if(event_after -> is_donor_found_or_annotation & 64 ) known_junction_supp ++;
					}
					to_be_supported[to_be_supported_count++] = event_after;
				}
				strcat(tmp_cigar, piece_cigar);
				if(strlen(tmp_cigar) > CORE_MAX_CIGAR_STR_LEN - 14){
					is_cigar_overflow=1;
					break;
				}
			}

			int mismatch_bases = 0;

			//#warning ">>>>>>>>>>>>>>>> COMMENT NEXT LINE <<<<<<<<<<<<<<<<<<<<<<<"
			//SUBREADprintf("ReadDebug:%s\t%s\n", explain_context -> read_name , tmp_cigar);
			if(is_cigar_overflow) sprintf(tmp_cigar, "%dM",  explain_context -> full_read_len);

			unsigned int final_position;

//			#warning "'0 &&' is because there could be indels in the high-confident region but this indel is finally disused."
			if( 0 && explain_context -> result_back_junction_numbers[back_i] + explain_context -> result_front_junction_numbers[front_i] <= 2) final_position = result -> selected_position;
			else final_position = explain_context -> result_back_junctions[back_i][0].abs_offset_for_start;

			if(0 && FIXLENstrcmp("simulated.11420793", explain_context->read_name)==0)SUBREADprintf("FFFINAL %s : POS=%u, ABS=%u\n", explain_context->read_name, final_position, explain_context -> result_back_junctions[back_i][0].abs_offset_for_start);

			int is_exonic_read_fraction_OK = 1;

			if( global_context -> config.minimum_exonic_subread_fraction > 0.0000001 && (!is_junction_read) && result -> used_subreads_in_vote>0)
			{
				int min_subreads = global_context -> config.minimum_exonic_subread_fraction * result-> used_subreads_in_vote; 
				if( result -> selected_votes < min_subreads )
					is_exonic_read_fraction_OK = 0 ;
			}



			int final_qual = 0, applied_mismatch = 0, non_clipped_length = 0, total_indel_length = 0, final_MATCH = 0, chromosomal_length = 0, full_section_clipped = 0;

			if(is_exonic_read_fraction_OK)
			{
				final_qual  = final_CIGAR_quality(global_context, thread_context, explain_context -> full_read_text, explain_context -> full_qual_text, explain_context -> full_read_len , tmp_cigar, final_position, is_first_section_negative != ((result->result_flags & CORE_IS_NEGATIVE_STRAND)?1:0), &mismatch_bases, result -> confident_coverage_start, result -> confident_coverage_end,  explain_context -> read_name, &non_clipped_length, &total_indel_length, & final_MATCH, & chromosomal_length, & full_section_clipped);
				//#warning ">>>>>>> COMMENT THIS <<<<<<<"
				//printf("OCT27-STEP2-%s:%d-POS%u-VOT%d-CIG-%s [ %d ]-INDELs=%llu; M/MM=%d,%d\n", explain_context -> read_name, explain_context  -> is_second_read + 1,  result -> selected_position, result -> selected_votes, tmp_cigar, is_cigar_overflow, ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_entry_table -> numOfElements, final_MATCH, mismatch_bases);


				applied_mismatch = is_junction_read? global_context->config.max_mismatch_junction_reads:global_context->config.max_mismatch_exonic_reads ;
				if(explain_context->full_read_len > EXON_LONG_READ_LENGTH)
					applied_mismatch = ((((explain_context->full_read_len+1)<<16) / 100) * applied_mismatch)>>16;

				if(global_context -> config.space_type == GENE_SPACE_COLOR) applied_mismatch += to_be_supported_count*2;
			}


			//#warning " ========== COMMENT THIS LINE !! ========="
			if(0 && FIXLENstrcmp("HWI-ST945:119:D0J2JACXX:1:1303:17374:199067", explain_context -> read_name) ==0){
				char outpos1[100];
				absoffset_to_posstr(global_context, final_position, outpos1);
				SUBREADprintf("FINALQUAL %s : FINAL_POS=%s ( %u )\tCIGAR=%s\tMM=%d / MAPLEN=%d > %d?\tVOTE=%d > %0.2f x %d ?  MASK=%d\tQUAL=%d\tBRNO=%d\nKNOWN_JUNCS=%d PENALTY=%d\n\n", explain_context -> read_name, outpos1 , final_position , tmp_cigar, mismatch_bases, non_clipped_length, applied_mismatch,  result -> selected_votes, global_context -> config.minimum_exonic_subread_fraction,result-> used_subreads_in_vote, result->result_flags, final_qual, explain_context -> best_read_id, known_junction_supp, explain_context -> best_indel_penalty);
				//exit(0);
			}


			if(mismatch_bases <= applied_mismatch && is_exonic_read_fraction_OK && fusions_in_read < 2 ){// && (0 == full_section_clipped || 0 == global_context -> config.do_breakpoint_detection)) {
				realignment_result_t * realign_res = final_realignments+final_alignment_number;
				final_alignment_number ++;

				realign_res -> realign_flags = result->result_flags;
				realign_res -> first_base_is_jumpped = 0;
				realign_res -> mapping_result = result;
				realign_res -> chromosomal_length = chromosomal_length;
				realign_res -> known_junction_supp = known_junction_supp;
				realign_res -> final_penalty = explain_context -> best_indel_penalty;

				if(mismatch_bases >  applied_mismatch ) realign_res -> realign_flags |= CORE_TOO_MANY_MISMATCHES;
				else realign_res -> realign_flags &= ~CORE_TOO_MANY_MISMATCHES;

				if(((result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0) != is_first_section_negative)
				{
					assert((global_context-> config.do_fusion_detection || global_context-> config.do_long_del_detection));
					realign_res -> first_base_is_jumpped = 1;
				}
				strcpy(realign_res -> cigar_string, tmp_cigar);

				if(1)
				{
					int is_RNA_from_positive = -1;
					unsigned long long read_id = 2llu * explain_context ->  pair_number + explain_context->is_second_read; 

					for(xk1= 0; xk1 < to_be_supported_count; xk1++)
					{
						if(xk1 >= MAX_EVENTS_IN_READ) break;

						if(to_be_supported [xk1] -> event_type !=CHRO_EVENT_TYPE_INDEL && is_junction_read){
							if(to_be_supported [xk1] -> event_type == CHRO_EVENT_TYPE_JUNCTION && to_be_supported [xk1] -> is_donor_found_or_annotation && is_RNA_from_positive == -1)
								is_RNA_from_positive = !(to_be_supported [xk1] -> is_negative_strand);
						}

						//final counts are added in function "add_realignment_event_support" in core.c
						
						realign_res -> supporting_chromosome_events[xk1] = to_be_supported[xk1];
						realign_res -> flanking_size_left[xk1] = flanking_size_left[xk1];
						realign_res -> flanking_size_right[xk1] = flanking_size_right[xk1];
						realign_res -> crirical_support[xk1] += (read_id == to_be_supported [xk1] -> critical_read_id);
					}
					if(to_be_supported_count < MAX_EVENTS_IN_READ ) 
						realign_res -> supporting_chromosome_events[to_be_supported_count] = NULL;
					
					result -> result_flags |= CORE_IS_FULLY_EXPLAINED;
					result -> read_length = explain_context->full_read_len;

					if(is_RNA_from_positive == -1)
					{
						realign_res -> realign_flags |= CORE_NOTFOUND_DONORS ;
						realign_res -> realign_flags &= ~(CORE_IS_GT_AG_DONORS);
					}
					else
					{
						realign_res -> realign_flags &= ~ (CORE_NOTFOUND_DONORS | CORE_IS_GT_AG_DONORS);

						if(is_RNA_from_positive)
							realign_res -> realign_flags |= CORE_IS_GT_AG_DONORS;
					}
				}

				realign_res -> first_base_position = final_position;
				realign_res -> final_quality = final_qual;
				realign_res -> final_mismatched_bases = mismatch_bases;
				realign_res -> final_matched_bases = (unsigned short)final_MATCH;
				realign_res -> best_second_diff_bases = (9<explain_context -> best_second_match_diff)?-1:explain_context -> best_second_match_diff; 

			}
		}
	}

	return final_alignment_number;
}




#define ceq(c,t) ((c)[0]==(t)[0] && (c)[1]==(t)[1])
#define c2eq(ch1, ch2, tg1, tg2) ((ceq(ch1, tg1) && ceq(ch2, tg2)) || (ceq(ch1, tg2) && ceq(ch2, tg1)) )

int paired_chars_full_core(char * ch1, char * ch2, int is_reverse)
{
	if (c2eq(ch1, ch2, "GT", "AG") || c2eq(ch1, ch2, "CT", "AC"))
	{
		if (is_reverse) if (ceq(ch1, "AG") || ceq(ch1, "AC")) return 2;
		if (!is_reverse) if (ceq(ch1, "CT") || ceq(ch1, "GT")) return 2;
	}
	else if ( c2eq(ch1, ch2,"GC","AG") || c2eq(ch1, ch2,"GC","CT") || c2eq(ch1, ch2,"AT","AC") || c2eq(ch1, ch2,"GT","AT"))
	{
		if (is_reverse) if (ceq(ch1, "GC") || ceq(ch1, "AT")  || ceq(ch1, "AG") || ceq(ch1, "AC")) return 1;
		if (!is_reverse) if (ceq(ch1, "GC") || ceq(ch1, "AT") ||ceq(ch1, "GT") || ceq(ch1, "CT")) return 1;
	}
	return 0;
}

int paired_chars_part_core(char * ch1, char * ch2, int is_reverse)
{
	if (c2eq(ch1, ch2, "GT", "AG") || c2eq(ch1, ch2, "CT", "AC")) {
		if (is_reverse){
			if (ceq(ch1, "AG") || ceq(ch1, "AC")) return 1;
		} else {
			if (ceq(ch1, "CT") || ceq(ch1, "GT")) return 1;
		}
	}
	return 0;
}

#define is_donor_chars_full(cc) (((cc)[0]=='G' && (cc)[1]=='T') || \
			    ((cc)[0]=='A' && (cc)[1]=='G') || \
			    ((cc)[0]=='A' && (cc)[1]=='C') || \
			    ((cc)[0]=='C' && (cc)[1]=='T') || \
			    ((cc)[0]=='G' && (cc)[1]=='C') || \
			    ((cc)[0]=='A' && (cc)[1]=='T') || \
			    ((cc)[0]=='A' && (cc)[1]=='C') ) 


#define is_donor_chars_part(cc) (((cc)[0]=='G' && (cc)[1]=='T') || \
			    ((cc)[0]=='A' && (cc)[1]=='G') || \
			    ((cc)[0]=='A' && (cc)[1]=='C') || \
			    ((cc)[0]=='C' && (cc)[1]=='T')) 

//#warning "=============== NO DONOR-RECEPTOR NEEDED =============="
//#define is_donor_chars(x) 1
//#define paired_chars(x,y,z) 1

#define is_donor_chars is_donor_chars_part
#define  paired_chars paired_chars_part_core




void print_big_margin(global_context_t * global_context, subread_read_number_t pair_number, int is_second_read){
	unsigned short * big_margin_record = _global_retrieve_big_margin_ptr(global_context,pair_number, is_second_read);
	int x1;

	SUBREADprintf("\n  >>> READ_NO=%u,  SECOND=%d, MEM=%p <<< \n", (unsigned int)pair_number, is_second_read, big_margin_record);
	for(x1 = 0; x1 < global_context->config.big_margin_record_size/3 ; x1++)
	{
		SUBREADprintf("%d %d~%d   ", big_margin_record[x1*3] , big_margin_record[x1*3+1] , big_margin_record[x1*3+2]);
	}
	SUBREADputs("");
}

#define ABGIGUOUS_TOLERANCE 3

int is_ambiguous_voting(global_context_t * global_context, subread_read_number_t pair_number, int is_second_read, int selected_vote, int max_start,int max_end, int read_len, int is_negative)
{
//	#warning "=========== THE NEXT LINE IS ONLY FOR COMPARING WITH STAR!! ============== "
//	return 0;
	if( global_context->config.big_margin_record_size<3) return 0;
	int xk1;
	int encounter = 0;

	if(is_negative)
	{
		int tmp = max_start;
		max_start = read_len - max_end;
		max_end = read_len - tmp;
	}

	unsigned short * big_margin_record = _global_retrieve_big_margin_ptr(global_context,pair_number, is_second_read);

	for(xk1 = 0; xk1 < global_context->config.big_margin_record_size/3 ; xk1++)
	{
		if(!big_margin_record[xk1*3])break;

		if(big_margin_record[xk1*3] >= selected_vote - 1)	// actually, max-1
		{
			if(0) {
				if ( max_start >= big_margin_record[xk1*3+1] - ABGIGUOUS_TOLERANCE && max_end <= big_margin_record[xk1*3+2] + ABGIGUOUS_TOLERANCE )
					encounter++;
				else if ( big_margin_record[xk1*3+1] >= max_start - ABGIGUOUS_TOLERANCE && big_margin_record[xk1*3+2] <= max_end + ABGIGUOUS_TOLERANCE )
					encounter++;

			} else {
			// 4 and 4 are the best setting for indel and fusion simulation.
				if(selected_vote >= big_margin_record[xk1*3]) {
					if(big_margin_record[xk1*3+1] >= max_start - 4 && big_margin_record[xk1*3+2] <= max_end + 4)
						encounter++;
				} else {
					if(big_margin_record[xk1*3+1] <= max_start + 4 && big_margin_record[xk1*3+2] >= max_end - 4)
						encounter++;
				}
			}
		}

	}

	if(encounter>1) return encounter;
	return 0;
}

#define JUNCTION_CONFIRM_WINDOW 17
// This function implements the same function of donor_score, except that the two halves are from different strands.
// Both halves are forced to positive strand and the split point is found.
// Note that the donor/receptor sides are still expected for distinguishing between Fusion Breaks and Fusion Junctions.

// Note that the read_text is on reversed mode. The guess points are on reversed mode too.
// "Left" and "Right" means the left/right half in the "reversed" read.
int donor_jumped_score(global_context_t * global_context, thread_context_t * thread_context, unsigned int small_virtualHead_abs_offset, unsigned int large_virtualHead_abs_offset, int guess_start, int guess_end,  char * read_text, int read_len, int is_small_half_negative, int is_large_half_negative, int small_half_on_left_reversed, int * final_split_point, int * is_GT_AG_strand, int * is_donor_found_or_annotation, int * small_side_increasing_coordinate, int * large_side_increasing_coordinate)
{
	gene_value_index_t * value_index = thread_context?thread_context->current_value_index:global_context->current_value_index ;
	// guess_end is the index of the first UNWANTED BASE.
	int most_likely_point_as_reversed = (guess_start+guess_end)/2;

	int selected_real_split_point = -1, selected_junction_strand = -1;
	//char donor_left[2], donor_right[2];
 
	int best_score = -111111;

	int real_split_point_i;
	int real_split_point_numbers = guess_end - guess_start;

	char positive_read[MAX_READ_LENGTH+1];
	strcpy(positive_read, read_text) ;
	reverse_read(positive_read, read_len, global_context->config.space_type);

	//printf("TEST_JUMPED: %u - %u\n", small_virtualHead_abs_offset, large_virtualHead_abs_offset);
	

	(*small_side_increasing_coordinate) = (small_half_on_left_reversed != is_small_half_negative);
	(*large_side_increasing_coordinate) = (small_half_on_left_reversed == is_large_half_negative);


	for(real_split_point_i = 0 ; real_split_point_i < real_split_point_numbers; real_split_point_i++)
	{
		int left_should_match, right_should_match;
		int left_should_not_match, right_should_not_match;
		int real_split_point_as_reversed = (real_split_point_i % 2)?-((real_split_point_i+1)/2):((1+real_split_point_i)/2);
		real_split_point_as_reversed += most_likely_point_as_reversed;

		if(real_split_point_as_reversed > read_len-JUNCTION_CONFIRM_WINDOW)continue;
		if(real_split_point_as_reversed < JUNCTION_CONFIRM_WINDOW)continue;

		int is_donor_test_ok=0;

		if(small_half_on_left_reversed)
		{
			unsigned int small_pos_test_begin = small_virtualHead_abs_offset + (is_small_half_negative?real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW:(read_len - real_split_point_as_reversed)); 
			char * small_pos_read_begin = (is_small_half_negative?read_text:positive_read) + (is_small_half_negative?
						(real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW)           :
						(read_len - real_split_point_as_reversed)
  						);

			unsigned int large_pos_test_begin = large_virtualHead_abs_offset + (is_large_half_negative?real_split_point_as_reversed:(read_len - real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW));
			char * large_pos_read_begin = (is_large_half_negative?read_text:positive_read) + (is_large_half_negative?
						(real_split_point_as_reversed)     :
						(read_len - real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW));

			left_should_match = match_chro(small_pos_read_begin , value_index , small_pos_test_begin , JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);
			right_should_match = match_chro(large_pos_read_begin , value_index , large_pos_test_begin , JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);
			left_should_not_match = right_should_not_match = 0;
		//match_chro(read_text + real_split_point - JUNCTION_CONFIRM_WINDOW, value_index, small_virtualHead_abs_offset + real_split_point - JUNCTION_CONFIRM_WINDOW , JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);

		}
		else
		{
			unsigned int small_pos_test_begin = small_virtualHead_abs_offset + (is_small_half_negative?real_split_point_as_reversed:(read_len - real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW)); 
			char * small_pos_read_begin = (is_small_half_negative?read_text:positive_read) + (is_small_half_negative?
							(real_split_point_as_reversed):(read_len - real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW));

			unsigned int large_pos_test_begin = large_virtualHead_abs_offset + (is_large_half_negative?(real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW):(read_len - real_split_point_as_reversed));
			char * large_pos_read_begin = (is_large_half_negative?read_text:positive_read) + (is_large_half_negative?
							  (real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW):(read_len - real_split_point_as_reversed));

			left_should_match = match_chro(small_pos_read_begin , value_index , small_pos_test_begin , JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);
			right_should_match = match_chro(large_pos_read_begin , value_index , large_pos_test_begin , JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);
			left_should_not_match = right_should_not_match = 0;

		}

		//#warning "============ REMOVE THE TWO '+ 1' FROM THE NEXT LINE ================="
		//#warning "============ ADD THE TWO '+ 1's IN THE BLANKETS FOR SVs GRANT APP ================="
		int mismatch_in_between_allowd = (global_context -> config.more_accurate_fusions)?(0):(1);
		if(left_should_match + right_should_match  >= JUNCTION_CONFIRM_WINDOW*2 - mismatch_in_between_allowd  &&
			left_should_not_match <= JUNCTION_CONFIRM_WINDOW -3 && right_should_not_match <= JUNCTION_CONFIRM_WINDOW -3)
		{
			int test_score = is_donor_test_ok*500+left_should_match + right_should_match - left_should_not_match - right_should_not_match;
			if(test_score > best_score)
			{
				selected_real_split_point = real_split_point_as_reversed;
				best_score = test_score;
			}
		}
	}

	if(best_score>0)
	{
		//printf("TEST_JUMPED: BSCORE=%d  SPLT=%d\n", best_score , selected_real_split_point);
		*final_split_point = selected_real_split_point;
		*is_donor_found_or_annotation = best_score>=500;
		*is_GT_AG_strand = selected_junction_strand;
		return best_score;
	}
	return 0;
}


int donor_score(global_context_t * global_context, thread_context_t * thread_context, unsigned int left_virtualHead_abs_offset, unsigned int right_virtualHead_abs_offset, int left_indel_offset, int right_indel_offset, int normally_arranged, int guess_start, int guess_end,  char * read_text, int read_len, int * final_split_point, int * is_GT_AG_strand, int * is_donor_found_or_annotation, int * final_inserted_bases, int * small_side_increasing_coordinate, int * large_side_increasing_coordinate, char * read_name)
{


	gene_value_index_t * value_index = thread_context?thread_context->current_value_index:global_context->current_value_index;
	int need_donor_test = global_context->config.do_breakpoint_detection && global_context -> config.check_donor_at_junctions && (!(  global_context-> config.do_fusion_detection ||   global_context-> config.do_long_del_detection));

	(*small_side_increasing_coordinate)=!normally_arranged;
	(*large_side_increasing_coordinate)= normally_arranged;
	
	// guess_end is the index of the first UNWANTED BASE.
	int most_likely_point = (guess_start+guess_end)/2;
	
	// "split_point" is the first base NOT IN piece 1; it is also the first base IN piece 2. 
	int selected_real_split_point = -1, selected_junction_strand = -1 , selected_inserted_bases = 0;
	char donor_left[3], donor_right[3];
	
 
	int best_score = -111111;
	int non_insertion_preferred = 0;

	int real_split_point_i;
	int real_split_point_numbers = guess_end - guess_start;

	if(0 && FIXLENstrcmp("R006856515", read_name) == 0) 
		SUBREADprintf("TESTDON: LR=%d; RR=%d\n", left_indel_offset, right_indel_offset);
	
	for(real_split_point_i = 0 ; real_split_point_i < real_split_point_numbers; real_split_point_i++)
	{
		int left_should_match, right_should_match = 0;
		int left_should_not_match = 0, right_should_not_match = 0;
		int real_split_point = (real_split_point_i % 2)?-((real_split_point_i+1)/2):((1+real_split_point_i)/2);
		real_split_point += most_likely_point;
		int is_donor_test_ok = 0;

		if(real_split_point > read_len-JUNCTION_CONFIRM_WINDOW)continue;
		if(real_split_point < JUNCTION_CONFIRM_WINDOW)continue;

		if(global_context->config.prefer_donor_receptor_junctions)
		{
			if(normally_arranged)
			{
				gvindex_get_string (donor_left, value_index, left_virtualHead_abs_offset + real_split_point + left_indel_offset, 2, 0);
				if(is_donor_chars(donor_left))
				{
					gvindex_get_string (donor_right, value_index, right_virtualHead_abs_offset + real_split_point + right_indel_offset - 2, 2, 0);
					if(is_donor_chars(donor_right))
					{
						is_donor_test_ok = paired_chars(donor_left, donor_right,0);
					}
				}
			}
			else
			{
				gvindex_get_string (donor_left, value_index, right_virtualHead_abs_offset + real_split_point + left_indel_offset, 2, 0);
				gvindex_get_string (donor_right, value_index, left_virtualHead_abs_offset + real_split_point + right_indel_offset - 2, 2, 0);
				is_donor_test_ok = is_donor_chars(donor_left) && is_donor_chars(donor_right) && paired_chars(donor_left, donor_right,0);
			}
		}

	//	donor_left[2]=0; donor_right[2]=0;

		if(0 && FIXLENstrcmp("R006856515", read_name) == 0) 
		{
			donor_left[2]=0;
			donor_right[2]=0;
			SUBREADprintf("TESTDON: %s %s; OFFSET=%d; DON_OK=%d; NORMAL=%d; LEFT_OFF=%d; RIGHT_OFF=%d\n", donor_left, donor_right, real_split_point_i, is_donor_test_ok, normally_arranged, left_indel_offset, right_indel_offset);
		}

		//#warning "============ REMOVE THE TWO '+ 1' FROM THE NEXT LINE ================="
		//#warning "============ ADD TWO '+ 1' IN THE BLANKETS FOR SVs GRANT APP ================="
		int mismatch_in_between_allowd = (global_context -> config.more_accurate_fusions)?(0) : (1);
		if(is_donor_test_ok || !need_donor_test)
		{
			if(normally_arranged)
			{
				int inserted_bases=0;

				left_should_match = match_chro(read_text + real_split_point - JUNCTION_CONFIRM_WINDOW, value_index, left_virtualHead_abs_offset + real_split_point - JUNCTION_CONFIRM_WINDOW + left_indel_offset , JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);	
				//printf("INS=%d; LM=%d\t\tLOL=%u, LOR=%u, SP=%d\n", inserted_bases, left_should_match, left_virtualHead_abs_offset, right_virtualHead_abs_offset, real_split_point);
				if(left_should_match > JUNCTION_CONFIRM_WINDOW- (global_context->config.max_insertion_at_junctions?5:2))
				{
					for(inserted_bases = 0; inserted_bases <= global_context->config.max_insertion_at_junctions; inserted_bases++)
					{

						right_should_match = match_chro(read_text + real_split_point + inserted_bases, value_index, right_virtualHead_abs_offset + real_split_point + right_indel_offset + inserted_bases, JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);	
					//	printf("INS=%d; LM=%d; RM=%d\t\tLOL=%u, LOR=%u, SP=%d\n", inserted_bases, left_should_match, right_should_match, left_virtualHead_abs_offset, right_virtualHead_abs_offset, real_split_point);
						if(right_should_match >= 2*JUNCTION_CONFIRM_WINDOW - left_should_match - mismatch_in_between_allowd)
						{
							left_should_not_match = match_chro(read_text + real_split_point + inserted_bases, value_index, left_virtualHead_abs_offset + real_split_point + left_indel_offset, JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);	
							right_should_not_match = match_chro(read_text + real_split_point - JUNCTION_CONFIRM_WINDOW, value_index, right_virtualHead_abs_offset  + real_split_point + right_indel_offset - JUNCTION_CONFIRM_WINDOW + inserted_bases, JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);	


							if(left_should_not_match <= JUNCTION_CONFIRM_WINDOW -5 && right_should_not_match <= JUNCTION_CONFIRM_WINDOW -5)
							{
								int test_score ; 
								if(global_context->config.max_insertion_at_junctions)
									test_score = 100*(is_donor_test_ok*3000+left_should_match + right_should_match) - (left_should_not_match + right_should_not_match) - 20*inserted_bases;
								else
									test_score = 100*(is_donor_test_ok*3000+left_should_match + right_should_match - left_should_not_match - right_should_not_match);

								if(test_score > best_score)
								{
									//if(left_virtualHead_abs_offset > 2729745284 - 200 && left_virtualHead_abs_offset< 2729745284 + 200)
									//	SUBREADprintf("INS=%d; BSS=%d; TSC=%d\n%s\n\n", inserted_bases , best_score, test_score, read_text);
									selected_junction_strand = (donor_left[0]=='G' || donor_right[1]=='G');
									selected_inserted_bases = inserted_bases;
									selected_real_split_point = real_split_point;	
									best_score = test_score;
								}
							}

						}
						if(global_context->config.max_insertion_at_junctions && 0 == inserted_bases && right_should_match >= 2*JUNCTION_CONFIRM_WINDOW - left_should_match - 5)
							non_insertion_preferred = 1;

					}
				}
			}
			else
			{
				right_should_match = match_chro(read_text + real_split_point - JUNCTION_CONFIRM_WINDOW, value_index, right_virtualHead_abs_offset + right_indel_offset + real_split_point - JUNCTION_CONFIRM_WINDOW , JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);
				left_should_match = match_chro(read_text + real_split_point, value_index, left_virtualHead_abs_offset + real_split_point + left_indel_offset, JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);	

				right_should_not_match = match_chro(read_text + real_split_point, value_index, right_virtualHead_abs_offset + real_split_point + right_indel_offset, JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);	
				left_should_not_match = match_chro(read_text + real_split_point - JUNCTION_CONFIRM_WINDOW, value_index, left_virtualHead_abs_offset + left_indel_offset + real_split_point - JUNCTION_CONFIRM_WINDOW, JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);	
	
				//printf("LEFT:MA=%d UMA=%d     RIGHT:MA=%d UMA=%d\n", left_should_match, left_should_not_match, right_should_match, right_should_not_match);

				if(left_should_match +right_should_match >= 2*JUNCTION_CONFIRM_WINDOW - mismatch_in_between_allowd && 
					left_should_not_match <= JUNCTION_CONFIRM_WINDOW -5 && right_should_not_match <= JUNCTION_CONFIRM_WINDOW -5)
				{
					
					int test_score;

					test_score = 100*(is_donor_test_ok*3000+left_should_match + right_should_match - left_should_not_match - right_should_not_match);
					if(test_score > best_score)
					{
						selected_junction_strand = (donor_left[0]=='G' || donor_right[1]=='G');
						selected_real_split_point = real_split_point;	
						best_score = test_score;
					}
				}
			}
		}
	}
	if(best_score>0 && (0==non_insertion_preferred || 0==selected_inserted_bases))
	{
		*final_split_point = selected_real_split_point;
		*is_donor_found_or_annotation = best_score>=290000;
		*is_GT_AG_strand = selected_junction_strand;
		*final_inserted_bases = selected_inserted_bases;

		if(0 && FIXLENstrcmp("R006856515", read_name)==0)
			SUBREADprintf("FINAL_INS_LEN=%d; BEST_SCORE=%d  %s\n", selected_inserted_bases, best_score, read_name);
		return (1+best_score)/100;
	}
	return 0;

}

void find_new_junctions(global_context_t * global_context, thread_context_t * thread_context, subread_read_number_t pair_number, char * read_name, char * read_text, char * qual_text, int read_len, int is_second_read, int best_read_id)
{
	mapping_result_t * result =_global_retrieve_alignment_ptr(global_context, pair_number, is_second_read, best_read_id);
	subjunc_result_t * subjunc_result =_global_retrieve_subjunc_ptr(global_context, pair_number, is_second_read, best_read_id);

	if(read_len > EXON_LONG_READ_LENGTH)
	{
		assert(result -> selected_position <= 0xffff0000);
		core_search_short_exons(global_context, thread_context,  read_text, qual_text, read_len, result -> selected_position, (subjunc_result -> minor_votes < 1)? result -> selected_position:subjunc_result -> minor_position, result -> confident_coverage_start, result -> confident_coverage_end);
	}

	int selected_real_split_point = subjunc_result->split_point;

	//#warning " =============== remove "+ 2" FROM THE NEXT LINE (FOR A HIGHER ACCURACY FROM SubFusion on 19 JAN 2015)  =================="
	if((global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection) && subjunc_result -> minor_votes < 1)return;
	if((!(global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection)) && subjunc_result -> minor_votes < 1)return;

	//if(result -> selected_votes < global_context->config.minimum_subread_for_first_read)return;

	if(global_context->config.do_big_margin_filtering_for_junctions)
		if(is_ambiguous_voting(global_context, pair_number, is_second_read, result->selected_votes, result -> confident_coverage_start, result -> confident_coverage_end, read_len, (result->result_flags & CORE_IS_NEGATIVE_STRAND)?1:0))return;

	unsigned int left_virtualHead_abs_offset = min(result -> selected_position, subjunc_result -> minor_position);
	unsigned int right_virtualHead_abs_offset = max(result -> selected_position, subjunc_result -> minor_position);

	int is_GT_AG_donors = result->result_flags & 0x3;
	int is_donor_found_or_annotation = is_GT_AG_donors<3;
	int is_strand_jumped = (result->result_flags & CORE_IS_STRAND_JUMPED)?1:0;

	if(selected_real_split_point>0)
	{
		unsigned int left_edge_wanted, right_edge_wanted;

		if(is_strand_jumped)
		{
			if(0){
		
				// note that splicing point and the coverage coordinates are "major negative" view.
				// recover the "negative view" splicing point location
				int S = (result->result_flags & CORE_IS_NEGATIVE_STRAND) ? selected_real_split_point : (read_len - selected_real_split_point);
				int Sbar = read_len - S;

				int is_abnormal_as_reversed = (subjunc_result->minor_coverage_start > result->confident_coverage_start) + (subjunc_result -> minor_position >  result -> selected_position) == 1;
				if(!(result->result_flags & CORE_IS_NEGATIVE_STRAND)) is_abnormal_as_reversed = !is_abnormal_as_reversed;
				int is_small_half_negative = ((result->result_flags & CORE_IS_NEGATIVE_STRAND)?1:0) + (subjunc_result->minor_position < result->selected_position) ==1;

				if(is_abnormal_as_reversed && is_small_half_negative)
				{
					left_edge_wanted = left_virtualHead_abs_offset + S;
					right_edge_wanted = right_virtualHead_abs_offset + Sbar;
				}
				else if(is_abnormal_as_reversed && !is_small_half_negative)
				{
					left_edge_wanted = left_virtualHead_abs_offset + Sbar - 1;
					right_edge_wanted = right_virtualHead_abs_offset + S - 1;
				}
				else if(!is_abnormal_as_reversed && is_small_half_negative)
				{
					left_edge_wanted = left_virtualHead_abs_offset + S - 1;
					right_edge_wanted = right_virtualHead_abs_offset + Sbar - 1;
				}
				else // if(!is_abnormal_as_reversed && !is_small_half_negative)
				{
					left_edge_wanted = left_virtualHead_abs_offset + Sbar;
					right_edge_wanted = right_virtualHead_abs_offset + S;
				}

				if(left_edge_wanted >= right_edge_wanted){
					SUBREADprintf("REVERSED NEW JUNC: %u ~ %u : ABN_REV=%d , SMALL_NEG=%d, LEFT_VH=%u, RIGHT_VH=%u, S/~S=%d/%d\n", left_edge_wanted, right_edge_wanted, is_abnormal_as_reversed, is_small_half_negative, left_virtualHead_abs_offset, right_virtualHead_abs_offset, S, Sbar);
				}

			}else{
				unsigned int major_half_smallest_coordinate, minor_half_smallest_coordinate;
				major_half_smallest_coordinate = result -> selected_position + selected_real_split_point;
				minor_half_smallest_coordinate = subjunc_result->minor_position + read_len - selected_real_split_point;
				left_edge_wanted = min(major_half_smallest_coordinate, minor_half_smallest_coordinate);
				right_edge_wanted = max(major_half_smallest_coordinate, minor_half_smallest_coordinate);
				int is_abnormal_as_reversed = (subjunc_result->minor_coverage_start > result->confident_coverage_start) + (minor_half_smallest_coordinate > major_half_smallest_coordinate) == 1;
				int is_small_half_negative = ((result->result_flags & CORE_IS_NEGATIVE_STRAND)?1:0) + (minor_half_smallest_coordinate < major_half_smallest_coordinate) ==1;
				if(!(result->result_flags & CORE_IS_NEGATIVE_STRAND)) is_abnormal_as_reversed = !is_abnormal_as_reversed;
				if(is_small_half_negative != is_abnormal_as_reversed)
				{
					left_edge_wanted -=1;	
					right_edge_wanted -=1;
				}
			}
		}
		else
		{
			int selected_real_split_point_for_left = selected_real_split_point;
			int selected_real_split_point_for_right = selected_real_split_point;
			if((subjunc_result->minor_coverage_start > result->confident_coverage_start) + (subjunc_result -> minor_position >  result -> selected_position) == 1) //abnormally arranged halves
				selected_real_split_point_for_right --;
			else	// normally arranged halves
				selected_real_split_point_for_left --;


			
			int minor_indel_offset = (subjunc_result->double_indel_offset & 0xf);
			int major_indel_offset = (subjunc_result->double_indel_offset >> 4) & 0xf;
			if(major_indel_offset>=8)major_indel_offset=-(16-major_indel_offset);
			//assert(minor_indel_offset==0);
			//assert(major_indel_offset==0);

			left_edge_wanted = left_virtualHead_abs_offset + selected_real_split_point_for_left + ((result -> selected_position > subjunc_result -> minor_position)?minor_indel_offset: major_indel_offset);
			right_edge_wanted = right_virtualHead_abs_offset + selected_real_split_point_for_right;
		}

		char * chro_name_left, *chro_name_right;
		int chro_pos_left,chro_pos_right;
			
		locate_gene_position( left_edge_wanted , &global_context -> chromosome_table, &chro_name_left, &chro_pos_left);
		locate_gene_position( right_edge_wanted , &global_context -> chromosome_table, &chro_name_right, &chro_pos_right);
		if((!( global_context-> config.do_fusion_detection ||  global_context-> config.do_long_del_detection) ) && chro_name_right!=chro_name_left) return;

		//insert event
		HashTable * event_table = NULL;
		chromosome_event_t * event_space = NULL;
		if(thread_context)
		{
			event_table = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
			event_space = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
		}
		else
		{
			event_table = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
			event_space = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
		}

		// note that selected_real_split_point is the first UNWANTED base after left half.
	
		//if(abs(left_edge_wanted-27286396) < 250 || abs(right_edge_wanted - 27286396)<250)
		if(0 && FIXLENstrcmp("R003738400", read_name) == 0) 
		{
			char leftpos[100], rightpos[100];
			absoffset_to_posstr(global_context, left_edge_wanted, leftpos);
			absoffset_to_posstr(global_context, right_edge_wanted, rightpos);
			SUBREADprintf("READ=%s, LEFT=%s, RIGHT=%s\n", read_name, leftpos, rightpos);
		}

		chromosome_event_t * found = NULL;
		chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];
		int found_events = search_event(global_context, event_table, event_space, left_edge_wanted , EVENT_SEARCH_BY_SMALL_SIDE,  CHRO_EVENT_TYPE_INDEL | CHRO_EVENT_TYPE_JUNCTION | CHRO_EVENT_TYPE_FUSION, search_return);

		mark_gapped_read(result);
		if(found_events)
		{
			int kx1; 
			for(kx1 = 0; kx1 < found_events ; kx1++)
			{
				if(search_return[kx1] -> event_large_side == right_edge_wanted)
				{
					found = search_return[kx1];	
					break;
				}
			}
		}

		//if( 1018082 == pair_number)
		//		SUBREADprintf("NEW_CHIMERISM_HERE [%u:%d: R_%d] : %s , %s , %u , %u, %c ; INC=%d %d\n", pair_number, best_read_id, is_second_read+1, chro_name_left, chro_name_right, chro_pos_left, chro_pos_right, is_strand_jumped?'X':'=', subjunc_result -> small_side_increasing_coordinate, subjunc_result -> large_side_increasing_coordinate);


		int is_key_fusion = 0;
		if(0)if(
			( 9566365 + 1210 - 200 <= left_edge_wanted && 9566365 + 1210 + 200 >= left_edge_wanted) && 
		    ( 36859887 + 1210 - 200  <= right_edge_wanted &&  36859887 + 1210 + 200 >= right_edge_wanted) 
		){
			SUBREADprintf("Read = %s, FOUND = %p in %d , %s:%u , %s:%u, INCs= %d, %d, JUMP=%d\n", read_name, found, found_events, chro_name_left, chro_pos_left, chro_name_right, chro_pos_right, subjunc_result -> small_side_increasing_coordinate, subjunc_result -> large_side_increasing_coordinate, is_strand_jumped);
			is_key_fusion = 1;
		}

		if(found) found -> supporting_reads ++;
		else
		{
			int event_no;


			if(thread_context)
				event_no = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> total_events ++;
			else
				event_no = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) ->  total_events ++;


			event_space = reallocate_event_space(global_context, thread_context, event_no);

			chromosome_event_t * new_event = event_space+event_no; 
			memset(new_event,0,sizeof(chromosome_event_t));
			new_event -> event_small_side = left_edge_wanted;
			new_event -> event_large_side = right_edge_wanted + subjunc_result->indel_at_junction;
			new_event -> critical_read_id = 2llu * pair_number + is_second_read;

			int new_event_type =(((global_context -> config.entry_program_name == CORE_PROGRAM_SUBJUNC && global_context ->  config.do_fusion_detection)||(global_context -> config.entry_program_name == CORE_PROGRAM_SUBJUNC && global_context ->  config.do_long_del_detection))&& !global_context -> config.prefer_donor_receptor_junctions)?CHRO_EVENT_TYPE_FUSION:CHRO_EVENT_TYPE_JUNCTION;

			//#warning "=========================== DELETE NEXT LINE !!! =================================="
			//new_event_type = CHRO_EVENT_TYPE_REMOVED;

			if(is_strand_jumped) new_event_type = CHRO_EVENT_TYPE_FUSION;
			if((subjunc_result->minor_coverage_start > result->confident_coverage_start) + (subjunc_result -> minor_position >  result -> selected_position) ==1)
				new_event_type = CHRO_EVENT_TYPE_FUSION;
			if(chro_name_right!=chro_name_left)
				new_event_type = CHRO_EVENT_TYPE_FUSION;
			if(right_edge_wanted - left_edge_wanted > global_context -> config.maximum_intron_length)
				if(!(global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection))
					new_event_type = CHRO_EVENT_TYPE_REMOVED;


			if(1)
			{
				unsigned int dist = new_event -> event_large_side -  new_event -> event_small_side;
				int origin_type = new_event_type;
				int fusion_cover_len = -1;

				if(dist > MAX_INSERTION_LENGTH && new_event_type == CHRO_EVENT_TYPE_FUSION)
				{
					int cov_end, cover_start, major_cov;
					cov_end = max(subjunc_result->minor_coverage_end, result->confident_coverage_end );
					cover_start = min(subjunc_result->minor_coverage_start, result->confident_coverage_start);

					major_cov =  result->confident_coverage_end  -  result->confident_coverage_start;

					fusion_cover_len = cov_end - cover_start ;

					if(fusion_cover_len < read_len - 15 || major_cov > read_len - 15)
						new_event_type = CHRO_EVENT_TYPE_REMOVED;
				}

				if(dist > MAX_INSERTION_LENGTH && new_event_type == CHRO_EVENT_TYPE_FUSION && subjunc_result -> minor_votes < 2)
					new_event_type = CHRO_EVENT_TYPE_REMOVED;
				else if(new_event_type == CHRO_EVENT_TYPE_FUSION && subjunc_result -> minor_votes < 1)
					new_event_type = CHRO_EVENT_TYPE_REMOVED;


				if(dist > MAX_INSERTION_LENGTH && new_event_type == CHRO_EVENT_TYPE_FUSION && result -> selected_votes < 2)
					new_event_type = CHRO_EVENT_TYPE_REMOVED;
				else if(new_event_type == CHRO_EVENT_TYPE_FUSION && result -> selected_votes < 1)
					new_event_type = CHRO_EVENT_TYPE_REMOVED;

				if(0 && origin_type == CHRO_EVENT_TYPE_FUSION)
				{
					char leftpos[100], rightpos[100];
					absoffset_to_posstr(global_context, new_event -> event_small_side, leftpos);
					absoffset_to_posstr(global_context, new_event -> event_large_side, rightpos);

					if(new_event_type == CHRO_EVENT_TYPE_REMOVED)
						SUBREADprintf("NEW_FUSION REMOVED %s SUGGEST %s ~ %s MAJOR COV=%d ~ %d, MINOR COV=%d ~ %d, RLEN=%d, COVED=%d, VOTES=%d, %d, %s, SPLIT=%d\n", read_name, leftpos, rightpos, result->confident_coverage_start, result->confident_coverage_end, subjunc_result->minor_coverage_start, subjunc_result->minor_coverage_end, read_len, fusion_cover_len, result -> selected_votes, subjunc_result -> minor_votes, is_strand_jumped?"JUMPED":"======", selected_real_split_point);
					else
						SUBREADprintf("NEW_FUSION WANTED %s SUGGEST %s ~ %s  MAJOR COV=%d ~ %d, MINOR COV=%d ~ %d, RLEN=%d, COVED=%d, VOTES=%d, %d, %s, SPLIT=%d\n", read_name, leftpos, rightpos, result->confident_coverage_start, result->confident_coverage_end, subjunc_result->minor_coverage_start, subjunc_result->minor_coverage_end, read_len, fusion_cover_len, result -> selected_votes, subjunc_result -> minor_votes, is_strand_jumped?"JUMPED":"======", selected_real_split_point);
				}

				if(dist > MAX_INSERTION_LENGTH && new_event_type == CHRO_EVENT_TYPE_FUSION && (selected_real_split_point < read_len * 0.2 || selected_real_split_point >= read_len *0.8000) )
					new_event_type = CHRO_EVENT_TYPE_REMOVED;
			}
		//if(pair_number == 13)
		//printf("MMMMX %d %u -- %u : TYPE %d\n" , event_no, left_edge_wanted, right_edge_wanted, new_event_type);


//			if((is_donor_found_or_annotation || !global_context -> config.check_donor_at_junctions) &&(!is_strand_jumped) && right_edge_wanted - left_edge_wanted <= global_context -> config.maximum_intron_length
//				&& (subjunc_result->minor_coverage_start > result->confident_coverage_start) + (subjunc_result -> minor_position >  result -> selected_position) !=1)


			if(is_key_fusion) SUBREADprintf("   INSERTED AS %d ( in %d or %d )\n", new_event_type , CHRO_EVENT_TYPE_JUNCTION, CHRO_EVENT_TYPE_FUSION);

			if(new_event_type == CHRO_EVENT_TYPE_JUNCTION)
			{
				new_event -> is_negative_strand= !is_GT_AG_donors;
				new_event -> event_type = CHRO_EVENT_TYPE_JUNCTION;

				new_event -> supporting_reads = 1;
				new_event -> indel_length = 0;
				new_event -> indel_at_junction = subjunc_result->indel_at_junction;
				new_event -> is_donor_found_or_annotation = is_donor_found_or_annotation; 

				new_event -> small_side_increasing_coordinate = subjunc_result -> small_side_increasing_coordinate;
				new_event -> large_side_increasing_coordinate = subjunc_result -> large_side_increasing_coordinate;
				
				put_new_event(event_table, new_event , event_no);
				
				if(0 && FIXLENstrcmp("R000000052", read_name) == 0) 
					SUBREADprintf("NEW_JUNCTION_HERE : %s , %u , %u  (%u, %u)\n", chro_name_right, chro_pos_left, chro_pos_right,  new_event -> event_small_side, new_event -> event_large_side);
			}
			else if(new_event_type == CHRO_EVENT_TYPE_FUSION)
			{
				if((global_context ->  config.do_fusion_detection || global_context ->  config.do_long_del_detection))
				{
					new_event -> event_type = CHRO_EVENT_TYPE_FUSION;
					new_event -> is_strand_jumped = is_strand_jumped;


					new_event -> supporting_reads = 1;
					new_event -> indel_length = 0;

					new_event -> small_side_increasing_coordinate = subjunc_result -> small_side_increasing_coordinate;
					new_event -> large_side_increasing_coordinate = subjunc_result -> large_side_increasing_coordinate;
					
					put_new_event(event_table, new_event , event_no);
					//if( 1018082 == pair_number)
					//	SUBREADprintf("NEW_CHIMERISM_HERE_FULL [%u:%d: R_%d] : %s , %s , %u , %u, %c ; INC=%d %d\n", pair_number, best_read_id, is_second_read+1, chro_name_left, chro_name_right, chro_pos_left, chro_pos_right, is_strand_jumped?'X':'=', new_event -> small_side_increasing_coordinate, new_event -> large_side_increasing_coordinate);
				}
			}
		}
	}
}

void write_translocation_results_final(void * key, void * buckv, HashTable * tab);
void write_inversion_results_final(void * key, void * buckv, HashTable * tab);

int write_fusion_final_results(global_context_t * global_context)
{
	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 
	char fn2 [MAX_FILE_NAME_LENGTH+30];

	sprintf(fn2,"%s.breakpoints.vcf", global_context->config.output_prefix);
	FILE * ofp = f_subr_open(fn2, "wb");
	fprintf(ofp,"##fileformat=VCFv4.1\n");
	fprintf(ofp,"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
	fprintf(ofp,"##INFO=<ID=MATEID,Number=1,Type=String,Description=\"Paired breakend id\">\n");
	fprintf(ofp,"##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Supporting read number\">\n");
	fprintf(ofp,"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n");

	int xk1, disk_is_full = 0;
	unsigned int all_junctions = 0;
	int no_sup_juncs = 0;
	int all_juncs = 0;

	for(xk1 = 0; xk1 < indel_context -> total_events ; xk1++)
	{ 
		char * chro_name_left,* chro_name_right;
		int chro_pos_left, chro_pos_right; 
		chromosome_event_t * event_body = indel_context -> event_space_dynamic +xk1;
		if(event_body -> event_type != CHRO_EVENT_TYPE_FUSION && (global_context->config.entry_program_name != CORE_PROGRAM_SUBREAD || event_body -> event_type != CHRO_EVENT_TYPE_JUNCTION))
			continue;

		all_juncs++;

		if(event_body->final_counted_reads<1|| event_body->critical_supporting_reads < 1 - 1)
		{
			no_sup_juncs++;
			continue;
		}
		locate_gene_position( event_body -> event_small_side , &global_context -> chromosome_table, &chro_name_left, &chro_pos_left);
		locate_gene_position( event_body -> event_large_side , &global_context -> chromosome_table, &chro_name_right, &chro_pos_right);

		chro_pos_left+=1;
		chro_pos_right+=1;
		all_junctions ++;

		int wlen;
		char alt_base[500];
		char ref_base;
		char bkt = event_body -> large_side_increasing_coordinate?'[':']';

		gene_value_index_t * current_index = find_current_value_index(global_context , event_body -> event_small_side , 1);
		ref_base = gvindex_get( current_index, event_body -> event_small_side);
		if(event_body -> small_side_increasing_coordinate)
			sprintf(alt_base,"%c%s:%u%c%c", bkt, chro_name_right, chro_pos_right, bkt, ref_base);
		else
			sprintf(alt_base,"%c%c%s:%u%c", ref_base, bkt, chro_name_right, chro_pos_right, bkt);

		wlen = fprintf(ofp,"%s\t%u\tbnd_%d\t%c\t%s\t.\tPASS\tSVTYPE=BND;MATEID=bnd_%d;SR=%d\n", chro_name_left, chro_pos_left, all_junctions *2 -1, ref_base, alt_base, all_junctions*2, event_body -> final_counted_reads);

		current_index = find_current_value_index(global_context , event_body -> event_large_side , 1);
		ref_base = gvindex_get( current_index, event_body -> event_large_side );
		bkt = event_body -> small_side_increasing_coordinate?'[':']';
		if(event_body -> large_side_increasing_coordinate)
			sprintf(alt_base,"%c%s:%u%c%c", bkt, chro_name_left, chro_pos_left, bkt, ref_base);
		else
			sprintf(alt_base,"%c%c%s:%u%c", ref_base, bkt, chro_name_left, chro_pos_left, bkt);

		wlen += fprintf(ofp,"%s\t%u\tbnd_%d\t%c\t%s\t.\tPASS\tSVTYPE=BND;MATEID=bnd_%d;SR=%d\n", chro_name_right, chro_pos_right, all_junctions *2, ref_base, alt_base, all_junctions*2 -1, event_body -> final_counted_reads);
		if(wlen <18) disk_is_full = 1;
// fprintf(ofp, "%s\t%u\t%s\t%u\t%s\t%d\t%s\t%s\n", chro_name_left, chro_pos_left, chro_name_right, chro_pos_right+1, event_body -> is_strand_jumped?"No":"Yes", event_body -> final_counted_reads, event_body -> small_side_increasing_coordinate?"Yes":"No", event_body -> large_side_increasing_coordinate?"Yes":"No");
	}

	global_context -> all_fusions = all_junctions;

	if(global_context->config.do_structural_variance_detection){
		global_context -> translocation_result_table.entry_table -> appendix1 = ofp;
		global_context -> translocation_result_table.entry_table -> appendix2 = global_context;
		HashTableIteration(global_context -> translocation_result_table.entry_table, write_translocation_results_final);
		global_context -> inversion_result_table.entry_table -> appendix1 = ofp;
		global_context -> inversion_result_table.entry_table -> appendix2 = global_context;
		HashTableIteration(global_context -> inversion_result_table.entry_table, write_inversion_results_final);
	}

	fclose(ofp);

	if(disk_is_full){
		unlink(fn2);
		SUBREADprintf("ERROR: disk is full. No fusion table is generated.\n");
	}
	return 0;
}

void write_inversion_results_final(void * bukey, void * buckv, HashTable * tab){
	int x1;
	bucketed_table_bucket_t * buck = buckv;

	FILE * ofp = (FILE *)tab -> appendix1;
	global_context_t * global_context = (global_context_t * )tab -> appendix2;
	for(x1 = 0; x1 < buck -> items; x1++)
	{
		if(buck->positions[x1] - buck->positions[x1] % buck -> maximum_interval_length == buck -> keyed_bucket)
		{
			inversion_result_t * inv_res = buck -> details[x1];

			char * src_chr;
			int src_pos;

			locate_gene_position(inv_res -> small_side,  &global_context -> chromosome_table, &src_chr , &src_pos);
			fprintf(ofp, "INV\t%s\t%d\t%s\t%u\t%s\n",  src_chr, src_pos + 1, src_chr, src_pos + 1 + inv_res -> length,  inv_res -> is_precisely_called ? "PRECISE":"IMPRECISE");
			fprintf(ofp, "INV\t%s\t%d\t%s\t%u\t%s\n",  src_chr, src_pos + 2, src_chr, src_pos + inv_res -> length,  inv_res -> is_precisely_called ? "PRECISE":"IMPRECISE");

			//fprintf(ofp, "INVERSION\t%s\t%u\t%u\t%u\t%u\n", src_chr, src_pos, inv_res -> length, inv_res -> all_sup_D , inv_res -> max_sup_E);
		}
	}
	
}

void write_translocation_results_final(void * bukey, void * buckv, HashTable * tab){
	int x1;
	bucketed_table_bucket_t * buck = buckv;

	FILE * ofp = (FILE *)tab -> appendix1;
	global_context_t * global_context = (global_context_t * )tab -> appendix2;
	for(x1 = 0; x1 < buck -> items; x1++)
	{
		if(buck->positions[x1] - buck->positions[x1] % buck -> maximum_interval_length == buck -> keyed_bucket)
		{
			char * src_chr, *targ_chr;
			int src_pos, targ_pos;

			translocation_result_t * trans_res = buck -> details[x1];

			locate_gene_position(trans_res -> source_left_side,  &global_context -> chromosome_table, &src_chr , &src_pos);
			locate_gene_position(trans_res -> target_left_side,  &global_context -> chromosome_table, &targ_chr , &targ_pos);

			//fprintf(ofp, "TRANSLOCATION\t%s\t%u\t%u\t%s\t%u\t%s\t%u\t%u\n", src_chr, src_pos, trans_res -> length, targ_chr, targ_pos, trans_res -> is_inv?"INV":"STR", trans_res -> all_sup_P , trans_res -> max_sup_QR);
			/*
			SUBREADprintf("ABS=%u, %u, PRECISE=%d\n", trans_res -> source_left_side, trans_res -> target_left_side, trans_res -> is_precisely_called);
			SUBREADprintf("%u, %u\n", src_pos, targ_pos);
			SUBREADprintf("%s, %s\n", src_chr, targ_chr);
			*/
			fprintf(ofp, "%s\t%s\t%u\t%s\t%d\t%s\t%s\n", src_chr == targ_chr?"ITX":"CTX", src_chr, src_pos + 1, targ_chr, targ_pos + 1, trans_res -> is_inv?"X":"=",  trans_res -> is_precisely_called ? "PRECISE":"IMPRECISE");
			fprintf(ofp, "%s\t%s\t%u\t%s\t%d\t%s\t%s\n", src_chr == targ_chr?"ITX":"CTX", src_chr, src_pos + trans_res -> length + 1, targ_chr, targ_pos + 1, trans_res -> is_inv?"X":"=", trans_res -> is_precisely_called ? "PRECISE":"IMPRECISE");
			fprintf(ofp, "DEL\t%s\t%d\t%u\t%s\n", src_chr, src_pos + 1, trans_res -> length ,  trans_res -> is_precisely_called ? "PRECISE":"IMPRECISE");
		}
	}

}

int write_junction_final_results(global_context_t * global_context)
{

	int no_sup_juncs = 0, disk_is_full = 0;

	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 
	char fn2 [MAX_FILE_NAME_LENGTH+30];

	sprintf(fn2,"%s.junction.bed", global_context->config.output_prefix);
	FILE * ofp = f_subr_open(fn2, "wb");

	fprintf(ofp, "#Chr, StartLeftBlock, EndRightBlock, Junction_Name, nSupport, Strand, StartLeftBlock, EndRightBlock, Color, nBlocks, BlockSizes, BlockStarts\n");

	int xk1;
	unsigned int all_junctions = 0;

	for(xk1 = 0; xk1 < indel_context -> total_events ; xk1++)
	{ 
		char * chro_name_left,* chro_name_right, indel_sect[10];
		int chro_pos_left, chro_pos_right; 
		chromosome_event_t * event_body = indel_context -> event_space_dynamic +xk1;
		if(event_body -> event_type != CHRO_EVENT_TYPE_JUNCTION)
			continue;

		//#warning "  ================================== remove '- 1' from the next line!!! ================================="
		if(event_body->final_counted_reads <  1 || ( event_body->critical_supporting_reads < 1 - 1&& event_body->indel_at_junction))
		{
			no_sup_juncs++;
			continue;
		}

		locate_gene_position( event_body -> event_small_side , &global_context -> chromosome_table, &chro_name_left, &chro_pos_left);
		locate_gene_position( event_body -> event_large_side , &global_context -> chromosome_table, &chro_name_right, &chro_pos_right);

		chro_pos_left++;


		unsigned int feature_start = chro_pos_left - event_body -> junction_flanking_left;
		if(chro_pos_left <= event_body -> junction_flanking_left){
			feature_start = 1;
			event_body -> junction_flanking_left = chro_pos_left - 1;
		}

		unsigned int feature_end = chro_pos_right + event_body -> junction_flanking_right;

		all_junctions ++;

		indel_sect[0]=0;
		if(event_body->indel_at_junction)
			sprintf(indel_sect,"INS%d", event_body->indel_at_junction);
		if(event_body-> is_donor_found_or_annotation &64)strcat(indel_sect,"ANNO");
		//else if(event_body->critical_supporting_reads < 1)
		//	strcpy(indel_sect, "NOCRT");


		int wlen = fprintf(ofp,"%s\t%u\t%u\tJUNC%08u%s\t%d\t%c\t%u\t%u\t%d,%d,%d\t2\t%d,%d\t0,%u\n", chro_name_left, feature_start,  feature_end,
												all_junctions, indel_sect,  event_body -> final_counted_reads, event_body->is_negative_strand?'-':'+',
												feature_start,  feature_end, event_body->is_negative_strand?0:255, /*event_body -> anti_supporting_reads*/ event_body->is_negative_strand?255:0, event_body->is_negative_strand?255:0,
												 event_body -> junction_flanking_left, event_body -> junction_flanking_right, feature_end-feature_start-event_body -> junction_flanking_right);
		if(wlen < 10) disk_is_full = 1;
	}

	fclose(ofp);
	if(disk_is_full){
		unlink(fn2);
		SUBREADprintf("ERROR: disk is full; no junction table is created.\n");
	}
	global_context -> all_junctions = all_junctions;
	//printf("Non-support juncs=%d;  Final juncs = %d\n", no_sup_juncs, all_junctions);
	return 0;
}



void get_chro_2base(char *buf, gene_value_index_t * index, unsigned int pos, int is_negative_strand)
{
	gvindex_get_string (buf, index, pos, 2, is_negative_strand);
}


int paired_chars_part(char * ch1, char * ch2, int is_reverse)
{
	if (c2eq(ch1, ch2, "GT", "AG") || c2eq(ch1, ch2, "CT", "AC"))
	{
		if (is_reverse) if (ceq(ch1, "AG") || ceq(ch1, "AC")) return 1;
		if (!is_reverse) if (ceq(ch1, "CT") || ceq(ch1, "GT")) return 1;
	}
	return 0;
}
#define is_donar_chars_part(cc) (((cc)[0]=='G' && (cc)[1]=='T') || \
			    ((cc)[0]=='A' && (cc)[1]=='G') || \
			    ((cc)[0]=='A' && (cc)[1]=='C') || \
			    ((cc)[0]=='C' && (cc)[1]=='T')) 


#define SHORT_EXON_MIN_LENGTH 18
#define EXON_EXTENDING_SCAN 0
#define SHORT_EXON_WINDOW 6 
#define SHORT_EXON_EXTEND 5000

void core_search_short_exons(global_context_t * global_context, thread_context_t * thread_context, char * read_text, char * qualityb0, int rl, unsigned int P1_Pos, unsigned int P2_Pos, short read_coverage_start, short read_coverage_end)
{
	char inb[MAX_READ_LENGTH], qualityb[MAX_READ_LENGTH];
	if ( (rl <= EXON_LONG_READ_LENGTH ) && (!EXON_EXTENDING_SCAN)) return;
	//return;
	gene_value_index_t * base_index = thread_context?thread_context->current_value_index:global_context->current_value_index ;
	//insert event
	HashTable * event_table = NULL;
	chromosome_event_t * event_space = NULL;
	if(thread_context)
	{
		event_table = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}
	else
	{
		event_table = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}

	strcpy(inb, read_text);
	strcpy(qualityb, qualityb0);

	unsigned int pos_small=min(P1_Pos, P2_Pos), pos_big = max(P1_Pos, P2_Pos);

	int max_score , test_score;
	unsigned int best_j1_edge=0 , best_j2_edge=0;
	int need_to_test = 0;

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
// SCAN TO THE HEAD  /////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

	if (read_coverage_start  > SHORT_EXON_MIN_LENGTH)
	{
		max_score = -1;

		int need_check2 = 1;
		if(qualityb[0])
		{
			float head_quality = read_quality_score(qualityb , SHORT_EXON_MIN_LENGTH , global_context->config.phred_score_format); 
			if(head_quality < 6 )
				need_check2 = 0;
		}


		if(need_check2)
			if(SHORT_EXON_MIN_LENGTH *0.6 < match_chro(inb, base_index, pos_small, SHORT_EXON_MIN_LENGTH , 0, global_context->config.space_type))
				need_check2 = 0; 


		if(need_check2)
		{

			int delta_pos, is_indel = 0;
			for(delta_pos=-3; delta_pos <=3; delta_pos ++)
			{
				if(match_chro(inb, base_index, pos_small + delta_pos, SHORT_EXON_MIN_LENGTH , 0, global_context->config.space_type) >= SHORT_EXON_MIN_LENGTH*.7)
				{
					is_indel = 1;
					break;
				}
			}
			// The head of the read is incorrect. Do we need to search a long way?
			// See if there is a donor in the head area.
			int test_donor_pos;
			char cc[3];
			cc[2]=0;

			if(!is_indel)
				for(test_donor_pos = SHORT_EXON_MIN_LENGTH ; test_donor_pos < read_coverage_start ; test_donor_pos ++)
				{
					get_chro_2base(cc, base_index, pos_small + test_donor_pos, 0);
					if(is_donar_chars_part(cc))
					{
						need_to_test = 1;
						break;
					}
				}
		}
	}

	max_score = -999;
	int max_is_GTAG = 0;

	if(need_to_test && pos_small >= SHORT_EXON_MIN_LENGTH)
	{
		unsigned int test_end = pos_small - SHORT_EXON_EXTEND;
		if(SHORT_EXON_EXTEND > pos_small) test_end = 0;

		unsigned int new_pos = pos_small-SHORT_EXON_MIN_LENGTH;
		while(1)
		{
			new_pos = match_chro_range(inb,  base_index, new_pos, 7 , new_pos - test_end , SEARCH_BACK);
			if(new_pos==0xffffffff) break;
			// There is an exact match. See if the donor/receptors are matched.
			// new_pos is the new head position of the read.
			int splice_point;
			for(splice_point = SHORT_EXON_MIN_LENGTH; splice_point < read_coverage_start ; splice_point ++)
			{
				char cc[3];
				cc[2]=0;
				char cc2[3];
				cc2[2]=0;

				get_chro_2base(cc, base_index, pos_small + splice_point -2, 0);
				if(is_donar_chars_part(cc))
				{
					// <<< EXON---|CC2---INTRON---CC|---EXON
					get_chro_2base(cc2, base_index, new_pos + splice_point, 0);
					if(is_donar_chars_part(cc2) && paired_chars_part(cc2 , cc, 0)) 
					{
						int matched_in_exon_old = match_chro(inb + splice_point, base_index, pos_small + splice_point , SHORT_EXON_WINDOW , 0, global_context->config.space_type);
						int matched_in_exon_new = match_chro(inb, base_index, new_pos , splice_point, 0, global_context->config.space_type);

						
						test_score = 1000000+ (matched_in_exon_new )*10000  + matched_in_exon_old * 1000 + new_pos - test_end;
						if(test_score <= max_score) continue;
						max_score = test_score;

						if(matched_in_exon_new < splice_point || matched_in_exon_old < SHORT_EXON_WINDOW ) 
							continue;

						max_is_GTAG = (cc2[0]=='G' || cc2[1]=='G');
						//printf("EX CC=%s\tCC2=%s\tis_GTAG=%d\n",cc,cc2,max_is_GTAG);
						best_j1_edge = new_pos + splice_point - 1;
						best_j2_edge = pos_small + splice_point;
					}
				}
			}
		}
	}


	if(best_j1_edge>0)
	{
		int event_no;
		chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];
		chromosome_event_t * found = NULL;

		int found_events = search_event(global_context, event_table, event_space, best_j1_edge , EVENT_SEARCH_BY_SMALL_SIDE, CHRO_EVENT_TYPE_JUNCTION|CHRO_EVENT_TYPE_FUSION, search_return);

		if(found_events)
		{
			int kx1; 
			for(kx1 = 0; kx1 < found_events ; kx1++)
			{
				if(search_return[kx1] -> event_large_side == best_j2_edge)
				{
					found = search_return[kx1];	
					break;
				}
			}
		}

		if(found) found -> supporting_reads ++;
		else
		{
			if(thread_context)
				event_no = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> total_events ++;
			else
				event_no = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) ->  total_events ++;

			event_space = reallocate_event_space(global_context, thread_context, event_no);

			chromosome_event_t * new_event = event_space+event_no; 
			memset(new_event,0,sizeof(chromosome_event_t));
			new_event -> event_small_side = best_j1_edge;
			new_event -> event_large_side = best_j2_edge;

			new_event -> is_negative_strand= !max_is_GTAG;
			new_event -> event_type = CHRO_EVENT_TYPE_JUNCTION;

			new_event -> supporting_reads = 1;
			new_event -> indel_length = 0;

			put_new_event(event_table, new_event , event_no);
		}
		//printf("FOUND NEW JUNCTION HEAD: %u - %u\n", best_j1_edge, best_j2_edge);
	}


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
// SCAN TO THE TAIL  /////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

	need_to_test = 0;
	max_score = -999;


	if (read_coverage_end< rl - SHORT_EXON_MIN_LENGTH)
	{
		int need_check2 = 1;
		if(qualityb[0])
		{
			float head_quality = read_quality_score(qualityb + rl - SHORT_EXON_MIN_LENGTH , SHORT_EXON_MIN_LENGTH , global_context->config.phred_score_format); 
			if(head_quality < 6 )
				need_check2 = 0;
		}


		if(SHORT_EXON_MIN_LENGTH *0.6 < match_chro(inb + rl - SHORT_EXON_MIN_LENGTH, base_index, pos_big + rl - SHORT_EXON_MIN_LENGTH , SHORT_EXON_MIN_LENGTH , 0, global_context->config.space_type))
			need_check2 = 0; 
		if(need_check2)
		{
			int delta_pos, is_indel = 0;
			for(delta_pos=-3; delta_pos <=3; delta_pos ++)
			{
				if(match_chro(inb + rl - SHORT_EXON_MIN_LENGTH, base_index, pos_big + rl - SHORT_EXON_MIN_LENGTH + delta_pos, SHORT_EXON_MIN_LENGTH , 0, global_context->config.space_type) >= SHORT_EXON_MIN_LENGTH*.7)
				{
					is_indel = 1;
					break;
				}
			}
			// The head of the read is incorrect. Do we need to search a long way?
			// See if there is a donor in the head area.
			int test_donor_pos;
			char cc[3];
			cc[2]=0;

			if(!is_indel)
				for(test_donor_pos = read_coverage_end  ; test_donor_pos < rl ; test_donor_pos ++)
				{
					get_chro_2base(cc, base_index, pos_big + test_donor_pos, 0);
					if(is_donar_chars_part(cc))
					{
						need_to_test = 1;
						break;
					}
				}
		}
	}

	best_j1_edge = 0;
	max_is_GTAG = 0;

	if(need_to_test)
	{
		unsigned int test_end = pos_big + SHORT_EXON_EXTEND;
		if(test_end > base_index -> length + base_index -> start_point) test_end = base_index -> length + base_index -> start_point;

		unsigned int new_pos = pos_big +rl - SHORT_EXON_MIN_LENGTH +16;

		while(1)
		{
			if(new_pos +  test_end - new_pos < base_index-> start_base_offset + base_index->length)
			{
				assert(new_pos<0xffff0000);
				new_pos = match_chro_range(inb + rl - SHORT_EXON_MIN_LENGTH,  base_index, new_pos, 7 , test_end - new_pos , SEARCH_FRONT);
			}
			else break;

			if(new_pos==0xffffffff) break;
			// There is an exact match. See if the donor/receptors are matched.
			// (new_pos + SHORT_EXON_MIN_LENGTH -rl + splice_point) is the new exon start.

			int splice_point;
			for(splice_point = read_coverage_end ; splice_point < rl -  SHORT_EXON_MIN_LENGTH; splice_point ++)
			{
				char cc[3];
				cc[2]=0;
				char cc2[3];
				cc2[2]=0;

				unsigned int new_pos_tail = (new_pos + SHORT_EXON_MIN_LENGTH -rl + splice_point);

				get_chro_2base(cc, base_index, pos_big + splice_point, 0);
				if(is_donar_chars_part(cc))
				{
					get_chro_2base(cc2, base_index, new_pos_tail -2, 0);
					if(is_donar_chars_part(cc2) && paired_chars_part(cc , cc2, 0)) 
					{
						int matched_in_exon_new = match_chro(inb + splice_point, base_index, new_pos_tail , rl - splice_point , 0, global_context->config.space_type);
						int matched_in_exon_old = match_chro(inb + splice_point - SHORT_EXON_WINDOW , base_index, pos_big + splice_point - SHORT_EXON_WINDOW , SHORT_EXON_WINDOW, 0, global_context->config.space_type);

						test_score = 1000000+ (matched_in_exon_new)*10000 + matched_in_exon_old * 1000  + test_end - new_pos;
						if(test_score <= max_score) continue;
						max_score = test_score;

						if(matched_in_exon_new < (rl - splice_point) || matched_in_exon_old < SHORT_EXON_WINDOW)
							continue;

						// EXON ---|CC---INTRON---CC2|--- EXON >>>
						max_is_GTAG = (cc[0]=='G'|| cc[1]=='G');
						best_j1_edge = pos_big + splice_point - 1;
						best_j2_edge = new_pos_tail;
					}
				}
			}

		}
	}


	if(best_j1_edge>0)
	{
		int event_no;
		chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];
		chromosome_event_t * found = NULL;

		int found_events = search_event(global_context, event_table, event_space, best_j1_edge , EVENT_SEARCH_BY_SMALL_SIDE, CHRO_EVENT_TYPE_JUNCTION|CHRO_EVENT_TYPE_FUSION, search_return);

		if(found_events)
		{
			int kx1; 
			for(kx1 = 0; kx1 < found_events ; kx1++)
			{
				if(search_return[kx1] -> event_large_side == best_j2_edge)
				{
					found = search_return[kx1];	
					break;
				}
			}
		}

		if(found) found -> supporting_reads ++;
		else
		{
			if(thread_context)
				event_no = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> total_events ++;
			else
				event_no = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) ->  total_events ++;


			event_space = reallocate_event_space(global_context, thread_context, event_no);

			chromosome_event_t * new_event = event_space+event_no; 
			memset(new_event,0,sizeof(chromosome_event_t));
			new_event -> event_small_side = best_j1_edge;
			new_event -> event_large_side = best_j2_edge;

			new_event -> is_negative_strand= !max_is_GTAG;
			new_event -> event_type = CHRO_EVENT_TYPE_JUNCTION;

			new_event -> supporting_reads = 1;
			new_event -> indel_length = 0;

			put_new_event(event_table, new_event , event_no);
			//printf("FOUND NEW JUNCTION TAIL: %u - %u\n", best_j1_edge, best_j2_edge);
		}
	}
}









int core_select_best_matching_halves_maxone(global_context_t * global_context, gene_vote_t * vote, unsigned int * best_pos1, unsigned int * best_pos2, int * best_vote1, int * best_vote2, char * is_abnormal, short * half_marks, int * is_reversed_halves, float accept_rate, int read_len, long long int hint_pos, int tolerable_bases, short * read_coverage_start, short * read_coverage_end, gene_vote_number_t * indel_in_p1, gene_vote_number_t * indel_in_p2, gehash_data_t max_pos, gene_vote_number_t max_votes, short max_start, short max_end, short max_mask, gene_vote_number_t * max_indel_recorder, int* best_select_max_votes, int rl)
{
	int best_splicing_point = -1, i,j;
	char * best_chro_name, is_reversed;
	int best_chro_pos;
	int selected_max_votes = -1;


	is_reversed = (max_mask & IS_NEGATIVE_STRAND)?1:0;
	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			char * chro_name;
			char is_partner_reversed;
			int chro_pos;

			int overlapped_len, overlap_start, overlap_end;
			// All logical conditions

			//if( (vote->votes[i][j] < vote-> coverage_start[i][j]) < 12 && (vote-> coverage_end[i][j] > rl - 12 )) continue;

			is_partner_reversed = (vote->masks [i][j] & IS_NEGATIVE_STRAND) ? 1:0;
			overlap_start = max(max_start , vote->coverage_start[i][j]);
			overlap_end   = min(max_end , vote->coverage_end[i][j]);
			overlapped_len =overlap_end - overlap_start;

			int coverage_len = max_end - max_start + vote->coverage_end[i][j] - vote->coverage_start[i][j];
			if (overlapped_len >0)coverage_len -= overlapped_len;
			//SUBREADprintf("MAX: %d-%d   OTHER %d-%d    COV=%d   OVLP=%d\n", max_start, max_end, vote->coverage_start[i][j], vote->coverage_end[i][j], coverage_len, overlapped_len);



			if(overlapped_len >=14)
				continue;

			long long int dist = vote->pos[i][j];
			dist -= max_pos;

			//SUBREADprintf ("D=%lld\n", abs(dist));
			if (abs(dist)<6)
				continue;

			int support_r1 = 1; 
			int support_r2 = 1;

			if (max_votes < support_r1 || vote->votes[i][j]<support_r2)
				continue;

			// Same chromosome
			if ((vote->coverage_start[i][j] < max_start) + is_reversed == 1)
			{
				locate_gene_position(max_pos + read_len, &(global_context -> chromosome_table) , &best_chro_name, &best_chro_pos);
				locate_gene_position(vote->pos[i][j] , &(global_context -> chromosome_table), &chro_name, &chro_pos);
			}else
			{
				locate_gene_position(max_pos , &(global_context -> chromosome_table), &best_chro_name, &best_chro_pos);
				locate_gene_position(vote->pos[i][j] +read_len, &(global_context -> chromosome_table), &chro_name, &chro_pos);
			}

			if (chro_name != best_chro_name)	// The pointers can be compared because they can be the same.
				continue;

			int is_fusion = 0;

			if(is_reversed != is_partner_reversed) is_fusion = 1; 

			if( is_reversed && ((max_pos > vote->pos[i][j]) + (vote->coverage_start[i][j] < max_start) != 1))is_fusion = 1;
			if((! is_reversed) && ((max_pos > vote->pos[i][j]) + (vote->coverage_start[i][j] > max_start) != 1)) is_fusion = 1;

			if(abs(dist) > 500000 || chro_name != best_chro_name) continue;

			int test_vote_value ;
			test_vote_value = 8888888 +  vote->votes[i][j]* 1000000 - abs(dist);
			if (hint_pos>=0)
			{
				long long int hint_dist = hint_pos;
				hint_dist -= vote->pos[i][j];
				if (abs (hint_dist) < 100000)
					test_vote_value += 100;
				if (abs (hint_dist) < 5000)
					test_vote_value += 100;
			}

			if (test_vote_value<selected_max_votes)continue;
			// Conditions of order of R3 and R5
			*half_marks &= ~IS_REVERSED_HALVES;
			if (vote->coverage_start[i][j] < max_start && (((max_pos < vote->pos[i][j]) && !is_reversed) || ((max_pos > vote->pos[i][j]) && is_reversed) ) )
				*half_marks |= IS_REVERSED_HALVES;
			if (vote->coverage_start[i][j] >= max_end  &&  (((max_pos > vote->pos[i][j]) && !is_reversed) || ((max_pos < vote->pos[i][j]) && is_reversed) ) )
				*half_marks |= IS_REVERSED_HALVES;

			if (vote->coverage_start[i][j] < max_start)
			{
				(*half_marks) = (*half_marks) & ~IS_R1_CLOSE_TO_5;
			}
			else
			{
				(*half_marks) |= IS_R1_CLOSE_TO_5;
			}

			if(max_mask & IS_NEGATIVE_STRAND)
				*half_marks = (*half_marks) |   IS_NEGATIVE_STRAND_R1;
			else
				*half_marks = (*half_marks) &  ~IS_NEGATIVE_STRAND_R1;

			if(vote->masks[i][j] & IS_NEGATIVE_STRAND)
				*half_marks = (*half_marks) |   IS_NEGATIVE_STRAND_R2;
			else
				*half_marks = (*half_marks) &  ~IS_NEGATIVE_STRAND_R2;
	

			
			best_splicing_point = ((vote->coverage_start[i][j] < max_start)? (vote->coverage_end[i][j]):(max_end)) + ((vote->coverage_start[i][j] < max_start)? (max_start):(vote->coverage_start[i][j]));


			best_splicing_point /=2;

			* best_pos1 = max_pos ;
			* best_pos2 = vote->pos[i][j] ;
			* best_vote1 = max_votes ;
			* best_vote2 = vote->votes[i][j] ;
			* read_coverage_start = min(vote->coverage_start[i][j] , max_start);
			* read_coverage_end = max(vote->coverage_end[i][j] , max_end);

			* read_coverage_start = max_start;
			* read_coverage_end = max_end;
			
			int k;
			for(k=0; k<MAX_INDEL_TOLERANCE ; k+=3)
				if(!max_indel_recorder[k+3])break;
			* indel_in_p1 = max_indel_recorder[k+2];

			for(k=0; k<MAX_INDEL_TOLERANCE ; k+=3)
				if(!vote->indel_recorder[i][j][k+3])break;
			* indel_in_p2 = vote->indel_recorder[i][j][k+2];


			* is_reversed_halves = is_reversed;

			if (test_vote_value >=100)
				*half_marks = (*half_marks) | IS_PAIRED_HINTED;
			else
				*half_marks = (*half_marks) & ~(IS_PAIRED_HINTED);

			if (is_fusion)
				*half_marks = (*half_marks)    | IS_FUSION;
			else
				*half_marks = (*half_marks) & ~( IS_FUSION);
	

			selected_max_votes = test_vote_value; 

		}
	*best_select_max_votes = selected_max_votes ;
	return best_splicing_point;
}



int core_select_best_matching_halves(global_context_t * global_context , gene_vote_t * vote, unsigned int * best_pos1, unsigned int * best_pos2, int * best_vote1, int * best_vote2, char * is_abnormal, short * half_marks, int * is_reversed_halves, float accept_rate, int read_len, long long int hint_pos, int tolerable_bases, short * read_coverage_start, short * read_coverage_end, char * indel_in_p1, char * indel_in_p2 , int * max_cover_start, int * max_cover_end, int rl, int repeated_pos_base, int is_negative, char * repeat_record, unsigned int index_valid_range)
{
	unsigned int tmp_best_pos1=0, tmp_best_pos2=0;
	int tmp_best_vote1=0, tmp_best_vote2=0, tmp_is_reversed_halves=0;
	char tmp_is_abnormal=0;
	gene_vote_number_t tmp_indel_in_p1=0, tmp_indel_in_p2=0;
	short tmp_half_marks=0, tmp_read_coverage_start=0, tmp_read_coverage_end=0;
	int ret = 0, best_ret = 0;	

	int i,j;
	int test_select_votes=-1, best_select_votes = 1000000;
	//int max_minor = 0;

	/*
	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			if(vote->votes[i][j] < vote->max_vote)continue;
			int ii,jj;
			for (ii=0; ii<GENE_VOTE_TABLE_SIZE;ii++)
				for(jj=0; jj< vote->items[ii]; jj++)
				{
					if(max_minor >= vote->votes[ii][jj]) continue;
					if(ii==i && jj==j)continue;
					long long int dist =  vote->pos[ii][jj];
					dist =abs(dist - vote->pos[i][j]);
					if(dist > 500000)
						continue;
					max_minor = vote->votes[ii][jj];
				}

		}

	int encountered = 0;


	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			if(vote->votes[i][j] < vote->max_vote)continue;
			int ii,jj;
			for (ii=0; ii<GENE_VOTE_TABLE_SIZE;ii++)
				for(jj=0; jj< vote->items[ii]; jj++)
				{
					if(max_minor != vote->votes[ii][jj]) continue;
					if(ii==i && jj==j)continue;
					long long int dist =  vote->pos[ii][jj];
					dist =abs(dist - vote->pos[i][j]);
					if(dist > 500000)
						continue;
					encountered++;
				}

		}
	*/

	int repeated_pos = repeated_pos_base;
	int offset_shifting = (rl > 220)?4:0;
	//int encounter = 0;

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			/*if((vote->votes[i][j] >=  vote->max_vote -1) && (vote->max_coverage_start >= vote-> coverage_start[i][j] - EXON_MAX_BIGMARGIN_OVERLAPPING ) &&  (vote->max_coverage_end <= vote-> coverage_end[i][j] + EXON_MAX_BIGMARGIN_OVERLAPPING))
				encounter++;*/
			if(repeated_pos_base>=0 && vote->pos[i][j]<=index_valid_range)
				if(vote->votes[i][j] >=  vote->max_vote && repeated_pos < repeated_pos_base+12)
				{
					repeat_record[repeated_pos] = (vote-> coverage_start[i][j] >> offset_shifting);
					repeat_record[repeated_pos+1] = (vote-> coverage_end[i][j] >> offset_shifting);
					repeat_record[repeated_pos+2] = (is_negative?0x80:0) | (vote->votes[i][j]&0x7f);
					repeated_pos+=3;
				}
		}
	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			if(repeated_pos_base>=0 && vote->pos[i][j]<=index_valid_range)
				if(vote->votes[i][j] ==  vote->max_vote -1 && repeated_pos < repeated_pos_base+12)
				{
					repeat_record[repeated_pos] = (vote-> coverage_start[i][j] >> offset_shifting);
					repeat_record[repeated_pos+1] = (vote-> coverage_end[i][j] >> offset_shifting);
					repeat_record[repeated_pos+2] = (is_negative?0x80:0) | (vote->votes[i][j]&0x7f);
					repeated_pos+=3;
				}
		}


	/*
	if(encounter>=2)
		return 0;
	*/

	ret = core_select_best_matching_halves_maxone(global_context, vote, &tmp_best_pos1, &tmp_best_pos2, &tmp_best_vote1, &tmp_best_vote2,  &tmp_is_abnormal,&tmp_half_marks, &tmp_is_reversed_halves, accept_rate, read_len, hint_pos,  tolerable_bases, &tmp_read_coverage_start, &tmp_read_coverage_end, &tmp_indel_in_p1, &tmp_indel_in_p2, vote -> max_position,  vote->max_vote, vote-> max_coverage_start, vote-> max_coverage_end,  vote-> max_mask, vote->max_indel_recorder, &test_select_votes, rl);
	test_select_votes += vote->max_vote*1000000;
			//SUBREADprintf("TSV=%d\n",test_select_votes);

	if(test_select_votes > best_select_votes)
	{
		best_select_votes = test_select_votes;
		*best_pos1 = tmp_best_pos1;
		*best_pos2 = tmp_best_pos2;
		*is_reversed_halves= tmp_is_reversed_halves;
		
		*best_vote1 = tmp_best_vote1;
		*best_vote2 = tmp_best_vote2;
		*is_abnormal = tmp_is_abnormal;
		*indel_in_p1 = tmp_indel_in_p1;
		*indel_in_p2 = tmp_indel_in_p2;
				
		*half_marks = tmp_half_marks;
		*read_coverage_start = tmp_read_coverage_start;
		*read_coverage_end = tmp_read_coverage_end;

		* max_cover_start = vote-> max_coverage_start;
		* max_cover_end = vote-> max_coverage_end;
		best_ret = ret;
	}		
	return best_ret;
}



#define EXON_DONOR_TEST_WINDOW 17


// pos1 must be small than pos2.
int core13_test_donor(char *read, int read_len, unsigned int pos1, unsigned int pos2, int guess_break_point, char negative_strand, int test_range, char is_soft_condition, int EXON_INDEL_TOLERANCE, int* real_break_point, gene_value_index_t * my_value_array_index, int indel_offset1, int indel_offset2, int is_reversed, int space_type, int * best_donor_score, int * is_GTAG)
{
	int bps_pos_x;
	int search_start = guess_break_point - test_range ;
	int search_end   = guess_break_point + test_range ;
	char h1_2ch[3], h2_2ch[3];

	h1_2ch[2] = h2_2ch[2]=0;
	search_start=max(10, search_start);
	search_end = min(read_len-10, search_end);
	int best_break = -1;
	int min_x = -9099;

	for (bps_pos_x = search_start; bps_pos_x < search_end ; bps_pos_x ++)
	{
		int paired_score = 0;
		get_chro_2base(h1_2ch, my_value_array_index, pos1 - indel_offset1+ bps_pos_x , is_reversed);
		get_chro_2base(h2_2ch, my_value_array_index, pos2 - 2 - indel_offset2 + bps_pos_x, is_reversed);


		//if(!is_reversed)
		//SUBREADprintf("C1=%s @%u, C2=%s @%u\n",h1_2ch, pos1 + bps_pos_x, h2_2ch,pos2 - 2 + indel_offset + bps_pos_x);
		if(h1_2ch[0]==h2_2ch[0] && h1_2ch[1]==h2_2ch[1]) continue;

		if(is_donar_chars_part(h1_2ch) && is_donar_chars_part(h2_2ch))
		{

			paired_score = paired_chars_part(h1_2ch, h2_2ch, is_reversed);

			if(paired_score)
			{
				int m1, m2, x1, x2;
				int break_point_half = is_reversed?(read_len - bps_pos_x):bps_pos_x;
				int first_exon_end,second_half_start;
				int donar_conf_len = 0;

				donar_conf_len = min(break_point_half , EXON_DONOR_TEST_WINDOW);
				donar_conf_len = min(read_len - break_point_half, donar_conf_len);
				//SUBREADprintf("DONOR_CONF_LEN=%d\n", donar_conf_len);

				if (is_reversed)
				{
					first_exon_end = pos2 + bps_pos_x - indel_offset2;
					second_half_start = pos1 + bps_pos_x- indel_offset1;

					m1 = match_chro(read + break_point_half - donar_conf_len , my_value_array_index, first_exon_end, donar_conf_len, is_reversed, space_type);
					m2 = match_chro(read + break_point_half , my_value_array_index, second_half_start-donar_conf_len , donar_conf_len, is_reversed, space_type);

					x1 = match_chro(read + break_point_half ,  my_value_array_index, first_exon_end - donar_conf_len, donar_conf_len , is_reversed, space_type);
					x2 = match_chro(read + break_point_half - donar_conf_len ,  my_value_array_index, second_half_start , donar_conf_len, is_reversed, space_type);
				}
				else
				{
					first_exon_end = pos1 + bps_pos_x - indel_offset1;
					second_half_start = pos2 + bps_pos_x - indel_offset2;

					m1 = match_chro(read + break_point_half - donar_conf_len, my_value_array_index, first_exon_end-donar_conf_len , donar_conf_len, is_reversed, space_type);
					m2 = match_chro(read + break_point_half , my_value_array_index, second_half_start, donar_conf_len, is_reversed, space_type);

					x1 = match_chro(read + break_point_half ,  my_value_array_index, first_exon_end, donar_conf_len , is_reversed,space_type);
					x2 = match_chro(read + break_point_half - donar_conf_len,  my_value_array_index, second_half_start - donar_conf_len, donar_conf_len , is_reversed,space_type);
				}

				#ifdef TEST_TARGET
				if(memcmp(read, TEST_TARGET, 15)==0)
				{
					SUBREADprintf("DONOR TEST STR=%s, %s ; pos=%d    %d %d ; M=%d %d ; X=%d %d\n", h1_2ch, h2_2ch, bps_pos_x, indel_offset1, indel_offset2, m1, m2, x1, x2);
				}
				#endif
	
				int threshold = 3;
				if (paired_score == 1)
					threshold = 3;

				#ifdef QUALITY_KILL
				if (m1 >= donar_conf_len-1    && m2>=donar_conf_len-1 )
					if(x1<donar_conf_len - threshold  && x2<donar_conf_len- threshold )
				#else
				if (m1 >= donar_conf_len-1    && m2>=donar_conf_len -1)
					if(x1<donar_conf_len - threshold  && x2<donar_conf_len - threshold)
				#endif
					{
						int score =  3000-(x1 + x2) + (m1+ m2) ;
						if (min_x < score)
						{
							min_x = score;
							best_break = bps_pos_x;
							*is_GTAG = 1==((is_reversed) + (h1_2ch[0]=='G' || h1_2ch[1]=='G'));	//"GT" or "AG"
							//printf("FL CC=%s\tCC2=%s\tis_GTAG=%d\tREV=%d\n",h1_2ch,h2_2ch,*is_GTAG, is_reversed);
							*best_donor_score = score;
						}
					}
			}
		}
	}

	if (best_break>0)
	{
				#ifdef TEST_TARGET
				if(memcmp(read, TEST_TARGET, 15)==0)
					SUBREADprintf("SELECRED!!!_BREAKPOINT=%d, RAW POS=%u,%u, R=%s\n",  best_break, pos1 , pos2, read);
				#endif
		//SUBREADprintf ("FINAL BREAK: %d   ; REV = %d\n ", best_break, is_reversed);
		*real_break_point = best_break;
		return 1;
	}
	else
	{
				#ifdef TEST_TARGET
				if(memcmp(read, TEST_TARGET, 15)==0)
					SUBREADprintf("KILLED!!!_BREAKPOINT=%d, R=%s\n",  best_break+ pos1, read);
				#endif
	}
	return 0;
}






#define EXON_LARGE_WINDOW 60 
#define ACCEPTED_SUPPORT_RATE 0.3

void core_fragile_junction_voting(global_context_t * global_context, thread_context_t * thread_context, char * rname, char * read, char * qual, unsigned int full_rl, int negative_strand, int color_space, unsigned int low_border, unsigned int high_border, gene_vote_t *vote_p1)
{
	int windows = full_rl / EXON_LARGE_WINDOW +1;
	float overlap = (1.0*windows * EXON_LARGE_WINDOW - full_rl) / (windows-1);

	int ww;
	int window_cursor = 0;

	HashTable * event_table = NULL;
	chromosome_event_t * event_space = NULL;
	if(thread_context)
	{
		event_table = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}
	else
	{
		event_table = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}

	int GENE_SLIDING_STEP = global_context->current_index -> index_gap;


	for(ww=0; ww<windows;ww++)
	{
		window_cursor = (int)(ww * EXON_LARGE_WINDOW - ww * overlap);
		int read_len = EXON_LARGE_WINDOW;
		if(ww == windows-1)
			read_len = full_rl -window_cursor;

		float subread_step = 3.00001;
		int i;
		int subread_no;
		char * InBuff;
		InBuff = read + window_cursor;
		char tmp_char = InBuff[read_len];
		InBuff[read_len] = 0;
		
		init_gene_vote(vote_p1);
		for(subread_no=0; ; subread_no++)
		{
			int subread_offset1 = (int)(subread_step * (subread_no+1));
			subread_offset1 -= subread_offset1%GENE_SLIDING_STEP;
			subread_offset1 += GENE_SLIDING_STEP-1;

			for(i=0; i<GENE_SLIDING_STEP ; i++)
			{
				int subread_offset = (int)(subread_step * subread_no); 
				subread_offset -= subread_offset%GENE_SLIDING_STEP -i;

				char * subread_string = InBuff + subread_offset;
				gehash_key_t subread_integer = genekey2int(subread_string, color_space);

				gehash_go_q(global_context->current_index, subread_integer , subread_offset, read_len,negative_strand, vote_p1, 5, subread_no,  low_border, high_border - read_len);
			}
			if(subread_offset1 >= read_len -16)
				break;
		}

		int ii, jj, kk;
		for(ii = 0; ii < GENE_VOTE_TABLE_SIZE; ii++) {
			for(jj = 0; jj < vote_p1 -> items[ii] ; jj++) {
				if(vote_p1 -> votes[ii][jj] < vote_p1 -> max_vote) continue;

				gene_vote_number_t * indel_recorder = vote_p1 -> indel_recorder[ii][jj];
				unsigned int voting_position =  vote_p1 -> pos[ii][jj];
				int last_indel = 0, last_correct_subread=0;

				for(kk =0; indel_recorder[kk]  && (kk < MAX_INDEL_SECTIONS); kk+=3){
					char movement_buffer[MAX_READ_LENGTH * 10 / 7];
					//chromosome_event_t * last_event = NULL;
					int last_event_id = -1;

					int indels = indel_recorder[kk+2] - last_indel;
					if(indels==0) continue;

					int next_correct_subread = indel_recorder[kk] -1;

					int last_correct_base = find_subread_end(read_len, global_context->config.total_subreads , last_correct_subread) - 9;
					int first_correct_base = find_subread_end(read_len, global_context->config.total_subreads , next_correct_subread) - 16 + 9;
					first_correct_base = min(first_correct_base+10, read_len);
					last_correct_base = max(0, last_correct_base);
					last_correct_base = min(read_len-1, last_correct_base);

					int x1, dyna_steps;

					dyna_steps = core_dynamic_align(global_context, thread_context, InBuff + last_correct_base, first_correct_base - last_correct_base, voting_position + last_correct_base + last_indel, movement_buffer, indels, rname);

					movement_buffer[dyna_steps]=0;

					if(0 && strcmp("MISEQ:13:000000000-A1H1M:1:1112:12194:5511", rname) == 0)
					{
						SUBREADprintf("IR= %d  %d~%d\n", dyna_steps, last_correct_base, first_correct_base);

						for(x1=0; x1<dyna_steps;x1++)
						{
							int mc, mv=movement_buffer[x1];
							if(mv==0)mc='=';
							else if(mv==1)mc='D';
							else if(mv==2)mc='I';
							else mc='X';
							SUBREADprintf("%c",mc);
						}
						SUBREADputs("");
					}
					unsigned int cursor_on_chromosome = voting_position + last_correct_base + last_indel, cursor_on_read = last_correct_base;
					int last_mv = 0;
					unsigned int indel_left_boundary = 0;
					int is_in_indel = 0, current_indel_len = 0, total_mismatch = 0;

					for(x1=0; x1<dyna_steps;x1++)
					{
						int mv=movement_buffer[x1];
						if(mv==3) total_mismatch++;
					}

					if(total_mismatch<2 || (global_context->config.maximise_sensitivity_indel && total_mismatch <= 2 ))
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
									gene_value_index_t * current_value_index = thread_context?thread_context->current_value_index:global_context->current_value_index; 
									int ambiguous_i, ambiguous_count=0;
									int best_matched_bases = match_chro(InBuff + cursor_on_read - 6, current_value_index, indel_left_boundary - 6, 6, 0, global_context->config.space_type)  +
												 match_chro(InBuff + cursor_on_read - min(current_indel_len,0), current_value_index, indel_left_boundary + max(0, current_indel_len), 6, 0, global_context->config.space_type);
									for(ambiguous_i=-5; ambiguous_i<=5; ambiguous_i++)
									{
										int left_match = match_chro(InBuff + cursor_on_read - 6, current_value_index, indel_left_boundary - 6, 6+ambiguous_i, 0, global_context->config.space_type);
										int right_match = match_chro(InBuff + cursor_on_read + ambiguous_i - min(current_indel_len,0), current_value_index, indel_left_boundary + ambiguous_i + max(0, current_indel_len), 6-ambiguous_i, 0,global_context->config.space_type);
										if(left_match+right_match == best_matched_bases) ambiguous_count ++;
									}

									if(0 && strcmp("MISEQ:13:000000000-A1H1M:1:1112:12194:5511", rname) == 0)
										SUBREADprintf("INDEL_DDADD: abs(I=%d); INDELS=%d; LOC=%u\n",i, current_indel_len, indel_left_boundary-1);
									if(abs(current_indel_len)<=global_context -> config.max_indel_length)
									{
										chromosome_event_t * new_event = local_add_indel_event(global_context, thread_context, event_table, InBuff + cursor_on_read + min(0,current_indel_len), indel_left_boundary - 1, current_indel_len, 1, ambiguous_count, 0, NULL);
										if(last_event_id >=0 && new_event){
											// the event space can be changed when the new event is added. the location is updated everytime.
											chromosome_event_t * event_space = NULL;
											if(thread_context)
												event_space = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
											else
												event_space = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
											chromosome_event_t * last_event = event_space + last_event_id;

											int dist = new_event -> event_small_side - last_event -> event_large_side +1;

											new_event -> connected_previous_event_distance = dist;
											last_event -> connected_next_event_distance = dist;
										}

										if (new_event)
											last_event_id = new_event -> global_event_id;
										else	last_event_id = -1;
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
					 last_correct_subread = indel_recorder[i+1]-1;
				}

			}
		}



		if(1)
		{
			finalise_vote(vote_p1);
			select_best_vote(vote_p1);
			//print_votes(vote_p1, global_context -> config.index_prefix);
			unsigned int best_pos1=0;
			unsigned int best_pos2=0;
			int best_vote1=0;
			int best_vote2=0;
			char is_abnormal=0;
			short half_marks=0;
			int is_reversed_halves=0, max_cover_start=0, max_cover_end=0;
			char indel_in_p1=0, indel_in_p2=0;
			short read_coverage_start =0, read_coverage_end=0;
			gene_value_index_t * base_index = thread_context?thread_context->current_value_index:global_context->current_value_index ;

			int splice_point = core_select_best_matching_halves(global_context, vote_p1, &best_pos1, &best_pos2, &best_vote1, &best_vote2, &is_abnormal ,&half_marks, &is_reversed_halves, ACCEPTED_SUPPORT_RATE, read_len, -1,  0, &read_coverage_start, &read_coverage_end, &indel_in_p1, &indel_in_p2, &max_cover_start, &max_cover_end, read_len, -1 , 0, NULL , 0xffffffff);

			//SUBREADprintf("RN=%s , WINDOW = %d ~ %d , SP=%d;  BV=%d;  BV2=%d\n", rname , window_cursor , window_cursor + read_len , splice_point, best_vote1, best_vote2);
			if (splice_point>0 && best_vote1 >= 1 && best_vote2>=1)
			{
				int test_real_break_point = -1, test_donor_score=-1;
				int is_GTAG = 0;
				int is_accepted = core13_test_donor(InBuff, read_len, min(best_pos1, best_pos2), max(best_pos1,best_pos2), splice_point, negative_strand, read_len/4, 0, 5, &test_real_break_point, base_index, 0, 0, negative_strand, color_space, &test_donor_score, &is_GTAG);

				if (is_accepted ){
					unsigned int pos_small = min(test_real_break_point+ best_pos1,  test_real_break_point+ best_pos2) - 1;
					unsigned int pos_big = max(test_real_break_point+ best_pos1,  test_real_break_point+ best_pos2);

					int event_no;
					chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];
					chromosome_event_t * found = NULL;

					int found_events = search_event(global_context, event_table, event_space, pos_small , EVENT_SEARCH_BY_SMALL_SIDE, CHRO_EVENT_TYPE_JUNCTION|CHRO_EVENT_TYPE_FUSION, search_return);

					if(found_events)
					{
						int kx1; 
						for(kx1 = 0; kx1 < found_events ; kx1++)
						{
							if(search_return[kx1] -> event_large_side == pos_big)
							{
								found = search_return[kx1];	
								break;
							}
						}
					}

					if(found) found -> supporting_reads ++;
					else
					{
						if(thread_context)
							event_no = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> total_events ++;
						else
							event_no = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) ->  total_events ++;

						event_space = reallocate_event_space(global_context, thread_context, event_no);

						chromosome_event_t * new_event = event_space+event_no; 
						memset(new_event,0,sizeof(chromosome_event_t));
						new_event -> event_small_side = pos_small;
						new_event -> event_large_side = pos_big;

						new_event -> is_negative_strand= !is_GTAG;
						new_event -> event_type = CHRO_EVENT_TYPE_JUNCTION;

						new_event -> supporting_reads = 1;
						new_event -> indel_length = 0;

						put_new_event(event_table, new_event , event_no);
			//			SUBREADprintf("ADD JUNCTION BY FRAGILE, %d-%d\n", pos_small, pos_big);
					}

				}

			}
		}
		InBuff[read_len] = tmp_char;
	}
}


void print_frags(global_context_t * global_context, fragment_list_t * fls){
	int x1;

	for(x1 =0; x1 < fls -> fragments; x1++){
		subread_read_number_t fno = fls -> fragment_numbers[x1] / 2;
		int f_is_B = fls -> fragment_numbers[x1] % 2;
		
		mapping_result_t * f_res = _global_retrieve_alignment_ptr(global_context, fno, f_is_B, 0);
		mapping_result_t * mate_res = _global_retrieve_alignment_ptr(global_context, fno, !f_is_B, 0);
		char outpos[100];
		char outposm[100];
		absoffset_to_posstr(global_context, f_res -> selected_position, outpos);
		absoffset_to_posstr(global_context, mate_res -> selected_position, outposm);

		int f_negative = (f_res -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;
		int mate_negative = (mate_res -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;

		if(f_is_B) f_negative=!f_negative;
		else mate_negative=!mate_negative;

		//SUBREADprintf("TRALOG: READ %09u %c AT %s (%c)  ;  MATE: %s (%c)\n", fno, f_is_B?'B':'A' , outpos, f_negative?'N':'P' , outposm, mate_negative?'N':'P');

	}
}

// fragnos_paired_B = B_fragment_no * 2 + is_mate_b   (is_mate_b points the mate that has the location in locations_mate_B)
// fragnos_paired_C = C_fragment_no * 2 + is_mate_c   (is_mate_c points the mate that has the location in locations_mate_C)
//
// locations_mate_B and locations_mate_C are the locations where the sequence is moved to. I.e., locations_mate_B and locations_mate_C are far far away from fragment A.
//
int find_translocation_BC_mates(global_context_t * global_context, mapping_result_t * res_A1, mapping_result_t * res_A2, fragment_list_t * listB, fragment_list_t * listC, int is_INV, unsigned long long * fragnos_paired_B, unsigned long long * fragnos_paired_C, unsigned int * locations_mate_B, unsigned int * locations_mate_C,unsigned int  * guessed_brkP_small_sum, unsigned int * guessed_moved_length_sum , unsigned int * guessed_brkQ_small_sum){

	int ret = 0, xk1, xk2;
	char * is_C_used = malloc(sizeof(char) * listC->fragments);
	memset(is_C_used, 0, sizeof(char) * listC->fragments);
	long long tmp_guessed_brkP_small_sum = 0, tmp_guessed_moved_length_sum = 0, tmp_guessed_brkQ_small_sum = 0;

	for(xk1 = 0; xk1 < listB->fragments; xk1++)
	{
		long long minimum_mate_distance = 0x7fffffff;
		int minimum_xk2 = -1;
		unsigned int mate_C_pos = 0;
		mapping_result_t * res_Ca = NULL, * res_Cc = NULL, * res_Ba = NULL, *res_Bb = NULL;
		mapping_result_t meta_C_res_body, res_Ca_body;
		res_Ca = &res_Ca_body;

		mapping_result_t * meta_C_res = &meta_C_res_body;

		subread_read_number_t B_read_no = listB->fragment_numbers[xk1]/2;
		int B_read_is_b = listB->fragment_numbers[xk1]%2;

		mapping_result_t meta_B_res_body, res_Ba_body;
		mapping_result_t * meta_B_res = &meta_B_res_body;
		res_Ba = &res_Ba_body;

		bigtable_readonly_result(global_context, NULL, B_read_no, 0, !B_read_is_b, meta_B_res, NULL);
		res_Bb = meta_B_res;

		bigtable_readonly_result(global_context, NULL, B_read_no, 0, B_read_is_b, res_Ba, NULL);

		for(xk2 = 0; xk2 < listC->fragments; xk2++)
		{
			if(is_C_used[xk2]) continue;

			subread_read_number_t C_read_no = listC->fragment_numbers[xk2]/2;
			int C_read_is_b = listC->fragment_numbers[xk2]%2;

			bigtable_readonly_result(global_context, NULL, C_read_no, 0, !C_read_is_b, meta_C_res, NULL);
			res_Cc = meta_C_res;

			bigtable_readonly_result(global_context, NULL, C_read_no, 0, C_read_is_b, res_Ca, NULL);

			int is_meta_B_negative = (meta_B_res  -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;
			if(!B_read_is_b) is_meta_B_negative = !is_meta_B_negative;

			int is_meta_C_negative = (meta_C_res  -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;
			if(!C_read_is_b) is_meta_C_negative = !is_meta_C_negative;

			//SUBREADprintf("TRALOG: MATES : B[%d] = %u (%c); C[%d] = %u (%c)\n", xk1, meta_B_res -> selected_position, is_meta_B_negative?'N':'P' , xk2, meta_C_res  -> selected_position, is_meta_C_negative?'N':'P');

			if(is_meta_B_negative != is_meta_C_negative &&
			   meta_B_res -> selected_position < meta_C_res -> selected_position &&
			   meta_C_res -> selected_position - meta_C_res -> selected_position < global_context -> config.maximum_translocation_length &&
			   meta_C_res -> selected_position - meta_B_res -> selected_position <  minimum_mate_distance)
			{
				minimum_mate_distance = meta_C_res -> selected_position - meta_B_res -> selected_position;
				minimum_xk2 = xk2;
				mate_C_pos =  meta_C_res -> selected_position;
			}
		}
		// read B has a mate of C[minimum xk2] if there is one.
		if(minimum_xk2>=0)
		{
			subread_read_number_t C_mate_fno = listC -> fragment_numbers[minimum_xk2] / 2;
			int C_mate_is_b = listC -> fragment_numbers[minimum_xk2] % 2;

			fragnos_paired_B[ret] = (B_read_no*2)+(!B_read_is_b);
			locations_mate_B[ret] = meta_B_res -> selected_position;
			
			fragnos_paired_C[ret] = (C_mate_fno*2)+(C_mate_is_b);
			locations_mate_C[ret] = mate_C_pos;

			is_C_used[minimum_xk2] = 1;


			int gapA, gapB, gapC;

			if(is_INV){
				gapA = res_Ca -> selected_position - res_A1 -> selected_position - res_A1 -> read_length;
				gapB = res_A2 -> selected_position - res_Ba -> selected_position - res_Ba -> read_length;
				gapC = res_Cc -> selected_position - res_Bb -> selected_position - res_Bb -> read_length;
			}else{
				gapA = res_Ba -> selected_position - res_A1 -> selected_position - res_A1 -> read_length;
				gapB = res_A2 -> selected_position - res_Ca -> selected_position - res_Ca -> read_length;
				gapC = res_Cc -> selected_position - res_Bb -> selected_position - res_Bb -> read_length;
			}

			tmp_guessed_brkP_small_sum += res_A1 -> selected_position + res_A1 -> read_length + gapA/2;
			tmp_guessed_moved_length_sum += res_A2 -> selected_position - res_A1 -> selected_position - res_A1 -> read_length - gapB/2 + gapA/2;
			tmp_guessed_brkQ_small_sum += res_Bb -> selected_position + res_Bb -> read_length + gapC/2;

			ret ++;
		}
	}

	free(is_C_used);

	if(ret>0){
		*guessed_brkP_small_sum= tmp_guessed_brkP_small_sum / ret;
		*guessed_moved_length_sum = tmp_guessed_moved_length_sum/ ret;
		*guessed_brkQ_small_sum = tmp_guessed_brkQ_small_sum / ret;
	}

	return ret;
}


// This function sees if all the mates of read B_x and C_y are at the same location.
// If mates of B_x and C_y spread on a large region, it is usually unreliable.
// posesB and posesB are linear absolute positions of the mate reads.
int find_translocation_BC_conformation(global_context_t * global_context, int PEmates, unsigned int  * posesB, unsigned int * posesC){

	unsigned int min_pos = 0xffffffff, max_pos = 0, xk1;
	if(PEmates<1) return 0;

	for(xk1 = 0; xk1 < PEmates; xk1++)
	{
		min_pos = min(min_pos, posesB[xk1]);
		min_pos = min(min_pos, posesC[xk1]);

		max_pos = max(max_pos, posesB[xk1]);
		max_pos = max(max_pos, posesC[xk1]);
	}

	if(max_pos - min_pos< 2*global_context -> config.maximum_pair_distance)return 1;
	return 0;
}


// fliB and fliB are : frag_[BC]_no * 2 + is_Read_b_close_to_BreakPoint_P
int breakpoint_PQR_supported(global_context_t * global_context , unsigned int brkPno , unsigned int brkQno, unsigned int brkRno, fragment_list_t * fliB, fragment_list_t * fliC, int isInv){
	int fli_i;
	int isFliB, nSupB=0, nSupC=0;

	for(isFliB = 0; isFliB < 2; isFliB++){
		fragment_list_t * fli = isFliB?fliB:fliC;
		int * nSup = isFliB?&nSupB:&nSupC;
		// fliB => support source_small ~ target_large if inv, or source_small ~ target_small if !inv
		// fliC => support source_large ~ target_small if inv, or source_large ~ target_large if !inv
		
		// the read that is close to BreakPoint_P should support source, the other read should support target
		for(fli_i = 0; fli_i < fli -> fragments; fli_i ++){
			subread_read_number_t frag_BC_no = fli -> fragment_numbers[fli_i]/2;
			int is_Read_b_close_to_BreakPoint_P = fli -> fragment_numbers[fli_i]%2;
			unsigned int source_small, source_large, target_smallQ, target_largeQ, target_smallR, target_largeR, target_large, target_small;

			get_event_two_coordinates(global_context, brkPno, NULL, NULL, &source_small, NULL, NULL, &source_large);
			get_event_two_coordinates(global_context, brkQno, NULL, NULL, &target_smallQ, NULL, NULL, &target_largeQ);
			get_event_two_coordinates(global_context, brkRno, NULL, NULL, &target_smallR, NULL, NULL, &target_largeR);


			if(target_smallQ <= target_smallR + BREAK_POINT_MAXIMUM_TOLERANCE && target_smallQ >= target_smallR - BREAK_POINT_MAXIMUM_TOLERANCE)
			{
				//target_smallQ is target, target_smallR is target
				target_large = target_smallR;
				target_small = target_smallQ;
			}else{

				//target_largeQ is target, target_largeR is target
				target_large = target_largeQ;
				target_small = target_largeR;
			}


			mapping_result_t res_BC_close_P_body, res_BC_close_Q_body;

			mapping_result_t * res_BC_close_P = &res_BC_close_P_body, * res_BC_close_Q = & res_BC_close_Q_body;

			bigtable_readonly_result(global_context, NULL, frag_BC_no, 0, is_Read_b_close_to_BreakPoint_P, res_BC_close_P, NULL);
			bigtable_readonly_result(global_context, NULL, frag_BC_no, 0, !is_Read_b_close_to_BreakPoint_P, res_BC_close_Q, NULL);

			unsigned int P_pos = isInv?( isFliB?source_large:source_small ):( isFliB?source_small:source_large );
			unsigned int Q_pos = isInv?( isFliB?target_large:target_small ):( isFliB?target_small:target_large );

			SUBREADprintf("TRALOG: PQR_TARGET P=%u~%u; Q=%u~%u, R=%u~%u ; Ppos=%u, Qpos=%u, Pread=%u, Qread=%u on %s\n", source_small, source_large, target_smallQ, target_largeQ, target_smallR, target_largeR,  P_pos, Q_pos, res_BC_close_P -> selected_position, res_BC_close_Q -> selected_position, isInv?"INV":"STR");

			long long dist;
			dist = res_BC_close_P -> selected_position;
			dist -= P_pos;
			if(abs(dist) < global_context -> config.maximum_pair_distance){
				dist = res_BC_close_Q -> selected_position;
				dist -= Q_pos;
				if(abs(dist) < global_context -> config.maximum_pair_distance)
					(*nSup)++;
			}
		}
	}
	//return nSupB + 1 >= fliB -> fragments/2 && nSupC + 1 >= fliC-> fragments/2 ;
	SUBREADprintf("TRALOG: PQR_NSUP: B=%d, C=%d on %s\n", nSupB, nSupC, isInv?"INV":"STR");
	return nSupB > 0 && nSupC > 0 && nSupB + 2 >= fliB->fragments / 2 && nSupC + 2 >= fliC->fragments / 2;
}

// fragnoD1_mates and fragnoD2_mates are poteltial E reads 1/2.
// D1: D's small read; D2: D's large read
// E2 ~ D2
// E1 ~ D1
// E2.start > Y.large
// E1.start > Y.small

int breakpoint_YZ_supported(global_context_t * global_context, unsigned int brkYno, unsigned int brkZno, unsigned long long * fragnoD1_mates, int fragnoD1len, unsigned long long * fragnoD2_mates, int fragnoD2len){
	int x1;
	int is_D2_mates;

	unsigned int inversion_small_edge, inversion_large_edge;
	get_event_two_coordinates(global_context, brkYno, NULL, NULL, &inversion_small_edge, NULL, NULL, &inversion_large_edge);


	int nSupD1mates = 0, nSupD2mates = 0;
	for(is_D2_mates = 0; is_D2_mates < 2; is_D2_mates ++){
		unsigned long long * fragno_Dmates = is_D2_mates?fragnoD2_mates:fragnoD1_mates;
		int fragno_Dno = is_D2_mates?fragnoD2len:fragnoD1len;
		int * nSupMates = is_D2_mates?&nSupD2mates:&nSupD1mates;
		for(x1 = 0; x1 < fragno_Dno; x1++){
			subread_read_number_t fragno_Dmate = fragno_Dmates[x1] / 2;
			int is_large_read_far_from_D  = fragno_Dmates[x1] % 2;

			mapping_result_t frag_D_mate_a_body, frag_D_mate_b_body;
			mapping_result_t * frag_D_mate_a = &frag_D_mate_a_body, * frag_D_mate_b = & frag_D_mate_b_body;

			bigtable_readonly_result(global_context, NULL, fragno_Dmate, 0, 0, frag_D_mate_a, NULL);
			bigtable_readonly_result(global_context, NULL, fragno_Dmate, 0, 1, frag_D_mate_b, NULL);

			mapping_result_t * frag_D_mate_1 = (frag_D_mate_a -> selected_position > frag_D_mate_b -> selected_position)?frag_D_mate_b:frag_D_mate_a;
			mapping_result_t * frag_D_mate_2 = (frag_D_mate_a -> selected_position <=frag_D_mate_b -> selected_position)?frag_D_mate_b:frag_D_mate_a;

			mapping_result_t * res_to_support_small_edge = (is_D2_mates ^ is_large_read_far_from_D)?frag_D_mate_2:frag_D_mate_1;
			mapping_result_t * res_to_support_large_edge = (is_D2_mates ^ is_large_read_far_from_D)?frag_D_mate_1:frag_D_mate_2;

			long long distsm;
			distsm = res_to_support_small_edge -> selected_position;
			distsm -= inversion_small_edge;

			long long distla;
			distla = res_to_support_large_edge -> selected_position;
			distla -=  inversion_large_edge;

			//SUBREADprintf("INVLOG: Dist_SM=%lld, Dist_LA=%lld\n", distsm, distla);

			if(distsm > -8 && distsm <  global_context -> config.maximum_pair_distance){
	
				if(distla > -8 && distla <  global_context -> config.maximum_pair_distance)
					(*nSupMates) ++;
			} 
			

		}
	}

	//SUBREADprintf("INVLOG: breakpoint_YZ_supported nSupD1=%d >= %d,  nSupD2=%d >= %d\n", nSupD1mates, fragnoD1len, nSupD2mates, fragnoD2len);
	return nSupD1mates > 0 && nSupD2mates > 0 && nSupD1mates + 2 >= fragnoD1len / 2 && nSupD2mates + 2 >= fragnoD2len / 2;
}

#define _PQR_LIST_SIZE 48 

int find_translocation_brk_PQR(global_context_t * global_context, mapping_result_t * resA1, mapping_result_t * resA2, fragment_list_t * fliB, fragment_list_t * fliC, unsigned int * brkPno,   unsigned int *  brkQno,  unsigned int *  brkRno, int isInv, unsigned int * is_cand_P_found)
{
	unsigned int event_pos_list_A1[_PQR_LIST_SIZE];
	void * event_ptr_list_A1[_PQR_LIST_SIZE];

	char * chroA=NULL;
	int posA1=0;
	
	locate_gene_position(resA1 -> selected_position,  &global_context -> chromosome_table, &chroA, &posA1);

	
	int candA1i, found_PQR = 0;
	int candA1Number = bktable_lookup(&global_context -> breakpoint_table_P, chroA, posA1, global_context -> config.maximum_pair_distance , event_pos_list_A1, event_ptr_list_A1, _PQR_LIST_SIZE);
	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 
	int candBrkPi , candBrkPNumber=0;

	//SUBREADprintf("A FOUND %d P ", candA1Number);

	for(candA1i = 0; candA1i < candA1Number ; candA1i++){
		int event_no = event_ptr_list_A1[candA1i] - NULL;
		chromosome_event_t * event_body = indel_context -> event_space_dynamic + event_no;

		long long small_dist = event_body -> event_small_side, large_dist = event_body -> event_large_side;
		small_dist -= resA1 -> selected_position;
		large_dist -= resA2 -> selected_position;
	}

	//SUBREADprintf(", (%d may be used)\n", candBrkPNumber);

	*is_cand_P_found = candBrkPNumber;

	for(candBrkPi = 0; candBrkPi < candBrkPNumber; candBrkPi++){
		unsigned int event_no_P = event_ptr_list_A1[candBrkPi] - NULL;
		chromosome_event_t * event_body_P = indel_context -> event_space_dynamic + event_no_P;

		unsigned int anchor_for_brkQ = isInv?event_body_P -> event_large_side:event_body_P -> event_small_side;
		unsigned int anchor_for_brkR = isInv?event_body_P -> event_small_side:event_body_P -> event_large_side;

		unsigned int event_pos_list_Q[_PQR_LIST_SIZE];
		void * event_ptr_list_Q[_PQR_LIST_SIZE];

		unsigned int event_pos_list_R[_PQR_LIST_SIZE];
		void * event_ptr_list_R[_PQR_LIST_SIZE];

		char * charAncQ = NULL, * charAncR = NULL;
		int posAncQ=0, posAncR = 0;
		locate_gene_position(anchor_for_brkQ, &global_context -> chromosome_table, &charAncQ, &posAncQ);
		locate_gene_position(anchor_for_brkR, &global_context -> chromosome_table, &charAncR, &posAncR);

		int candQi, candQnumber = bktable_lookup(&global_context -> breakpoint_table_QR, charAncQ, posAncQ - BREAK_POINT_MAXIMUM_TOLERANCE , 2* BREAK_POINT_MAXIMUM_TOLERANCE , event_pos_list_Q, event_ptr_list_Q, _PQR_LIST_SIZE);
		int candRi, candRnumber = bktable_lookup(&global_context -> breakpoint_table_QR, charAncR, posAncR - BREAK_POINT_MAXIMUM_TOLERANCE , 2* BREAK_POINT_MAXIMUM_TOLERANCE , event_pos_list_R, event_ptr_list_R, _PQR_LIST_SIZE);

		SUBREADprintf("P [%s] FOUND %d Q AT %s:%u and %d R AT %s:%u\n", isInv?"INV":"STR", candQnumber, charAncQ, posAncQ, candRnumber, charAncR, posAncR);

		for(candQi = 0 ; candQi < candQnumber ; candQi++){
			unsigned int event_no_Q = event_ptr_list_Q[candQi] - NULL;
			chromosome_event_t * event_body_Q = indel_context -> event_space_dynamic + event_no_Q;

			long long cand_Q_small_dist = event_body_Q -> event_small_side;
			cand_Q_small_dist -= isInv?event_body_P -> event_large_side:event_body_P -> event_small_side;

			int is_Q_small_side_close_to_P = abs(cand_Q_small_dist) <= BREAK_POINT_MAXIMUM_TOLERANCE;

			SUBREADprintf("Q: SMALL_CLOSE_P = %d, DIR = %c %c\n", is_Q_small_side_close_to_P,  event_body_Q -> small_side_increasing_coordinate?'>':'<', event_body_Q -> large_side_increasing_coordinate?'>':'<');

			if(  is_Q_small_side_close_to_P  && event_body_Q -> large_side_increasing_coordinate == 1) continue;  // the large side is the target location.
			if((!is_Q_small_side_close_to_P) && event_body_Q -> small_side_increasing_coordinate == 1) continue;  // the small side is the target location.


			if(  isInv  && event_body_Q -> large_side_increasing_coordinate != event_body_Q -> small_side_increasing_coordinate) continue;
			if((!isInv) && event_body_Q -> large_side_increasing_coordinate == event_body_Q -> small_side_increasing_coordinate) continue;

			for(candRi = 0 ; candRi < candRnumber ; candRi++){
				unsigned int event_no_R = event_ptr_list_R[candRi] - NULL;
				chromosome_event_t * event_body_R = indel_context -> event_space_dynamic + event_no_R;

				srInt_64 cand_R_dist_to_Q = is_Q_small_side_close_to_P?event_body_Q -> event_large_side:event_body_Q -> event_small_side;
				cand_R_dist_to_Q -= is_Q_small_side_close_to_P?event_body_R -> event_large_side:event_body_R-> event_small_side;

				SUBREADprintf("R: candDist=%lld, DIR = %c %c\n", cand_R_dist_to_Q,  event_body_Q -> small_side_increasing_coordinate?'>':'<', event_body_Q -> large_side_increasing_coordinate?'>':'<');

				if(abs(cand_R_dist_to_Q) > BREAK_POINT_MAXIMUM_TOLERANCE) continue;
				int is_R_small_side_close_to_P = is_Q_small_side_close_to_P;

				if(  is_R_small_side_close_to_P  && !event_body_R -> large_side_increasing_coordinate) continue;
				if(!(is_R_small_side_close_to_P) && !event_body_R -> small_side_increasing_coordinate) continue;

				if(  isInv  && event_body_R -> large_side_increasing_coordinate != event_body_R -> small_side_increasing_coordinate) continue;
				if(!(isInv) && event_body_R -> large_side_increasing_coordinate == event_body_R -> small_side_increasing_coordinate) continue;
				(*brkPno) = event_no_P;
				(*brkQno) = event_no_Q;
				(*brkRno) = event_no_R;
				found_PQR++;
				return 1;
			}
		}
	}

	return found_PQR;
}


void get_event_two_coordinates(global_context_t * global_context, unsigned int event_no, char ** small_chro, int * small_pos, unsigned int * small_abs, char ** large_chro, int * large_pos, unsigned int * large_abs){

	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 
	chromosome_event_t * event_body = indel_context -> event_space_dynamic + event_no;

	if(small_abs)(*small_abs) = event_body -> event_small_side;
	if(large_abs)(*large_abs) = event_body -> event_large_side;

	if(small_chro && small_pos)
		locate_gene_position(event_body -> event_small_side,  &global_context -> chromosome_table, small_chro, small_pos);
	if(large_chro && large_pos)
		locate_gene_position(event_body -> event_large_side,  &global_context -> chromosome_table, large_chro, large_pos);
}


void create_or_update_translocation_imprecise_result(global_context_t * global_context , unsigned int guessed_P_small, unsigned int guessed_tra_len, unsigned int guessed_Q_small , int paired_BC_reads, int isInv){

	char * brkPchr;
	int brkPsmall;
	void * trans_old_ptrs [_PQR_LIST_SIZE];
	unsigned int trans_old_poses [_PQR_LIST_SIZE];

	locate_gene_position(guessed_P_small,  &global_context -> chromosome_table, &brkPchr, &brkPsmall);

	int is_trans_found = 0, old_res_i, old_res_number = bktable_lookup(&global_context -> translocation_result_table, brkPchr, brkPsmall - BREAK_POINT_MAXIMUM_TOLERANCE, 2*BREAK_POINT_MAXIMUM_TOLERANCE, trans_old_poses, trans_old_ptrs, _PQR_LIST_SIZE);
	for(old_res_i = 0; old_res_i < old_res_number; old_res_i++){
		translocation_result_t * old_res = (translocation_result_t * )trans_old_ptrs[old_res_i];

		long long target_dist = old_res -> target_left_side;
		target_dist -= guessed_Q_small;

		if(abs(target_dist) < BREAK_POINT_MAXIMUM_TOLERANCE && isInv == old_res -> is_inv){
			target_dist = old_res -> length;
			target_dist -= guessed_tra_len;
			if(abs(target_dist) < BREAK_POINT_MAXIMUM_TOLERANCE){
				old_res -> all_sup_P ++;
				old_res -> max_sup_QR = max(old_res -> max_sup_QR , paired_BC_reads);
				is_trans_found = 1;
				break;
			}
		} 
	}

	if(0 == is_trans_found){
		translocation_result_t * new_res = malloc(sizeof(translocation_result_t));
		memset(new_res, 0, sizeof(translocation_result_t));
		new_res -> target_left_side = guessed_Q_small;
		new_res -> length = guessed_tra_len;
		new_res -> source_left_side = guessed_P_small;
		new_res -> is_precisely_called = 0;
		new_res -> all_sup_P = 1;
		new_res -> max_sup_QR = paired_BC_reads;
		new_res -> is_inv = isInv;

		bktable_append(&global_context -> translocation_result_table,brkPchr, brkPsmall, new_res);
	}

}

void create_or_update_translocation_result(global_context_t * global_context , unsigned int brkPno, unsigned int brkQno, unsigned int brkRno , int paired_BC_reads, int isInv){

	char *brkPchr, *brkQchr, *tmpchr;
	int brkPsmall, brkPlarge, brkQsmall, tmpint;
	unsigned int brkPabs_small, brkQabs_small, brkRabs_small, brkRabs_large, brkQabs_large;

	SUBREADprintf("\nTRALOG: FINALLY_CONFIRMED: %s ; %d PE_MATES\n", isInv?"INV":"STR", paired_BC_reads);

	get_event_two_coordinates(global_context, brkPno, &brkPchr, &brkPsmall, &brkPabs_small,  &tmpchr, &brkPlarge, NULL);
	get_event_two_coordinates(global_context, brkQno, &brkQchr, &brkQsmall, &brkQabs_small,  &tmpchr, &tmpint, &brkQabs_large);
	get_event_two_coordinates(global_context, brkRno, NULL, NULL, &brkRabs_small,  NULL, NULL, &brkRabs_large);

	SUBREADprintf("TRARES: %s:%u (len=%d) => %s:%u   (Coor: last_base_before)\n", brkPchr, brkPsmall, brkPlarge - brkPsmall - 1, brkQchr, brkQsmall);

	void * trans_old_ptrs [_PQR_LIST_SIZE];
	unsigned int trans_old_poses [_PQR_LIST_SIZE];

	unsigned int new_target_left_side, new_length;


	if(brkQabs_small >= brkRabs_small - BREAK_POINT_MAXIMUM_TOLERANCE && brkQabs_small <= brkRabs_small + BREAK_POINT_MAXIMUM_TOLERANCE)
	{
		// Q small and R large are target 
		new_target_left_side = brkQabs_small;
	} else{
		// Q large and R small are target 
		new_target_left_side = brkQabs_large;
	}
	
	new_length = brkPlarge - brkPsmall - 1;

	int is_trans_found = 0, old_res_i, old_res_number = bktable_lookup(&global_context -> translocation_result_table, brkPchr, brkPsmall - BREAK_POINT_MAXIMUM_TOLERANCE, 2*BREAK_POINT_MAXIMUM_TOLERANCE, trans_old_poses, trans_old_ptrs, _PQR_LIST_SIZE);
	for(old_res_i = 0; old_res_i < old_res_number; old_res_i++){
		translocation_result_t * old_res = (translocation_result_t * )trans_old_ptrs[old_res_i];

		long long target_dist = old_res -> target_left_side;
		target_dist -= new_target_left_side;

		if(abs(target_dist) < BREAK_POINT_MAXIMUM_TOLERANCE && isInv == old_res -> is_inv){
			target_dist = old_res -> length;
			target_dist -= new_length;
			if(abs(target_dist) < BREAK_POINT_MAXIMUM_TOLERANCE){
				old_res -> all_sup_P ++;
				old_res -> max_sup_QR = max(old_res -> max_sup_QR , paired_BC_reads);
				is_trans_found = 1;
				break;
			}
		} 
	}

	if(0 == is_trans_found){

		translocation_result_t * new_res = malloc(sizeof(translocation_result_t));
		memset(new_res, 0, sizeof(translocation_result_t));
		new_res -> target_left_side = new_target_left_side;
		new_res -> length = new_length;
		new_res -> source_left_side = brkPabs_small;
		new_res -> is_precisely_called = 1;
		new_res -> event_P_number = brkPno;
		new_res -> event_Q_number = brkQno;
		new_res -> event_R_number = brkRno;
		new_res -> all_sup_P = 1;
		new_res -> max_sup_QR = paired_BC_reads;
		new_res -> is_inv = isInv;

		bktable_append(&global_context -> translocation_result_table,brkPchr, brkPsmall, new_res);
	}
}


void finalise_translocations(global_context_t * global_context){

	void ** s1_ptrs, **s2_ptrs;
	unsigned int * s1_poses, * s2_poses;

	s1_ptrs = malloc(sizeof(void *) * S12_LIST_CAPACITY);
	s2_ptrs = malloc(sizeof(void *) * S12_LIST_CAPACITY);

	s1_poses = malloc(sizeof(int) * S12_LIST_CAPACITY);
	s2_poses = malloc(sizeof(int) * S12_LIST_CAPACITY);

	unsigned long long * s1_selected_list = malloc(sizeof(long long) * S12_LIST_CAPACITY);	// fragment_no * 2 + is_second_read
	unsigned long long * s2_selected_list = malloc(sizeof(long long) * S12_LIST_CAPACITY);

	mapping_result_t ** s1_result_ptr_list =  malloc(sizeof(mapping_result_t *) * S12_LIST_CAPACITY);
	mapping_result_t ** s2_result_ptr_list =  malloc(sizeof(mapping_result_t *) * S12_LIST_CAPACITY);

	int frag_Q_larger_read;
	subread_read_number_t frag_A_i;

	for(frag_A_i = 0; frag_A_i < global_context -> funky_list_A.fragments; frag_A_i ++){
		fragment_list_t fli_STR_B, fli_STR_C, fli_INV_B, fli_INV_C;

		fraglist_init(&fli_STR_B);
		fraglist_init(&fli_STR_C);
		fraglist_init(&fli_INV_B);
		fraglist_init(&fli_INV_C);

		subread_read_number_t frag_A_no = global_context -> funky_list_A.fragment_numbers[frag_A_i];

		mapping_result_t q_res_A_body, q_res_B_body;

		mapping_result_t * q_res_A = &q_res_A_body;
		mapping_result_t * q_res_B = &q_res_B_body;

		bigtable_readonly_result(global_context, NULL, frag_A_no, 0, 0, q_res_A, NULL);
		bigtable_readonly_result(global_context, NULL, frag_A_no, 0, 1, q_res_B, NULL);

		mapping_result_t * q_res_1 = q_res_A -> selected_position >  q_res_B -> selected_position?q_res_B:q_res_A;
		mapping_result_t * q_res_2 = q_res_A -> selected_position <= q_res_B -> selected_position?q_res_B:q_res_A;

		/***************************************************************************************************
 		 * 
 		 *  is_q1_negative and is_q2_negative describes the strandness of the original FASTQ read sequence.
 		 *
 		 *  For the very normal mappings, is_q1_negative must be 0 and is_q2_negative must be 1.
 		 *
 		 *  If is_q1_negative != is_q2_negative, then there is a strand-jumpping fusion between the two reads.
 		 */

		int is_q1_negative = (q_res_1 -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;
		int is_q2_negative = (q_res_2 -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;

		if(q_res_B == q_res_1)is_q1_negative=!is_q1_negative;
		if(q_res_B == q_res_2)is_q2_negative=!is_q2_negative;

		long long dist = q_res_A ->selected_position;
		dist -= q_res_B->selected_position;

		if( abs(dist) < 1000 && !(is_q1_negative == 0 && is_q2_negative == 1))
		{
			SUBREADprintf("TRALOG: STRANDNESS_BUG %08llu\n", frag_A_no);
		}


		for(frag_Q_larger_read = 0; frag_Q_larger_read < 2; frag_Q_larger_read++){
			void ** s_ptrs = frag_Q_larger_read?s2_ptrs:s1_ptrs;
			unsigned int * s_poses = frag_Q_larger_read?s2_poses:s1_poses;
			int q_res_offset = 0;
			mapping_result_t * q_res = frag_Q_larger_read?q_res_2:q_res_1;

			char * q_res_chro = NULL;
			locate_gene_position(q_res -> selected_position,  &global_context -> chromosome_table, &q_res_chro, &q_res_offset);
			q_res_offset +=1 ; // all tables are one-based.

			unsigned int q_search_start = q_res_offset;
			if(q_search_start > FUNKY_COLOCATION_TOLERANCE) q_search_start -= FUNKY_COLOCATION_TOLERANCE;
			else q_search_start = 0;

			int cand_i, canidate_s_items = bktable_lookup(&global_context -> funky_table_BC, q_res_chro, q_search_start, 2*FUNKY_COLOCATION_TOLERANCE, s_poses, s_ptrs, S12_LIST_CAPACITY);

			if(0 && frag_A_no == 143736){
				SUBREADprintf("TRALOG: SEARCH CLOSE TO %s READ: %s:%u ; HAD %d HITS\n", frag_Q_larger_read?"LARGE":"SMALL", q_res_chro, q_search_start, canidate_s_items);
			}

			// scan if candidate is reversed.
			// s_ptrs - NULL is the fragment no.
			for(cand_i = 0; cand_i < canidate_s_items; cand_i ++){
				subread_read_number_t frag_S_no = (s_ptrs[cand_i] - NULL)/ 2;
				int frag_S_is_read_B = (s_ptrs[cand_i] - NULL) % 2;

				mapping_result_t read_S_res_body, mate_S_res_body;
				mapping_result_t * read_S_res = &read_S_res_body;
				mapping_result_t * mate_S_res = &mate_S_res_body;

				bigtable_readonly_result(global_context, NULL, frag_S_no, 0, frag_S_is_read_B, read_S_res, NULL);
				bigtable_readonly_result(global_context, NULL, frag_S_no, 0, !frag_S_is_read_B, mate_S_res, NULL);

				int is_read_S_negative = (read_S_res -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;
				int is_mate_S_negative = (mate_S_res -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;
				if(frag_S_is_read_B) is_read_S_negative = !is_read_S_negative;
				else is_mate_S_negative = !is_mate_S_negative;

				int is_INV_TRA = is_mate_S_negative == is_read_S_negative;

				if(is_INV_TRA && is_read_S_negative == !frag_Q_larger_read){
					if(frag_Q_larger_read)
						fraglist_append(&fli_INV_B, frag_S_no * 2 + frag_S_is_read_B);
					else	
						fraglist_append(&fli_INV_C, frag_S_no * 2 + frag_S_is_read_B);
				}

				if((!is_INV_TRA) && is_read_S_negative == !frag_Q_larger_read){
					if(frag_Q_larger_read)
						fraglist_append(&fli_STR_C, frag_S_no * 2 + frag_S_is_read_B);
					else	
						fraglist_append(&fli_STR_B, frag_S_no * 2 + frag_S_is_read_B);
				}
			}
		}

		unsigned int guesed_p_small, guessed_tra_length, guessed_q_small, is_brkP_cand_found = 0;

		if(fli_INV_B.fragments >= 1 && fli_INV_C.fragments >= 1){
			int PEmates = find_translocation_BC_mates(global_context, q_res_1, q_res_2, &fli_INV_B, &fli_INV_C, 1, s1_selected_list, s2_selected_list, s1_poses, s2_poses, &guesed_p_small, &guessed_tra_length, &guessed_q_small);
			int ConformPE = find_translocation_BC_conformation(global_context, PEmates, s1_poses, s2_poses);
			int brkPQR_are_found = 0;
			unsigned int brkPno, brkQno, brkRno;

			char out1pos[100], out2pos[100];
			absoffset_to_posstr(global_context, q_res_1 -> selected_position, out1pos);
			absoffset_to_posstr(global_context, q_res_2 -> selected_position, out2pos);
			SUBREADprintf("TRALOG: A_READ: %09llu: INV : %s ~ %s ; %d PE_MATES (%s)\n", frag_A_no, out1pos, out2pos, PEmates, ConformPE?"CONFORMABLE":"INCONSISTENT");

			//SUBREADputs("TRALOG: INV_C:");
			//print_frags(global_context,&fli_INV_C);
			//SUBREADputs("TRALOG: INV_B:");
			//print_frags(global_context,&fli_INV_B);
			if(PEmates)
				brkPQR_are_found = find_translocation_brk_PQR(global_context, q_res_1, q_res_2, &fli_INV_B, &fli_INV_C, &brkPno, &brkQno, &brkRno, 1, &is_brkP_cand_found);

			if(brkPQR_are_found){
				brkPQR_are_found = breakpoint_PQR_supported(global_context , brkPno , brkQno, brkRno, &fli_INV_B, &fli_INV_C, 1);
				SUBREADprintf("TRALOG: A_READ: INV BRK_PQR_SUPPED=%d\n", brkPQR_are_found);
			}
			if(brkPQR_are_found)
				create_or_update_translocation_result( global_context , brkPno, brkQno, brkRno , PEmates, 1);
			else if(ConformPE && fli_INV_B.fragments > 2 && fli_INV_C.fragments > 2 && is_brkP_cand_found)
				create_or_update_translocation_imprecise_result(global_context, guesed_p_small, guessed_tra_length, guessed_q_small, PEmates, 1);
		}

		if(fli_STR_B.fragments >= 1 && fli_STR_C.fragments >= 1){
			int PEmates = find_translocation_BC_mates(global_context, q_res_1, q_res_2, &fli_STR_B, &fli_STR_C, 0, s1_selected_list, s2_selected_list, s1_poses, s2_poses, &guesed_p_small, &guessed_tra_length, &guessed_q_small);
			int ConformPE = find_translocation_BC_conformation(global_context, PEmates, s1_poses, s2_poses);

			char out1pos[100], out2pos[100];
			absoffset_to_posstr(global_context, q_res_1 -> selected_position, out1pos);
			absoffset_to_posstr(global_context, q_res_2 -> selected_position, out2pos);

			SUBREADprintf("TRALOG: A_READ: %09llu: TRA : %s ~ %s ; %d PE_MATES (%s)\n", frag_A_no, out1pos, out2pos, PEmates, ConformPE?"CONFORMABLE":"INCONSISTENT");

			//SUBREADputs("TRALOG: STR_B:");
			//print_frags(global_context,&fli_STR_B);
			//SUBREADputs("TRALOG: STR_C:");
			//print_frags(global_context,&fli_STR_C);

			int brkPQR_are_found = 0;
			unsigned int brkPno, brkQno, brkRno;

			if(PEmates)
				brkPQR_are_found = find_translocation_brk_PQR(global_context, q_res_1, q_res_2, &fli_STR_B, &fli_STR_C, &brkPno, &brkQno, &brkRno, 0, &is_brkP_cand_found);

			if(brkPQR_are_found){
				brkPQR_are_found = breakpoint_PQR_supported(global_context , brkPno , brkQno, brkRno, &fli_STR_B, &fli_STR_C, 0);
			}

			if(brkPQR_are_found)
				create_or_update_translocation_result( global_context , brkPno, brkQno, brkRno , PEmates, 0);
			else if(ConformPE && fli_INV_B.fragments > 2 && fli_INV_C.fragments > 2 && is_brkP_cand_found)
				create_or_update_translocation_imprecise_result(global_context, guesed_p_small, guessed_tra_length, guessed_q_small, PEmates, 0);
		}

		fraglist_destroy(&fli_STR_B);
		fraglist_destroy(&fli_STR_C);
		fraglist_destroy(&fli_INV_B);
		fraglist_destroy(&fli_INV_C);
	}
	
	free(s1_result_ptr_list);
	free(s2_result_ptr_list);
	free(s1_ptrs);
	free(s2_ptrs);
	free(s1_poses);
	free(s2_poses);
	free(s1_selected_list);
	free(s2_selected_list);

}

void finalise_inversions(global_context_t * global_context){
	subread_read_number_t frag_A_i;
	void ** s1_ptrs, **s2_ptrs;
	unsigned int * s1_poses, * s2_poses;

	s1_ptrs = malloc(sizeof(void *) * S12_LIST_CAPACITY);
	s2_ptrs = malloc(sizeof(void *) * S12_LIST_CAPACITY);

	s1_poses = malloc(sizeof(int) * S12_LIST_CAPACITY);
	s2_poses = malloc(sizeof(int) * S12_LIST_CAPACITY);

	unsigned long long * s1_selected_list = malloc(sizeof(long long) * S12_LIST_CAPACITY);	// fragment_no * 2 + is_second_read
	unsigned long long * s2_selected_list = malloc(sizeof(long long) * S12_LIST_CAPACITY);

	mapping_result_t ** s1_result_ptr_list =  malloc(sizeof(mapping_result_t *) * S12_LIST_CAPACITY);
	mapping_result_t ** s2_result_ptr_list =  malloc(sizeof(mapping_result_t *) * S12_LIST_CAPACITY);

	int frag_Q_larger_read, xk1, xk2;

	for(frag_A_i = 0; frag_A_i < global_context -> funky_list_DE.fragments; frag_A_i ++){
		int s1_list_items = 0, s2_list_items = 0;

		subread_read_number_t frag_A_no = global_context -> funky_list_DE.fragment_numbers[frag_A_i];

		mapping_result_t q_res_A_body, q_res_B_body;

		mapping_result_t * q_res_A = &q_res_A_body, * q_res_B = &q_res_B_body;

		bigtable_readonly_result(global_context, NULL, frag_A_no, 0, 0, q_res_A, NULL);
		bigtable_readonly_result(global_context, NULL, frag_A_no, 0, 1, q_res_B, NULL);

		mapping_result_t * q_res_1 = q_res_A -> selected_position >  q_res_B -> selected_position?q_res_B:q_res_A;
		mapping_result_t * q_res_2 = q_res_A -> selected_position <= q_res_B -> selected_position?q_res_B:q_res_A;


		/***************************************************************************************************
 		 * 
 		 *  is_q1_negative and is_q2_negative describes the strandness of the original FASTQ read sequence.
 		 *
 		 *  For the very normal mappings, is_q1_negative must be 0 and is_q2_negative must be 1.
 		 *
 		 *  If is_q1_negative != is_q2_negative, then there is a strand-jumpping fusion between the two reads.
 		 */

		int is_q1_negative = (q_res_1 -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;
		int is_q2_negative = (q_res_2 -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;

		if(q_res_B == q_res_1)is_q1_negative=!is_q1_negative;
		if(q_res_B == q_res_2)is_q2_negative=!is_q2_negative;

		if(is_q1_negative == 0 && is_q2_negative == 0)	// D READ
		{
			for(frag_Q_larger_read = 0; frag_Q_larger_read < 2; frag_Q_larger_read++){
				int * s_list_items = frag_Q_larger_read?&s2_list_items:&s1_list_items;
				void ** s_ptrs = frag_Q_larger_read?s2_ptrs:s1_ptrs;
				unsigned int * s_poses = frag_Q_larger_read?s2_poses:s1_poses;
				int q_res_offset = 0;
				mapping_result_t * q_res = frag_Q_larger_read?q_res_2:q_res_1;
				unsigned long long * s_selected_list = frag_Q_larger_read?s2_selected_list:s1_selected_list;
				mapping_result_t ** s_result_ptr_list = frag_Q_larger_read?s2_result_ptr_list:s1_result_ptr_list;


				char * q_res_chro = NULL;
				locate_gene_position(q_res -> selected_position,  &global_context -> chromosome_table, &q_res_chro, &q_res_offset);
				q_res_offset +=1 ; // all tables are one-based.

				unsigned int q_search_start = q_res_offset;
				if(q_search_start > FUNKY_COLOCATION_TOLERANCE) q_search_start -= FUNKY_COLOCATION_TOLERANCE;
				else q_search_start = 0;

				int cand_i, canidate_s_items = bktable_lookup(&global_context -> funky_table_DE, q_res_chro, q_search_start, 2*FUNKY_COLOCATION_TOLERANCE, s_poses, s_ptrs, S12_LIST_CAPACITY);
				// scan if candidate is reversed.
				// s_ptrs - NULL is the fragment no.
				for(cand_i = 0; cand_i < canidate_s_items; cand_i ++){
					subread_read_number_t frag_S_no = (s_ptrs[cand_i] - NULL)/2;
					int frag_S_larger_read = (s_ptrs[cand_i] - NULL)%2;

					if(frag_S_no == frag_A_no) continue;

					if(frag_S_larger_read == frag_Q_larger_read){

						mapping_result_t res_S_A_body, res_S_B_body;
						mapping_result_t * res_S_A = &res_S_A_body , * res_S_B = &res_S_B_body;

						bigtable_readonly_result(global_context, NULL, frag_S_no, 0, 0, res_S_A, NULL);
						bigtable_readonly_result(global_context, NULL, frag_S_no, 0, 1, res_S_B, NULL);

						mapping_result_t * res_S_1 = res_S_A -> selected_position >  res_S_B -> selected_position?res_S_B:res_S_A;
						mapping_result_t * res_S_2 = res_S_A -> selected_position <= res_S_B -> selected_position?res_S_B:res_S_A;

						mapping_result_t * co_locatted_S_res = frag_S_larger_read?res_S_2:res_S_1;

						int is_s1_negative = (res_S_1 -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;
						int is_s2_negative = (res_S_2 -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;

						if(res_S_B == res_S_1) is_s1_negative = !is_s1_negative;
						if(res_S_B == res_S_2) is_s2_negative = !is_s2_negative;


						if( is_s1_negative != 0 && is_s2_negative != 0 ){	// E READ
							s_selected_list[*s_list_items] = frag_S_no * 2 + frag_S_larger_read;
							s_result_ptr_list[*s_list_items] = co_locatted_S_res;
							(*s_list_items)++;
						}
					}
				}
			} 
		}

		int found_INV_frags = 0;
		srInt_64 guessed_Z_large_abs_sum = 0, guessed_Y_small_abs_sum = 0;

		for(xk1 = 0; xk1 < s1_list_items; xk1++){
			for(xk2 = 0; xk2 < s2_list_items ; xk2 ++){
				if(s1_selected_list[xk1]/2 == s2_selected_list[xk2]/2)
				{
					found_INV_frags ++;
					// now there is only one D fragment. here we found the E fragment for it (E fragment is in s1[xk1] and s2[xk2])
					// s1 is the E read that is close to D_1;  s2 is the E read that is close to D_2;   D_1 is the D read with smaller coordinate.
					// res_E1 is the read that is close to D_2; mapping location of E_1 should be larger than D_2

					mapping_result_t * res_D1 = q_res_1;
					mapping_result_t * res_D2 = q_res_2;

					mapping_result_t * res_E1 = s2_result_ptr_list[xk2];
					mapping_result_t * res_E2 = s1_result_ptr_list[xk1];

					int Gap_a_length = res_E2 -> selected_position - res_D1 -> selected_position - res_D1 -> read_length;
					int Gap_b_length = res_E1 -> selected_position - res_D2 -> selected_position - res_D2 -> read_length;
					int average_gap_len = (Gap_b_length + Gap_a_length)/2;
					guessed_Y_small_abs_sum += res_D1 -> selected_position + res_D1 -> read_length - average_gap_len / 2;
					guessed_Z_large_abs_sum += res_E1 -> selected_position - average_gap_len / 2;
					SUBREADprintf("INVLOG: GUESSED_LEN = %d + %d / 2 = %d\n", Gap_a_length, Gap_b_length, average_gap_len);
				}
			}
		}
		
		unsigned int brkYno=0xffffffff, brkZno=0xffffffff;
		int cand_YZ_breakpoints = 0;
		if(found_INV_frags > 0)
		{
			char * q_small_chro = NULL;
			int q_small_pos = 0;

			guessed_Y_small_abs_sum /= found_INV_frags;
			guessed_Z_large_abs_sum /= found_INV_frags;
			SUBREADprintf("INVLOG: GUESSED_YZ=%lld, %lld\n", guessed_Y_small_abs_sum, guessed_Z_large_abs_sum);

			locate_gene_position(q_res_1 -> selected_position,  &global_context -> chromosome_table, &q_small_chro, &q_small_pos);
			int cand_Y, cand_Z;
			cand_YZ_breakpoints = bktable_lookup(&global_context -> breakpoint_table_YZ, q_small_chro, q_small_pos, global_context -> config.maximum_pair_distance , s1_poses, s1_ptrs, S12_LIST_CAPACITY);

			//SUBREADprintf("INVLOG: %09u FOUND %d CANDIDATE BKs AT %s:%u\n", frag_A_no, cand_YZ_breakpoints, q_small_chro, q_small_pos);

			indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 

			for(cand_Y = 0; cand_Y < cand_YZ_breakpoints ; cand_Y ++){
				if(brkYno < 0xffffffff) break;

				int event_no_Y = s1_ptrs[cand_Y] - NULL;
				chromosome_event_t * event_body_Y = indel_context -> event_space_dynamic + event_no_Y;

				if(event_body_Y -> small_side_increasing_coordinate) continue; 
				if(event_body_Y -> small_side_increasing_coordinate != event_body_Y -> large_side_increasing_coordinate)
					assert(0);

				if(abs(event_body_Y -> event_large_side - q_res_2 -> selected_position) < global_context -> config.maximum_pair_distance){

					for(cand_Z = 0; cand_Z < cand_YZ_breakpoints ; cand_Z ++){
						int event_no_Z = s1_ptrs[cand_Z] - NULL;
						chromosome_event_t * event_body_Z = indel_context -> event_space_dynamic + event_no_Z;

						if(!event_body_Z -> small_side_increasing_coordinate) continue; 
						if(event_body_Z -> small_side_increasing_coordinate != event_body_Z -> large_side_increasing_coordinate)
							assert(0);

						long long dist_small = event_body_Z -> event_small_side , dist_large = event_body_Z -> event_large_side;
						dist_small -= event_body_Y -> event_small_side;
						dist_large -= event_body_Y -> event_large_side;

						long long dist_small_large_diff = dist_small;
						dist_small_large_diff -= dist_large;

						if(abs(dist_small_large_diff)  <= BREAK_POINT_MAXIMUM_TOLERANCE && abs(dist_large) <= BREAK_POINT_MAXIMUM_TOLERANCE && event_body_Z -> small_side_increasing_coordinate != event_body_Y -> small_side_increasing_coordinate){

							brkYno = event_no_Y;
							brkZno = event_no_Z;

							break;
						}
					}


					if(1)
					{
						char outpos1[100], outpos2[100];
						absoffset_to_posstr(global_context, event_body_Y -> event_small_side, outpos1);
						absoffset_to_posstr(global_context, event_body_Y -> event_large_side, outpos2);

						SUBREADprintf("INVLOG: %09llu FOUND BREAKPOINT YZ: %s ~ %s, INC_COR: %c %c , nSUP=%d\n", frag_A_no, outpos1, outpos2, event_body_Y -> small_side_increasing_coordinate?'>':'<', event_body_Y -> large_side_increasing_coordinate?'>':'<' , event_body_Y -> final_counted_reads);

					}

				}
			}
		}

		
		char *brkYchr = "NULL";
		unsigned int brkYabs_small = 0, brkYabs_large = 0;
		int brkYsmall = 0, brkYlarge = 0;
		int is_precisely_called = 0, is_roughly_called = 0;
		if(brkYno < 0xffffffff){
			// s1_selected_list : 2 * fragment_S_no + frag_S_larger_read
			int is_passed_YZ = breakpoint_YZ_supported(global_context, brkYno, brkZno, s1_selected_list, s1_list_items, s2_selected_list, s2_list_items); 
			if(is_passed_YZ)
			{
				is_precisely_called = 1;

				get_event_two_coordinates(global_context, brkYno, &brkYchr, &brkYsmall, &brkYabs_small, &brkYchr, &brkYlarge, &brkYabs_large);

			}
			else is_roughly_called = 1;
			//SUBREADprintf("\nINVLOG: FINALLY_%sCONFIRMED: %09u  %s:%u (len=%d) INVERSED.\n", is_passed_YZ?"":"NOT ", frag_A_no, brkYchr, brkYsmall, brkYlarge - brkYsmall);
		}

		//SUBREADprintf("\nINVLOG: FINALLY_GUESSED: %09u  found_INV_frags=%d, s1_list_items=%d, s2_list_items=%d, cand_YZ_breakpoints=%d\n", frag_A_no, found_INV_frags, s1_list_items, s2_list_items, cand_YZ_breakpoints);

		//for(xk1 = 0; xk1 < s1_list_items; xk1++) SUBREADprintf("INVLOG: %09d S_1 MATES: %09llu\n" , frag_A_no , s1_selected_list[xk1]/2);
		//for(xk1 = 0; xk1 < s2_list_items; xk1++) SUBREADprintf("INVLOG: %09d S_2 MATES: %09llu\n" , frag_A_no , s2_selected_list[xk1]/2);



		/*
		if(found_INV_frags >= min(s1_list_items , s2_list_items) - 2 && found_INV_frags > 1 && !is_precisely_called && cand_YZ_breakpoints>0){
			// guess brkYlarge, brkYsmall, brkZlarge, brkZsmall, brkYabsLarge, brkZabsLarge...
			locate_gene_position(guessed_Y_small_abs_sum,  &global_context -> chromosome_table, &brkYchr, &brkYsmall);
			locate_gene_position(guessed_Z_large_abs_sum,  &global_context -> chromosome_table, &brkYchr, &brkYlarge);
			//SUBREADprintf("\nINVLOG: FINALLY_GUESSED: %09u  %s:%u (len=%llu) INVERSED.\n", frag_A_no, brkYchr, brkYsmall, guessed_Z_large_abs_sum - guessed_Y_small_abs_sum);
			is_roughly_called = 1;
		}*/

		if( is_precisely_called || is_roughly_called )
		{
			void * old_ptrs[_PQR_LIST_SIZE];
			unsigned int old_poses[_PQR_LIST_SIZE];
			int old_found = 0, old_i, old_inversions = bktable_lookup(&global_context -> inversion_result_table, brkYchr, brkYsmall - BREAK_POINT_MAXIMUM_TOLERANCE, 2*BREAK_POINT_MAXIMUM_TOLERANCE, old_poses, old_ptrs, _PQR_LIST_SIZE);
			for(old_i = 0; old_i < old_inversions; old_i ++){
				inversion_result_t * inv_res_old = (inversion_result_t *) old_ptrs[old_i];
				long long old_dist = inv_res_old -> length;
				old_dist -= brkYlarge - brkYsmall;	// the difference on inversion length.
				if(abs(old_dist) < BREAK_POINT_MAXIMUM_TOLERANCE){
					inv_res_old -> all_sup_D ++;
					inv_res_old -> max_sup_E = max(inv_res_old -> max_sup_E , found_INV_frags);
					old_found = 1;
					break;
				}
			}

			if(0 == old_found){
				inversion_result_t * inv_res_new = malloc(sizeof(chromosome_event_t));
				memset(inv_res_new, 0 , sizeof(chromosome_event_t));

				inv_res_new -> length = brkYlarge - brkYsmall;
				inv_res_new -> is_precisely_called = is_precisely_called;
				if(is_precisely_called){
					inv_res_new -> event_Y_number = brkYno;
					inv_res_new -> event_Z_number = brkZno;
					inv_res_new -> small_side = brkYabs_small;
				}else{
					inv_res_new -> event_Y_rough_small_abs = guessed_Y_small_abs_sum;
					inv_res_new -> event_Z_rough_large_abs = guessed_Z_large_abs_sum;
					inv_res_new -> small_side = guessed_Y_small_abs_sum;
				}
				inv_res_new -> all_sup_D = 1;
				inv_res_new -> max_sup_E = found_INV_frags;

				bktable_append(&global_context -> inversion_result_table, brkYchr, brkYsmall, inv_res_new);
			}
		}
	}

	free(s1_result_ptr_list);
	free(s2_result_ptr_list);
	free(s1_ptrs);
	free(s2_ptrs);
	free(s1_poses);
	free(s2_poses);
	free(s1_selected_list);
	free(s2_selected_list);
}

void build_breakpoint_tables(global_context_t  * global_context){

	int xk1;
	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 

	for(xk1 = 0; xk1 < indel_context -> total_events ; xk1++)
	{ 
		char * chro_name_left= NULL,* chro_name_right = NULL;
		int chro_pos_left= 0, chro_pos_right = 0;

		chromosome_event_t * event_body = indel_context -> event_space_dynamic + xk1;

		if(event_body -> event_type != CHRO_EVENT_TYPE_FUSION && event_body -> event_type != CHRO_EVENT_TYPE_JUNCTION)
			continue;

		locate_gene_position(event_body -> event_small_side,  &global_context -> chromosome_table, &chro_name_left, &chro_pos_left);
		locate_gene_position(event_body -> event_large_side,  &global_context -> chromosome_table, &chro_name_right, &chro_pos_right);

		long long dist = chro_pos_left;
		dist -= chro_pos_right;
		if(dist<0)dist=-dist;

		int breakpoint_group = -1;

		if(event_body -> is_strand_jumped){
			// breakpoint QR or YZ
			if(chro_name_left != chro_name_right || dist > global_context -> config.maximum_translocation_length)
				breakpoint_group = 2;	// QR
			else
				breakpoint_group = 3;	// YZ
		}else{
			// breakpoint QR or P
			if(chro_name_left != chro_name_right || dist > global_context -> config.maximum_translocation_length)
				breakpoint_group = 2;	// QR
			else
				breakpoint_group = 1;	// P 
		}


		bucketed_table_t * index_table = breakpoint_group == 1?
							&global_context -> breakpoint_table_P :
							(breakpoint_group == 2?
								&global_context -> breakpoint_table_QR:
								(breakpoint_group == 3?
									&global_context -> breakpoint_table_YZ:
									NULL
								)
							);

		//SUBREADprintf("BPLOG: %s:%u ~ %s:%u (%c) GRP=%d (%p)\n", chro_name_left, chro_pos_left, chro_name_right, chro_pos_right, event_body -> is_strand_jumped?'X':'=', breakpoint_group, index_table);

		if(index_table)	bktable_append(index_table, chro_name_left, chro_pos_left, NULL + xk1);
		if(index_table)	bktable_append(index_table, chro_name_right, chro_pos_right, NULL + xk1);
	}
}

void finalise_structural_variances(global_context_t * global_context){
	SUBREADprintf("Funky Tables: A:%u, BC:%u, DE:%u\n", (unsigned int) global_context -> funky_list_A.fragments, (unsigned int)global_context -> funky_table_BC.fragments / 2, (unsigned int)global_context -> funky_list_DE.fragments);

	build_breakpoint_tables(global_context);
	SUBREADprintf("Breakpoint Tables: P:%u, QR:%u, YZ:%u\n", (unsigned int)global_context -> breakpoint_table_P.fragments, (unsigned int)global_context -> breakpoint_table_QR.fragments, (unsigned int)global_context -> breakpoint_table_YZ.fragments);
	finalise_translocations(global_context);
	finalise_inversions(global_context);
}
