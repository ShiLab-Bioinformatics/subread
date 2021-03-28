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
#include <math.h>
#include <signal.h>
#include <dirent.h> 
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include "hashtable.h"
#include "HelperFunctions.h"
#include "gene-value-index.h"
#include "gene-algorithms.h"
#include "input-files.h"
#include "core.h"
#include "core-indel.h"
#include "core-junction.h"
#include "core-bigtable.h"
#include "sublog.h"


//#define indel_debug
//#define is_target_window_X(x) 0

int DIGIN_COUNTER=0;
int FRESH_COUNTER=0;
extern unsigned int BASE_BLOCK_LENGTH;

void is_bad_read(char * r, int len, char * w)
{
}

void merge_event_sides(void * arr, int start, int items, int items2)
{
	void ** arrl  =  (void **)arr;
	chromosome_event_t * event_space = (chromosome_event_t*) arrl[1];
	int * event_id_list = (int *)arrl[0];
	int is_small_side = (arrl[2]!=NULL);

	int * merge_tmp = malloc(sizeof(int)* (items+items2));
	int x1, first_cursor = start, second_cursor = start+items;

	for(x1=0;x1<items+items2;x1++)
	{
		int choose_2 = 0;
		if(first_cursor >= start+items)
			choose_2 = 1;
		else if(second_cursor < start+items+items2)
		{

			chromosome_event_t * levent = event_space + event_id_list[first_cursor];
			chromosome_event_t * revent = event_space + event_id_list[second_cursor];

			if(is_small_side){
				if(levent->event_small_side > revent ->event_small_side)
					choose_2 = 1;
			}else
			{
				if(levent->event_large_side > revent ->event_large_side)
					choose_2 = 1;
			}
		
		}

		if(choose_2)
			merge_tmp[x1] = event_id_list[second_cursor++];
		else
			merge_tmp[x1] = event_id_list[first_cursor++];
	
	}

	memcpy(event_id_list + start, merge_tmp, sizeof(int)* (items+items2));

	free(merge_tmp);
}

void exchange_event_sides(void * arr, int l, int r)
{
	void ** arrl  =  (void **)arr;
	int * event_id_list = (int *)arrl[0];

	int tmp ;
	tmp = event_id_list[l];
	event_id_list[l] = event_id_list[r];
	event_id_list[r] = tmp;
	
}

int compare_event_inner(chromosome_event_t * event_space, int l, int r, int is_small_side)
{

	chromosome_event_t * revent = event_space +r;
	chromosome_event_t * levent = event_space +l;
	if(is_small_side){
		if(levent->event_small_side <revent ->event_small_side) return -1;
		else if(levent->event_small_side ==revent ->event_small_side) return 0;
		else return 1;
	}else
	{
		if(levent->event_large_side <revent ->event_large_side) return -1;
		else if(levent->event_large_side ==revent ->event_large_side) return 0;
		else return 1;
	}


}

int compare_event_sides(void * arr, int l, int r)
{
	void ** arrl  =  (void **)arr;
	chromosome_event_t * event_space = (chromosome_event_t*) arrl[1];
	int * event_id_list = (int *)arrl[0];
	int is_small_side = (arrl[2]!=NULL);

	return compare_event_inner(event_space, event_id_list[l], event_id_list[r], is_small_side);
}

int BINsearch_event(chromosome_event_t * event_space, int * event_ids, int is_small_side, unsigned int pos, int events)
{
	int imax, imin;
	assert(events > 0);
	imax = events - 1; imin = 0;

	while (1)	
	{
		int imid = (imin+ imax)/2;
		unsigned tested_mid_pos;
		//assert(imax<events);
		//assert(imin<events);
		//assert(imid<events);
		assert(event_ids[imid]<events);
		chromosome_event_t * mid_event = event_space+event_ids[imid];
		tested_mid_pos = is_small_side?mid_event->event_small_side:mid_event->event_large_side;

		if(tested_mid_pos == pos)
			return imid;
		else if(tested_mid_pos < pos)
			imin = imid +1;
		else	imax = imid -1;

		if(imin > imax)
			return imax;
	}
}


typedef struct{
	int thread_id;
	int block_start;
	int block_end;
	HashTable * result_tab;
	int * small_side_ordered_event_ids, * large_side_ordered_event_ids;
	chromosome_event_t * event_space;
	global_context_t * global_context;
} AT_context_t;


#define ANTI_SUPPORTING_READ_LIMIT 100

void * anti_support_thread_run(void * v_cont){
	AT_context_t * atcont = v_cont;

	int * cancelled_event_list = malloc(sizeof(int)*ANTI_SUPPORTING_READ_LIMIT), x2;
	int * small_side_ordered_event_ids = atcont -> small_side_ordered_event_ids, * large_side_ordered_event_ids = atcont -> large_side_ordered_event_ids;
	global_context_t * global_context = atcont -> global_context;
	chromosome_event_t * event_space = atcont -> event_space;
	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 
	if(indel_context->total_events<1)return NULL;

	int current_read_number, cancelled_events, x1;

	for(current_read_number = atcont -> block_start; current_read_number < atcont -> block_end; current_read_number++)
	{
		int is_second_read;
		for (is_second_read = 0; is_second_read < 1 + global_context -> input_reads.is_paired_end_reads; is_second_read ++)
		{
			int best_read_id;
			for(best_read_id = 0; best_read_id < global_context -> config.multi_best_reads; best_read_id++)
			{
				mapping_result_t *current_result = _global_retrieve_alignment_ptr(global_context, current_read_number, is_second_read, best_read_id);
				if(current_result -> selected_votes<1) break;
				if(!global_context->config.report_multi_mapping_reads)if(current_result -> result_flags & CORE_IS_BREAKEVEN) continue;
				//if(current_result -> result_flags & CORE_IS_GAPPED_READ) continue;
				if(current_result->selected_votes < global_context->config.minimum_subread_for_first_read)
					continue;

				cancelled_events=0;
				unsigned int mapping_head = current_result -> selected_position;
				unsigned int coverage_start = mapping_head + current_result -> confident_coverage_start;
				unsigned int coverage_end = mapping_head + current_result -> confident_coverage_end;

				int small_side_left_event = BINsearch_event(event_space, small_side_ordered_event_ids,1, coverage_start - 1, indel_context->total_events) + 1;
				int large_side_left_event = BINsearch_event(event_space, large_side_ordered_event_ids,0, coverage_start - 1, indel_context->total_events) + 1;

				int small_side_right_event = BINsearch_event(event_space, small_side_ordered_event_ids,1, coverage_end, indel_context->total_events)+20;
				int large_side_right_event = BINsearch_event(event_space, large_side_ordered_event_ids,0, coverage_end, indel_context->total_events)+20;

				for(x1 = small_side_left_event ; x1 <= small_side_right_event ; x1++)
				{
					if(x1 >= indel_context->total_events) break;
					if(ANTI_SUPPORTING_READ_LIMIT <= cancelled_events) break;
					chromosome_event_t * event_body = event_space+ small_side_ordered_event_ids[x1];
					if(event_body -> event_small_side <= coverage_start + 5) continue;
					if(event_body -> event_small_side >= coverage_end - 5) continue;

					long long cur_count = HashTableGet(atcont -> result_tab , NULL+small_side_ordered_event_ids[x1]+1 ) - NULL;
					cur_count++;
					HashTablePut( atcont -> result_tab , NULL+small_side_ordered_event_ids[x1]+1  , NULL+cur_count );

					cancelled_event_list[cancelled_events++] = small_side_ordered_event_ids[x1];
				}

				for(x1 = large_side_left_event ; x1 <= large_side_right_event ; x1++)
				{
					if(x1 >= indel_context->total_events) break;
					chromosome_event_t * event_body = event_space+ large_side_ordered_event_ids[x1];
					if(event_body -> event_large_side <= coverage_start + 5) continue;
					if(event_body -> event_large_side >= coverage_end- 5) continue;

					int to_be_add = 1;
					for(x2=0; x2<cancelled_events; x2++)
						if(cancelled_event_list[x2] == large_side_ordered_event_ids[x1])
						{
							to_be_add = 0;
							break;
						}

					if(to_be_add){
					//	printf("OCT27-ANTISUP-READ @ %u has LARGE @ %u~%u, INDELS = %d, ASUP = %d\n", coverage_end, event_body -> event_small_side,  event_body -> event_large_side, event_body -> indel_length, event_body -> anti_supporting_reads);
						long long cur_count = HashTableGet(atcont -> result_tab , NULL+large_side_ordered_event_ids[x1]+1 ) - NULL;
						cur_count++;
						HashTablePut( atcont -> result_tab , NULL+large_side_ordered_event_ids[x1]+1  , NULL+cur_count );
					}
				}
			}
		}
	} 

	free(cancelled_event_list);
	return NULL;
}


void anti_support_add_count(void * ky, void * va, HashTable * tab){
	chromosome_event_t * event_space = tab -> appendix1 ;
	int eno = ky-NULL-1;
	chromosome_event_t * event_body = event_space + eno;
	event_body -> anti_supporting_reads+= (va - NULL);
}

int anti_supporting_read_scan(global_context_t * global_context)
{

	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 
	if(indel_context->total_events<1)return 0;

	chromosome_event_t * event_space = indel_context -> event_space_dynamic;

	int x1, * small_side_ordered_event_ids, * large_side_ordered_event_ids;

	small_side_ordered_event_ids = malloc(sizeof(int)*indel_context -> total_events);
	large_side_ordered_event_ids = malloc(sizeof(int)*indel_context -> total_events);

	for(x1=0; x1<indel_context->total_events; x1++) {
		small_side_ordered_event_ids[x1]=x1;
		large_side_ordered_event_ids[x1]=x1;
	}

	void * sort_data[3];
	sort_data[0] = small_side_ordered_event_ids;
	sort_data[1] = event_space;
	sort_data[2] = (void *)0xffff;	// small_side_ordering

	merge_sort(sort_data, indel_context->total_events, compare_event_sides, exchange_event_sides, merge_event_sides);

	sort_data[0] = large_side_ordered_event_ids;
	sort_data[2] = NULL;	// large_side_ordering

	merge_sort(sort_data, indel_context->total_events, compare_event_sides, exchange_event_sides, merge_event_sides);


	AT_context_t ATconts[64];
	pthread_t AThreads[64];
	
	int thread_no, last_end = 0;
	for(thread_no=0; thread_no< global_context -> config.all_threads; thread_no++){
		AT_context_t * atc = ATconts + thread_no;
		atc -> thread_id = thread_no;
		atc -> block_start = last_end;
		atc -> block_end = last_end = global_context -> processed_reads_in_chunk / global_context -> config.all_threads * thread_no;
		if(thread_no == global_context -> config.all_threads-1)  atc -> block_end = global_context -> processed_reads_in_chunk;

		atc -> global_context = global_context;
		atc -> result_tab = HashTableCreate(200000);
		atc -> small_side_ordered_event_ids = small_side_ordered_event_ids;
		atc -> large_side_ordered_event_ids = large_side_ordered_event_ids;
		atc -> event_space = event_space;

		pthread_create(AThreads + thread_no , NULL, anti_support_thread_run, atc);
	}

	for(thread_no=0; thread_no< global_context -> config.all_threads; thread_no++){
		pthread_join(AThreads[thread_no], NULL);
		ATconts[thread_no].result_tab -> appendix1 = event_space;
		HashTableIteration( ATconts[thread_no].result_tab, anti_support_add_count );
		HashTableDestroy(ATconts[thread_no].result_tab);
	}

	free(small_side_ordered_event_ids);
	free(large_side_ordered_event_ids);

	return 0;
}

chromosome_event_t * reallocate_event_space( global_context_t* global_context,thread_context_t * thread_context,int event_no)
{
	int max_event_no;

	if(thread_context) {
		max_event_no = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> current_max_event_number;
		if(max_event_no<=event_no)
		{
			((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> current_max_event_number *= 1.6;
			 ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic = 
				realloc(((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic, sizeof(chromosome_event_t) * 
					((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> current_max_event_number);
		}
		return ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	} else {
		max_event_no = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> current_max_event_number;
		if(max_event_no<=event_no)
		{
			((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> current_max_event_number *= 1.6;
			((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic =
				realloc(((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic, sizeof(chromosome_event_t) *
					((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> current_max_event_number); 
		}
		return ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}

}

int is_ambiguous_indel_score(chromosome_event_t * e)
{
	return 0;
	//if(e -> indel_length == 1 && e->is_ambiguous>4)return 1;
	//if(e -> indel_length <= 3 && e->is_ambiguous>3 && (e-> event_small_side % 7 < 3))return 1;
	//if(e -> indel_length == 4) return e->is_ambiguous>4;
	//return 0;
}

int event_neighbour_sort_compare(void * arr, int l, int r){
	unsigned int ** sort_data = (unsigned int **) arr;
	if(sort_data[1][l] > sort_data[1][r]) return 1;
	if(sort_data[1][l] < sort_data[1][r]) return -1;
	return 0;
}

void event_neighbour_sort_exchange(void * arr, int l, int r){
	unsigned int ** sort_data = (unsigned int **) arr;
	unsigned int tmpi;

	tmpi = sort_data[0][l];
	sort_data[0][l] = sort_data[0][r];
	sort_data[0][r] = tmpi;

	tmpi = sort_data[1][l];
	sort_data[1][l] = sort_data[1][r];
	sort_data[1][r] = tmpi;
}

// smallFirst
void event_neighbour_sort_merge(void * arr, int start, int items, int items2){
	unsigned int ** sort_data = (unsigned int **) arr;

	unsigned int * tmp_global_id_list = malloc(sizeof(int) * (items+items2));
	unsigned int * tmp_offset_list = malloc(sizeof(int) * (items+items2));

	int i1_cursor = start, i2_cursor = items + start;
	int tmp_cursor = 0;

	while(1){
		if(i1_cursor == items + start && i2_cursor == items + items2 + start )break;
		int select_items_1 = (i1_cursor < items + start && event_neighbour_sort_compare(arr, i1_cursor, i2_cursor) <= 0) || (i2_cursor == start + items + items2);
		if(select_items_1){
			tmp_global_id_list[tmp_cursor] = sort_data[0][i1_cursor];
			tmp_offset_list[tmp_cursor ++ ] = sort_data[1][i1_cursor ++];
		}else{
			tmp_global_id_list[tmp_cursor] = sort_data[0][i2_cursor];
			tmp_offset_list[tmp_cursor ++ ] = sort_data[1][i2_cursor ++];
		}
	}

	memcpy( sort_data[0] + start, tmp_global_id_list, sizeof(int) * (items+items2) );
	memcpy( sort_data[1] + start, tmp_offset_list, sizeof(int) * (items+items2) );
	free(tmp_global_id_list);
	free(tmp_offset_list);
}

int test_redundant_event( global_context_t * global_context, chromosome_event_t* testee, chromosome_event_t* neighbour ){
	int length_delta = 3, is_same = 0;
	if(testee -> event_type == CHRO_EVENT_TYPE_INDEL) length_delta = 0;

	long long dist = testee -> event_large_side;
	dist -= neighbour -> event_large_side;

	if(abs(dist) <= length_delta){
		if(testee -> event_type == CHRO_EVENT_TYPE_INDEL){
			int indel_diff = testee -> indel_length;
			indel_diff -= neighbour -> indel_length;

			if(abs(indel_diff) <= length_delta)
				is_same = 1;
		}else is_same = 1;
	}

	int is_redund = 0;
	if(is_same){
		if(testee -> indel_at_junction > neighbour -> indel_at_junction) is_redund = 1;
		else if(testee -> indel_at_junction == neighbour -> indel_at_junction){
			if(testee -> supporting_reads < neighbour -> supporting_reads) is_redund = 1;
			else if( testee -> supporting_reads == neighbour -> supporting_reads ){
				if(testee -> event_small_side < neighbour -> event_small_side) is_redund = 1; 
			}
		}
	}

	return is_redund;
}


#define reallocate_to_be_removed_ids		if(to_be_removed_number >= remove_id_size-1){\
							remove_id_size = remove_id_size*3/2;\
							to_be_removed_ids = realloc(to_be_removed_ids, remove_id_size);\
						}

#define MIN_EVENT_SUPPORT_NO 2

void remove_neighbour(global_context_t * global_context)
{
	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 

	HashTable * event_table = indel_context -> event_entry_table;
	chromosome_event_t * event_space = indel_context -> event_space_dynamic;
	int xk1,xk2,xk3;
	int * to_be_removed_ids;
	int to_be_removed_number = 0, all_junctions = 0;
	int remove_id_size = 999999;

	to_be_removed_ids = malloc(sizeof(int) * (1+remove_id_size));
	for(xk1=0; xk1<indel_context->total_events; xk1++)
	{
		chromosome_event_t * event_body = &event_space[xk1];
		if(CHRO_EVENT_TYPE_REMOVED == event_body->event_type) continue;
		if(global_context -> config.do_remove_neighbour_for_scRNA && global_context -> config.scRNA_input_mode && event_body -> supporting_reads < MIN_EVENT_SUPPORT_NO){
			reallocate_to_be_removed_ids;
			to_be_removed_ids[to_be_removed_number++] = event_body -> global_event_id;
			continue;
		}
		all_junctions++;
		int event_type;
		for(event_type = 0; event_type < 2; event_type++)
		{
			if(event_type) // for INDELs
			{
				int neighbour_range = 3;
				if(event_body->event_type != CHRO_EVENT_TYPE_INDEL) continue;

				reallocate_to_be_removed_ids;
				if(is_ambiguous_indel_score(event_body))to_be_removed_ids[to_be_removed_number++] = event_body -> global_event_id;
				else for(xk2=-neighbour_range; xk2<=neighbour_range; xk2++)
				{
					//if(!xk2)continue;
					unsigned int test_pos_small = event_body -> event_small_side + xk2;
					chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];
					
					int found_events = search_event(global_context, event_table, event_space, test_pos_small  , EVENT_SEARCH_BY_SMALL_SIDE, CHRO_EVENT_TYPE_INDEL, search_return);
					//if(test_pos_small > 284136 && test_pos_small < 284146) printf("POS=%u\tTRES=%d\n", test_pos_small, found_events);

					for(xk3 = 0; xk3<found_events; xk3++)
					{
						reallocate_to_be_removed_ids;
						chromosome_event_t * tested_neighbour = search_return[xk3];
						long long int length_diff = tested_neighbour -> indel_length;
						length_diff -=  event_body -> indel_length;
						if(length_diff == 0 && xk2==0) continue;
						if(is_ambiguous_indel_score(tested_neighbour))
							to_be_removed_ids[to_be_removed_number++] = event_body -> global_event_id;
						else
							if(length_diff==0)
							{
								if(event_body -> event_small_side - event_body->connected_previous_event_distance + 1 == tested_neighbour-> event_large_side) continue;
								if(tested_neighbour -> event_small_side - tested_neighbour->connected_previous_event_distance + 1 == event_body-> event_large_side) continue;

								if(event_body -> event_large_side + event_body->connected_next_event_distance - 1 == tested_neighbour-> event_small_side) continue;
								if(tested_neighbour -> event_large_side + tested_neighbour->connected_next_event_distance - 1 == event_body-> event_small_side) continue;

								if(event_body->event_quality < tested_neighbour -> event_quality)
									to_be_removed_ids[to_be_removed_number++] = event_body -> global_event_id;
								else if(event_body->event_quality ==  tested_neighbour -> event_quality &&
									(event_body -> supporting_reads < tested_neighbour -> supporting_reads||(event_body -> supporting_reads == tested_neighbour -> supporting_reads && xk2<0)))
								{
									to_be_removed_ids[to_be_removed_number++] = event_body -> global_event_id;
								}
							}	
					}
				}
			}
			else // for junctions/fusions
			{
				int neighbour_range = 11;
				int indel_range = 4 , delta_small = 0;// max_range = 0;
				if(event_body -> event_type == CHRO_EVENT_TYPE_INDEL) continue;
				if(event_body -> is_donor_found_or_annotation & 64) continue; // do not remove known events
				//if(! global_context->config.check_donor_at_junctions )max_range = 10;

				//for(delta_small = - max_range; delta_small <= max_range ; delta_small ++)
					for(xk2=-neighbour_range; xk2<=neighbour_range; xk2++)
					{
						if(!xk2)continue;
						unsigned int test_pos_small = event_body -> event_small_side + xk2;
						chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];
						
						int found_events = search_event(global_context,event_table, event_space, test_pos_small + delta_small , EVENT_SEARCH_BY_SMALL_SIDE, CHRO_EVENT_TYPE_JUNCTION|CHRO_EVENT_TYPE_FUSION, search_return);
						for(xk3 = 0; xk3<found_events; xk3++)
						{

							reallocate_to_be_removed_ids;
							chromosome_event_t * tested_neighbour = search_return[xk3];

							if(tested_neighbour -> indel_at_junction > event_body -> indel_at_junction) continue;
							if(event_body -> indel_at_junction > tested_neighbour -> indel_at_junction && abs(tested_neighbour -> event_large_side - tested_neighbour -> event_small_side + tested_neighbour -> indel_at_junction - event_body -> event_large_side + event_body -> event_small_side - event_body -> indel_at_junction) <= 16)
								to_be_removed_ids[to_be_removed_number++] = event_body -> global_event_id;
							else if(tested_neighbour -> event_large_side >=  event_body -> event_large_side -indel_range + xk2 && tested_neighbour -> event_large_side <= event_body -> event_large_side + indel_range + xk2 && (event_body -> supporting_reads < tested_neighbour -> supporting_reads || (event_body -> supporting_reads ==  tested_neighbour -> supporting_reads  && xk2<0)))
								to_be_removed_ids[to_be_removed_number++] = event_body -> global_event_id;
						}
					}

				if((global_context-> config.do_fusion_detection || global_context-> config.do_long_del_detection) && event_body -> event_type == CHRO_EVENT_TYPE_FUSION)
				{
					for(xk2=-10 ; xk2 < 10 ; xk2++)
					{
						if(!xk2)continue;

						reallocate_to_be_removed_ids;
						unsigned int test_pos_small = event_body -> event_small_side + xk2;
						chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];
						int found_events = search_event(global_context,event_table, event_space, test_pos_small + delta_small , EVENT_SEARCH_BY_BOTH_SIDES, CHRO_EVENT_TYPE_JUNCTION|CHRO_EVENT_TYPE_FUSION, search_return); 

						for(xk3 = 0; xk3<found_events; xk3++)
						{
							chromosome_event_t * tested_neighbour = search_return[xk3];
							if(found_events && event_body -> supporting_reads < tested_neighbour -> supporting_reads)
							{
								to_be_removed_ids[to_be_removed_number++] = event_body -> global_event_id;
								break;
							}
						}
					}
				}
			}
		}
	}

	for(xk1=0; xk1<to_be_removed_number; xk1++) {
		int event_no = to_be_removed_ids[xk1];
		chromosome_event_t * deleted_event =  &event_space[event_no];
		
		for(xk2=0;xk2<2;xk2++){
			unsigned int pos = xk2?deleted_event->event_large_side:deleted_event->event_small_side;
			unsigned int * res = HashTableGet(event_table, NULL+pos);
			if(res){
				int current_size = res[0]&0x0fffffff;
				int wrt_ptr = 1;
				for(xk3 = 1 ; xk3< current_size+1 ; xk3++){
					if(!res[xk3])break;
					if(res[xk3] -1 == event_no) continue;
					if(wrt_ptr != xk3) res[wrt_ptr]=res[xk3];
					wrt_ptr++;
				}
				if(0 && wrt_ptr==1){
					HashTableRemove(event_table, NULL+pos);
					free(res);
				}else if(wrt_ptr < current_size +1)res[wrt_ptr]=0;
			}
		}
		if(deleted_event -> event_type == CHRO_EVENT_TYPE_INDEL && deleted_event -> inserted_bases)
			free(deleted_event -> inserted_bases);
		deleted_event -> event_type = CHRO_EVENT_TYPE_REMOVED;
	}

	free(to_be_removed_ids);
}



int localPointerCmp_forEventEntry(const void *pointer1, const void *pointer2)
{
	unsigned int v1 = (unsigned int)(pointer1-NULL);
	unsigned int v2 = (unsigned int)(pointer2-NULL);
	return v1!=v2;
}

srUInt_64 localPointerHashFunction_forEventEntry(const void *pointer)
{
	unsigned int v = (pointer-NULL);
	return (unsigned long) (v << 24) ^ (v>>8);
}


int init_indel_tables(global_context_t * context)
{
	indel_context_t * indel_context = (indel_context_t*)malloc(sizeof(indel_context_t));

	// event_entry_table is a table from an absolute chromosome location to a list of integers.
	// each integer is the index of event in the event array + 1.
	// unused values are set to 0.

	indel_context -> event_entry_table = NULL;

	context -> module_contexts[MODULE_INDEL_ID] = indel_context;

	indel_context -> total_events = 0;
	indel_context -> current_max_event_number = 0;
	//indel_context -> current_max_event_number = context->config.init_max_event_number;
	//indel_context -> event_space_dynamic = (chromosome_event_t *)malloc(sizeof(chromosome_event_t)*indel_context -> current_max_event_number);
	
	indel_context -> event_space_dynamic = NULL;

	if(context -> config.all_threads < 2) {
		indel_context -> event_entry_table = HashTableCreate(399997);

		indel_context -> event_entry_table -> appendix1=malloc(1024 * 1024 * 64);
		indel_context -> event_entry_table -> appendix2=malloc(1024 * 1024 * 64);
		memset(indel_context -> event_entry_table -> appendix1, 0, 1024 * 1024 * 64);
		memset(indel_context -> event_entry_table -> appendix2, 0, 1024 * 1024 * 64);

		HashTableSetKeyComparisonFunction(indel_context->event_entry_table, localPointerCmp_forEventEntry);
		HashTableSetHashFunction(indel_context->event_entry_table, localPointerHashFunction_forEventEntry);

		indel_context -> total_events = 0;
		indel_context -> current_max_event_number = context->config.init_max_event_number;
		indel_context -> event_space_dynamic = malloc(sizeof(chromosome_event_t)*indel_context -> current_max_event_number);
		if(!indel_context -> event_space_dynamic)
		{
			sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_FATAL, "Cannot allocate memory for threads. Please try to reduce the thread number.");
			return 1;
		}
	}
	if(context->config.is_third_iteration_running) {
		char * fns = malloc(200);
		fns[0]=0;
		exec_cmd("ulimit -n", fns, 200);
		int max_open_file = atoi(fns);
		//SUBREADprintf("SYS FILE LIMIT=%d\n", max_open_file);
		free(fns);
		max_open_file = max(100, max_open_file);
		max_open_file = min(3000, max_open_file);


		indel_context -> local_reassembly_pileup_files = HashTableCreate(100);
		indel_context -> local_reassembly_pileup_files -> appendix1 = NULL + max_open_file  * 2/ 3;
		HashTableSetDeallocationFunctions(indel_context -> local_reassembly_pileup_files, NULL, NULL);
		HashTableSetKeyComparisonFunction(indel_context -> local_reassembly_pileup_files, my_strcmp);
		HashTableSetHashFunction(indel_context -> local_reassembly_pileup_files ,HashTableStringHashFunction);
	}

	indel_context -> dynamic_align_table = malloc(sizeof(short*)*MAX_READ_LENGTH);
	indel_context -> dynamic_align_table_mask = malloc(sizeof(char *)*MAX_READ_LENGTH);

	int xk1;
	for(xk1=0;xk1<MAX_READ_LENGTH; xk1++) {
		indel_context -> dynamic_align_table[xk1] = malloc(sizeof(short)*MAX_READ_LENGTH);
		indel_context -> dynamic_align_table_mask[xk1] = malloc(sizeof(char)*MAX_READ_LENGTH);
	}

	for(xk1=0; xk1<EVENT_BODY_LOCK_BUCKETS; xk1++)
		subread_init_lock(indel_context -> event_body_locks+xk1);

	return 0;
}

int init_indel_thread_contexts(global_context_t * global_context, thread_context_t * thread_context, int task)
{
	indel_thread_context_t * indel_thread_context = (indel_thread_context_t*)malloc(sizeof(indel_thread_context_t));
	indel_context_t * indel_context = (indel_context_t *) global_context -> module_contexts[MODULE_INDEL_ID];
	
	if(task == STEP_VOTING) {
		indel_thread_context -> event_entry_table = HashTableCreate(399997 * ( global_context -> config.scRNA_input_mode?18:2 ));
		indel_thread_context -> event_entry_table -> appendix1=NULL;//indel_context -> event_entry_table-> appendix1;
		indel_thread_context -> event_entry_table -> appendix2=NULL;//indel_context -> event_entry_table-> appendix2;
		HashTableSetKeyComparisonFunction(indel_thread_context->event_entry_table, localPointerCmp_forEventEntry);
		HashTableSetHashFunction(indel_thread_context->event_entry_table, localPointerHashFunction_forEventEntry);

		indel_thread_context -> total_events = 0;
		indel_thread_context -> current_max_event_number = global_context->config.init_max_event_number;
		indel_thread_context -> event_space_dynamic = malloc(sizeof(chromosome_event_t)*indel_thread_context -> current_max_event_number);
		if(!indel_thread_context -> event_space_dynamic)
		{
			sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_FATAL, "Cannot allocate memory for threads. Please try to reduce the thread number.");
			return 1;
		}

		indel_thread_context -> dynamic_align_table = malloc(sizeof(short*)*MAX_READ_LENGTH);
		indel_thread_context -> dynamic_align_table_mask = malloc(sizeof(char *)*MAX_READ_LENGTH);

		int xk1;
		for(xk1=0;xk1<MAX_READ_LENGTH; xk1++)
		{
			indel_thread_context -> dynamic_align_table[xk1] = malloc(sizeof(short)*MAX_READ_LENGTH);
			indel_thread_context -> dynamic_align_table_mask[xk1] = malloc(sizeof(char)*MAX_READ_LENGTH);
		}
	}
	else if(task == STEP_ITERATION_TWO)
	{
		indel_thread_context -> event_entry_table = indel_context -> event_entry_table; 
		indel_thread_context -> total_events = indel_context -> total_events;
		indel_thread_context -> event_space_dynamic = indel_context -> event_space_dynamic;
		indel_thread_context -> final_reads_mismatches_array = (unsigned short*)malloc(sizeof(unsigned short)*indel_context -> total_events);
		indel_thread_context -> final_counted_reads_array = (unsigned short*)malloc(sizeof(unsigned short)*indel_context -> total_events);
		assert(indel_thread_context -> final_counted_reads_array);

		memset(indel_thread_context -> final_reads_mismatches_array , 0, sizeof(unsigned short)*indel_context -> total_events);
		memset(indel_thread_context -> final_counted_reads_array , 0, sizeof(unsigned short)*indel_context -> total_events);

		subread_init_lock(&thread_context -> output_lock);

	}

	thread_context -> module_thread_contexts[MODULE_INDEL_ID] = indel_thread_context;
	return 0;
}


void destory_event_entry_table(HashTable * old)
{
	int bucket;
	KeyValuePair * cursor;
	for(bucket=0; bucket<old -> numOfBuckets; bucket++)
	{
		cursor = old -> bucketArray[bucket];
		while(1)
		{
			if (!cursor) break;
			int * id_list = (int *) cursor ->value;
//#warning "=============== UNCOMMENT THE NEXT LINE IN THE RELEASE!!! THIS IS FOR REMOVING A SEGFAULT ERROR =================="
			free(id_list);

			cursor = cursor->next;
		}
	}
}

typedef struct {
	unsigned int scanning_positons;
	unsigned int thread_bodytable_number;
} scanning_events_record_t; 



#define _get_global_body(iii) ( indel_context ->  event_space_dynamic + records[iii].thread_bodytable_number)

int scanning_events_compare(void * arr, int l, int r){
	void ** arrr = (void **) arr;
	indel_context_t * indel_context = arrr[0];
	scanning_events_record_t * records = arrr[1];
	chromosome_event_t * body_l = _get_global_body(l);
	chromosome_event_t * body_r = _get_global_body(r);

	if(records[l].scanning_positons > records[r].scanning_positons)return 1;
	if(records[l].scanning_positons < records[r].scanning_positons)return -1;

	if((body_l -> is_donor_found_or_annotation & 64)!=0 && (body_r -> is_donor_found_or_annotation & 64)==0) return 1; // prefer known truth.
	if((body_l -> is_donor_found_or_annotation & 64)==0 && (body_r -> is_donor_found_or_annotation & 64)!=0) return -1;  
	if(body_l -> supporting_reads > body_r -> supporting_reads) return -1;
	if(body_l -> supporting_reads < body_r -> supporting_reads) return 1;
	if(abs(body_l -> indel_length) < abs(body_r -> indel_length)) return 1;
	if(abs(body_l -> indel_length) > abs(body_r -> indel_length)) return -1;
	if(body_l -> indel_length > body_r -> indel_length) return -1; // same length, but L is del and R is ins -- prefer del than ins
	if(body_l -> indel_length < body_r -> indel_length) return 1;
	
	if(body_l -> event_small_side > body_r -> event_small_side)return 1;
	if(body_l -> event_small_side < body_r -> event_small_side)return -1;
	if(body_l -> event_large_side > body_r -> event_large_side)return 1;
	return -1;
}

void null_scanning_events_merge(void * arr,  int start, int items, int items2){}
void scanning_events_merge(void * arr,  int start, int items, int items2){
	void ** arrr = (void **) arr;
	scanning_events_record_t * records = arrr[1];

	int read_1_ptr = start, read_2_ptr = start+items, write_ptr;
	scanning_events_record_t * merged_records = malloc(sizeof(scanning_events_record_t) * (items+items2));

	for(write_ptr=0; write_ptr<items+items2; write_ptr++){
		if((read_1_ptr >= start+items)||(read_2_ptr < start+items+items2 && scanning_events_compare(arr, read_1_ptr, read_2_ptr) > 0))
			memcpy(merged_records+write_ptr, records+(read_2_ptr++), sizeof(scanning_events_record_t));
		else
			memcpy(merged_records+write_ptr, records+(read_1_ptr++), sizeof(scanning_events_record_t));
	}	
	memcpy(records + start, merged_records, sizeof(scanning_events_record_t) * (items+items2));
	free(merged_records);
}

void scanning_events_exchange(void * arr, int l, int r){
	void ** arrr = (void **) arr;
	scanning_events_record_t * records = arrr[1];

	unsigned long tmpi;

	tmpi = records[l].scanning_positons;
	records[l].scanning_positons = records[r].scanning_positons;
	records[r].scanning_positons = tmpi;

	tmpi = records[l].thread_bodytable_number;
	records[l].thread_bodytable_number = records[r].thread_bodytable_number;
	records[r].thread_bodytable_number = tmpi;
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


int sort_junction_entry_table(global_context_t * global_context){
	indel_context_t * indel_context = (indel_context_t*)global_context -> module_contexts[MODULE_INDEL_ID];
	chromosome_event_t * event_space = indel_context -> event_space_dynamic;

	if(indel_context -> event_entry_table){
		if(indel_context -> event_entry_table->appendix1)
		{
			free(indel_context -> event_entry_table -> appendix1);
			free(indel_context -> event_entry_table -> appendix2);
		}

		// free the entry table: the pointers.
		destory_event_entry_table(indel_context -> event_entry_table);
		HashTableDestroy(indel_context -> event_entry_table);
	}

	indel_context -> event_entry_table = HashTableCreate(399997);
	HashTableSetKeyComparisonFunction(indel_context->event_entry_table, localPointerCmp_forEventEntry);
	HashTableSetHashFunction(indel_context->event_entry_table, localPointerHashFunction_forEventEntry);

	if(global_context -> config.use_bitmap_event_table) {
		indel_context -> event_entry_table -> appendix1=malloc(1024 * 1024 * 64);
		indel_context -> event_entry_table -> appendix2=malloc(1024 * 1024 * 64);
		memset(indel_context -> event_entry_table -> appendix1, 0, 1024 * 1024 * 64);
		memset(indel_context -> event_entry_table -> appendix2, 0, 1024 * 1024 * 64);
	} else {
		indel_context -> event_entry_table -> appendix1=NULL;
		indel_context -> event_entry_table -> appendix2=NULL;
	}

	int xx1, current_record_number=0, current_record_size=10000;

	scanning_events_record_t * records = malloc(sizeof(scanning_events_record_t)*current_record_size);

	for(xx1 = 0; xx1 < indel_context -> total_events; xx1++){
		chromosome_event_t * body = event_space + xx1;
		_test_record_size;
		_add_record;
	}

	void * sort_arr[2];
	sort_arr[0] = indel_context;
	sort_arr[1] = records;

	// many repeated elements 
	// do not use quick sort.
	merge_sort(sort_arr, current_record_number, scanning_events_compare, scanning_events_exchange, scanning_events_merge);
	unsigned int last_scannping_pos = records[0].scanning_positons;
	int merge_start = 0;
	HashTable * event_table = indel_context -> event_entry_table;

	for(xx1 =0; xx1 <= current_record_number; xx1++){
		scanning_events_record_t * this_record = NULL;
		if(xx1<current_record_number) this_record = records+xx1;

		if(xx1>0){
			if(xx1 == current_record_number || last_scannping_pos!= this_record -> scanning_positons){
				int merge_end = xx1, merge_i;
				if( merge_end-merge_start > MAX_EVENT_ENTRIES_PER_SITE )
					merge_end = merge_start + MAX_EVENT_ENTRIES_PER_SITE;
				unsigned int * id_list = malloc(sizeof(int) * (1+merge_end-merge_start));
				assert(id_list);
				id_list[0] = merge_end-merge_start;
				for(merge_i = merge_start; merge_i < merge_end; merge_i++){
					chromosome_event_t * body = _get_global_body(merge_i);
					id_list[merge_i - merge_start + 1] = records[merge_i].thread_bodytable_number + 1;

					mark_event_bitmap(event_table->appendix1, body -> event_small_side);
					mark_event_bitmap(event_table->appendix2, body -> event_large_side);
				}
				merge_start = xx1;

				//#warning "======= UNCOMMENT NEXT ======="
				HashTablePut(indel_context -> event_entry_table, NULL + last_scannping_pos, id_list);
			}
		}
		if(xx1 == current_record_number) break;
		last_scannping_pos = this_record -> scanning_positons;
	}

	free(records);
	return 0;
}

#define _get_event_body(iii) (records[iii].thread_no<0?( indel_context -> event_space_dynamic +  records[iii].thread_bodytable_number ):(((indel_thread_context_t*)(thread_contexts[records[iii].thread_no].module_thread_contexts[MODULE_INDEL_ID]))-> event_space_dynamic + records[iii].thread_bodytable_number ))

typedef struct {
	unsigned int thread_bodytable_number;
	short thread_no;
} concatinating_events_record_t; 

#define _test_conc_size if(conc_rec_items >= conc_rec_size - 1){\
			conc_rec_size *= 1.5;\
			records = realloc(records, conc_rec_size * sizeof(concatinating_events_record_t));\
			if(NULL == records) return -1;\
		}

int conc_sort_compare(void * arr, int l, int r){
	void ** arrr = (void **)arr;
	concatinating_events_record_t * records = arrr[0];
	indel_context_t * indel_context = arrr[1]; 
	thread_context_t * thread_contexts = arrr[2];

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


void conc_sort_merge(void * arr,  int start, int items, int items2){
	void ** arrr = (void **) arr;
	concatinating_events_record_t * records = arrr[0];

	int read_1_ptr = start, read_2_ptr = start+items, write_ptr;
	concatinating_events_record_t * merged_records = malloc(sizeof(concatinating_events_record_t) * (items+items2));

	for(write_ptr=0; write_ptr<items+items2; write_ptr++){
		if((read_1_ptr >= start+items)||(read_2_ptr < start+items+items2 && conc_sort_compare(arr, read_1_ptr, read_2_ptr) > 0))
			memcpy(merged_records+write_ptr, records+(read_2_ptr++), sizeof(concatinating_events_record_t));
		else
			memcpy(merged_records+write_ptr, records+(read_1_ptr++), sizeof(concatinating_events_record_t));
	}	
	memcpy(records + start, merged_records, sizeof(concatinating_events_record_t) * (items+items2));
	free(merged_records);
}



void conc_sort_exchange(void * arr, int l, int r){

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

// sort and merge events from all threads and the global event space.
int sort_global_event_table(global_context_t * global_context){
	return finalise_indel_and_junction_thread(global_context, NULL, STEP_VOTING);
}

int finalise_indel_and_junction_thread(global_context_t * global_context, thread_context_t * thread_contexts, int task)
{
	indel_context_t * indel_context = (indel_context_t*)global_context -> module_contexts[MODULE_INDEL_ID];
	if(task == STEP_VOTING) {
		concatinating_events_record_t * records;
		int conc_rec_size = 10000, conc_rec_items = 0;
		records = malloc(sizeof(concatinating_events_record_t) * conc_rec_size);

		int xk1, thn;
		if(thread_contexts)
			for(thn =0; thn < global_context->config.all_threads; thn++){
				thread_context_t * thread_context = thread_contexts+thn;
				indel_thread_context_t * indel_thread_context = (indel_thread_context_t*)thread_context -> module_thread_contexts[MODULE_INDEL_ID];

				for(xk1 = 0; xk1 < indel_thread_context -> total_events; xk1++){
					chromosome_event_t * old_body = indel_thread_context  ->  event_space_dynamic + xk1;
					if(old_body -> event_type == CHRO_EVENT_TYPE_REMOVED) continue;

					_test_conc_size;
					records[conc_rec_items].thread_no=thn;
					records[conc_rec_items++].thread_bodytable_number=xk1;
				}
			}

		for(xk1 = 0; xk1 < indel_context -> total_events; xk1++){
			chromosome_event_t * old_body = indel_context ->  event_space_dynamic + xk1;
			if(old_body -> event_type == CHRO_EVENT_TYPE_REMOVED) continue;
			_test_conc_size;
			records[conc_rec_items].thread_no=-1;
			records[conc_rec_items++].thread_bodytable_number=xk1;
		}

		void * sort_arr[3];
		sort_arr[0] = records;
		sort_arr[1] = indel_context;
		sort_arr[2] = thread_contexts;

		// many repeated elements -- do not use quick sort.
		merge_sort(sort_arr, conc_rec_items, conc_sort_compare, conc_sort_exchange, conc_sort_merge);

		chromosome_event_t * prev_env = NULL;
		int merge_target_size = 10000;
		int merge_target_items = 0;
		chromosome_event_t * merge_target = malloc(sizeof(chromosome_event_t) * merge_target_size);

		int merge_start = 0;
		for(xk1 = 0; xk1 <= conc_rec_items; xk1++){
			chromosome_event_t * this_event = (xk1 == conc_rec_items)?NULL:_get_event_body(xk1);
			if(xk1 > 0){
				int compret = 0;
				if(xk1 < conc_rec_items) compret = conc_sort_compare(sort_arr, xk1-1, xk1);
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
						else    merged_body -> connected_previous_event_distance = max(merged_body -> connected_previous_event_distance, old_body -> connected_previous_event_distance);

						if(old_body->inserted_bases && old_body -> event_type == CHRO_EVENT_TYPE_INDEL  && old_body -> indel_length<0){
							//printf("OCT27-FREEMEM-INS [%d] : thread-%d + %d %u-%p\n", xk_merge, records[xk_merge].thread_no, records[xk_merge].thread_bodytable_number, old_body->event_small_side, old_body->inserted_bases );
							if(merged_body -> inserted_bases && merged_body -> inserted_bases != old_body -> inserted_bases){
					//			SUBREADprintf("FREE PTR=%p\n", old_body -> inserted_bases);
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

		if(thread_contexts)
			for(thn =0; thn < global_context->config.all_threads; thn++){
				thread_context_t * thread_context = thread_contexts+thn;
				indel_thread_context_t * indel_thread_context = (indel_thread_context_t*)thread_context -> module_thread_contexts[MODULE_INDEL_ID];

				destory_event_entry_table(indel_thread_context -> event_entry_table);
				HashTableDestroy(indel_thread_context -> event_entry_table);

				free(indel_thread_context -> event_space_dynamic);

				for(xk1=0;xk1<MAX_READ_LENGTH; xk1++)
				{
					free(indel_thread_context -> dynamic_align_table[xk1]);
					free(indel_thread_context -> dynamic_align_table_mask[xk1]);
				}

				free(indel_thread_context->dynamic_align_table);
				free(indel_thread_context->dynamic_align_table_mask);
				free(indel_thread_context);
			}

		if(indel_context -> event_space_dynamic ) free(indel_context -> event_space_dynamic);
		indel_context -> event_space_dynamic = NULL;

		indel_context -> event_space_dynamic = merge_target; 
		indel_context -> current_max_event_number = merge_target_size;
		indel_context -> total_events = merge_target_items;
	} else if(task == STEP_ITERATION_TWO) {
		int xk1, thn;
		if(thread_contexts)
			for(thn =0; thn < global_context->config.all_threads; thn++){
				thread_context_t * thread_context = thread_contexts+thn;
				indel_thread_context_t * indel_thread_context = (indel_thread_context_t*)thread_context -> module_thread_contexts[MODULE_INDEL_ID];

				for(xk1 = 0; xk1 < indel_thread_context -> total_events; xk1++)
				{
					indel_context -> event_space_dynamic [xk1] . final_counted_reads +=indel_thread_context -> final_counted_reads_array[xk1];
					indel_context -> event_space_dynamic [xk1] . final_reads_mismatches +=indel_thread_context -> final_reads_mismatches_array[xk1];
				}
				free(indel_thread_context -> final_counted_reads_array);
				free(indel_thread_context -> final_reads_mismatches_array);
				free(indel_thread_context);
			}
	}
	return 0;
}


void add_annotation_to_junctions(void * key, void * val, HashTable * tab){
	global_context_t * global_context = tab->appendix1;
	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 
	char * gene_chr = key;
	int * pos_strands = val;
	int current_lagre_end = -1;
	char chroname[MAX_CHROMOSOME_NAME_LEN];

	int x1, ww=0, wi=0;
	for(x1 =0 ; ; x1++){
		if(ww){
			chroname[wi++]=gene_chr[x1];
			chroname[wi]=0;
		}
		if(0==gene_chr[x1]) break;
		if(gene_chr[x1]==':')
			ww=1;
	}

	//#warning ">>>>>>> COMMENt NEXT <<<<<<<<<<<<<<"
	//if(1) SUBREADprintf("SEE LOCATION: %s has %d\n",  gene_chr ,  pos_strands[0]+1 );

	for(x1 = 1; x1 < pos_strands[0]+1; x1+=3){
		if(pos_strands[x1]==0)break;
		//#warning ">>>>>>> COMMENt NEXT <<<<<<<<<<<<<<"
		//SUBREADprintf("SEE LOCATION ___ : %d >0 && %d < %d %d %d\n", current_lagre_end , current_lagre_end,  pos_strands[x1],  pos_strands[x1+1],  pos_strands[x1+2]);
		if(current_lagre_end > 0 && current_lagre_end < pos_strands[x1]){
			// add a junction between current_lagre_end and pos_strands[x1]
			if(indel_context -> total_events >= indel_context -> current_max_event_number -1){
				indel_context -> current_max_event_number*=1.2;
				indel_context -> event_space_dynamic = realloc(indel_context -> event_space_dynamic , sizeof( chromosome_event_t )*indel_context -> current_max_event_number);
			}
			chromosome_event_t * newbody = indel_context -> event_space_dynamic  + indel_context -> total_events ;
			memset(newbody, 0, sizeof(chromosome_event_t));
			newbody -> global_event_id = indel_context -> total_events++;
			newbody -> event_type = CHRO_EVENT_TYPE_JUNCTION;
			
			newbody -> event_small_side = linear_gene_position(&global_context->chromosome_table , chroname, current_lagre_end);
			newbody -> event_large_side = linear_gene_position(&global_context->chromosome_table , chroname, pos_strands[x1]);

			//#warning ">>>>>>> COMMENt NEXT <<<<<<<<<<<<<<"
			if(0){
				if ( 1 || newbody -> event_small_side == 625683723){
					SUBREADprintf("Loaded event :");
					debug_show_event(global_context, newbody);
				}
			}
			newbody -> is_negative_strand = pos_strands[x1+2];
			newbody -> is_donor_found_or_annotation = 127;
		//	printf("09NOV: add_annotation_to_junctions: %u -> %u '%s', %d,%d\n", newbody -> event_small_side ,  newbody -> event_large_side, chroname, current_lagre_end, pos_strands[x1]);
		}

		if(current_lagre_end < pos_strands[x1+1])
			current_lagre_end = pos_strands[x1+1];
	}
}


typedef struct {
	global_context_t * global_context;
	HashTable * feature_sorting_table;
} do_load_juncs_context_t;

int do_juncs_add_feature(char * gene_name, char * transcript_id, char * chro_name, unsigned int feature_start, unsigned int feature_end, int is_negative_strand, void * context){
	//#warning ">>>>>>>>>>>>>>>>> COMMENT NEXT <<<<<<<<<<<<<<<<<<<<"
	//return 0;
	//#warning ">>>>>>>>>>>>>>>>> COMMENT NEXT <<<<<<<<<<<<<<<<<<<<"
	//SUBREADprintf("INJ LOCS: %s : %u, %u\n", chro_name, feature_start, feature_end);
	do_load_juncs_context_t * do_load_juncs_context = context;
	HashTable * feature_sorting_table = do_load_juncs_context -> feature_sorting_table;

	global_context_t * global_context = do_load_juncs_context -> global_context;
	char tmp_chro_name[MAX_CHROMOSOME_NAME_LEN];
	if(global_context -> sam_chro_to_anno_chr_alias){
		char * sam_chro = get_sam_chro_name_from_alias(global_context -> sam_chro_to_anno_chr_alias, chro_name);
		if(sam_chro!=NULL) chro_name = sam_chro;
	}
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

	char sort_key[ MAX_CHROMOSOME_NAME_LEN * 3 + 2 ];
	sprintf(sort_key, "%s:%s", gene_name, chro_name);
	int * old_features = HashTableGet(feature_sorting_table, sort_key), x1;
	int written_space = -1, written_last_item = 0;
	if(old_features){
		int has_space = 0;
		for(x1= 1; x1 < old_features[0] + 1 ; x1+=3){
			if(old_features[x1]==0){
				has_space = 1;
				break;
			}
		}
		if(!has_space){
			int old_size = old_features[0];
			old_features[0] *= 1.5;
			old_features[0] -= old_features[0]%3;

			// I used malloc but not realloc here, because the old memory block will be freed by HashTablePutReplace.
			int * old_features2 = malloc( sizeof(int) * (old_features[0] +1));
			memcpy(old_features2, old_features, (old_size + 1) * sizeof(int));
			old_features2[old_size + 1] = 0;

			HashTablePutReplace(feature_sorting_table, sort_key, old_features2, 0);
			old_features = old_features2;
		}
		for(x1= 1; x1 < old_features[0] + 1 ; x1+=3){
			if(old_features[x1]> feature_start || old_features[x1]==0){
				written_space = x1;
				if(old_features[x1]> feature_start){
					int x2;
					for(x2 = x1+3; x2 <  old_features[0] + 1;x2+=3){
						if(old_features[x2]==0)break;
					}
					x2 += 2;
					assert(x2 < old_features[0] + 1);
					if(x2 + 1 <  old_features[0] + 1) old_features[x2+1]=0;

					for(; x2 >= x1 + 3; x2 --)
						old_features[x2] = old_features[x2-3];
				} else written_last_item = 1;
				break;
			}
		}
	}else{
		int init_space_sort = 15;
		old_features = malloc((init_space_sort + 1) * sizeof(int));
		old_features[0]=init_space_sort;
		char * mm_sort_key = malloc(strlen(sort_key)+1);
		strcpy(mm_sort_key, sort_key);
		HashTablePut(feature_sorting_table, mm_sort_key, old_features);
		written_space = 1;
		written_last_item = 1;
	}
	assert(written_space >0);
	old_features[written_space]=feature_start - 1;
	old_features[written_space+1]=feature_end - 1;
	old_features[written_space+2]=is_negative_strand;
	if(written_last_item && written_space + 3 < old_features[0]+1) old_features[written_space+3]=0;

	return 0;
}

int load_known_junctions(global_context_t * global_context){

	HashTable * feature_sorting_table = HashTableCreate(90239);
	HashTableSetKeyComparisonFunction(feature_sorting_table, my_strcmp);
	HashTableSetHashFunction(feature_sorting_table, HashTableStringHashFunction);
	HashTableSetDeallocationFunctions(feature_sorting_table, free, free);

	do_load_juncs_context_t do_load_juncs_context;
	memset(&do_load_juncs_context, 0, sizeof(do_load_juncs_context));

	do_load_juncs_context.global_context = global_context;
	do_load_juncs_context.feature_sorting_table = feature_sorting_table;

	int features = load_features_annotation(global_context->config.exon_annotation_file , global_context->config.exon_annotation_file_type, global_context->config.exon_annotation_gene_id_column, NULL, global_context->config.exon_annotation_feature_name_column, &do_load_juncs_context, do_juncs_add_feature);

	feature_sorting_table -> appendix1 = global_context;
	HashTableIteration(feature_sorting_table, add_annotation_to_junctions);
	HashTableDestroy(feature_sorting_table);

	if(features <0)return -1;
	else{
		return 0;
	}
}


int there_are_events_in_range(char * bitmap, unsigned int pos, int sec_len)
{
	int ret = 0;
	if(!bitmap) return 1;
	unsigned int offset = pos >> 3;
	unsigned int offset_end = (pos + sec_len) >> 3;
	unsigned int offset_byte = (offset >> 3)&0x3ffffff;
	unsigned int offset_byte_end =1+((offset_end >> 3)&0x3ffffff);
	for(; offset_byte < offset_byte_end ; offset_byte++)
	{
		if(bitmap[offset_byte]){
			ret = 1;
			break;
		}
	}
//	printf("TEST_JUMP: THERE ARE %d at %u\n" , ret , pos);
	return ret;
	
}
// bitmap has 64M bytes.
int check_event_bitmap(unsigned char * bitmap, unsigned int pos)
{
	unsigned int offset = pos >> 3;
	unsigned int offset_byte = (offset >> 3)&0x3ffffff;
	char keybyte = bitmap[offset_byte];
	if(!keybyte) return 0;

	int offset_bit = offset & 7;
	return (keybyte & (1<<offset_bit))!=0;
}

void mark_event_bitmap(unsigned char * bitmap, unsigned int pos)
{

	// bitmap is as large as 64MBytes
	// pos_all_bit = pos / 8
	// pos_byte = pos_all_bit / 8
	// pos_bit  = pos_all_bit % 8

	unsigned int offset = pos >> 3;
	unsigned int offset_byte = (offset >> 3)&0x3ffffff;
	int offset_bit = offset & 7;
	bitmap[offset_byte] |= (1<<offset_bit);
}


void put_new_event(HashTable * event_table, chromosome_event_t * new_event , int event_no)
{
	unsigned int sides[2];
	sides[0] = new_event->event_small_side;
	sides[1] = new_event->event_large_side;

	new_event -> global_event_id = event_no;
	int xk1, xk2;

	for(xk1=0; xk1<2; xk1++)
	{
		if(0 == sides[xk1])continue ;
		unsigned int * id_list = HashTableGet(event_table, NULL+sides[xk1]);
		if(!id_list)
		{
			id_list = malloc(sizeof(unsigned int)*(1+EVENT_ENTRIES_INIT_SIZE));
			id_list[0]=EVENT_ENTRIES_INIT_SIZE;
			id_list[1]=0;
			HashTablePut(event_table , NULL+sides[xk1] , id_list);
		}

		for(xk2=1;xk2< id_list[0] ; xk2++)
			if(id_list[xk2]==0) break;

		if(xk2 <  id_list[0])
			id_list[xk2] = event_no+1;
		if(xk2 < id_list[0]) id_list[xk2+1] = 0;
	}

	if(event_table->appendix1)
	{
		mark_event_bitmap(event_table->appendix1, sides[0]);
		mark_event_bitmap(event_table->appendix2, sides[1]);
	}
}
int search_event(global_context_t * global_context, HashTable * event_table, chromosome_event_t * event_space, unsigned int pos, int search_type, unsigned char event_type, chromosome_event_t ** return_buffer)
{
	int ret = 0;

	if(pos < 1 || pos > 0xffff0000)return 0;

	if(event_table->appendix1)
	{
		if(search_type == EVENT_SEARCH_BY_SMALL_SIDE)
			if(! check_event_bitmap(event_table -> appendix1, pos))return 0;

		if(search_type == EVENT_SEARCH_BY_LARGE_SIDE)
			if(! check_event_bitmap(event_table -> appendix2, pos))return 0;

		if(search_type == EVENT_SEARCH_BY_BOTH_SIDES)
			if(!(check_event_bitmap(event_table -> appendix1, pos) || check_event_bitmap(event_table -> appendix2, pos)))
				return 0;
	}

	unsigned int * res = HashTableGet(event_table, NULL+pos); 
	if(res)
	{
		int xk2;
		int current_size = res[0]&0x0fffffff;
		for(xk2=1; xk2< current_size+1 ; xk2++)
		{
			if(!res[xk2])break;
			//if(res[xk2] - 1>= ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> current_max_event_number ) { SUBREADprintf("FATAL ERROR: Event id out-of-boundary: %u > %u.\n", res[xk2],  ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> current_max_event_number ); continue;}
			chromosome_event_t * event_body = &event_space[res[xk2]-1]; 

			if((event_body -> event_type & event_type) == 0)continue;
			//SUBREADprintf("VB1:%u\n", res[xk2]);
			if(search_type == EVENT_SEARCH_BY_SMALL_SIDE && event_body -> event_small_side != pos)continue;
			if(search_type == EVENT_SEARCH_BY_LARGE_SIDE && event_body -> event_large_side != pos)continue;
			if(search_type == EVENT_SEARCH_BY_BOTH_SIDES && event_body -> event_small_side != pos && event_body -> event_large_side != pos)continue;

			return_buffer[ret++] = event_body;
		}
	}

	return ret;
}


void get_insertion_sequence(global_context_t * global_context, thread_context_t * thread_context , char * binary_bases , char * read_text , int insertions)
{
	int xk1;
	read_text[0]=0;

	for(xk1=0; xk1<insertions; xk1++)
	{
		int byte_no = xk1/4;
		int bit_no = 2*(xk1%4);
		
		read_text[xk1] = int2base((3&(binary_bases[byte_no]>> bit_no)));
	}
	read_text[insertions]=0;
}

void set_insertion_sequence(global_context_t * global_context, thread_context_t * thread_context , char ** binary_bases , char * read_text , int insertions)
{
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

chromosome_event_t * local_add_indel_event(global_context_t * global_context, thread_context_t * thread_context, HashTable * event_table, char * read_text, unsigned int left_edge, int indels, int score_supporting_read_added, int is_ambiguous, int mismatched_bases,int * old_event_id)
{
		chromosome_event_t * found = NULL;
		chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];

		chromosome_event_t * event_space = NULL;
		if(thread_context)
			event_space = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
		else
			event_space = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;


		int found_events = search_event(global_context, event_table, event_space, left_edge, EVENT_SEARCH_BY_SMALL_SIDE, CHRO_EVENT_TYPE_INDEL|CHRO_EVENT_TYPE_LONG_INDEL, search_return);

		//#warning ">>>>>>>>>>>>>>>>>> COMMENt THIS <<<<<<<<<<<<<<<<<<<<<<<"
		//printf("OCT27-STEPRR-TST LEFT=%u, INDEL=%d; FOUND=%d\n", left_edge , indels , found_events);


		if(found_events)
		{
			int kx1; 
			for(kx1 = 0; kx1 < found_events ; kx1++)
			{
				if(search_return[kx1] -> indel_length == indels)
				{
					found = search_return[kx1];	
					break;
				}
			}
		}

		if(found){
			if(old_event_id) (*old_event_id)=found->global_event_id;
			found -> supporting_reads += score_supporting_read_added;
			//found -> is_ambiguous = max(is_ambiguous , found -> is_ambiguous );
			return NULL;
		}
		else
		{
			int event_no;
			//event_no =(thread_context)?((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> total_events: ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) ->  total_events;

			if(thread_context)
				event_no = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> total_events ++;
			else
				event_no = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) ->  total_events ++;

			event_space = reallocate_event_space(global_context, thread_context, event_no);


			chromosome_event_t * new_event = event_space+event_no; 
			memset(new_event,0,sizeof(chromosome_event_t));
			if(indels <0)	// if it is an insertion:
				set_insertion_sequence(global_context, thread_context , &new_event -> inserted_bases , read_text, -indels);

			new_event -> event_small_side = left_edge;
			new_event -> event_large_side = left_edge + 1 + max(0, indels);
			//assert(new_event -> event_small_side && new_event -> event_large_side);
			new_event -> event_type = CHRO_EVENT_TYPE_INDEL;
			new_event -> indel_length = indels;
			new_event -> supporting_reads = 1;
			new_event -> event_quality = 1;//pow(0.5 , 3*mismatched_bases);
			//new_event -> is_ambiguous = is_ambiguous;
			
			//#warning ">>>>>>>>>>>>>>>>>> COMMENt THIS <<<<<<<<<<<<<<<<<<<<<<<"
			//printf("OCT27-STEPRR-NEW LEFT=%u RIGHT=%u, INDEL=%d\n", new_event -> event_small_side ,  new_event -> event_large_side, new_event -> indel_length );

			put_new_event(event_table, new_event , event_no);
			return new_event;
		}
	
}

float EXON_RECOVER_MATCHING_RATE = 0.9;
float EXON_INDEL_MATCHING_RATE_TAIL = 0.8;
float EXON_INDEL_MATCHING_RATE_HEAD = 0.8;




int core_extend_covered_region_15(global_context_t * global_context, gene_value_index_t *array_index, unsigned int read_start_pos, char * read, int read_len, int cover_start, int cover_end, int window_size, int req_match_5end , int req_match_3end, int indel_tolerance, int space_type, int tail_indel, short * head_indel_pos, int * head_indel_movement, short * tail_indel_pos, int * tail_indel_movement, int is_head_high_quality, char * qual_txt, int qual_format, float head_matching_rate, float tail_matching_rate){
	int is_head;
	
	if(0){
		char posout[100];
		absoffset_to_posstr(global_context, read_start_pos, posout);
		SUBREADprintf("RTXT=%s, MAPP=%s\n", read, posout);
	}

	for(is_head = 0; is_head < 2; is_head ++){
		int move_i;
		int best_match_n = -1, best_movement = 0;
		for(move_i = 0; move_i < 2* indel_tolerance-1; move_i ++){
			int indel_move = (move_i+1)/2 * (move_i %2?1:-1);
			int indel_movement = indel_move + (is_head?0:tail_indel);
			int this_match_n;

			if(is_head)
				this_match_n = match_chro(read, array_index, read_start_pos - indel_movement, window_size, 0, space_type );
			else
				this_match_n = match_chro(read + read_len - window_size, array_index, read_start_pos + read_len - window_size + indel_movement, window_size, 0, space_type );
			
			// indel_movement : negativ = insertion, positive = deletion before or after the covered region

			//SUBREADprintf("MOVE: HEAD=%d, MATCH=%d, MOV=%d\n", is_head, this_match_n, indel_movement);
			if(this_match_n > best_match_n){
				best_match_n = this_match_n;
				best_movement = indel_movement;
			}
		}

		int max_score=-1, best_splice_on_read=0;
		if(best_match_n >0 && best_movement!=0){
			int test_start = is_head?window_size:(cover_end);
			int test_end = is_head?cover_start - max(0, -best_movement):(read_len - window_size - max(0, -best_movement));
			int indel_first_piece_end_on_read, indel_second_piece_start_on_read;
			
			for(indel_first_piece_end_on_read = test_start ; indel_first_piece_end_on_read < test_end; indel_first_piece_end_on_read++){
				unsigned int indel_first_piece_end_on_chro = read_start_pos + indel_first_piece_end_on_read + (is_head? -best_movement :tail_indel);
				indel_second_piece_start_on_read = indel_first_piece_end_on_read - min(0, best_movement);
				unsigned int indel_second_piece_start_on_chro = indel_first_piece_end_on_chro + max(0,best_movement);
				int first_piece_match = match_chro(read + indel_first_piece_end_on_read - window_size , array_index , indel_first_piece_end_on_chro - window_size, window_size, 0, space_type);
				int second_piece_match = match_chro(read + indel_second_piece_start_on_read , array_index , indel_second_piece_start_on_chro, window_size, 0, space_type);
				int score = first_piece_match + second_piece_match;
				if(score > max_score){
					max_score = score;
					best_splice_on_read = indel_first_piece_end_on_read;
				}
				if(score == 2*window_size)break;
				//SUBREADprintf("FIND_SPLICE_POINT: IS_HEAD=%d, FIRST_PIECE_END=%d-%d, CHRO=%u-%u MATCHED=%d,%d\n", is_head, indel_first_piece_end_on_read, indel_second_piece_start_on_read, indel_first_piece_end_on_chro, indel_second_piece_start_on_chro, first_piece_match, second_piece_match);
			}
		}

		if(max_score >=  2*window_size - 1){
			if(is_head){
				(*head_indel_pos) = best_splice_on_read;
				(*head_indel_movement) = best_movement;
			} else {
				(*tail_indel_pos) = best_splice_on_read;
				(*tail_indel_movement) = best_movement;
			}
		}
	}

	return 0;

}



int core_extend_covered_region_13(gene_value_index_t *array_index, unsigned int read_start_pos, char * read, int read_len, int cover_start, int cover_end, int window_size, int req_match_5end , int req_match_3end, int indel_tolerance, int space_type, int tail_indel, short * head_indel_pos, int * head_indel_movement, short * tail_indel_pos, int * tail_indel_movement, int is_head_high_quality, char * qual_txt, int qual_format, float head_matching_rate, float tail_matching_rate)
{
	int ret = 0;
	*head_indel_pos = -1;
	*tail_indel_pos = -1;
	if (cover_start >= window_size && head_matching_rate < 1.0001) 
	{
		int head_test_len =  cover_start;
		int roughly_mapped = match_chro(read, array_index, read_start_pos, head_test_len , 0, space_type);
		if (roughly_mapped >= head_test_len  * EXON_RECOVER_MATCHING_RATE - 0.0001)
		{
			ret |= 1;
		}
		else
		{
			int window_end_pos = cover_start + window_size-1;
			int not_too_bad = 1;
			int right_match_number=0;

			while (1)
			{
				int matched_bases_in_window = match_chro_wronglen(read+ window_end_pos - window_size, array_index, read_start_pos + window_end_pos - window_size, window_size, space_type, NULL, &right_match_number);
				int best_indel_movement = -1000;
				int best_indel_pos = -1;

				if (matched_bases_in_window >= req_match_5end) // this window is not bad enough so that an indel is considered
					window_end_pos--;
				else
				{
					roughly_mapped = match_chro(read, array_index, read_start_pos, window_end_pos - right_match_number , 0, space_type);
					if (roughly_mapped < (int)(0.5+ (window_end_pos - right_match_number)  * EXON_RECOVER_MATCHING_RATE))
					{
						// the quality of this window is too low (at most 1 base is matched). I must consider if it is an indel.
						int indel_movement_i ;

						int best_matched_bases_after_indel = -1;

						not_too_bad = 0;
						for (indel_movement_i = 0; indel_movement_i < 2* indel_tolerance-1 ; indel_movement_i ++)
						{
							int indel_movement = (indel_movement_i+1)/2 * (indel_movement_i %2?1:-1);
							int test_length = window_end_pos /*- 1*/ - max(0, indel_movement) -  right_match_number;

							if (test_length < window_size) continue;
							//if (test_length <= 1+ abs(indel_movement)/4) continue;
							//test_length = min(10, test_length);
							if (abs(indel_movement) > indel_tolerance) continue;

							int test_start =0;// window_end_pos - max(0, indel_movement) -  right_match_number - test_length;

							float matched_bases_after_indel = match_chro_support(read +test_start, array_index, read_start_pos + indel_movement +test_start, test_length,0, space_type, qual_txt, qual_format);
							SUBREADprintf("HEAD : MATCHED_AFTER_INDEL = %f ; MVMT=%d ; WINDOW_END=%d\n", matched_bases_after_indel, indel_movement, window_end_pos);

							float test_rate = head_matching_rate;
							
							if(test_length < 3) test_rate = 1;

							if(best_matched_bases_after_indel <matched_bases_after_indel  && matched_bases_after_indel >= (int)(0.5+test_length * test_rate))
							{
								not_too_bad = 1;

								best_matched_bases_after_indel = matched_bases_after_indel;
								best_indel_movement = indel_movement;
								best_indel_pos = window_end_pos - right_match_number;


								
								*head_indel_pos = best_indel_pos;
								*head_indel_movement = best_indel_movement;
							}
						}
						if(best_indel_pos<0) *head_indel_pos =  window_end_pos - right_match_number;
					//	break;
					}
					window_end_pos--;
				}
				if (window_end_pos - window_size <= 0) break;
			}
			if(not_too_bad)
				ret |=1;
			//else
				//*head_indel_pos = -1;
		}
	}
	else ret |= 1;



	if (cover_end <= read_len - window_size && tail_matching_rate < 1.0001) 
	{
		int tail_test_len = read_len - cover_end;
		int roughly_mapped = match_chro(read + cover_end, array_index, read_start_pos + cover_end + tail_indel, tail_test_len , 0, space_type);
		if (roughly_mapped >= tail_test_len  * EXON_RECOVER_MATCHING_RATE - 0.0001)
		{
			ret |= 2;
		}
		else
		{
			int window_start_pos = cover_end - window_size +1;
			int not_too_bad = 1;

			while (1)
			{
				int left_match_number = 0;
				int matched_bases_in_window = match_chro_wronglen(read+ window_start_pos, array_index, read_start_pos + window_start_pos + tail_indel, window_size, space_type, &left_match_number, NULL);
				int best_indel_movement = -1000;
				int best_indel_pos = -1;

				//printf("MOVE 4-bp WINDOW : %d ; matched = %d >=? %d\n", window_start_pos, matched_bases_in_window, req_match_3end);

				if (matched_bases_in_window >= req_match_3end) // this window is not bad enough so that an indel is considered
					window_start_pos++;
				else
				{
					roughly_mapped = match_chro(read+window_start_pos + left_match_number, array_index, read_start_pos + window_start_pos + tail_indel + left_match_number, read_len - window_start_pos - left_match_number , 0, space_type);

					if (roughly_mapped < (int)(0.5+(read_len - window_start_pos -left_match_number )  * EXON_RECOVER_MATCHING_RATE))
					{
						// the quality of this window is too low (at most 1 base is matched). I must consider if it is an indel.
						int indel_movement_i ;
						int best_matched_bases_after_indel = -1;
						not_too_bad = 0;
						for (indel_movement_i = 0; indel_movement_i < 2* indel_tolerance; indel_movement_i ++)
						{
							
							int indel_adjustment = (indel_movement_i+1)/2 * (indel_movement_i %2?1:-1);
							int indel_movement = tail_indel + indel_adjustment;
							int test_length = read_len - window_start_pos  - left_match_number + min(0, indel_adjustment);
							//test_length = min(10, test_length);


							if (test_length < window_size) continue;
							//if (test_length <= 1 + abs(indel_movement)/4) continue;
							if (abs(indel_movement) > indel_tolerance) continue;

							char * qual_txt_moved = qual_txt;
							if(qual_txt[0]) qual_txt_moved+=window_start_pos - min(0, indel_movement) + left_match_number;

							float matched_bases_after_indel = match_chro_support(read + window_start_pos - min(0, indel_movement) + left_match_number, array_index, read_start_pos + window_start_pos  + max(0,indel_movement) +left_match_number , test_length,0, space_type, qual_txt_moved , qual_format);
							SUBREADprintf("TAIL : MATCHED_AFTER_INDEL = %f ; MVMT=%d ; WINDOW_END=%d\n", matched_bases_after_indel, indel_movement, window_start_pos - min(0, indel_movement) + left_match_number);

							float test_rate = tail_matching_rate;
							
							if(test_length < 3) test_rate = 1;

							//printf("%d <? %d && %d >=? %d\n", best_matched_bases_after_indel, matched_bases_after_indel, matched_bases_after_indel, (int)(0.5+test_length * test_rate));
							if(best_matched_bases_after_indel <matched_bases_after_indel  && matched_bases_after_indel >= (int)(0.5+test_length * test_rate))
							{
								not_too_bad = 1;
								best_matched_bases_after_indel = matched_bases_after_indel;
								best_indel_movement = indel_movement;
								best_indel_pos = window_start_pos + left_match_number ;//-1;
								*tail_indel_movement = best_indel_movement ;
							}
						}
					
						if(best_indel_pos<0)
							*tail_indel_pos =  window_start_pos + left_match_number ;
						else
							*tail_indel_pos = best_indel_pos;
					}
					window_start_pos++;
				}
				if (window_start_pos + window_size >= read_len) break;
			}
			if(not_too_bad)
				ret |=2;
			//else
			//	*tail_indel_pos =-1;
		}
	}
	else ret |= 2;

	return ret;
}






int core_dynamic_align(global_context_t * global_context, thread_context_t * thread_context, char * read, int read_len, unsigned int begin_position, char * movement_buffer, int expected_offset, char * read_name);

int find_new_indels(global_context_t * global_context, thread_context_t * thread_context, int pair_number, char * read_name, char * read_text, char * qual_text, int read_len, int is_second_read, int best_read_id)
{
	int i,j;

	int last_correct_subread = 0;
	int last_indel = 0;

	gene_vote_number_t * indel_recorder;
	unsigned int voting_position;

	HashTable * event_table = NULL;

	if(thread_context)
		event_table = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
	else
		event_table = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_entry_table; 


	mapping_result_t *current_result = _global_retrieve_alignment_ptr(global_context, pair_number, is_second_read, best_read_id);

	if(global_context->config.do_big_margin_filtering_for_reads)
		if(is_ambiguous_voting(global_context, pair_number, is_second_read , current_result->selected_votes, current_result->confident_coverage_start, current_result->confident_coverage_end, read_len, (current_result->result_flags & CORE_IS_NEGATIVE_STRAND)?1:0))return 0;



	indel_recorder = current_result -> selected_indel_record; 
	voting_position = current_result -> selected_position; 


	if(!indel_recorder[0])
		return 0;

	gene_value_index_t * current_value_index = thread_context?thread_context->current_value_index:global_context->current_value_index; 

	if(0){ // test "shifting indels" -- not needed when "NEW" go_Q used
		int last_offset = 0;
		for(i=0; indel_recorder[i]  && (i<MAX_INDEL_SECTIONS); i+=3){
			last_offset = indel_recorder[i+2];
		}

		if(0==last_offset) return 0;
	}

	for(i=0; indel_recorder[i]  && (i<MAX_INDEL_SECTIONS); i+=3)
	{
		int indels = indel_recorder[i+2] - last_indel;
			// -1 : 1 insert; 1: 1 delete
		int next_correct_subread = indel_recorder[i] -1;
		//int next_correct_subread_last = indel_recorder[i+1] -1;

		if(indels)
		{
			int last_correct_base = find_subread_end(read_len, global_context->config.total_subreads , last_correct_subread) - 9;
			int first_correct_base = find_subread_end(read_len, global_context->config.total_subreads , next_correct_subread) - 16 + 9;

			last_correct_base = max(0, last_correct_base);
			last_correct_base = min(read_len-1, last_correct_base);
			first_correct_base = min(first_correct_base, read_len-1);
			first_correct_base = max(0, first_correct_base);
			first_correct_base = max(first_correct_base, last_correct_base);

			if(first_correct_base < last_correct_base || first_correct_base > last_correct_base + 3000)
				SUBREADprintf("WRONG ORDER: F=%u, L=%d\n", first_correct_base , last_correct_base);

			//printf("OCT27-STEPRR %s POS=%u, I=%d; INDELS=%d; COVER=%d -- %d (vote_no : %d - %d)\n", read_name, current_result->selected_position, i, indels, last_correct_base, first_correct_base, last_correct_subread, next_correct_subread);

			if(global_context -> config.use_dynamic_programming_indel || read_len > EXON_LONG_READ_LENGTH)
			{
				char movement_buffer[MAX_READ_LENGTH * 10 / 7];
				//chromosome_event_t * last_event = NULL;
				int last_event_id = -1;

				first_correct_base = min(first_correct_base+10, read_len);
				int x1, dyna_steps;

				dyna_steps = core_dynamic_align(global_context, thread_context, read_text + last_correct_base, first_correct_base - last_correct_base, voting_position + last_correct_base + last_indel, movement_buffer, indels, read_name);

				movement_buffer[dyna_steps]=0;

				if(0 && FIXLENstrcmp("simulated.2467286", read_name) == 0){
					SUBREADprintf("DYNAMIC: indels=%d, %s\n", indels, movement_buffer);
				}

				//#warning ">>>>>>> COMMENT NEXT BLOCK IN RELEASE <<<<<<<<"
				if(0) {
					char outstr[1000];
					sprintf(outstr, "OCT27-STEPDD-IR %s     %d  %d~%d ", read_name, dyna_steps, last_correct_base, first_correct_base);

					for(x1=0; x1<dyna_steps;x1++)
					{
						int mc, mv=movement_buffer[x1];
						if(mv==0)mc='=';
						else if(mv==1)mc='D';
						else if(mv==2)mc='I';
						else mc='X';
						sprintf(outstr+strlen(outstr),"%c",mc);
					}
					puts(outstr);
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

				if(total_mismatch<=2 || (global_context->config.maximise_sensitivity_indel && total_mismatch <= 2 ))
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

								int ambiguous_count=0;
								if(0){
									int ambiguous_i, best_matched_bases = match_chro(read_text + cursor_on_read - 6, current_value_index, indel_left_boundary - 6, 6, 0, global_context->config.space_type)  +
												 match_chro(read_text + cursor_on_read - min(current_indel_len,0), current_value_index, indel_left_boundary + max(0, current_indel_len), 6, 0, global_context->config.space_type);
									for(ambiguous_i=-5; ambiguous_i<=5; ambiguous_i++)
									{
										int left_match = match_chro(read_text + cursor_on_read - 6, current_value_index, indel_left_boundary - 6, 6+ambiguous_i, 0, global_context->config.space_type);
										int right_match = match_chro(read_text + cursor_on_read + ambiguous_i - min(current_indel_len,0), current_value_index, indel_left_boundary + ambiguous_i + max(0, current_indel_len), 6-ambiguous_i, 0,global_context->config.space_type);
										if(left_match+right_match == best_matched_bases) ambiguous_count ++;
									}
								}

								if(abs(current_indel_len)<=global_context -> config.max_indel_length)
								{
									int  old_event_id = -1;

									if(0)if(indel_left_boundary >= 46481340 + 1210 - 200 && indel_left_boundary  <= 46481340 + 1210 + 200 && abs(current_indel_len )==4){
										SUBREADprintf("ADD AN INDEL: %s : %u ; len = %d\n", read_name, indel_left_boundary, current_indel_len);
									}
									chromosome_event_t * new_event = local_add_indel_event(global_context, thread_context, event_table, read_text + cursor_on_read + min(0,current_indel_len), indel_left_boundary - 1, current_indel_len, 1, ambiguous_count, 0, &old_event_id);
									mark_gapped_read(current_result);

									chromosome_event_t * event_space = NULL;
									if(thread_context)
										event_space = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
									else
										event_space = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;

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
			else
			{
				if(first_correct_base - last_correct_base  < -min(0,indels)) continue;

				int best_score=-99999, is_ambiguous_indel= 0;
				unsigned int best_pos = 0;
				int cutoff = first_correct_base-last_correct_base  + min(0, indels) - global_context -> config.flanking_subread_indel_mismatch;
				int mid_j = (first_correct_base-last_correct_base )/2, min_j=99990;
				int mismatched_bases = 0;

				if(first_correct_base > last_correct_base + 4)
					for(j=0; j<(first_correct_base-last_correct_base + min(0, indels)); j++)
					{
						// j is the first UNWANTED base after the first matched part.
						int score , s2, s1;

						s1 = match_chro(read_text + last_correct_base, current_value_index, voting_position + last_correct_base + last_indel, j, 0, global_context->config.space_type);
						if(s1 >= j - global_context -> config.flanking_subread_indel_mismatch)
						{
							s2 = match_chro(read_text + last_correct_base + j - min(0, indels) , current_value_index, voting_position + last_correct_base + j + max(0,indels) + last_indel, first_correct_base-last_correct_base - j + min(0, indels), 0, global_context->config.space_type);

							score = s2+s1;
							if(score >= cutoff)
							{
								if(score == best_score)
									is_ambiguous_indel ++;

								if(score>best_score || (abs(mid_j-j)>min_j && best_score == score))
								{
									mismatched_bases = first_correct_base-last_correct_base  + min(0, indels) - score;
									best_score = score;
									// best_pos here is the absolute offset of the FIRST UNWANTED BASE.
									best_pos = voting_position  + last_correct_base + j + last_indel;
									min_j = abs(mid_j-j);
									if(score>best_score)
										is_ambiguous_indel = 0;
								}
							}
						}
					}


				
				/*
				SUBREADprintf("\n%s\nRANGE:%d - %d ; BSCORE:%d @  %d   ; AMB=%d ; INDELS=%d ; MM=%d\n", read_name, first_correct_base, last_correct_base, best_score, last_indel + last_correct_base + best_j, is_ambiguous_indel, indels, mismatched_bases);

				SUBREADprintf("%s\n", read_text);
				for(j=0;j<  last_indel + last_correct_base + best_j; j++)
					putc(' ' , stderr);
				putc('\\' , stderr);
				for(j=0;j<  indels; j++)
					putc(gvindex_get( current_value_index, best_pos + j)  , stderr);
				putc('\n', stderr);

				*/

				if(best_score >0) {
					//#warning ">>>>>>>>>>>>>>>>>> COMMEN THIS <<<<<<<<<<<<<<<<<<<<<<"
					//printf("OCT27-STEPAD-%s POS=%u INDEL=%d\n", read_name, best_pos , indels);
					local_add_indel_event(global_context, thread_context, event_table, read_text + last_correct_base + last_indel, best_pos - 1, indels, 1, is_ambiguous_indel, mismatched_bases, NULL);
					mark_gapped_read(current_result);
				}
			}
		}

		last_correct_subread = indel_recorder[i+1]-1;
		last_indel = indel_recorder[i+2];
	}



	if(global_context -> config.extending_search_indels && global_context -> config.max_indel_length>0)
	{
		// USING EXTENDING TO FIND INDELS
		// THIS IS UNSTABLE!
		int s_head = match_chro(read_text, current_value_index, voting_position, 8, 0, global_context->config.space_type);
		int s_tail = match_chro(read_text+read_len - 8, current_value_index, voting_position + read_len + current_result -> indels_in_confident_coverage - 8, 8, 0, global_context->config.space_type);

		//SUBREADprintf("EXT_START %d, %d\n", s_head, s_tail);

		//#warning "============ REMOVE THE FIRST '1' FROM THE NEXT LINE! =========="
		if(s_head<=6 || s_tail<=6)
		{
			int is_head_high_quality = (current_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?0:1;
			int k;
			int cover_start = find_subread_end(read_len, global_context -> config.total_subreads, current_result -> selected_indel_record[0]-1)- 15;

			for (k=0; current_result -> selected_indel_record[k]; k+=3);
			int cover_end = find_subread_end(read_len, global_context -> config.total_subreads, current_result -> selected_indel_record[k-2]-1) + max(0,-current_result -> selected_indel_record[k-1]);

			int head_qual = is_head_high_quality?0x7fffffff:0, tail_qual = is_head_high_quality?0:0x7fffffff;


			if(qual_text && qual_text[0])
			{
				int j;
				head_qual = 0;
				tail_qual = 0;
				for (j = 0; j<cover_start; j++)
				{
					if(FASTQ_PHRED64 == global_context -> config.phred_score_format)
					{
						head_qual += (1000000-get_base_error_prob64i(qual_text[j]));
					}
					else
					{
						head_qual += (1000000-get_base_error_prob33i(qual_text[j]));
					}
				}
				for (j = 0; j<read_len - cover_end; j++)
				{
					if(FASTQ_PHRED64 == global_context -> config.phred_score_format)
					{
						tail_qual += (1000000-get_base_error_prob64i(qual_text[read_len - j -1]));
					}
					else
					{
						tail_qual += (1000000-get_base_error_prob33i(qual_text[read_len - j -1]));
					}
				}
			}

			float exon_tail_matching_rate=0.9, exon_head_matching_rate=0.9;

			int head_must_correct = 4;
			if(cover_start > 0)
			{
				if (head_qual / cover_start < 850000)
				{
					exon_head_matching_rate = 0.75;
					head_must_correct =2;
				}
				else if (head_qual / cover_start < 950000 )
				{
					exon_head_matching_rate = 0.85;
					head_must_correct =3;
				}
				else exon_head_matching_rate = 0.92;
			}else	exon_head_matching_rate = 9999;


			int tail_must_correct = 4;
			if(read_len-cover_end > 0)
			{
				if (tail_qual / (read_len-cover_end) < 850000)
				{
					exon_tail_matching_rate = 0.75;
					tail_must_correct =2;
				}
				else if (tail_qual / (read_len-cover_end) < 950000)
				{
					exon_tail_matching_rate = 0.85;
					tail_must_correct =3;
				}
				else exon_tail_matching_rate = 0.92;
			}else	exon_tail_matching_rate = 9999;

			short head_indel_pos=-1 , tail_indel_pos=-1;
			int head_indel_movement=0, tail_indel_movement=0;

			if(0)SUBREADprintf("HQ=%.4f; TQ=%.4f; HM=%d; TM=%d; COVG = %d ~ %d\n", exon_head_matching_rate, exon_tail_matching_rate, head_must_correct, tail_must_correct, cover_start, cover_end);
			/*
			if(exon_head_matching_rate<1)
				exon_head_matching_rate = 0.97;
			if(exon_tail_matching_rate<1)
				exon_tail_matching_rate = 0.97;
 			*/

			core_extend_covered_region_15(global_context, current_value_index, current_result->selected_position, read_text, read_len, cover_start, cover_end, 7, head_must_correct,tail_must_correct, global_context -> config.max_indel_length + 1, global_context -> config.space_type, current_result -> selected_indel_record[k-1], & head_indel_pos, & head_indel_movement, & tail_indel_pos, & tail_indel_movement, is_head_high_quality, qual_text,  global_context -> config.phred_score_format, exon_head_matching_rate, exon_tail_matching_rate);

			//head_indel_movement = -head_indel_movement;

			if(0)SUBREADprintf("HMV=%d; TMV=%d; TPOS=%d\n", head_indel_movement, tail_indel_movement, tail_indel_pos);
				
			if(head_indel_movement && s_head < 7)
			{
				unsigned int head_indel_left_edge = head_indel_pos + current_result->selected_position - 1;
				//head_indel_left_edge -= max(0, head_indel_movement);
				if(abs(head_indel_movement)<=global_context -> config.max_indel_length)
				{
					local_add_indel_event(global_context, thread_context, event_table, read_text + head_indel_pos, head_indel_left_edge, head_indel_movement, 1, 1, 0,NULL);
					mark_gapped_read(current_result);
				}
			}
			if(tail_indel_movement && s_tail < 7)
			{
				unsigned int tail_indel_left_edge = tail_indel_pos + current_result->selected_position + current_result -> indels_in_confident_coverage - 1;
				
				if(abs(tail_indel_movement)<=global_context -> config.max_indel_length)
				{
					local_add_indel_event(global_context, thread_context, event_table, read_text + tail_indel_pos, tail_indel_left_edge, tail_indel_movement, 1, 1, 0,NULL);
					mark_gapped_read(current_result);
				}

			}
		}

	}
	return 0;
}

void print_indel_table(global_context_t * global_context){
}

int write_indel_final_results(global_context_t * global_context)
{

	int xk1, disk_is_full = 0;
	char * inserted_bases = NULL;
	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 
	char * fn2, *ref_bases, *alt_bases;
	FILE * ofp = NULL;

	fn2 = malloc(MAX_FILE_NAME_LENGTH+30);
	sprintf(fn2, "%s.indel.vcf", global_context->config.output_prefix);

	ofp = f_subr_open(fn2, "wb");

	//if(!ofp)
	//	printf("HOW??? %s\n", fn2);
	inserted_bases = malloc(MAX_INSERTION_LENGTH + 2); 
	ref_bases = malloc(1000);
	alt_bases = malloc(1000);

	fputs("##fileformat=VCFv4.0\n##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n##INFO=<ID=SR,Number=1,Type=String,Description=\"Number of supporting reads for variants\">\n", ofp);
	fputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", ofp);

	for(xk1 = 0; xk1 <  indel_context -> total_events ; xk1++)
	{ 
		char * chro_name;
		int chro_pos; 

		chromosome_event_t * event_body = indel_context -> event_space_dynamic +xk1;


//		#warning " ================= REMOVE '- 1' from the next LINE!!! ========================="
		if((event_body -> event_type != CHRO_EVENT_TYPE_INDEL && event_body->event_type != CHRO_EVENT_TYPE_LONG_INDEL  && event_body -> event_type != CHRO_EVENT_TYPE_POTENTIAL_INDEL)|| (event_body -> final_counted_reads < 1  && event_body -> event_type == CHRO_EVENT_TYPE_INDEL) )
			continue;

		//assert((-event_body -> indel_length) < MAX_INSERTION_LENGTH);

		if(event_body -> indel_length<0)
		{
			get_insertion_sequence(global_context, NULL , event_body -> inserted_bases , inserted_bases , -event_body -> indel_length);
			free(event_body -> inserted_bases);
		}

		//assert(event_body -> event_small_side < event_body -> event_large_side );
		locate_gene_position( event_body -> event_small_side , &global_context -> chromosome_table, &chro_name, &chro_pos);

		// VCF format: chr21  1000001  .   AACC  AAGGGGGCC  29  .  INDEL;DP=20
		if(event_body -> event_type == CHRO_EVENT_TYPE_INDEL || event_body -> event_type == CHRO_EVENT_TYPE_LONG_INDEL)
		{
			ref_bases[0]=0;
			alt_bases[0]=0;

			gene_value_index_t * current_index = find_current_value_index(global_context , event_body -> event_small_side-1 , max(0, event_body -> indel_length) + 2);
			if(current_index)
			{
				int rlen = max(0, event_body -> indel_length) + 2;
				assert(rlen<900);
				gvindex_get_string(ref_bases, current_index, event_body -> event_small_side-1 , rlen, 0);
				ref_bases[rlen] = 0;

				if(event_body -> indel_length > 0)  // deletion
				{
					alt_bases[0] = ref_bases[0];
					alt_bases[1] = ref_bases[event_body -> indel_length+1];
					alt_bases[2] = 0;
				}
				else // insertion 
				{
					alt_bases[0] = ref_bases[0];
					strcpy(alt_bases+1,inserted_bases);
					strcat(alt_bases,ref_bases+1);
				}
			}

			if(CHRO_EVENT_TYPE_INDEL == event_body -> event_type)
			{
				if(event_body -> final_counted_reads<30)
					event_body -> event_quality = pow(0.5,(30-event_body -> final_counted_reads));
				else	event_body -> event_quality = 1;

			}
			int write_len = fprintf(ofp, "%s\t%u\t.\t%s\t%s\t%d\t.\tINDEL;DP=%d;SR=%d\n", chro_name, chro_pos, ref_bases, alt_bases, (int)(max(1, 250 + 10*log(event_body -> event_quality)/log(10))), event_body -> final_counted_reads + event_body -> anti_supporting_reads, event_body -> final_counted_reads);
			if(write_len < 10) disk_is_full = 1;
		}

		global_context->all_indels++;
	}

	fclose(ofp);
	free(ref_bases);
	free(alt_bases);
	free(inserted_bases);
	if(disk_is_full){
		unlink(fn2);	
		SUBREADprintf("ERROR: disk is full. Unable to write into the indel list.\n");
	}
	free(fn2);

	return 0;
}

int destroy_indel_module(global_context_t * global_context)
{
	int xk1;
	indel_context_t * indel_context = global_context -> module_contexts[MODULE_INDEL_ID];
	if(indel_context -> event_entry_table->appendix1)
	{
		free(indel_context -> event_entry_table -> appendix1);
		free(indel_context -> event_entry_table -> appendix2);
	}

	// free the entry table: the pointers.
	destory_event_entry_table(indel_context -> event_entry_table);
	HashTableDestroy(indel_context -> event_entry_table);

	free(indel_context -> event_space_dynamic);

	for(xk1=0;xk1<MAX_READ_LENGTH; xk1++)
	{
		free(indel_context -> dynamic_align_table[xk1]);
		free(indel_context -> dynamic_align_table_mask[xk1]);
	}

	for(xk1=0; xk1<EVENT_BODY_LOCK_BUCKETS; xk1++)
		subread_destroy_lock(indel_context -> event_body_locks+xk1);

	free(indel_context -> dynamic_align_table);
	free(indel_context -> dynamic_align_table_mask);
	return 0;
}



int trim_read(global_context_t *global_context , thread_context_t *thread_context, char * read_text, char * qual_text, int read_len, int * delta_position)
{
	int xk1, bad_bases = 0;
	if(qual_text[0]==0) return read_len;


	int tail_trim_pos, last_ok = read_len /2;

	bad_bases=0;
	for(xk1=read_len/2 ; xk1< read_len ; xk1++)
	{
		if(qual_text[xk1] <= '#'+2)
			bad_bases++;
		else
			last_ok = xk1;
		if(bad_bases>=2)
			break;
	}
	tail_trim_pos = last_ok;

	int head_trim_pos;

	bad_bases=0;
	last_ok = read_len /2;
	for(xk1=read_len/2 ; xk1>=0 ; xk1--)
	{
		if(qual_text[xk1] <= '#'+2)
			bad_bases++;
		else
			last_ok = xk1;
		if(bad_bases>=2)
			break;
	}
	if(bad_bases>=2)
		head_trim_pos = last_ok;
	else	head_trim_pos = 0;

	if(3*(tail_trim_pos - head_trim_pos )<read_len)return -1;

	for(xk1=0; xk1<tail_trim_pos - head_trim_pos; xk1++)
	{
		read_text[xk1] = read_text[xk1+head_trim_pos];
		qual_text[xk1] = qual_text[xk1+head_trim_pos];
	}

	(* delta_position) = head_trim_pos;
	read_text[xk1]=0;
	qual_text[xk1]=0;

	return tail_trim_pos - head_trim_pos;

}

#define READ_HEAD_TAIL_LEN 22 
int write_local_reassembly(global_context_t *global_context, HashTable *pileup_fp_table, unsigned int anchor_pos, char * read_name , char * read_text ,char * qual_text , int read_len, int is_anchor_certain)
{

	char * chro_name;
	int chro_offset;
	int delta_pos = 0;
	int new_read_len = trim_read(global_context, NULL, read_text, qual_text, read_len, &delta_pos);


	anchor_pos += delta_pos;
	if(new_read_len * 3<read_len * 2)
		return 0;

	read_len = new_read_len;

	gene_value_index_t * base_index = &global_context->all_value_indexes[0] ;

	if(is_anchor_certain)
	{
		int head_match = match_chro(read_text, base_index , anchor_pos, READ_HEAD_TAIL_LEN, 0, global_context -> config.space_type);
		int tail_match = match_chro(read_text + read_len - READ_HEAD_TAIL_LEN, base_index , anchor_pos + read_len - READ_HEAD_TAIL_LEN, READ_HEAD_TAIL_LEN, 0, global_context -> config.space_type);
		if(head_match > READ_HEAD_TAIL_LEN-3 && tail_match> READ_HEAD_TAIL_LEN-3)
			return 0;
	}




	if(0 == locate_gene_position_max(anchor_pos, &global_context -> chromosome_table, &chro_name, &chro_offset, NULL, NULL, read_len))
	{
		char temp_file_name[MAX_FILE_NAME_LENGTH+40];
		int close_now = 0;

		sprintf(temp_file_name,"%s@%s-%04u.bin", global_context -> config.temp_file_prefix, chro_name , chro_offset / BASE_BLOCK_LENGTH );

		FILE * pileup_fp = get_temp_file_pointer(temp_file_name, pileup_fp_table, &close_now); 
		//assert(read_len == strlen(read_text) && read_len > 90);
		write_read_block_file(pileup_fp , 0, read_name, 0, chro_name , chro_offset, NULL, 0, read_text , qual_text, read_len , 1 , is_anchor_certain , anchor_pos, read_len, 0);
		if(close_now) fclose(pileup_fp);
	}


	return 0;
}

int build_local_reassembly(global_context_t *global_context , thread_context_t *thread_context, int pair_number , char * read_name , char * read_text ,char * qual_text , int read_len , int mate_read_len , int is_second_read, int best_read_id, int use_mate_pos, mapping_result_t *current_result, mapping_result_t *mate_result)
{
	unsigned int read_anchor_position;
	if(!read_text) return 0;
	int is_position_certain = 1;

	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 
	//mapping_result_t *current_result = _global_retrieve_alignment_ptr(global_context, pair_number, is_second_read, best_read_id);
	//mapping_result_t * mate_result = _global_retrieve_alignment_ptr(global_context, pair_number, !is_second_read, best_read_id);



	// No complex cigar strings are allowed.
	if(use_mate_pos)
	{
		read_anchor_position = mate_result -> selected_position;

		if(is_second_read + ((mate_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0) == 1)
			// I'm the second read so I'm after the first one.
			read_anchor_position += global_context -> config.expected_pair_distance;
		else
			// I'm the first read so I'm ahead the second one.
			read_anchor_position -= global_context -> config.expected_pair_distance;

		is_position_certain = 0;

	}
	else
	{
		read_anchor_position = current_result -> selected_position;
	}

	//gene_alue_index_t * value_index = &global_context->all_value_indexes[0] ;

	if(qual_text[0] &&global_context -> config.phred_score_format == FASTQ_PHRED64)
		fastq_64_to_33(qual_text);

	//int window_no = 2*(read_anchor_position / REASSEMBLY_WINDOW_LENGTH);


	write_local_reassembly(global_context, indel_context -> local_reassembly_pileup_files, read_anchor_position, read_name , read_text , qual_text , read_len, is_position_certain);

	return 0;
}

int insert_pileup_read(global_context_t * global_context , reassembly_block_context_t * block_context , unsigned int first_base_pos ,char * read_text , char * qual_text , int read_len, int is_read_head_aligned)
{
	int window_i,xk1;

	for(window_i= -1; window_i<1; window_i++)
	{
		int de_brujin_grp_id = 2*(first_base_pos / REASSEMBLY_WINDOW_LENGTH) + window_i;
		int read_pos;
		if(de_brujin_grp_id<0)continue;

		HashTable * current_table = block_context -> de_bruijn_graphs [de_brujin_grp_id]; 
		for(read_pos = 0; read_pos <= read_len - global_context -> config.k_mer_length; read_pos +=1)
		{
			unsigned long long int node_key = 0;
			for(xk1=0; xk1<global_context ->config.k_mer_length; xk1++)
			{
				node_key = (node_key << 2) | (base2int(read_text[xk1 + read_pos]));
			}

			
			unsigned long long int next_key_prefix = 0x8000000000000000ull | (node_key >> 2);

			if(!read_pos && is_read_head_aligned)
			{
				int start_offset = (first_base_pos - de_brujin_grp_id * REASSEMBLY_WINDOW_LENGTH/2);
				if(start_offset < block_context -> start_offsets [de_brujin_grp_id])
				{
					block_context -> start_offsets [de_brujin_grp_id] = start_offset;
					block_context -> start_keys[1+10*de_brujin_grp_id]=0;
					block_context -> start_keys[10*de_brujin_grp_id]=next_key_prefix;
				}
				else if(start_offset == block_context -> start_offsets [de_brujin_grp_id])
				{
					int start_key_id;
					for(start_key_id=0; start_key_id<10&&block_context -> start_keys[de_brujin_grp_id * 10 + start_key_id]; start_key_id++);
					if(start_key_id<10) block_context -> start_keys[start_key_id + 10 * de_brujin_grp_id] = next_key_prefix;
				}
			}

			int next_node_type = node_key&0x03;
			unsigned int next_jump_array = HashTableGet(current_table ,NULL+next_key_prefix) - NULL;
			if(next_jump_array) next_jump_array-=1;

			unsigned int current_read_count = 0xff&(next_jump_array >> (next_node_type*8));
			if(current_read_count<254) current_read_count ++;

			next_jump_array = 1+((next_jump_array & ~(0xff << (next_node_type*8))) | ((current_read_count) << (next_node_type*8)));

			HashTablePut(current_table, NULL + next_key_prefix , NULL + next_jump_array);
		}
	}
	return 0;
}


int finalise_db_graph(global_context_t * global_context, reassembly_block_context_t * block_context, int de_brujin_graph_id, HashTable * de_brujin_graph, unsigned long long int start_key, short start_offset)
{
	
	unsigned long long int current_key = start_key;
	char inserted_sequence[MAX_INSERTION_LENGTH+REASSEMBLY_WINDOW_LENGTH]; 
	unsigned long long int used_keys[MAX_INSERTION_LENGTH+REASSEMBLY_WINDOW_LENGTH];
	int inserted_length = 0;
	int is_looped = 0;

	int xk1;

	for(xk1=0; xk1<global_context -> config.k_mer_length-1; xk1++)
	{
		int intbase = (start_key >> (2*(global_context -> config.k_mer_length-xk1-2))) & 3;
		inserted_sequence[xk1] = "AGCT"[intbase];
	}
	inserted_length = global_context -> config.k_mer_length-1;

	while(1)
	{
		unsigned int next_moves;
		int max_move=-1;
		next_moves = HashTableGet(de_brujin_graph , NULL + current_key) - NULL;
		if(!next_moves) break;
		else
		{
			int test_next;
			int all_reads = 0;
			int ignored_moves=0;
			int max_move_reads=0;
			unsigned long long int previous_key = current_key<<2;

			next_moves-=1;
			for(test_next=0; test_next<4; test_next++)
				all_reads += 0xff&(next_moves>>(8*test_next));

			while(1)
			{
				max_move = -1;
				max_move_reads = -1;
				for(test_next=0; test_next<4; test_next++)
				{
					if(ignored_moves & (1<<test_next)) continue;
					int test_move_reads=0xff&(next_moves>>(8*test_next));
					if(test_move_reads > max_move_reads)
					{
						max_move_reads = test_move_reads;
						max_move = test_next;
					}
				}

				if(max_move < 0)
					break;

				inserted_sequence[inserted_length] = "AGCT"[max_move];
				int xk2, is_used_key = 0;
				current_key = (previous_key | max_move);

				for(xk2=0; xk2<inserted_length- global_context -> config.k_mer_length + 1; xk2++)
				{
					if(current_key == used_keys[xk2])
					{
						is_used_key = 1;
						is_looped = 1;
						ignored_moves |= (1<<max_move);
						break;
					}
				}
				if(!is_used_key || is_looped) break;
			}
			used_keys[inserted_length - global_context -> config.k_mer_length + 1] = current_key;
			current_key &= (0x3fffffffffffffffull >> (2*(32 - global_context -> config.k_mer_length)));
			current_key |= 0x8000000000000000ull;
			inserted_length++;
		}
		if(inserted_length > MAX_INSERTION_LENGTH + REASSEMBLY_WINDOW_LENGTH - 1 || max_move<0|| is_looped)
		{
			is_looped = 1;
			inserted_sequence[inserted_length]=0;

			break;
		}
	}

	if(!is_looped)
	{
		inserted_sequence[inserted_length]=0;


		assert(global_context -> index_block_number == 1);	
		gene_value_index_t * base_index = &global_context->all_value_indexes[0] ;
		unsigned int window_start = block_context -> block_start_linear_pos + de_brujin_graph_id * REASSEMBLY_WINDOW_LENGTH / 2; 
		int best_matched_bases = -1, best_end_pos = -1, best_start_pos = start_offset;

		#define END_POS_CONFIRM_LENGTH 24 

		for(xk1 = inserted_length + start_offset - END_POS_CONFIRM_LENGTH ; xk1>= start_offset+15 ; xk1--)	//find the end of insertion
		{
			if(inserted_length - xk1 - END_POS_CONFIRM_LENGTH + best_start_pos >= MAX_INSERTION_LENGTH) break;
			int matched_bases = match_chro(inserted_sequence + inserted_length - END_POS_CONFIRM_LENGTH, base_index , window_start + xk1, END_POS_CONFIRM_LENGTH, 0, global_context -> config.space_type );
			if(matched_bases>best_matched_bases)
			{
				best_matched_bases = matched_bases;
				best_end_pos = xk1 + END_POS_CONFIRM_LENGTH;
			}
		}


		int new_base_length = inserted_length - best_end_pos + best_start_pos;
		if(new_base_length>0 && best_matched_bases >= END_POS_CONFIRM_LENGTH - 1)
		{
			best_matched_bases = -1;
			int best_insertion_point = -1;
			for(xk1 = 15; xk1<best_end_pos - best_start_pos; xk1++) 
			{
				int insertion_matched = 
					match_chro(inserted_sequence + xk1 - 15, base_index , window_start + best_start_pos + xk1 - 15, 15, 0, global_context -> config.space_type) + 
					match_chro(inserted_sequence + xk1 + new_base_length , base_index , window_start + best_start_pos + xk1 , 15, 0 , global_context -> config.space_type );
				if(insertion_matched>best_matched_bases && insertion_matched>=29)
				{
					best_insertion_point = xk1;
					best_matched_bases = insertion_matched;
				}
			}
			if(best_insertion_point>=0)
			{
				unsigned int best_pos = window_start + best_start_pos + best_insertion_point;
				int indels = -new_base_length;

				//inserted_sequence[best_insertion_point + new_base_length] = 0;


				indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 
				HashTable * event_table = indel_context -> event_entry_table;
				chromosome_event_t * event_space = indel_context -> event_space_dynamic;

				chromosome_event_t * found = NULL;
				chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];
				int found_events = search_event(global_context, event_table, event_space, best_pos - 1, EVENT_SEARCH_BY_SMALL_SIDE, CHRO_EVENT_TYPE_INDEL, search_return);

				if(found_events)
				{
					int kx1; 
					for(kx1 = 0; kx1 < found_events ; kx1++)
					{
						if(search_return[kx1] -> indel_length == indels)
						{
							found = search_return[kx1];	
							break;
						}
					}
				}

				if(found){
					found -> supporting_reads ++;
				}
				else 
				{
					int event_no;
					event_no = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) ->  total_events ++;

					event_space = reallocate_event_space(global_context, NULL, event_no);

					chromosome_event_t * new_event = event_space+event_no; 
					memset(new_event,0,sizeof(chromosome_event_t));

					set_insertion_sequence(global_context, NULL , &new_event -> inserted_bases , inserted_sequence+best_insertion_point, -indels);

					new_event -> event_small_side = best_pos - 1;
					new_event -> event_large_side = best_pos + max(0, indels);
					assert(new_event -> event_small_side && new_event -> event_large_side);
					new_event -> event_type = CHRO_EVENT_TYPE_INDEL;
					new_event -> indel_length = indels;
					new_event -> supporting_reads = 1;
					
					put_new_event(event_table, new_event , event_no);

				}
	
				fprintf(global_context -> long_insertion_FASTA_fp, ">BEGIN%u_%dM%dI%dM\n", window_start + best_start_pos, best_insertion_point , -indels , inserted_length + indels - best_insertion_point);
				int write_cursor;
				for(write_cursor = 0 ; write_cursor < inserted_length ; write_cursor += 70)
				{
					fwrite(inserted_sequence + write_cursor , min(70, inserted_length - write_cursor), 1, global_context -> long_insertion_FASTA_fp);
					fputs("\n",global_context -> long_insertion_FASTA_fp);
				}
			}
		}

	}
	return 0;
}


#define  finalise_pileup_file_by_voting  finalise_pileup_file

#define REASSEMBLY_SUBREAD_GAP 1

int search_window_once(global_context_t * global_context, reassembly_by_voting_block_context_t * block_context, int window_no, int start_read_no, int search_direction, int anti_loop, float allele_quality);
 // search_direction = 1 : from 5' end to 3' end
 // search_direction = -1 : from 3'end to 5' end
void search_window(global_context_t * global_context, reassembly_by_voting_block_context_t * block_context, int window_no, int start_read_no, int search_direction)
{

	block_context->rebuilt_window[0] = 0;
	block_context->rebuilt_window_size = 0;

	memset(block_context->final_alleles, 0, sizeof(struct reassmebly_window_allele ) * global_context -> config.reassembly_window_alleles);

	block_context -> search_cost = 0;
	block_context -> used_read_number = 0;
	block_context -> total_matched_bases = 0;
	block_context -> max_matched_bases = 0;

	/*int ret = */search_window_once(global_context, block_context, window_no, start_read_no, search_direction, 0, 1);

}


// this function returns if the read is not rectifiable (i.e. likely to be misalignment)
int rectify_read_text(global_context_t * global_context, reassembly_by_voting_block_context_t * block_context, int window_no, gehash_t * index, char * next_read_txt, int really_make_change)
{
//	return 1;
	int read_offset,i,j;
	int read_len = strlen(next_read_txt);

	//int is_test_target = (strcmp("AATGGAACAACCGCCTCAAATGGAATGGAATGGCCTCCAAAGGAAGGGAGTGACATTTAATGAATGCATTCTAATGGAATGGCATGGAATGGACTGGAA", next_read_txt)==0);
	int is_test_target = 0;

	init_gene_vote(block_context -> vote_list_rectify);
	for(read_offset = 0; read_offset < read_len - global_context -> config.reassembly_subread_length+1; read_offset += 1)
	{
		unsigned int subread_int = 0, xk2;
		for(xk2=0; xk2< global_context -> config.reassembly_subread_length; xk2++)
			subread_int = (subread_int << 2) | base2int(next_read_txt[read_offset + xk2]);
		gehash_go_q(index, subread_int , read_offset, read_len, 0, block_context -> vote_list_rectify, 22, read_offset,  0, 0x7fffffff);
	}

	memset(block_context -> read_rectify_space, 0, MAX_READ_LENGTH * 4 * sizeof(short));

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< block_context -> vote_list_rectify -> items[i]; j++)
		{
			int xk1;

			int neighbour_read_no =  (block_context -> vote_list_rectify -> pos[i][j] + MAX_READ_LENGTH+1)/10000;
			int neighbour_read_mapped_start, this_read_mapping_start, overlapped_length;
			unsigned long long int neighbour_read_key = (window_no * 1llu +1)<<32 ;
			neighbour_read_key |= (neighbour_read_no& 0xfffff)<<12; 
			char * neighbour_read_text = HashTableGet(block_context->read_sequence_table,NULL+neighbour_read_key);
			assert(neighbour_read_text);
			int neighbour_read_len = strlen(neighbour_read_text);
	
			// neighbour_read_mapped_start is the location on the neigibour read of the first base of the current read.
			neighbour_read_mapped_start = - neighbour_read_no*10000 + block_context -> vote_list_rectify -> pos[i][j];// + block_context -> vote_list -> coverage_start[i][j];

			if(neighbour_read_mapped_start>=0)
			{
				this_read_mapping_start = 0;
				overlapped_length = min(neighbour_read_len - neighbour_read_mapped_start, read_len);
			}
			else
			{
				// this_read_mapping_start is the location on this read of the first base of the neighbour read.
				this_read_mapping_start = -neighbour_read_mapped_start;
				neighbour_read_mapped_start = 0;
				overlapped_length = min(read_len - this_read_mapping_start, neighbour_read_len);
				if(is_test_target)
				{
					for(xk1 = 0; xk1 < 3+this_read_mapping_start; xk1++) fputc(' ',stdout );
					printf("%s\n", neighbour_read_text);
				}
			}
			

			for(xk1 = 0; xk1 < overlapped_length; xk1++)
			{
				int this_read_pos = this_read_mapping_start + xk1;
				int nb_read_pos = neighbour_read_mapped_start + xk1;
				block_context -> read_rectify_space[this_read_pos * 4 +base2int(neighbour_read_text[nb_read_pos])] ++;
			}
		}

	if(is_test_target)
		printf("   %s\n", next_read_txt);
	char ori[2000];
	strcpy(ori, next_read_txt);
	is_bad_read(ori, read_len, "RECT");
	int rectified = 0;


	if(is_test_target)
		printf("R0=%s\n", next_read_txt);
	char changed[200];
//	memset(changed, ' ', 200);
	for(i=0; i<read_len; i++)
	{
		int max_support = block_context -> read_rectify_space[i*4+base2int(next_read_txt[i])];
		int base_sel = -1;
		for(j=0; j<4; j++)
		{
			if(j!=base2int(next_read_txt[i]) && block_context -> read_rectify_space[i*4+j] > max_support)
			{
				base_sel = int2base(j);
				max_support = block_context -> read_rectify_space[i*4+j];
			}
		}

		if(base_sel>0) 
		{
			if(really_make_change)
			{
				next_read_txt[i]=base_sel;
			//	changed[i] = '~';
			}
			rectified ++;
		}
	}

	changed[i]=0;

	if(is_test_target)
		printf("R1=%s\n   %s\n\n", next_read_txt, changed);
	is_bad_read(next_read_txt, read_len, "REC2");
	return rectified<=2;
}

int search_window_once(global_context_t * global_context, reassembly_by_voting_block_context_t * block_context, int window_no, int start_read_no, int search_direction, int anti_loop, float allele_quality)
{

	block_context -> search_cost++;
	if(block_context -> search_cost > 1000) return 0;
	

	int is_start_certain = 1, is_digged_in = 0;
	unsigned long long int starter_read_key = (window_no * 1llu +1)<<32 ;
	starter_read_key |= (start_read_no& 0xfffff)<<12; 

	char * start_read_txt = HashTableGet(block_context->read_sequence_table,NULL+starter_read_key);
	assert(start_read_txt);


	int start_read_pos = HashTableGet(block_context->read_position_table, NULL+starter_read_key) - NULL ;
	if(start_read_pos == 0) is_start_certain = 0;
	else start_read_pos --;
	

	int read_len = strlen(start_read_txt), read_offset;
	init_gene_vote(block_context -> vote_list);

	int minimum_overlap_length = read_len  /4;
	//if(!is_start_certain) minimum_overlap_length *= 1.5;
	//minimum_overlap_length = max(1,16-global_context ->config.reassembly_subread_length);

	
	for(read_offset = 0; read_offset < read_len - global_context -> config.reassembly_subread_length+1; read_offset += 1)
	{
		unsigned int subread_int = 0, xk2;
		for(xk2=0; xk2< global_context -> config.reassembly_subread_length; xk2++)
			subread_int = (subread_int << 2) | base2int(start_read_txt[read_offset + xk2]);


		if(global_context -> config.reassembly_tolerable_voting)
			gehash_go_q_tolerable(&block_context->voting_indexes[window_no], subread_int , read_offset, read_len, 0, block_context -> vote_list, 1, 22, 24, 0, read_offset, global_context -> config.reassembly_tolerable_voting, global_context -> config.reassembly_subread_length,  0, 0x7fffffff);
		else
			gehash_go_q(&block_context->voting_indexes[window_no], subread_int , read_offset, read_len, 0, block_context -> vote_list,0,read_offset,  0, 0x7fffffff);
	}

	int WINDOW_SEARCH_TREE_WIDTH = 1;
	int i,j;

	int max_read_scoring [8];
	unsigned int max_read_ids [8];
	int max_read_nonoverlap [8];
	float max_allele_quality [8];

	memset(max_read_scoring, 0, sizeof(int)* 8);
	memset(max_allele_quality, 0, sizeof(float)* 8);
	memset(max_read_ids, 0, sizeof(int)* 8);
	memset(max_read_nonoverlap, 0, sizeof(int)* 8);

	//intif(is_target_window_X(window_no))
	//if(block_context -> window_start_pos > 1124357402 && block_context -> window_start_pos < 1124359402)
	//	print_votes(block_context -> vote_list, global_context -> config.index_prefix);

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< block_context -> vote_list -> items[i]; j++)
		{
			//if(is_digged_in)break;
			//if((block_context -> vote_list -> coverage_end[i][j] - block_context -> vote_list -> coverage_start[i][j]) * 6 > read_len ) 
//			if(block_context -> window_start_pos > 1124357402 && block_context -> window_start_pos < 1124359402)
//				SUBREADprintf("LLLX PRINT_VOTE @%u : V=%d\n", block_context -> vote_list -> pos[i][j],  block_context -> vote_list -> votes[i][j] );

			if( block_context -> vote_list -> votes[i][j] > 4) 
			{
				int next_read_no = (block_context -> vote_list -> pos[i][j] + MAX_READ_LENGTH +1 )/10000;
				if(next_read_no == start_read_no) continue;
				int next_read_nonoverlap = next_read_no * 10000 - block_context -> vote_list -> pos[i][j];
				float next_read_allele_quality = allele_quality;
				int anchor_to_index = 0;

				if(abs(next_read_nonoverlap) < 2)
					anchor_to_index = 0;	// PARALLEL or VALUELESS
				else if(next_read_nonoverlap < 0)
					anchor_to_index = -1;	// INDEX -> ANCHOR
				else	anchor_to_index = 1;    // ANCHOR -> INDEX


				unsigned long long int test_next_read_key = (window_no * 1llu +1)<<32 ;
				test_next_read_key |=  ( next_read_no& 0xfffff)<<12; 
				int is_next_certain = 1;
				int next_read_pos = HashTableGet(block_context->read_position_table, NULL+test_next_read_key) - NULL;
				if(next_read_pos==0)
					is_next_certain =0;
				else next_read_pos--;
				
				/*
				if(!is_start_certain || !is_next_certain)
				{
					if(anchor_to_index == 1 && (block_context -> vote_list -> coverage_start[i][j] - next_read_nonoverlap) > 5)
						anchor_to_index = 0;	// PARALLEL or VALUELESS

					if(anchor_to_index == -1 && block_context -> vote_list -> coverage_start[i][j] > 5)
						anchor_to_index = 0;	// PARALLEL or VALUELESS
				}*/

				char * next_read_txt = HashTableGet(block_context->read_sequence_table,NULL+test_next_read_key);
				assert(next_read_txt);
				int next_read_len = strlen(next_read_txt);
				int two_read_overlapped = next_read_len - next_read_nonoverlap;
				if(next_read_nonoverlap<0)
					two_read_overlapped = read_len + next_read_nonoverlap;
				if(two_read_overlapped< minimum_overlap_length)
					anchor_to_index = 0;
	
				if(search_direction>0)
				{
					if( next_read_nonoverlap + strlen(next_read_txt) <= read_len)
						anchor_to_index = 0;
				}

				if(anchor_to_index * search_direction > 0)
				{
	
					long long int pos_delta = next_read_pos;
					pos_delta -= start_read_pos;

					int xk3, next_read_scoring = 0, next_read_matching = 0;

					int order_reward = (is_next_certain && is_start_certain && search_direction * pos_delta >0)? 25000:0;
					if(is_start_certain * is_next_certain == 0) order_reward = 10000;
					int wrong_bases=0;

					if(anchor_to_index > 0)
					{
						if(next_read_len+ next_read_nonoverlap< read_len)continue;// no new bases are added by the new read.

						for(xk3 = next_read_nonoverlap ; xk3 < read_len ; xk3++)
							if(start_read_txt[xk3] == next_read_txt[xk3 - next_read_nonoverlap]) next_read_matching++;
						wrong_bases = read_len - next_read_nonoverlap - next_read_matching;
						if(wrong_bases>=2)
						{
							next_read_scoring = -1;
						}
						else
						{
							block_context -> total_matched_bases += next_read_matching;

							next_read_scoring = two_read_overlapped * 10000 + order_reward
   						    	+(unsigned int)(HashTableGet(block_context -> read_quality_table, NULL+test_next_read_key)-NULL)/10 - wrong_bases * 20000;
						}
					}
					else
					{
						if(0< next_read_nonoverlap) continue; // no new bases are added by the new read.  
						int should_match = min(next_read_len, read_len - next_read_nonoverlap);
						if(next_read_len >  read_len - next_read_nonoverlap) continue; // the head of the next read is untested: bad!
						for(xk3 = - next_read_nonoverlap ; xk3 <  should_match ; xk3++)
							if(next_read_txt[xk3] == start_read_txt[xk3 + next_read_nonoverlap]) next_read_matching++;

						wrong_bases = should_match + next_read_nonoverlap - next_read_matching;

						if(wrong_bases>=2)
						{
							next_read_scoring = -1;
						}
						else
						{
							block_context -> total_matched_bases += next_read_matching;

							next_read_scoring = two_read_overlapped * 10000 +order_reward
						 	+ (unsigned int)(HashTableGet(block_context -> read_quality_table, NULL+test_next_read_key)-NULL)/10 - wrong_bases * 20000;
						}

					}


					if(wrong_bases)
						next_read_allele_quality *= 0.5;

					int is_rectifiable = rectify_read_text(global_context, block_context, window_no, &block_context->voting_indexes[window_no], next_read_txt , 0);


/*
					if(block_context -> window_start_pos > 1124357402 && block_context -> window_start_pos < 1124359402)
					{
						SUBREADprintf("LLLX: NREAD #%d HAS QUAL %d ; is_rectifiable=%d ; WB=%d; ANCHOR_TO_I=%d\n", next_read_no, next_read_scoring, is_rectifiable, wrong_bases, anchor_to_index);
						SUBREADprintf("LLLX: NREAD OVL=%d, RLEN=%d, NRLEN=%d, NRMATCH=%d\n", next_read_nonoverlap, read_len, next_read_len, next_read_matching);
						SUBREADprintf("LLLX: THIS_READ=%s\n", start_read_txt);
						SUBREADprintf("LLLX: NEXT_READ=%s\n", next_read_txt);
					}
*/
					if((next_read_scoring > max_read_scoring[ WINDOW_SEARCH_TREE_WIDTH-1]) && is_rectifiable)
					{
						int is_unique = 1, is_replace = 0;
						for(xk3=0; xk3<WINDOW_SEARCH_TREE_WIDTH; xk3++)
						{

							if(global_context -> config.reassembly_window_alleles>1&&max_read_scoring[xk3]>0)
							{
								int position_diff = max_read_nonoverlap[xk3] - next_read_nonoverlap;
								int xk4;
								unsigned long long int test_prev_read_key = (window_no * 1llu +1)<<32 ;
								test_prev_read_key |=  (max_read_ids[xk3] & 0xfffff)<<12; 
								char  *prev_read_txt = HashTableGet(block_context->read_sequence_table,NULL+test_prev_read_key);
								int mismatch=0;	

								for(xk4 = abs(position_diff) ; ; xk4++)
								{
									int prev_read_pos = xk4-(position_diff>0?position_diff:0);
									int next_read_pos = xk4+(position_diff<0?position_diff:0);
									if(prev_read_txt[prev_read_pos]==0 || next_read_txt[next_read_pos]==0) break;
									mismatch += (prev_read_txt[prev_read_pos]!=next_read_txt[next_read_pos]);
								}


								if(mismatch<2)
								{
									is_unique = 0;
									if(position_diff > 0)
										is_replace = xk3+1;
								}
							}

						}
						if(is_replace)
						{
							max_read_scoring[is_replace-1] = next_read_scoring;
							max_read_ids[is_replace-1] = next_read_no;
							max_read_nonoverlap[is_replace-1] = next_read_nonoverlap;
							max_allele_quality[is_replace-1] = next_read_allele_quality;
						}
						else
						{
							if(is_unique) for(xk3=0; xk3<WINDOW_SEARCH_TREE_WIDTH; xk3++)
							{
								if(next_read_scoring > max_read_scoring[xk3])
								{
									int xk4;
									for(xk4 = WINDOW_SEARCH_TREE_WIDTH-1; xk4>xk3; xk4--)
									{
										max_read_scoring[xk4] = max_read_scoring[xk4-1];
										max_read_ids[xk4] = max_read_ids[xk4-1];
										max_read_nonoverlap[xk4] = max_read_nonoverlap[xk4-1];
										max_allele_quality[xk4] = max_allele_quality[xk4-1];
									}
									max_read_scoring[xk3] = next_read_scoring;
									max_read_ids[xk3] = next_read_no;
									max_read_nonoverlap[xk3] = next_read_nonoverlap;
									max_allele_quality[xk3] = next_read_allele_quality;
									//break;
								}
							}
						}
					}
				}

			}
		}


	for(i = 0; i<WINDOW_SEARCH_TREE_WIDTH; i++)
	{
		unsigned int next_read_no = max_read_ids[i];
		int is_fresh = 1;

//			if(block_context -> window_start_pos > 1124357402 && block_context -> window_start_pos < 1124359402)
//				SUBREADprintf("LLLX: TEST: %d\n", i);

		if(max_read_scoring[i]<=0) break;
		is_digged_in = 1;
		for(j=0; j< block_context -> used_read_number; j++)
		{
			if(block_context -> used_read_ids[j] == next_read_no)
			{
				is_fresh = 0;

				FRESH_COUNTER ++;
				break;
			}
		}

		if(!is_fresh)continue;

		block_context -> used_read_ids[block_context ->used_read_number++] = next_read_no;

		float next_allele_quality = max_allele_quality[i]; 
		int next_read_nonoverlap = max_read_nonoverlap[i];


		unsigned long long int test_next_read_key = (window_no * 1llu +1)<<32 ;
		test_next_read_key |=  ( next_read_no& 0xfffff)<<12; 


		char * next_read_txt = HashTableGet(block_context->read_sequence_table,NULL+test_next_read_key);
		assert(next_read_txt);

		int new_bases = (search_direction>0)?strlen(next_read_txt)-(read_len - next_read_nonoverlap):-next_read_nonoverlap;
		if (new_bases + block_context -> rebuilt_window_size > 1000){
			return 1;
		}

		int basic_length = 1;//max(read_len, strlen(next_read_txt))/3;
		int max_length = 99999;

		//if(anti_loop) max_length =  min(read_len, strlen(next_read_txt))/3;
		
		int current_rebuilt_size = block_context -> rebuilt_window_size;
		int current_matched_bases = block_context -> total_matched_bases;
		//int current_used_read_number = block_context -> used_read_number;
		if(new_bases>=basic_length && new_bases < max_length)
		{
			if(strlen(block_context -> rebuilt_window) != block_context ->rebuilt_window_size) return 1;

			//int is_rectifiable = rectify_read_text(global_context, block_context, window_no, &block_context->voting_indexes[window_no], next_read_txt , 1);
			//if(!is_rectifiable) continue;

			if(search_direction > 0)
			{
				int high_quality_offset = 0*min(2, strlen(next_read_txt) - new_bases);
				assert(high_quality_offset>=0);
				if(block_context ->rebuilt_window_size > high_quality_offset)
				{
					block_context -> rebuilt_window[ block_context ->rebuilt_window_size-high_quality_offset] =0;
					strcat(block_context -> rebuilt_window , next_read_txt + strlen(next_read_txt) - new_bases - high_quality_offset);
				}
				else
					strcat(block_context -> rebuilt_window , next_read_txt + strlen(next_read_txt) - new_bases);
			}
			else
			{
				int xk3;
				for(xk3= block_context -> rebuilt_window_size ; xk3>=0 ; xk3--)
					block_context -> rebuilt_window[new_bases + xk3] = block_context -> rebuilt_window[xk3];

				int high_quality_offset = 0*min(2, strlen(next_read_txt) - new_bases);
				//high_quality_offset = 0;
				assert(high_quality_offset>=0);
				memcpy(block_context -> rebuilt_window , next_read_txt , new_bases + high_quality_offset);
				block_context -> rebuilt_window [block_context -> rebuilt_window_size +new_bases] = 0;
			}

			block_context -> rebuilt_window_size += new_bases;
//			if(block_context -> window_start_pos > 1124357402 && block_context -> window_start_pos < 1124359402)
//				SUBREADprintf("LLLX: DIG IN: %d -> %d  (DIR=%d)\n", start_read_no, next_read_no, search_direction);

			search_window_once(global_context, block_context, window_no, next_read_no, search_direction, anti_loop, next_allele_quality);
			DIGIN_COUNTER ++;
			if(block_context -> search_cost > 1000) return 1;

		}
		if(search_direction > 0)
			// if the search is to the 3' end, the stack is simply dumpped.
			block_context -> rebuilt_window[current_rebuilt_size] = 0;
		else
		{
			// if the search is to the 5' end, we have to move the tail back to head.
			int xk4;
			int bytes_to_remove = block_context -> rebuilt_window_size - current_rebuilt_size;
			for(xk4 = bytes_to_remove; xk4 <  block_context -> rebuilt_window_size; xk4++) 
				block_context -> rebuilt_window[xk4-bytes_to_remove] = block_context -> rebuilt_window[xk4];
			block_context -> rebuilt_window[current_rebuilt_size] = 0;
		}

		block_context -> rebuilt_window_size = current_rebuilt_size;
		block_context -> total_matched_bases = current_matched_bases;
		//block_context -> used_read_number = current_used_read_number;
	}

//	if(block_context -> window_start_pos > 1124357402 && block_context -> window_start_pos < 1124359402)
//		SUBREADprintf("LLLX: NOT DIG IN: %d :: %d   DIR=%d\n", is_digged_in, block_context -> rebuilt_window_size, search_direction);

	if(!is_digged_in && block_context -> rebuilt_window_size > 0)
	{
		//if(block_context -> total_matched_bases * 10000 /block_context -> rebuilt_window_size > block_context -> max_matched_bases)
		int xk1, xk2;
		for(xk1 = 0; xk1 < global_context -> config.reassembly_window_alleles; xk1++)
		{
			//int is_same_allele = 0;

//			if(block_context -> window_start_pos > 1124357402 && block_context -> window_start_pos < 1124359402)
//				SUBREADprintf("LLLX [%d] CPY: %d<%d\n", xk1, block_context -> final_alleles[xk1].rebuilt_size , block_context -> rebuilt_window_size);

			if(block_context -> final_alleles[xk1].rebuilt_size < block_context -> rebuilt_window_size)
			{
				//int xk3=0;
				for(xk2= global_context -> config.reassembly_window_alleles - 1 ; xk2 > xk1; xk2--)
					memcpy(&block_context -> final_alleles[xk2], &block_context -> final_alleles[xk2 - 1], sizeof(struct reassmebly_window_allele));

				strcpy(block_context -> final_alleles[xk1].rebuilt_window , block_context -> rebuilt_window);
				block_context -> max_matched_bases = block_context -> total_matched_bases * 10000 / block_context -> rebuilt_window_size;
				block_context -> final_alleles[xk1].rebuilt_size = block_context -> rebuilt_window_size;
				block_context -> final_alleles[xk1].allele_quality = allele_quality;
				break;
			}
		}
	}

	return 0;
}


#define MAX_INDELS_IN_WINDOW 10  
#define PROBE_KEY_NUMBER 12
#define INDEL_FULL_ALIGN_MAX_MISMATCH 3
#define SINGLE_PROBE_MAX_MISMATCH 1

int full_indel_alignment(global_context_t * global_context, reassembly_by_voting_block_context_t * block_context, char * full_rebuilt_window, int full_rebuilt_window_size, gene_value_index_t * base_index, unsigned int window_start_pos, unsigned int * perfect_segment_start_pos, int * perfect_segment_lengths,int * indels_after_perfect_segments, short * indels_read_positions, float * indel_quality, unsigned int * contig_start_pos, unsigned int * contig_end_pos, int * head_removed_bases, int * tail_removed_bases, int * good_quality_indel)
{
	int xk1;
	int is_unreliable = 0;

	unsigned int xk2 = window_start_pos;
	unsigned int probe_poses [PROBE_KEY_NUMBER]; 
	int used_probe_in_rebuilt_window[PROBE_KEY_NUMBER];
	int used_probes = 0;

	if(window_start_pos>300) xk2-=300;
	if(full_rebuilt_window_size<101)return 0;
	unsigned int last_best_xk2 = xk2;

	memset(probe_poses ,0 , sizeof(unsigned int)*PROBE_KEY_NUMBER);
	(*contig_start_pos) = 0;
	(*contig_end_pos) = 0xffffffff;

	for(xk1=0; xk1<PROBE_KEY_NUMBER; xk1++)
	{
		char * probe_key = full_rebuilt_window + (full_rebuilt_window_size-global_context -> config.reassembly_key_length) * xk1 /(PROBE_KEY_NUMBER-1);

		int max_probe_match = -1;
		unsigned int next_best_xk2 = 0;


		// xk2 is the first base of the location matching this probe.
		for(xk2 = last_best_xk2; xk2 < REASSEMBLY_WINDOW_LENGTH + window_start_pos + 300 + MAX_INSERTION_LENGTH; xk2++)
		{
			int probe_match = match_chro(probe_key, base_index, xk2, global_context -> config.reassembly_key_length, 0,  global_context -> config.space_type);
			if(probe_match > max_probe_match )
			{
				next_best_xk2 = xk2;
				max_probe_match = probe_match;
			}
		}

		if(max_probe_match >= global_context -> config.reassembly_key_length - SINGLE_PROBE_MAX_MISMATCH)
		{
				// calculate the contig start;
			if((*contig_start_pos)==0)
			{
				(*contig_start_pos) = next_best_xk2 + 1 - (probe_key-full_rebuilt_window);
				(*head_removed_bases) = probe_key-full_rebuilt_window;
			}

			(*contig_end_pos) = next_best_xk2 + 1 - (probe_key-full_rebuilt_window) + full_rebuilt_window_size;
			(*tail_removed_bases) = full_rebuilt_window_size - (probe_key-full_rebuilt_window) - global_context -> config.reassembly_key_length;

			probe_poses[used_probes] = next_best_xk2+1;
			used_probe_in_rebuilt_window[used_probes] = (full_rebuilt_window_size-global_context -> config.reassembly_key_length) * xk1 /(PROBE_KEY_NUMBER-1);
			//if(max_probe_match<global_context -> config.reassembly_key_length - 1)
			//	is_unreliable += 2;

			used_probes +=1;
		}
		else
		{
			if(xk1 ==0 || xk1 == PROBE_KEY_NUMBER-1)
			{
				is_unreliable += 3;
			}
		}

	//	if(next_best_xk2>0)
	}


	for(xk1=0; xk1<used_probes - 1; xk1++)
	{
		if(probe_poses[xk1]>= probe_poses[xk1+1]) return 0;
	}

	// find indels between every two mapped probes.
	// and we know the indel length already.

	int ret = 0, total_unmatched = 0;

	//int H0_matched = match_chro(full_rebuilt_window ,  base_index, probe_poses[0]-1, full_rebuilt_window_size, 0, global_context -> config.space_type);

	int total_delta = 0;
	for(xk1=0; xk1 < used_probes-1; xk1++)
	{

		// indels_in_section  < 0 : insertion.

		// used_probe_in_rebuilt_window is the location of the probe on reads

		// probe_poses is the location of the probe on the chromosome.
		long long int indels_in_section = probe_poses[xk1+1];
		indels_in_section -= probe_poses[xk1];
		indels_in_section -= used_probe_in_rebuilt_window[xk1+1]; 
		indels_in_section += used_probe_in_rebuilt_window[xk1]; 

		total_delta += indels_in_section;
		if(indels_in_section)
		{
			// indels_in_section = -5: 5 insertion; indels_in_section = 3: 3 deletion
			// find the edge.
			int best_edge_score = -1;
			unsigned int section_best_edge = 0;

			// the smallest possible location of the second half is the location if the first half + deletion length
			xk2 = probe_poses[xk1] - 1 + max(0, indels_in_section);


			// find the exact indel point; xk2 is the tested point on the chromosome ( first WANTED base on the second half ).
			for(; xk2 < probe_poses[xk1+1] -1; xk2++)
			{
				int left_matched = match_chro(full_rebuilt_window + used_probe_in_rebuilt_window[xk1], base_index, probe_poses[xk1] - 1, (xk2 - probe_poses[xk1]+ 1 - max(0, indels_in_section)), 0,  global_context -> config.space_type);
				int right_matched = match_chro(full_rebuilt_window + used_probe_in_rebuilt_window[xk1] + (xk2 - probe_poses[xk1] + 1) -indels_in_section , base_index, xk2,  probe_poses[xk1+1] - 1 - xk2, 0,  global_context -> config.space_type);
				if(left_matched+right_matched > best_edge_score)
				{
					best_edge_score = left_matched+right_matched;
					section_best_edge = xk2;
				}
			}

			int section_mismatch = probe_poses[xk1+1]-probe_poses[xk1]-max(0, indels_in_section) - best_edge_score;
			is_unreliable += section_mismatch;


			perfect_segment_start_pos[ret] = probe_poses[xk1] - 1;
			// perfect_segment_lengths is : first WANTED base on the second half - first base pos on the first half - deletions 
			perfect_segment_lengths[ret] = section_best_edge - probe_poses[xk1] + 1 - max(0, indels_in_section);
			
			good_quality_indel[ret] = section_mismatch < INDEL_FULL_ALIGN_MAX_MISMATCH;
			indel_quality[ret] = pow(2,-is_unreliable);
			indels_after_perfect_segments[ret] = indels_in_section;
			indels_read_positions[ret] = used_probe_in_rebuilt_window[xk1] + perfect_segment_lengths[ret]; 
			ret++;
			total_unmatched += section_mismatch;
		}
		else
		{

			int matched = match_chro(full_rebuilt_window + used_probe_in_rebuilt_window[xk1], base_index, probe_poses[xk1] - 1, (probe_poses[xk1+1] - probe_poses[xk1]), 0,  global_context -> config.space_type);
			total_unmatched += ((probe_poses[xk1+1] - probe_poses[xk1]) - matched);
		}
	}

	if(ret>3|| total_unmatched>INDEL_FULL_ALIGN_MAX_MISMATCH || total_delta == 0)
	{
		return 0;
	}


	for(xk1=0; xk1<ret; xk1++)
	{
		indel_quality[xk1] *= pow(2.,-(ret * (1+total_unmatched)));
	}

	//printf("\nRET=%d\n\n", ret);

	return ret;

}



int match_str(char * c1, char * c2, int len)
{
	int xk1,ret=0;
	for(xk1=0; xk1<len; xk1++)
		ret += c1[xk1]==c2[xk1];
	return ret;
}


int find_best_location_for_probe(global_context_t * global_context ,gene_value_index_t * base_index, char * probe, unsigned int search_start, unsigned int search_len, unsigned int * ret_location)
{
	unsigned int xk2;
	unsigned int best_location = 0;
	int best_match = -1;
	for(xk2=search_start; xk2<search_start+search_len; xk2++)
	{
		int test_match = match_chro(probe, base_index, xk2, global_context -> config.reassembly_key_length, 0, global_context->config.space_type);
		if(test_match>best_match)
		{
			best_match = test_match;
			best_location = xk2;
		}
	}

	* ret_location = best_location;
	return best_match;

}

int find_potential_ultralong_indels(global_context_t * global_context , reassembly_by_voting_block_context_t * block_context , int window_no, char * contig_1, char * contig_2, unsigned int * potential_ins_left, unsigned int * potential_ins_len)
{
	int contig_len_1 = strlen(contig_1);
	int contig_len_2 = strlen(contig_2);
	unsigned int xk1, is_1_left;
	for(is_1_left = 0; is_1_left < 2; is_1_left ++)
	{
		int best_match = -1;

		char * test_key = is_1_left?contig_2:contig_1;
		char * testee = is_1_left?contig_1:contig_2;
		int testee_len = is_1_left?contig_len_1:contig_len_2;
		
		
		for(xk1=0; xk1< testee_len - global_context -> config.reassembly_key_length; xk1++)
		{
			int test_match = match_str(test_key, testee + xk1, global_context -> config.reassembly_key_length);
			if(test_match > best_match)
				best_match = test_match;
		}

		if(best_match >= global_context -> config.reassembly_key_length-1)
			return 0;
	}

	unsigned int window_left_pos = block_context->block_start_linear_pos + window_no * REASSEMBLY_WINDOW_LENGTH / 2 ;
	unsigned int search_start = window_left_pos;
	if(search_start>300)search_start-=300;
	else search_start=0;
	unsigned int search_end = 300 + global_context -> config.max_indel_length + window_left_pos;


	gene_value_index_t * base_index = &global_context->all_value_indexes[0] ;
	if(search_end > base_index -> length + base_index -> start_base_offset) search_end =  base_index -> length + base_index -> start_base_offset; 

	unsigned int contig_1_left_pos;
	int contig_1_left_match = find_best_location_for_probe(global_context, base_index, contig_1, search_start, search_end - search_start , &contig_1_left_pos);

	unsigned int contig_1_right_pos;
	int contig_1_right_match = find_best_location_for_probe(global_context, base_index, contig_1 + contig_len_1 - global_context -> config.reassembly_key_length, search_start, search_end - search_start , &contig_1_right_pos);



	unsigned int contig_2_left_pos;
	int contig_2_left_match = find_best_location_for_probe(global_context, base_index, contig_2, search_start, search_end - search_start , &contig_2_left_pos);

	unsigned int contig_2_right_pos;
	int contig_2_right_match = find_best_location_for_probe(global_context, base_index, contig_2 + contig_len_2 - global_context -> config.reassembly_key_length, search_start, search_end - search_start , &contig_2_right_pos);


	if(contig_1_left_match >= global_context -> config.reassembly_key_length -1 && contig_1_right_match >= global_context -> config.reassembly_key_length -1)
		return 0;

	if(contig_2_left_match >= global_context -> config.reassembly_key_length -1 && contig_2_right_match >= global_context -> config.reassembly_key_length -1)
		return 0;

	if(contig_1_left_match < global_context -> config.reassembly_key_length -1 && contig_1_right_match < global_context -> config.reassembly_key_length -1)
		return 0;

	if(contig_2_left_match < global_context -> config.reassembly_key_length -1 && contig_2_right_match < global_context -> config.reassembly_key_length -1)
		return 0;

	int ret = 0;
	unsigned int basic_ins_left = 0, basic_ins_right = 0;
	int is_1_at_left;

	if(contig_1_left_match >= global_context -> config.reassembly_key_length -1 && contig_2_right_match >= global_context -> config.reassembly_key_length -1  && contig_1_left_pos<contig_2_right_pos)
	{
		if(contig_len_1 + contig_len_2 - global_context -> config.reassembly_key_length > contig_2_right_pos - contig_1_left_pos)
		{
			*potential_ins_len = contig_len_1 + contig_len_2 - global_context -> config.reassembly_key_length - (contig_2_right_pos - contig_1_left_pos);
			*potential_ins_left = contig_1_left_pos;// + global_context -> config.reassembly_key_length;
			basic_ins_left = contig_1_left_pos;
			basic_ins_right = contig_2_right_pos;
			is_1_at_left = 1; 
			ret = 1;
		}
	}

	if(contig_1_right_match >= global_context -> config.reassembly_key_length -1  && contig_2_left_match >= global_context -> config.reassembly_key_length  -1 && contig_1_right_pos>contig_2_left_pos)
	{
		if( contig_len_1 + contig_len_2 - global_context -> config.reassembly_key_length> contig_1_right_pos - contig_2_left_pos)
		{
			*potential_ins_len = contig_len_1 + contig_len_2 - global_context -> config.reassembly_key_length - (contig_1_right_pos - contig_2_left_pos);
			*potential_ins_left = contig_2_left_pos;// + global_context -> config.reassembly_key_length;
			basic_ins_right = contig_1_right_pos;
			basic_ins_left = contig_2_left_pos;
			is_1_at_left = 0; 
			ret = 1;
		}
	}


	if(ret)
	{
		unsigned int x1=0,left_ins_edge = * potential_ins_left, best_edge = 0;
		char * left_contig = is_1_at_left ? contig_1: contig_2;
		char * right_contig = is_1_at_left ? contig_2: contig_1;
		char probe_window[global_context -> config.reassembly_key_length+1];
		probe_window[global_context -> config.reassembly_key_length]=0;

		for(; left_ins_edge < (*potential_ins_left) + global_context -> config.reassembly_key_length; left_ins_edge++)
			probe_window[x1++]=gvindex_get(base_index, left_ins_edge);

		for(; left_ins_edge < basic_ins_right; left_ins_edge++)
		{
			int delta = left_ins_edge - global_context -> config.reassembly_key_length - (*potential_ins_left);
			char * test_key = left_contig + delta;

			if(probe_window[global_context -> config.reassembly_key_length-1] == test_key[global_context -> config.reassembly_key_length-1]) best_edge = left_ins_edge ; 

			int matched = match_str(test_key, probe_window, global_context -> config.reassembly_key_length);
			if(matched< global_context -> config.reassembly_key_length-1) break;


			for(x1=0; x1<global_context -> config.reassembly_key_length - 1; x1++)
				probe_window[x1] = probe_window[x1+1];
			probe_window[global_context -> config.reassembly_key_length-1] = gvindex_get(base_index, left_ins_edge);
		}

		if(best_edge>0)
			*potential_ins_left = best_edge;

		if( ( basic_ins_right + global_context -> config.reassembly_key_length - best_edge) >strlen(right_contig) ||base_index -> length + base_index -> start_base_offset <= basic_ins_right + global_context -> config.reassembly_key_length || basic_ins_left <  base_index -> start_base_offset || best_edge >= base_index -> length + base_index -> start_base_offset) ret = 0;
		else 
		{
			int left_match = match_chro(left_contig, base_index, basic_ins_left, best_edge - basic_ins_left, 0, global_context->config.space_type);
			int right_match = match_chro(right_contig + strlen(right_contig)- ( basic_ins_right + global_context -> config.reassembly_key_length - best_edge), base_index, best_edge,  basic_ins_right + global_context -> config.reassembly_key_length - best_edge, 0, global_context->config.space_type);

			if(left_match+right_match< basic_ins_right+global_context -> config.reassembly_key_length-basic_ins_left - 2) ret=0;
		}
	}

	return ret;
}


int put_long_indel_event(global_context_t * global_context, unsigned int best_pos, int indels, float indel_quality, char * inserted_bases, int event_type)
{
	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 
	HashTable * event_table = indel_context -> event_entry_table;
	chromosome_event_t * event_space = indel_context -> event_space_dynamic;

	chromosome_event_t * found = NULL;
	chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];

	int pos_delta;
	for(pos_delta = -10; pos_delta <= 10; pos_delta++)
	{
		int found_events = search_event(global_context, event_table, event_space, best_pos - 1 + pos_delta, EVENT_SEARCH_BY_SMALL_SIDE, CHRO_EVENT_TYPE_LONG_INDEL|CHRO_EVENT_TYPE_INDEL, search_return);

		if(found_events)
		{
			int kx1; 
			for(kx1 = 0; kx1 < found_events ; kx1++)
			{
				if(search_return[kx1] -> indel_length == indels || event_type == CHRO_EVENT_TYPE_POTENTIAL_INDEL )
				{
					found = search_return[kx1];	
					break;
				}
			}
		}
		if(found)break;
	}

/*
				if(best_pos >= 1124357402&& best_pos < 1124359402)
				{
					fprintf(stderr, "LLLX-BB: %u :: FOUND=%p\n", best_pos, found);
				}

*/

	if(found)
	{
		found -> supporting_reads ++;
		found -> event_quality = max(found ->event_quality, indel_quality);
	}
	else 
	{
		int event_no;
		event_no = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) ->  total_events ++;
		event_space = reallocate_event_space(global_context , NULL, event_no);
		chromosome_event_t * new_event = event_space+event_no; 
		memset(new_event,0,sizeof(chromosome_event_t));

		if(indels<0 && inserted_bases)
			set_insertion_sequence(global_context, NULL , &new_event -> inserted_bases ,  inserted_bases , -indels);
		new_event -> event_small_side = best_pos - 1;
		new_event -> event_large_side = best_pos + max(0, indels);
		assert(new_event -> event_small_side && new_event -> event_large_side);
		new_event -> event_type = event_type;
		new_event -> indel_length = indels;
		new_event -> supporting_reads = 1;
		new_event -> event_quality = indel_quality;
		
		put_new_event(event_table, new_event , event_no);
	}

	return 0;
}


#define BASE_NEEDED_PER_WINDOW 60 

#define PILEUP_BOX_LEN_CERTAIN 12
#define MAX_PILEUP_BOX_COVERAGE_CERTAIN 12

#define PILEUP_BOX_LEN 50 
#define MAX_PILEUP_BOX_COVERAGE 50 

int finalise_pileup_file_by_voting(global_context_t * global_context , char * temp_file_name, char * chro_name, int block_no)
{
	int global_read_no = 0;
	char * read_text, * qual_text; 
	int xk1;
	int total_window_number =2 * (BASE_BLOCK_LENGTH / REASSEMBLY_WINDOW_LENGTH+1), window_adjust;

	unsigned int block_start_linear_pos = linear_gene_position(&global_context->chromosome_table, chro_name, block_no*BASE_BLOCK_LENGTH); 
	//if(block_start_linear_pos> 150000000) return 0;

	FILE * tmp_fp = f_subr_open(temp_file_name,"rb");
	if(!tmp_fp) 
		return 0;

	print_in_box(80,0,0,"Processing %s\n", temp_file_name);

	unsigned int * window_center_read_no = (unsigned int *) malloc(sizeof(unsigned int)* total_window_number * global_context -> config.reassembly_start_read_number);
	unsigned int * window_base_number = (unsigned int *) malloc(sizeof(unsigned int)* total_window_number);
	unsigned int * window_center_read_qual = (unsigned int *) malloc(sizeof(unsigned int)* total_window_number * global_context -> config.reassembly_start_read_number);
	unsigned int * window_center_read_masks = (unsigned int *) malloc(sizeof(unsigned int)* total_window_number);
	struct reassmebly_window_allele *first_half_alleles = (struct reassmebly_window_allele *) malloc(sizeof(struct reassmebly_window_allele ) * global_context -> config.reassembly_window_alleles);
	struct reassmebly_window_allele *second_half_alleles = (struct reassmebly_window_allele *) malloc(sizeof(struct reassmebly_window_allele ) * global_context -> config.reassembly_window_alleles);

	unsigned char * mapquality_lower_bounds_certain = calloc(1, BASE_BLOCK_LENGTH / PILEUP_BOX_LEN_CERTAIN * MAX_PILEUP_BOX_COVERAGE_CERTAIN);
	unsigned char * box_read_number_certain = calloc(1, BASE_BLOCK_LENGTH / PILEUP_BOX_LEN_CERTAIN);

	unsigned char * mapquality_lower_bounds = calloc(1, BASE_BLOCK_LENGTH / PILEUP_BOX_LEN * MAX_PILEUP_BOX_COVERAGE);
	unsigned char * box_read_number = calloc(1, BASE_BLOCK_LENGTH / PILEUP_BOX_LEN);

	char * full_rebuilt_window = (char *)malloc(3000);

	reassembly_by_voting_block_context_t block_context;

	block_context.voting_indexes = (gehash_t*)malloc(sizeof(gehash_t)* total_window_number);
	block_context.block_start_linear_pos = block_start_linear_pos;
	block_context.start_keys = (unsigned long long int *) calloc(sizeof(unsigned long long int),10 * total_window_number); 
	block_context.start_offsets = (short *)malloc(sizeof( short) * total_window_number);
	block_context.read_no_counter = (unsigned int *) calloc(sizeof(unsigned int),total_window_number); 
	block_context.read_rectify_space = (short *) malloc(sizeof(unsigned short) * 4 * MAX_READ_LENGTH);
	block_context.final_alleles = malloc(sizeof(struct reassmebly_window_allele ) * global_context -> config.reassembly_window_alleles);

	read_text = (char *) malloc(MAX_READ_LENGTH);
	qual_text = (char *) malloc(MAX_READ_LENGTH);

	HashTable *read_sequence_table=HashTableCreate(300000);
	HashTable *read_position_table=HashTableCreate(300000);
	HashTable *read_quality_table=HashTableCreate(300000);
	block_context.read_sequence_table = read_sequence_table;
	block_context.read_quality_table = read_quality_table;
	block_context.read_position_table = read_position_table;

	memset(block_context.final_alleles, 0, sizeof(struct reassmebly_window_allele ) * global_context -> config.reassembly_window_alleles);
	memset(window_center_read_no, 0, sizeof(unsigned int)* total_window_number * global_context -> config.reassembly_start_read_number);
	memset(window_center_read_qual, 0, sizeof(unsigned int)* total_window_number*global_context -> config.reassembly_start_read_number);
	memset(window_center_read_masks, 0, sizeof(unsigned int)* total_window_number);
	memset(window_base_number, 0, sizeof(unsigned int)* total_window_number);


	if(!tmp_fp) return 1;

	while(!feof(tmp_fp))
	{
		short read_len;
		long long int first_base_pos;
		//unsigned int qual_of_this_read = 1;
		base_block_temp_read_t read_rec;
		int rlen = fread(&read_rec, sizeof(read_rec), 1, tmp_fp);
		int is_position_certain = read_rec.strand;
		if(rlen<1) break;

		rlen = fread(&read_len, sizeof(short), 1, tmp_fp);
		if(rlen < 1){
			SUBREADprintf("ERROR: Unable to parse piles\n");
			return -1;
		}
		rlen = fread(read_text, sizeof(char), read_len, tmp_fp);
		if(rlen <read_len){
			SUBREADprintf("ERROR: Unable to parse piles\n");
			return -1;
		}
		rlen = fread(qual_text, sizeof(char), read_len, tmp_fp);
		if(rlen <read_len){
			SUBREADprintf("ERROR: Unable to parse piles\n");
			return -1;
		}

		first_base_pos = read_rec.pos;
		first_base_pos -= block_no * BASE_BLOCK_LENGTH;

		for(window_adjust = 0; window_adjust < global_context -> config.reassembly_window_multiplex ; window_adjust ++)
		{
			int window_no = 2*(first_base_pos / REASSEMBLY_WINDOW_LENGTH)+ window_adjust , xk2;
			if(window_no >=total_window_number)
				break;
			window_base_number[window_no]+=read_len;

			int qual_of_this_read = 127, read_offset;
			for(read_offset = 0; read_offset < read_len ; read_offset++)
			{
				int nch = qual_text[read_offset] - '!';
				if(nch < 14) qual_of_this_read -=4;
				else if(nch < 22) qual_of_this_read -= 2;
				else if(nch < 33) qual_of_this_read --;
			}
			qual_of_this_read=max(1,qual_of_this_read);

			
			if(is_position_certain)
			{
				int box_no = first_base_pos / PILEUP_BOX_LEN_CERTAIN; 
				for(xk1 = box_no * MAX_PILEUP_BOX_COVERAGE_CERTAIN ; xk1 < (box_no+1) * MAX_PILEUP_BOX_COVERAGE_CERTAIN; xk1++)
				{
					if(qual_of_this_read > mapquality_lower_bounds_certain[xk1])
					{
						for(xk2 = (box_no+1) * MAX_PILEUP_BOX_COVERAGE_CERTAIN - 1; xk2 > xk1; xk2--)
							mapquality_lower_bounds_certain[xk2] = mapquality_lower_bounds_certain[xk2-1];
						mapquality_lower_bounds_certain[xk1] = qual_of_this_read;
						break;
					}
				}
			}
			else
			{
				int box_no = first_base_pos / PILEUP_BOX_LEN; 
				for(xk1 = box_no * MAX_PILEUP_BOX_COVERAGE ; xk1 < (box_no+1) * MAX_PILEUP_BOX_COVERAGE; xk1++)
				{
					if(qual_of_this_read > mapquality_lower_bounds[xk1])
					{
						for(xk2 = (box_no+1) * MAX_PILEUP_BOX_COVERAGE - 1; xk2 > xk1; xk2--)
							mapquality_lower_bounds[xk2] = mapquality_lower_bounds[xk2-1];
						mapquality_lower_bounds[xk1] = qual_of_this_read;
						break;
					}
				}

			}
			if(is_position_certain && window_adjust == 1) break;
		}

	}
	fclose(tmp_fp);

	for(xk1=0; xk1< total_window_number; xk1++)
	{
		if(window_base_number[xk1] >= BASE_NEEDED_PER_WINDOW)
		{
			gehash_create(&block_context.voting_indexes[xk1], max(4000,window_base_number[xk1]*4), 1);
			block_context . start_offsets[xk1] = 0x7fff;
		}
	}


	tmp_fp = f_subr_open(temp_file_name,"rb");
	while(!feof(tmp_fp))
	{
		short read_len;
		long long int first_base_pos;
		base_block_temp_read_t read_rec;
		int rlen = fread(&read_rec, sizeof(read_rec), 1, tmp_fp);
		int is_position_certain = read_rec.strand;

		if(rlen<1) break;
		rlen = fread(&read_len, sizeof(short), 1, tmp_fp);
		if(rlen <1){
			SUBREADprintf("ERROR: cannot parse piles.\n");
			return -1;
		}

		rlen = fread(read_text, sizeof(char), read_len, tmp_fp);
		if(rlen <read_len){
			SUBREADprintf("ERROR: cannot parse piles.\n");
			return -1;
		}

		rlen = fread(qual_text, sizeof(char), read_len, tmp_fp);
		if(rlen <read_len){
			SUBREADprintf("ERROR: cannot parse piles.\n");
			return -1;
		}
		first_base_pos = read_rec.pos;
		first_base_pos -= block_no * BASE_BLOCK_LENGTH;
		assert(first_base_pos>=0);
		read_text [read_len]=0;
		qual_text [read_len]=0;

		is_bad_read(read_text, read_len, "LOAD");

		int read_offset;

		/*
		unsigned int qual_of_this_read = 1;
		for(read_offset = 0; read_offset < read_len ; read_offset++)
			qual_of_this_read += qual_text[read_offset];

		qual_of_this_read = qual_of_this_read * 100 / read_len + (is_position_certain?200:0);
		if(qual_of_this_read<1)qual_of_this_read=1;
		*/
			
		int qual_of_this_read = 127;
		for(read_offset = 0; read_offset < read_len ; read_offset++)
		{
			int nch = qual_text[read_offset] - '!';


			if(nch < 14) qual_of_this_read -=4;
			else if(nch < 22) qual_of_this_read -= 2;
			else if(nch < 33) qual_of_this_read --;
		}

		qual_of_this_read=max(1,qual_of_this_read);
	
		for(window_adjust = 0; window_adjust < global_context -> config.reassembly_window_multiplex; window_adjust++)
		{
			if(is_position_certain && window_adjust >= 2) break;
			int window_no = 2*(first_base_pos / REASSEMBLY_WINDOW_LENGTH) + window_adjust;

			if(window_no >= total_window_number) break;
			if(window_base_number[window_no]<BASE_NEEDED_PER_WINDOW) continue;


			if(is_position_certain)
			{
				int box_no = first_base_pos / PILEUP_BOX_LEN_CERTAIN; 
				if(qual_of_this_read < mapquality_lower_bounds_certain[box_no * MAX_PILEUP_BOX_COVERAGE_CERTAIN + MAX_PILEUP_BOX_COVERAGE_CERTAIN - 1])
					continue;
				if(box_read_number_certain[box_no]>MAX_PILEUP_BOX_COVERAGE_CERTAIN)continue;
				box_read_number_certain[box_no]++;
				//qual_of_this_read += 200;
			}
			else
			{
				int box_no = first_base_pos / PILEUP_BOX_LEN; 
	
				if(qual_of_this_read < mapquality_lower_bounds[box_no * MAX_PILEUP_BOX_COVERAGE + MAX_PILEUP_BOX_COVERAGE - 1])
					continue;
				if(box_read_number[box_no]>MAX_PILEUP_BOX_COVERAGE)continue;
				box_read_number[box_no]++;
			}

			int read_no = block_context.read_no_counter[window_no];
			block_context.read_no_counter[window_no]++;


			if(window_center_read_qual[window_no * global_context -> config.reassembly_start_read_number + global_context -> config.reassembly_start_read_number-1] < qual_of_this_read)//&& is_position_certain)
			{

				int xk44,xk3;
				for(xk44=0; xk44<global_context -> config.reassembly_start_read_number; xk44++)
					if(window_center_read_qual[window_no * global_context -> config.reassembly_start_read_number + xk44] < qual_of_this_read)
						break;
				for(xk3=global_context -> config.reassembly_start_read_number -2; xk3>=xk44 ; xk3--)
				{
					window_center_read_no[window_no * global_context -> config.reassembly_start_read_number + xk3 +1] =   window_center_read_no[window_no * global_context -> config.reassembly_start_read_number + xk3] ;
					window_center_read_qual[window_no * global_context -> config.reassembly_start_read_number + xk3 +1] = window_center_read_qual[window_no * global_context -> config.reassembly_start_read_number + xk3];
				}

				window_center_read_no[window_no * global_context -> config.reassembly_start_read_number + xk44] = read_no + 1;
				window_center_read_qual[window_no  * global_context -> config.reassembly_start_read_number + xk44] = qual_of_this_read; 

			}
			for(read_offset = 0; read_offset < read_len - global_context -> config.reassembly_subread_length+1; read_offset += REASSEMBLY_SUBREAD_GAP)
			{
				unsigned int subread_int = 0;
				int xk2;
				for(xk2=0; xk2< global_context -> config.reassembly_subread_length; xk2++)
					subread_int = (subread_int << 2) | base2int(read_text[read_offset + xk2]);

				gehash_insert_limited(&block_context.voting_indexes[window_no], subread_int , read_offset +10000 *read_no, 400, 10000);
			}



		/*
			if(strstr(temp_file_name, "chr1-0005"))
				if(window_no == 63763 || window_no == 63762)
				{
					SUBREADprintf("LLLX-INSERT_WIN: [%d] : %s  RN=%d  RL=%d  ITEMS=%llu\n",window_no, read_text, read_no, read_len, block_context.voting_indexes[window_no].current_items);

				}
		*/

			char * sequence_for_hash = (char *) malloc(read_len+1);
			assert(sequence_for_hash);
			strcpy(sequence_for_hash, read_text);
			unsigned long long int read_hashed_key = (window_no * 1llu +1)<<32 ;
			read_hashed_key |= (read_no & 0xfffff)<<12;

			if(is_position_certain)
				HashTablePut(read_position_table, NULL + read_hashed_key, NULL+first_base_pos+1);
			HashTablePut(read_sequence_table, NULL + read_hashed_key, sequence_for_hash);
			HashTablePut(read_quality_table, NULL + read_hashed_key, NULL + qual_of_this_read);
		}

		global_read_no ++;
	}
	fclose(tmp_fp);

	for(xk1=0; xk1< total_window_number; xk1++)
	{
		if(window_base_number[xk1]>=BASE_NEEDED_PER_WINDOW) 
			gehash_sort(&block_context.voting_indexes[xk1]);
	}


	gene_vote_t *vote_list;
	vote_list = (gene_vote_t *) malloc(sizeof(gene_vote_t));
	assert(vote_list);
	block_context.vote_list = vote_list;
	block_context.vote_list_rectify = (gene_vote_t *) malloc(sizeof(gene_vote_t));

	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 

	HashTable * event_table = indel_context -> event_entry_table;

	for(xk1 = 0; xk1 <  total_window_number; xk1++)
	{
		// xk1 is the number of the current window.

		if(window_base_number[xk1]<BASE_NEEDED_PER_WINDOW)
		{
			//printf("SKIPPED : %d : %d/%d\n", xk1, window_base_number[xk1], BASE_NEEDED_PER_WINDOW);
			continue;
		}
		// full_window = 5'end -> anchor_read -> 3'end



		//int IGNORED_BASES_IN_WINDOW  = 10;
		//for(IGNORED_BASES_IN_WINDOW = 0; IGNORED_BASES_IN_WINDOW <=60; IGNORED_BASES_IN_WINDOW+=30)
		{
			int start_read_i;
			for(start_read_i = 0; start_read_i < global_context -> config.reassembly_start_read_number; start_read_i ++)
			{

				// from middle read to 5' end

				if(!window_center_read_no[xk1 * global_context -> config.reassembly_start_read_number + start_read_i ]) continue;

				unsigned int this_start_read_no = window_center_read_no[xk1 * global_context -> config.reassembly_start_read_number + start_read_i ]-1 ;
				unsigned int window_start_pos = block_context.block_start_linear_pos + xk1 * REASSEMBLY_WINDOW_LENGTH / 2; 
				block_context.window_start_pos = window_start_pos;

				search_window(global_context, &block_context, xk1, this_start_read_no , -1);
				memcpy(first_half_alleles, block_context.final_alleles, global_context -> config.reassembly_window_alleles * sizeof(struct reassmebly_window_allele ));
				search_window(global_context, &block_context, xk1, this_start_read_no,  1);
				memcpy(second_half_alleles, block_context.final_alleles, global_context -> config.reassembly_window_alleles * sizeof(struct reassmebly_window_allele));

				/*
				if(window_start_pos >= 1124357402&& window_start_pos < 1124359402)
				{
					fprintf(stderr, "LLLX: %u [%d] :: START READ = %u ; FN=%s ; BASES=%u\n", window_start_pos, xk1, this_start_read_no, temp_file_name , window_base_number[xk1]);
					char fname [100];
					sprintf(fname, "DP-%d-%d", getpid(), xk1);
					gehash_dump(&block_context.voting_indexes[xk1], fname);
				}*/

				unsigned long long int starter_read_key = (xk1 * 1llu +1)<<32 ;
				starter_read_key |=  ( (this_start_read_no)& 0xfffff)<<12; 

				char * start_read_txt = HashTableGet(block_context.read_sequence_table,NULL+starter_read_key);
				assert(start_read_txt);

				int first_allele_no, second_allele_no;

				for(first_allele_no = 0; first_allele_no <  global_context -> config.reassembly_window_alleles && first_half_alleles[first_allele_no].rebuilt_size>0 ; first_allele_no++)
					for(second_allele_no = 0; second_allele_no <  global_context -> config.reassembly_window_alleles && second_half_alleles[second_allele_no].rebuilt_size>0  ; second_allele_no++)
					{
						int all_indels_in_window = 0;
						strcpy(full_rebuilt_window, first_half_alleles [first_allele_no]. rebuilt_window);
						strcat(full_rebuilt_window, start_read_txt);
						strcat(full_rebuilt_window, second_half_alleles[second_allele_no].rebuilt_window);

			
						int xk2, full_rebuilt_window_size = strlen(full_rebuilt_window);

						gene_value_index_t * base_index = &global_context->all_value_indexes[0] ;

						unsigned int perfect_segment_start_pos[MAX_INDELS_IN_WINDOW], contig_start_pos, contig_end_pos, all_fresh = 0;
						int perfect_segment_lengths[MAX_INDELS_IN_WINDOW], head_removed_bases=0, tail_removed_bases=0;
						int indels_after_perfect_segments[MAX_INDELS_IN_WINDOW];
						int good_quality_indels[MAX_INDELS_IN_WINDOW];
						float quality_of_indels_aln[MAX_INDELS_IN_WINDOW];
						short indels_read_positions[MAX_INDELS_IN_WINDOW];
						char contig_CIGAR[200];

						memset(good_quality_indels, 0, sizeof(int) * MAX_INDELS_IN_WINDOW);
						int indels_in_window = full_indel_alignment(global_context, &block_context, full_rebuilt_window, full_rebuilt_window_size, base_index, window_start_pos, perfect_segment_start_pos, perfect_segment_lengths, indels_after_perfect_segments, indels_read_positions, quality_of_indels_aln, &contig_start_pos, &contig_end_pos, &head_removed_bases, &tail_removed_bases, good_quality_indels);
						contig_CIGAR[0]=0;						
						int read_position_cursor = head_removed_bases;
						int is_indel_contig = 0;

						for(xk2 = 0; xk2 < indels_in_window; xk2++)
						{
							unsigned int best_pos = perfect_segment_start_pos [xk2] + perfect_segment_lengths[xk2];
							int indels = indels_after_perfect_segments[xk2], is_fresh = 1;
							float quality_of_this_indel = quality_of_indels_aln[xk2] * first_half_alleles [first_allele_no].allele_quality * second_half_alleles[second_allele_no].allele_quality;

							if(abs(indels) >= global_context -> config.max_indel_length){
								//#warning "====================== REMOVE THE debug output ==============="
								continue;
							}

							/*if(abs(indels) <= 16 && 0){
								continue;
							}*/

							is_indel_contig=1;
							sprintf(contig_CIGAR+strlen(contig_CIGAR), "%dM%d%c", indels_read_positions[xk2] - read_position_cursor, abs(indels), indels<0?'I':'D');  
							read_position_cursor = indels_read_positions[xk2];
							if(indels<0) read_position_cursor -= indels;

							all_indels_in_window++;

							if(indels >= -global_context -> config.max_indel_length && good_quality_indels[xk2])
							{
								int neighbour_delta;
								for(neighbour_delta = - 30; neighbour_delta < 30 ; neighbour_delta++)
								{
									chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];
									chromosome_event_t * event_space = indel_context -> event_space_dynamic;

									int xk3, found_events = search_event(global_context, event_table, event_space, best_pos + neighbour_delta , EVENT_SEARCH_BY_SMALL_SIDE, CHRO_EVENT_TYPE_LONG_INDEL|CHRO_EVENT_TYPE_INDEL, search_return);

									for(xk3 = 0; xk3<found_events; xk3++)
									{
										chromosome_event_t * tested_neighbour = search_return[xk3];

										if(tested_neighbour -> event_quality >= quality_of_this_indel)
											is_fresh = 0;
										else tested_neighbour -> event_type = CHRO_EVENT_TYPE_REMOVED; 
									}
					

								}


								if(is_fresh)
									put_long_indel_event(global_context, best_pos,  indels, quality_of_this_indel, full_rebuilt_window + indels_read_positions[xk2], CHRO_EVENT_TYPE_LONG_INDEL);
								all_fresh += is_fresh;
							}
						}
	

						//#warning "====================== Make sure ' is_indel_contig '  in the next line is harmless before RELEASE ==============="
						if( is_indel_contig || all_fresh)
						{
							int write_cursor;
							char * chro_begin;
							int chro_offset_start = 0;
							int chro_offset_end = 0;
							locate_gene_position_max(contig_start_pos + head_removed_bases ,& global_context -> chromosome_table, &chro_begin, &chro_offset_start, NULL, NULL, 0);
							locate_gene_position_max(contig_end_pos - tail_removed_bases ,& global_context -> chromosome_table, &chro_begin, &chro_offset_end, NULL, NULL, 0);
							if(full_rebuilt_window_size - read_position_cursor - tail_removed_bases)
								fprintf(global_context -> long_insertion_FASTA_fp, ">%s-%u-%u-%s%dM\n", chro_name, chro_offset_start, chro_offset_end - 1, contig_CIGAR, full_rebuilt_window_size - read_position_cursor - tail_removed_bases);
							else
								fprintf(global_context -> long_insertion_FASTA_fp, ">%s-%u-%u-%s0M\n", chro_name, chro_offset_start, chro_offset_end - 1, contig_CIGAR);
							full_rebuilt_window[full_rebuilt_window_size - tail_removed_bases] = 0;
							for(write_cursor = head_removed_bases ; write_cursor < full_rebuilt_window_size - tail_removed_bases ; write_cursor += 70)
							{
								fwrite(full_rebuilt_window + write_cursor , min(70, full_rebuilt_window_size - write_cursor), 1, global_context -> long_insertion_FASTA_fp);
								fputs("\n",global_context -> long_insertion_FASTA_fp);
							}
						}

					}
			}
		}

		
	}

	free(vote_list);
	free(block_context.vote_list_rectify);

	for(xk1=0; xk1< total_window_number; xk1++)
	{
		if(window_base_number[xk1]>=BASE_NEEDED_PER_WINDOW)
			gehash_destory(&block_context.voting_indexes[xk1]);
	}


	free_values_destroy(read_sequence_table);

	free(mapquality_lower_bounds_certain);
	free(mapquality_lower_bounds);
	free(box_read_number);
	free(box_read_number_certain);
	free(full_rebuilt_window);
	free(block_context.read_no_counter);
	free(block_context.voting_indexes);
	free(block_context.start_keys);
	free(block_context.start_offsets);
	free(block_context.read_rectify_space);
	free(block_context.final_alleles);
	free(first_half_alleles);
	free(second_half_alleles);
	free(window_base_number);
	free(window_center_read_no);
	free(window_center_read_qual);
	free(window_center_read_masks);
	free(read_text);
	free(qual_text);
	

	return 0;
}
int finalise_pileup_file_by_debrujin(global_context_t * global_context , char * temp_file_name, char * chro_name, int block_no)
{
	FILE * tmp_fp = f_subr_open(temp_file_name,"rb");
	if(!tmp_fp)
		return 0;
	char * read_text, * qual_text; 
	int xk1;
	unsigned int block_start_linear_pos = linear_gene_position(&global_context->chromosome_table, chro_name, block_no*BASE_BLOCK_LENGTH); 
	int total_window_number =2 * (BASE_BLOCK_LENGTH / REASSEMBLY_WINDOW_LENGTH+1);

	reassembly_block_context_t * block_context = (reassembly_block_context_t*)malloc(sizeof(reassembly_block_context_t));
	block_context -> block_start_linear_pos = block_start_linear_pos;
	block_context -> start_keys = (unsigned long long int *) calloc(sizeof(unsigned long long int),10 * total_window_number); 
	block_context -> start_offsets = (short *)malloc(sizeof( short) * total_window_number);

	assert(block_context -> start_offsets);
	assert(block_context -> start_keys);

	block_context -> de_bruijn_graphs = ( HashTable ** )malloc(sizeof( HashTable *) * total_window_number);
	for(xk1=0; xk1< total_window_number; xk1++)
	{
		block_context -> de_bruijn_graphs[xk1] = HashTableCreate(10* REASSEMBLY_WINDOW_LENGTH);
		block_context -> start_offsets[xk1] = 0x7fff;
	}

	read_text = (char *) malloc(MAX_READ_LENGTH);
	qual_text = (char *) malloc(MAX_READ_LENGTH);
	while(!feof(tmp_fp))
	{
		short read_len;
		unsigned int first_base_pos;
		base_block_temp_read_t read_rec;
		int rlen = fread(&read_rec, sizeof(read_rec), 1, tmp_fp);
		if(rlen<1) break;
		rlen = fread(&read_len, sizeof(short), 1, tmp_fp);
		if(rlen<1) return -1;
		rlen = fread(read_text, sizeof(char), read_len, tmp_fp);
		if(rlen<read_len) return -1;
		rlen = fread(qual_text, sizeof(char), read_len, tmp_fp);
		if(rlen<read_len) return -1;
		first_base_pos = read_rec.pos - block_no * BASE_BLOCK_LENGTH;

		insert_pileup_read(global_context , block_context , first_base_pos , read_text , qual_text , read_len, read_rec.read_pos>0);
	}
	
	for(xk1=0; xk1< total_window_number; xk1++)
	{
		if(block_context -> start_offsets[xk1]<0x7fff && block_context -> de_bruijn_graphs[xk1] -> numOfElements >3)
		{
			finalise_db_graph(global_context, block_context, xk1, block_context -> de_bruijn_graphs[xk1] , block_context -> start_keys [xk1*10], block_context -> start_offsets [xk1]);
		}
		HashTableDestroy(block_context -> de_bruijn_graphs[xk1]);
	}


	free(block_context -> de_bruijn_graphs);
	free(block_context -> start_keys);
	free(block_context -> start_offsets);
	free(block_context);
	free(read_text);
	free(qual_text);
	fclose(tmp_fp);
	return 0;
}

int finalise_long_insertions_by_hashtable(global_context_t * global_context)
{
	int chro_i;
	assert(global_context -> index_block_number == 1);
	unsigned int chro_start_pos = 0;
	char tmp_fname[MAX_FILE_NAME_LENGTH+40];

	sprintf(tmp_fname,"%s.reassembly.fa", global_context->config.output_prefix);
	global_context->long_insertion_FASTA_fp = f_subr_open(tmp_fname ,"wb");

	for(chro_i=0; chro_i < global_context -> chromosome_table.total_offsets ; chro_i++)
	{
		unsigned int chro_length = global_context -> chromosome_table.read_offsets[chro_i] - chro_start_pos;
		unsigned int chro_current_pos = 0;
		for(chro_current_pos = 0; chro_current_pos < chro_length ; chro_current_pos += BASE_BLOCK_LENGTH)
		{
			char temp_file_name[MAX_FILE_NAME_LENGTH + 50];
			int block_no = chro_current_pos / BASE_BLOCK_LENGTH;

			sprintf(temp_file_name,"%s@%s-%04u.bin", global_context -> config.temp_file_prefix, global_context -> chromosome_table.read_names+chro_i*MAX_CHROMOSOME_NAME_LEN , block_no );

			finalise_pileup_file(global_context , temp_file_name, global_context -> chromosome_table.read_names+MAX_CHROMOSOME_NAME_LEN*chro_i, block_no);
	
			unlink(temp_file_name);
		}
		chro_start_pos = global_context -> chromosome_table.read_offsets[chro_i];
	}

	fclose(global_context -> long_insertion_FASTA_fp);

	/*
	for(xk1=0; xk1 < global_context -> config.long_indel_iterations ; xk1++)
	{
		explorer_insertions_in_unmapped(global_context);
	}
	*/

	return 0;

}
void destroy_pileup_table(HashTable* local_reassembly_pileup_files)
{
	int bucket;
	KeyValuePair * cursor;

	for(bucket=0; bucket< local_reassembly_pileup_files -> numOfBuckets; bucket++)
	{
		cursor = local_reassembly_pileup_files -> bucketArray[bucket];
		while (1)
		{
			if (!cursor) break;
			FILE * fp = (FILE *)cursor->value;
			if(fp != NULL+1)
				fclose(fp);
			free((void *)cursor->key);
			cursor = cursor->next;
		}
	}
	
	HashTableDestroy(local_reassembly_pileup_files);


}

char * _COREMAIN_delete_temp_prefix = NULL;
void COREMAIN_SIGINT_hook(int param)
{
	#ifdef MAKE_STANDALONE
	int xk1, last_slash = -1;
	if(_COREMAIN_delete_temp_prefix != NULL)
	{
		char del2[MAX_FILE_NAME_LENGTH], del_suffix[MAX_FILE_NAME_LENGTH], del_name[MAX_FILE_NAME_LENGTH];
		SUBREADprintf("\n\nReceived a terminal signal. The temporary files were removed.\n");
		for(xk1=0; _COREMAIN_delete_temp_prefix[xk1]; xk1++)
		{
			if(_COREMAIN_delete_temp_prefix[xk1]=='/') last_slash = xk1;
			else if(_COREMAIN_delete_temp_prefix[xk1]=='\\')
			{
				SUBREADprintf("The file name is unknown.\n");
				return;
			}
		}
		if(last_slash>=0)
		{
			memcpy(del2, _COREMAIN_delete_temp_prefix, last_slash);
			del2[last_slash] = 0;
			strcpy(del_suffix , _COREMAIN_delete_temp_prefix + last_slash + 1);
		}
		else
		{
			strcpy(del2,".");
			strcpy(del_suffix , _COREMAIN_delete_temp_prefix);
		}
	
		if(strlen(del_suffix)>8)
		{
			DIR           *d;
			struct dirent *dir;

			d = opendir(del2);
			if (d)
			{
				while ((dir = readdir(d)) != NULL)
				{
					if(strstr(dir->d_name, del_suffix))
					{
						//printf("%s\n", dir->d_name);
						strcpy(del_name, del2);
						strcat(del_name, "/");
						strcat(del_name, dir->d_name);
						unlink(del_name);
					}
				}
				closedir(d);
			}
		}
			
	}

	exit(param);
	#endif
}



int finalise_long_insertions(global_context_t * global_context)
{
	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 
	HashTable * local_reassembly_pileup_files = indel_context -> local_reassembly_pileup_files;

	destroy_pileup_table(local_reassembly_pileup_files);

	return finalise_long_insertions_by_hashtable(global_context);
}

void init_global_context(global_context_t * context)
{
	memset(context->module_contexts, 0, 5*sizeof(void *));
	memset(&context->config, 0, sizeof(configuration_t));

	context->config.fast_run = 0;
	context->config.memory_use_multiplex = 1;
	context->config.report_sam_file = 1;
	context->config.do_breakpoint_detection = 0;
	context->config.do_fusion_detection = 0;
	context->config.do_long_del_detection = 0;
	context->config.do_structural_variance_detection = 0;
	context->config.more_accurate_fusions = 1;
	context->config.report_multi_mapping_reads = 0;

	//#warning "============= best values for the SVs application: 8; 5; 32 ==============="
	context->config.top_scores = 8 - 5;
	context->config.max_vote_combinations = 5 - 3;
	context->config.max_vote_simples = 5;
	context->config.max_vote_number_cutoff = 1;

	context->config.experiment_type = 0;
	context->config.prefer_donor_receptor_junctions = 1;
	context->config.maximum_translocation_length = 10000;
	context->config.maximum_colocating_distance = 500;
	context->config.do_big_margin_filtering_for_reads = 0;
	context->config.do_big_margin_filtering_for_junctions = 0;
	context->config.maximum_intron_length = 500000;
	context->config.use_hamming_distance_in_exon = 0;
	context->config.is_third_iteration_running = 0;
	context->input_reads.is_paired_end_reads = 0;
	context->will_remove_input_file = 0;
	context->config.ignore_unmapped_reads = 0;
	context->config.report_unmapped_using_mate_pos = 1;
	context->config.downscale_mapping_quality=0;
	context->config.ambiguous_mapping_tolerance = 39;
	context->config.use_hamming_distance_break_ties = 0;
	context->config.use_quality_score_break_ties = 0;
	context->config.extending_search_indels = 0;
	context->config.PE_predominant_weight = 0;

	//#warning "============= best values for the SVs application: 3 ====================="
	context->config.multi_best_reads = 1;
	context->config.reported_multi_best_reads = 1;
	context->config.is_SAM_file_input=0;
	context->config.use_dynamic_programming_indel=0;
	context->config.use_bitmap_event_table = 1;
	context->config.convert_color_to_base = 0;
	context->config.is_gzip_fastq = 0;

	context->config.is_BAM_output = 1;
	context->config.is_BAM_input = 0;
	context->config.is_input_read_order_required = 0;
	context->config.read_trim_5 = 0;
	context->config.read_trim_3 = 0;
	context->config.minimum_exonic_subread_fraction = -1.0;
	context->config.is_first_read_reversed = 0;
	context->config.is_second_read_reversed = 1;
	context->config.space_type = GENE_SPACE_BASE;
	context->config.minimum_pair_distance = 50;
	context->config.maximum_pair_distance = 600;
	context->config.expected_pair_distance = 300;
	context->config.restrected_read_order = 1;
	context->config.k_mer_length = 28;

	context->config.reassembly_start_read_number = 2;
	context->config.reassembly_key_length = 28;
	context->config.reassembly_subread_length = 12;
	context->config.reassembly_window_multiplex = 3;
	context->config.reassembly_tolerable_voting = 1;
	context->config.reassembly_window_alleles = 1;
	context->config.do_superlong_indel_detection = 0;
	context->config.flanking_subread_indel_mismatch = 1;

	context->config.total_subreads = 10;
	context->config.minimum_subread_for_first_read = 3;
	context->config.minimum_subread_for_second_read = 1;
	context->config.min_mapped_fraction = 0;
	context->config.max_mismatch_exonic_reads = 5;
	context->config.max_mismatch_junction_reads = 2;
	context->config.all_threads = 1;
	context->config.is_first_iteration_running = 1;
	context->config.is_second_iteration_running = 1;

	context->config.reads_per_chunk = 20*1024*1024;

//#warning "=========== 2*1024*1024 IS FOR TESTING BLOCKING AND SHOULD BE COMMENTED ==============="
//	context->config.reads_per_chunk = 2*1024*1024;

	context->config.use_memory_buffer = 1;
	context->config.is_methylation_reads = 0;
	context->config.report_no_unpaired_reads = 0;
	context->config.limited_tree_scan = 0;
	context->config.high_quality_base_threshold = 500000;
	context->config.report_multiple_best_in_pairs = 0;
	context->config.realignment_minimum_variant_distance = 16;

	context->config.init_max_event_number = 70000;
	context->config.show_soft_cliping = 1 ;
	context->config.big_margin_record_size = 9;

	context->config.read_group_id[0] = 0;
	context->config.read_group_txt[0] = 0;
	context->config.first_read_file[0] = 0;
	context->config.second_read_file[0] = 0;
	context->config.index_prefix[0] = 0;
	context->config.output_prefix[0] = 0;
	context->config.exon_annotation_file[0] = 0;
	context->config.SAM_extra_columns = 0;
	context->config.exon_annotation_file_type = FILE_TYPE_GTF;
	strcpy(context->config.exon_annotation_gene_id_column, "gene_id");
	strcpy(context->config.exon_annotation_feature_name_column, "exon");

	context->config.DP_penalty_create_gap = -1;
	context->config.DP_penalty_extend_gap = 0;
	context->config.DP_match_score = 2;
	context->config.DP_mismatch_penalty = 0;

	context->config.max_insertion_at_junctions=0;
	context->config.check_donor_at_junctions=1;
	memset(&context -> input_reads, 0, sizeof(read_input_t));

	signal (SIGTERM, COREMAIN_SIGINT_hook);
	signal (SIGINT, COREMAIN_SIGINT_hook);


	int seed_rand[2];
	double double_time = miltime();
	memcpy(seed_rand, &double_time, 2*sizeof(int));
	myrand_srand(seed_rand[0]^seed_rand[1]); // the seed is NOT used in R because myrand_srand will always takes four random numbers from R's RNG

	context->config.max_indel_length = 5;
	context->config.phred_score_format = FASTQ_PHRED33;
	context->start_time = miltime();

	context->timecost_load_index = 0;
	context->timecost_voting = 0;
	context->timecost_before_realign = 0;
	context->timecost_for_realign = 0;
}

void init_core_temp_path(global_context_t * context){
	int x1;
	char mac_rand[13];
	mac_or_rand_str(mac_rand);

	context->config.temp_file_prefix[0] = 0;
	if(context->config.output_prefix[0]){
		for(x1 = strlen(context->config.output_prefix); x1 >=0; x1--){
			if(context->config.output_prefix[x1]=='/'){
				memcpy(context->config.temp_file_prefix, context->config.output_prefix, x1);
				context->config.temp_file_prefix[x1]=0;
				break;
			}
		}
	}

	if(context->config.temp_file_prefix[0] == 0)strcpy(context->config.temp_file_prefix, "./");
	sprintf(context->config.temp_file_prefix+strlen(context->config.temp_file_prefix), "/core-temp-sum-%06u-%s", getpid(), mac_rand );
	_COREMAIN_delete_temp_prefix = context->config.temp_file_prefix;
}

#define INDEL_MASK_BY_INSERTION 1
#define INDEL_MASK_BY_DELETION 2
#define INDEL_MASK_BY_MATCH 0
#define INDEL_MASK_BY_MISMATCH 3


	int CORE_DPALIGN_CREATEGAP_PENALTY = -1;
	int CORE_DPALIGN_EXTENDGAP_PENALTY = 0;
	int CORE_DPALIGN_MATCH_SCORE = 2;
	int CORE_DPALIGN_MISMATCH_PENALTY = 0;


int core_dynamic_align(global_context_t * global_context, thread_context_t * thread_context, char * read, int read_len, unsigned int begin_position, char * movement_buffer, int expected_offset, char * read_name)
// read must be converted to the positive strand.
// movement buffer: 0:match, 1: read-insert, 2: gene-insert, 3:mismatch
// the size of the movement buffer must be equal to the length of the read plus max_indel * 3.
{
	int max_indel = min(16 , global_context->config.max_indel_length);
	int i,j;

	CORE_DPALIGN_CREATEGAP_PENALTY = global_context -> config.DP_penalty_create_gap;
	CORE_DPALIGN_EXTENDGAP_PENALTY = global_context -> config.DP_penalty_extend_gap;
	CORE_DPALIGN_MATCH_SCORE = global_context -> config.DP_match_score;
	CORE_DPALIGN_MISMATCH_PENALTY = global_context -> config.DP_mismatch_penalty;

	if(read_len < 3 || abs(expected_offset) > max_indel)
		return 0;
	if(expected_offset < 0 && read_len < (3-expected_offset))
		return 0;

	gene_value_index_t * current_value_index = thread_context?thread_context->current_value_index:global_context->current_value_index; 

	indel_context_t * indel_context = (indel_context_t*)global_context -> module_contexts[MODULE_INDEL_ID];

	//unsigned long long table_ptr = (unsigned long long) indel_context -> dynamic_align_table;

	short ** table = indel_context -> dynamic_align_table;
	char ** table_mask = indel_context -> dynamic_align_table_mask;
	if(thread_context)
	{
		indel_thread_context_t * indel_thread_context = (indel_thread_context_t*)thread_context -> module_thread_contexts[MODULE_INDEL_ID];
		table = indel_thread_context -> dynamic_align_table;
		table_mask = indel_thread_context -> dynamic_align_table_mask;
	}


	if(0 && strcmp(read_name, "MISEQ:13:000000000-A1H1M:1:1112:12194:5511") == 0)
	{
		int ii;
		for(ii = 0; ii<read_len - expected_offset; ii++)
			SUBREADprintf("%c",gvindex_get(current_value_index, begin_position + ii));

		SUBREADprintf ("\n%s\n", read);
	}


	// vertical move: deletion (1)
	// horizontal move: insertion (2)
	// cross move: match (0) or mismatch (3)
	// i: vertical move; j: horizontal move

	//SUBREADprintf("DM[%d]: %p %d,%d\n", thread_context -> thread_id, table_mask, read_len +  expected_offset,  read_len);
	for (i=0; i<read_len +  expected_offset; i++)
	{
		for(j=0; j<read_len; j++)
		{
			if(0&&(i >= read_len +  expected_offset || j >= read_len))
				SUBREADprintf("XXDM[%d]: %p %d,%d\n", thread_context -> thread_id, table_mask, read_len +  expected_offset,  read_len);
			table_mask[i][j]=0;

			if(0&&(i >= read_len +  expected_offset || j >= read_len))
				SUBREADprintf("YYDM[%d]: %p %d,%d\n", thread_context -> thread_id, table_mask, read_len +  expected_offset,  read_len);

			if (j < i - max_indel || j > max_indel + i)
			{
				table[i][j]=-9999;
				if(0 && strcmp(read_name, "MISEQ:13:000000000-A1H1M:1:1112:12194:5511") == 0)
				{
					putchar('\t');
				}
				continue;
			}

			short from_upper;

			if (i>0) from_upper = table[i-1][j] + (table_mask[i-1][j] == INDEL_MASK_BY_DELETION?CORE_DPALIGN_EXTENDGAP_PENALTY:CORE_DPALIGN_CREATEGAP_PENALTY);
			else     from_upper = -9999;

			short from_left;

			if (j>0) from_left = table[i][j-1] + (table_mask[i][j-1] == INDEL_MASK_BY_INSERTION?CORE_DPALIGN_EXTENDGAP_PENALTY:CORE_DPALIGN_CREATEGAP_PENALTY);
			else     from_left = -9999;

			char chromo_ch = gvindex_get(current_value_index, begin_position + i);
			char is_matched_ij = (chromo_ch == read[j])?CORE_DPALIGN_MATCH_SCORE:CORE_DPALIGN_MISMATCH_PENALTY;
			
			short from_upperleft;

			if (i>0 && j>0) from_upperleft = table[i-1][j-1] + is_matched_ij;
			else if(i==0 && j==0) from_upperleft = is_matched_ij;
			else	    from_upperleft = -9999; 

			if (from_upperleft == from_upper && from_upperleft > from_left)
			{
				table_mask[i][j]= INDEL_MASK_BY_DELETION;
				table[i][j] = from_upper;
			}
			else if(from_upperleft == from_left && from_upperleft > from_upper)
			{
				table_mask[i][j]= INDEL_MASK_BY_INSERTION;
				table[i][j] = from_left;
			}
			else if(from_upperleft > from_left && from_upperleft > from_upper)
			{
				table_mask[i][j]= (chromo_ch == read[j])?INDEL_MASK_BY_MATCH:INDEL_MASK_BY_MISMATCH;
				table[i][j] = from_upperleft;
			}
			else if(from_upperleft == from_left && from_upperleft == from_upper)
			{
				table_mask[i][j]= (chromo_ch == read[j])?INDEL_MASK_BY_MATCH:INDEL_MASK_BY_MISMATCH;
				table[i][j] = from_upperleft;
			}
			else if(from_left > from_upper)
			{
				table_mask[i][j]= INDEL_MASK_BY_INSERTION;
				table[i][j] = from_left;
			}
			else if(from_left <= from_upper)
			{
				table_mask[i][j]= INDEL_MASK_BY_DELETION;
				table[i][j] = from_upper;
			}

			if(0 && strcmp(read_name, "MISEQ:13:000000000-A1H1M:1:1112:12194:5511") == 0)
				SUBREADprintf("%c%c\t", chromo_ch, read[j]);
		}
		if(0 && strcmp(read_name, "MISEQ:13:000000000-A1H1M:1:1112:12194:5511") == 0)
			SUBREADputs("");

	}
	#ifdef indel_debug
	//SUBREADputs("");
	//SUBREADputs("");

	#endif

	short path_i = read_len + expected_offset - 1;
	int out_pos = 0, delta=0;
	j = read_len - 1;

	if(0 && strcmp(read_name, "MISEQ:13:000000000-A1H1M:1:1112:12194:5511") == 0)
	{
		int ii,jj;
		for(ii=0;ii< path_i+1; ii++)
		{
			SUBREADprintf("%d\t", ii);
			for(jj=0; jj<j+1; jj++)
				SUBREADprintf("% 6d",table[ii][jj]);
			SUBREADputs("");
		}
		SUBREADprintf("  \t");
		for(jj=0; jj<j+1; jj++)
			SUBREADprintf("#%4d ",jj);
		SUBREADputs("");
		SUBREADputs("");

		for(ii=0;ii< path_i+1; ii++)
		{
			for(jj=0; jj<j+1; jj++)
				SUBREADprintf("% 6d",table_mask[ii][jj]);
			SUBREADputs("");
		}
	}

	while(1)
	{
		if(table_mask[path_i][j] == INDEL_MASK_BY_INSERTION)
		{
			j--;
			delta --;
			movement_buffer[out_pos++] = 2;
		}
		else if(table_mask[path_i][j] == INDEL_MASK_BY_DELETION)
		{
			path_i--;
			delta ++;
			movement_buffer[out_pos++] = 1;
		}
		else if(table_mask[path_i][j] == INDEL_MASK_BY_MATCH || table_mask[path_i][j] == INDEL_MASK_BY_MISMATCH)
		{
			movement_buffer[out_pos++] = table_mask[path_i][j] == INDEL_MASK_BY_MATCH?0:3;
			path_i--;
			j--;
		}

		if(path_i == -1 && j == -1) break;
		if(j<0 || path_i<0) return 0;
	}
	//out_pos++;
	//movement_buffer[out_pos]= 0;

	if(expected_offset!=delta)return 0;
	for(i=0; i<out_pos/2; i++)
	{
		char tmp;
		tmp = movement_buffer[out_pos-1-i];
		movement_buffer[out_pos-1-i] = movement_buffer[i];
		movement_buffer[i] = tmp;
	}

	if(0 && strcmp(read_name, "MISEQ:13:000000000-A1H1M:1:1112:12194:5511") == 0)
	{
		for(i=0; i<out_pos; i++)
		{
			char tmp = movement_buffer[i];
			switch(tmp){
				case 0: putchar('=');break; 
				case 1: putchar('D');break; 
				case 2: putchar('I');break; 
				case 3: putchar('X');break; 
			}
		}
		putchar('\n');
	}
	return out_pos;

}


