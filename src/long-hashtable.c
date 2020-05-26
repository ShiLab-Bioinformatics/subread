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
  
  
#include<assert.h>
#include <stdlib.h>
#include <string.h>

#ifndef MACOS
#ifndef FREEBSD
#include <malloc.h>
#endif
#endif

#include<math.h>
#include"core.h"
#include"long-hashtable.h"

int lnhash_create(lnhash_t * the_table, unsigned int num_buckets)
{
	int x1;
	the_table -> num_elements = 0;
	the_table -> num_buckets = num_buckets;
	the_table -> is_sorted = 0;
	the_table -> subread_repeat_max = 25600;

	the_table -> buckets = malloc(num_buckets * sizeof(lnhash_buckct_t));
	if(sizeof(void*) == 4)
		the_table -> key_repeated_numbers = malloc(sizeof(short) * 0xfffffff);
	else
		the_table -> key_repeated_numbers = malloc(sizeof(short) * 0x100000000llu);
	for(x1=0; x1<num_buckets; x1++)
	{
		the_table -> buckets[x1].num_elements = 0;
		the_table -> buckets[x1].space_elements = 0;
		the_table -> buckets[x1].key_array = NULL;
		the_table -> buckets[x1].value_array = NULL;
	}

	return 0;
}

void lnhash_destroy(lnhash_t * the_table)
{
	int x1;
	free(the_table -> key_repeated_numbers);
	for(x1=0; x1<the_table -> num_buckets; x1++)
	{
		if(the_table -> buckets[x1].key_array)
			free(the_table -> buckets[x1].key_array);
		if(the_table -> buckets[x1].value_array)
			free(the_table -> buckets[x1].value_array);
	}

	free(the_table -> buckets);
}

unsigned int lnhash_locate_bucket(lnhash_t * the_table, gehash_key_t key)
{
	return key % the_table -> num_buckets;
}

void lnhash_resize_bucket(lnhash_buckct_t * buck, unsigned int new_size)
{
	if(buck -> space_elements < new_size)
	{
		if(buck -> space_elements>0)
		{
			buck -> space_elements *= LNHASH_BUCKET_SIZE_INCREMENT;
			buck -> key_array = realloc(buck -> key_array, sizeof(gehash_key_t) * buck -> space_elements);
			buck -> value_array = realloc(buck -> value_array, sizeof(lnhash_data_t) * buck -> space_elements);
		}
		else
		{
			buck -> space_elements = LNHASH_INIT_BUCKET_SIZE;
			buck -> key_array = malloc(sizeof(gehash_key_t) * buck -> space_elements);
			buck -> value_array = malloc(sizeof(lnhash_data_t) * buck -> space_elements);
		}
	}
}

int lnhash_insert(lnhash_t * the_table, gehash_key_t key, lnhash_data_t data)
{
	int bucket_no = lnhash_locate_bucket(the_table, key);
	lnhash_buckct_t * buck = &(the_table -> buckets[bucket_no]);

	if(the_table -> key_repeated_numbers[key] < the_table -> subread_repeat_max)
	{
		the_table -> key_repeated_numbers[key]++;
		
		lnhash_resize_bucket(buck, buck -> num_elements +1);
		buck -> key_array[buck -> num_elements] = key;
		buck -> value_array[buck -> num_elements] = data;
		buck -> num_elements++;

		return 0;
	}
	else if(the_table -> key_repeated_numbers[key] == the_table -> subread_repeat_max)
	{
		return 1;

		int read_c, write_c=0;
		for(read_c =0; read_c < buck -> num_elements; read_c++)
		{
			if(buck -> key_array[read_c] != key)
			{
				if(write_c!=read_c)
				{
					buck -> key_array[write_c] = buck -> key_array[read_c];
					buck -> value_array[write_c] = buck -> value_array[read_c];
				}
				write_c++;
			}
		}

		buck -> num_elements = write_c;
		the_table -> key_repeated_numbers[key] = 65535;
	}

	return 1;
}

//void merge_sort(void * arr, int arr_size, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2));
//
int lnhash_mergesort_compare(void * arr, int l, int r)
{
	lnhash_buckct_t  * buck = arr;
	gehash_key_t lk = buck->key_array[l];
	gehash_key_t rk = buck->key_array[r];
	if(lk==rk)return 0;
	if(lk> rk)return 1;
	return -1;
}

void lnhash_mergesort_exchange(void * arr, int l, int r)
{
	lnhash_buckct_t  * buck = arr;
	gehash_key_t tmpk;
	lnhash_data_t tmpd;
	tmpk = buck->key_array[l];
	tmpd = buck->value_array[l];

	buck->value_array[l] = buck->value_array[r];
	buck->key_array[l] = buck->key_array[r];

	buck->key_array[r] = tmpk;
	buck->value_array[r] = tmpd;
}

void lnhash_mergesort_merge(void * arr, int start, int items, int items2)
{
	lnhash_buckct_t  * buck = arr;
	int i1_cursor = start, i2_cursor = start + items, out_cursor = 0;
	int end_point = start + items + items2;

	gehash_key_t * tmpk_arr = malloc(sizeof(gehash_key_t) * (items + items2));
	lnhash_data_t * tmpd_arr = malloc(sizeof(lnhash_data_t) * (items + items2));

	while(1)
	{
		if(i1_cursor == start + items && i2_cursor == end_point) break;
		assert(out_cursor <items + items2);

		int choose_2 = 0;

		if(i1_cursor >= start+items) 
			choose_2 = 1;

		if(i2_cursor < end_point)
			if(buck->key_array[i1_cursor] > buck->key_array[i2_cursor])
				choose_2 = 1;

		if(choose_2)
		{
			tmpk_arr[out_cursor] = buck->key_array[i2_cursor];
			tmpd_arr[out_cursor] = buck->value_array[i2_cursor];
			i2_cursor++;
		}
		else
		{
			tmpk_arr[out_cursor] = buck->key_array[i1_cursor];
			tmpd_arr[out_cursor] = buck->value_array[i1_cursor];
			i1_cursor++;
		}
		out_cursor++;
	}

	memcpy(buck->key_array   + start, tmpk_arr, sizeof(gehash_key_t ) * (items + items2));
	memcpy(buck->value_array + start, tmpd_arr, sizeof(lnhash_data_t) * (items + items2));

	free(tmpk_arr);
	free(tmpd_arr);
}

void lnhash_sort_bucket(lnhash_buckct_t * buck)
{
	merge_sort(buck, buck->num_elements, lnhash_mergesort_compare, lnhash_mergesort_exchange, lnhash_mergesort_merge);
}

void lnhash_resort(lnhash_t * the_table)
{
	int x1;
	for(x1 = 0;x1 < the_table -> num_buckets; x1++)
	{
		lnhash_buckct_t * buck = &(the_table -> buckets[x1]);
		lnhash_sort_bucket(buck);
	}

	the_table -> is_sorted = 1;
}

int lnhash_get_record_column(lnhash_data_t pos)
{
	return (int)(pos % LNHASH_VOTE_ARRAY_HEIGHT);
}

int lnhash_update_votes(lnhash_vote_t * votes, int record_column, lnhash_data_t test_pos, short dist_to_head)
{
	int row_cursor;

	for(row_cursor=0; row_cursor<votes->vote_record_items[record_column]; row_cursor++)
	{
		if(votes -> vote_records[record_column][row_cursor].head_position == test_pos)
		{
			votes -> vote_records[record_column][row_cursor].num_votes++;
			votes -> vote_records[record_column][row_cursor].coverage_start = min(
				votes -> vote_records[record_column][row_cursor].coverage_start, dist_to_head);
			votes -> vote_records[record_column][row_cursor].coverage_end = max(
				votes -> vote_records[record_column][row_cursor].coverage_end, dist_to_head + 16);
			return 1;
		}
	}
	return 0;
}

void lnhash_add_votes(lnhash_vote_t * votes, int record_column, lnhash_data_t test_pos, short dist_to_head)
{
	if(votes->vote_record_items[record_column]<LNHASH_VOTE_ARRAY_WIDTH)
	{
		int row_cursor = votes->vote_record_items[record_column];
		votes -> vote_records[record_column][row_cursor].head_position = test_pos;
		votes -> vote_records[record_column][row_cursor].num_votes = 1;
		votes -> vote_records[record_column][row_cursor].coverage_start = dist_to_head;
		votes -> vote_records[record_column][row_cursor].coverage_end = dist_to_head + 16;
		votes->vote_record_items[record_column]++;
	}
}

// no indels are considered.
int lnhash_query(lnhash_t * the_table, lnhash_vote_t * votes, gehash_key_t queried_key, short dist_to_head)
{
//	assert(the_table -> is_sorted);
	unsigned int bucket_no = lnhash_locate_bucket(the_table , queried_key);
	int lower_cursor, higher_cursor, match_cursor = 0, hits=0;

	lnhash_buckct_t * buck = &(the_table -> buckets[bucket_no]);

	if(the_table -> is_sorted)
	{
		lower_cursor = 0;
		higher_cursor = buck -> num_elements;
		if(higher_cursor<1) return 0;
		
		while(1)
		{
			if(lower_cursor> higher_cursor)
				return 0;
			int middle_no = (lower_cursor+higher_cursor)/2;

			if(buck->key_array[middle_no] == queried_key)
			{
				match_cursor = middle_no;
				break;
			}
			else if(buck->key_array[middle_no] > queried_key)
				higher_cursor = middle_no-1;
			else	lower_cursor = middle_no+1;
		}

		while( match_cursor > 0 )
		{
			if(buck->key_array[match_cursor-1] != queried_key)
				break;
			match_cursor--;
		}
	}

	for(; match_cursor < buck -> num_elements; match_cursor++)
	{
		if(buck->key_array[match_cursor]!=queried_key)
		{
			if(the_table -> is_sorted) break;
			else continue;
		}

		lnhash_data_t test_pos = buck -> value_array[match_cursor];

		test_pos -= dist_to_head;
		//if((test_pos & 0x3fff) > 1300) continue;

		int record_column = lnhash_get_record_column(test_pos);

		int found = lnhash_update_votes(votes, record_column, test_pos, dist_to_head);
		if(!found)  lnhash_add_votes   (votes, record_column, test_pos, dist_to_head);
		hits ++;
	}
	return hits;
}

int compare_voting_items(void * arrV, int i, int j)
{
	lnhash_vote_record_t * arr = (lnhash_vote_record_t*) arrV;
	if(arr[i].head_position > arr[j].head_position) return 1;
	else if(arr[i].head_position < arr[j].head_position) return -1;
	return 0;
}

void exchange_voting_items(void * arrV, int i, int j)
{
	lnhash_vote_record_t * arr = (lnhash_vote_record_t*) arrV;
	lnhash_vote_record_t tmp;
	memcpy(&tmp, arr + i, sizeof(lnhash_vote_record_t));
	memcpy(arr + i, arr + j, sizeof(lnhash_vote_record_t));
	memcpy(arr + j, &tmp, sizeof(lnhash_vote_record_t));
}

void merge_vorting_items(void * arrV, int start, int items1, int items2)
{
	lnhash_vote_record_t * arr = (lnhash_vote_record_t*) arrV;
	lnhash_vote_record_t * arrtmp = malloc( sizeof(lnhash_vote_record_t) * (items1 + items2));

	int cursor_1 = 0, cursor_2 = 0, xk1;

	for(xk1=0; xk1<items1+items2; xk1++)
	{
		lnhash_vote_record_t * curr_pnt1_item = arr + start + cursor_1; 
		lnhash_vote_record_t * curr_pnt2_item = arr + start + items1+ cursor_2; 
		int use_item1 = cursor_1 < items1;

		if(use_item1 && cursor_2 < items2)
			if(curr_pnt1_item -> head_position >= curr_pnt2_item -> head_position)
				use_item1 = 0;
		if(use_item1){
			memcpy(arrtmp+xk1, arr+start+cursor_1, sizeof(lnhash_vote_record_t));
			cursor_1++;
		}
		else{
			memcpy(arrtmp+xk1, arr+start+items1+cursor_2, sizeof(lnhash_vote_record_t));
			cursor_2++;
		}
	}

	memcpy(arr + start, arrtmp, sizeof(lnhash_vote_record_t) * (items1 + items2));
	free(arrtmp);
}


int sorted_voting_table_EX(lnhash_vote_t * vote, lnhash_vote_record_t ** result_buffer, int minimum_votes, int offsetted) 
{
	int ret_size = LNHASH_VOTE_ARRAY_HEIGHT * 3;
	int ret_items = 0;
	lnhash_vote_record_t * ret = malloc(sizeof(lnhash_vote_record_t) * ret_size);

	int vote_X, vote_Y;
	for(vote_X = 0; vote_X < LNHASH_VOTE_ARRAY_HEIGHT; vote_X++)
		for(vote_Y = 0; vote_Y < vote -> vote_record_items[vote_X]; vote_Y++)
		{
			lnhash_vote_record_t * rec = &(vote -> vote_records[vote_X][vote_Y]);
			if(rec -> num_votes >= minimum_votes)
			{
				if(ret_size <= ret_items)
				{
					ret_size *= 1.3;
					ret = realloc(ret, sizeof(lnhash_vote_record_t) * ret_size);
				}

				if(offsetted) rec -> head_position += rec -> coverage_start;
				memcpy(ret + ret_items, rec, sizeof(lnhash_vote_record_t));

				ret_items ++;
			}
		}


	merge_sort(ret, ret_items, compare_voting_items , exchange_voting_items, merge_vorting_items);

	(*result_buffer) = ret;
	return ret_items;

}

int sorted_voting_table_offset(lnhash_vote_t * vote, lnhash_vote_record_t ** result_buffer, int minimum_votes){
	return sorted_voting_table_EX(vote, result_buffer, minimum_votes, 1);
}
int sorted_voting_table(lnhash_vote_t * vote, lnhash_vote_record_t ** result_buffer, int minimum_votes) 
{
	return sorted_voting_table_EX(vote, result_buffer, minimum_votes, 0);
}

