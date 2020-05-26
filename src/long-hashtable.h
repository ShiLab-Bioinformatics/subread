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
  
  
#ifndef _LONG_HASHTABLE_H_
#define _LONG_HASHTABLE_H_
#include <stdlib.h>
#include <stdio.h>
#include "subread.h"

#define LNHASH_BUCKET_SIZE_INCREMENT 1.4
#define LNHASH_INIT_BUCKET_SIZE 300

typedef unsigned long long lnhash_data_t;

typedef struct {
	unsigned int num_elements;
	unsigned int space_elements;

	gehash_key_t * key_array;
	lnhash_data_t * value_array;
} lnhash_buckct_t;



typedef struct {
	int is_sorted;
	unsigned long long int num_elements;
	unsigned int num_buckets;
	unsigned short * key_repeated_numbers;
	unsigned short subread_repeat_max;

	lnhash_buckct_t * buckets;
} lnhash_t;

#define LNHASH_VOTE_ARRAY_WIDTH 240 
#define LNHASH_VOTE_ARRAY_HEIGHT 233
//#define LNHASH_VOTE_ARRAY_WIDTH 48 
//#define LNHASH_VOTE_ARRAY_HEIGHT 31

typedef struct{
	lnhash_data_t head_position;
	short coverage_start;
	short coverage_end;
	short num_votes;
} lnhash_vote_record_t;

typedef struct{
	int max_votes;
	int vote_record_items [LNHASH_VOTE_ARRAY_HEIGHT];
	lnhash_vote_record_t vote_records [LNHASH_VOTE_ARRAY_HEIGHT][LNHASH_VOTE_ARRAY_WIDTH];
} lnhash_vote_t; 

#define lnhash_init_votes(v) memset(v->vote_record_items, 0, sizeof(int) * LNHASH_VOTE_ARRAY_HEIGHT)  

// The EX version creates a hashtable with the given version number
int lnhash_create(lnhash_t * the_table, unsigned int num_buckets);

// This function puts a data item into the table. If there is duplication, it insert another copy into the table but do not overlap on the old one.
int lnhash_insert(lnhash_t * the_table, gehash_key_t key, lnhash_data_t data);

// this function must be called before querying or dumpping.
void lnhash_resort(lnhash_t * the_table);

// this function returns the number of items being hit
int lnhash_query(lnhash_t * the_table, lnhash_vote_t * vote, gehash_key_t queried_key, short dist_to_head);

// Free all memory that is allocated for the table. Only the table structure itself is not freed.
void lnhash_destroy(lnhash_t * the_table);

// This function returns the number of votes that has at least "minimum_votes" votes; the vote records are sorted by their head_position.
int sorted_voting_table(lnhash_vote_t * vote, lnhash_vote_record_t ** result_buffer, int minimum_votes);

// This function adds rec -> coverage_start to rec -> head_position.
int sorted_voting_table_offset(lnhash_vote_t * vote, lnhash_vote_record_t ** result_buffer, int minimum_votes);

#endif
