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
  
  
#ifndef _SORTED_HASHTABLE_H_
#define _SORTED_HASHTABLE_H_
#include <stdlib.h>
#include <stdio.h>

#include "subread.h"
#include "gene-algorithms.h"
#include "gene-value-index.h"

#define GEHASH_DEFAULT_SIZE	2000000000
#define GEHASH_BUCKET_LENGTH 100

#define gehash_fast_t gehash_t
#define gehash_destory_fast gehash_destory


// This function creates a new hash table. The invoter may provide the expected size of the table or -1 for a default size (2 billions)
// This function returns 0 if success or an errno
// The version of the hash table created by this function must be SUBINDEX_VER0.
int gehash_create(gehash_t * the_table, size_t expected_size, char is_small_table);

// The EX version creates a hashtable with the given version number
int gehash_create_ex(gehash_t * the_table, size_t expected_size, char is_small_table, int version_number, int GENE_SLIDING_STEP, int padding);

// This function puts a data item into the table. If there is duplication, it insert another copy into the table but do not overlap on the old one.
int gehash_insert(gehash_t * the_table, gehash_key_t key, gehash_data_t data, unsigned int * bucket_sizes);
void gehash_try_insert_measure(unsigned int * bucket_sizes, int bucket_no, gehash_key_t key);

// This function does what gehash_insert does, but insert nothing if the key has occured max_key_occurance times.
int gehash_insert_limited(gehash_t * the_table, gehash_key_t key, gehash_data_t data, int max_key_occurance, int prob_replace);

// This function queries the table and put the matched data item into data_result.
// This function returns 0 if not found, or the number of matched items.
// The invoter is in charge of allocating memory for results.
size_t gehash_get(gehash_t * the_table, gehash_key_t key, gehash_data_t * data_result, size_t max_result_space);

// Test existance, disregarding numbers.
// Return 1 if exist, 0 if not.
int gehash_exist(gehash_t * the_table, gehash_key_t key);


size_t gehash_go_q_CtoT(gehash_t * the_table, gehash_key_t key, int offset, int read_len, int is_reversed, gene_vote_t * vote,gene_vote_number_t weight, int max_match_number, int indel_tolerance, int subread_number,int max_error_bases, unsigned int low_border, unsigned int high_border);

size_t gehash_go_q_tolerable(gehash_t * the_table, gehash_key_t key, int offset, int read_len, int is_reversed, gene_vote_t * vote, gene_vote_number_t weight, gene_quality_score_t quality, int max_match_number, int indel_tolerance, int subread_number,int max_error_bases, int subread_len, unsigned int low_border, unsigned int high_border);

size_t gehash_go_q(gehash_t * the_table, gehash_key_t key, int offset, int read_len, int is_reversed, gene_vote_t * vote,int indel_tolerance, int subread_number, unsigned int low_border, unsigned int high_border);
size_t gehash_go_X(gehash_t * the_table, gehash_key_t key, int offset, int read_len, int is_reversed, gene_vote_t * vote,int indel_tolerance, int subread_number, unsigned int low_border, unsigned int high_border, int run_round, unsigned int * shift_indel_locations, unsigned int * shift_indel_NO);

// This function performs the same functionality, but runs only on AMD-64 cpus, and the length of each key must be 4 bytes.
size_t gehash_get_hpc(gehash_t * the_table, gehash_key_t key, gehash_data_t * data_result, size_t max_result_space);

// This function removes all items under the key. It returns the number of items that has been removed in this call.
size_t gehash_remove(gehash_t * the_table, gehash_key_t key);

// Free all memory that is allocated for the table. Only the table structure itself is not freed.
void gehash_destory(gehash_t * the_table);

// This function conpletely dumps a table into a disk file.
// It returns 0 if success, otherwise -1.
int gehash_dump(gehash_t * the_table, const char fname []);

void finalise_vote(gene_vote_t * vote);

// This function loads a dumpped hash table.
// The invoker does not need to initialise the table; it will be initialised in the function.
// It returns 0 if success, otherwise -1.
int gehash_load(gehash_t * the_table, const char fname []);

void gehash_prealloc(gehash_t * the_table);

size_t gehash_update(gehash_t * the_table, gehash_key_t key, gehash_data_t data_new);

short indel_recorder_copy(gene_vote_number_t *dst, gene_vote_number_t* src);

void assign_best_vote(gene_vote_t * vote, int i, int j);

void select_best_vote(gene_vote_t * vote);
void gehash_sort(gehash_t * the_table);
int gehash_load_option(const char fname [], int option_no, int * result);

// calculate # of buckets for estimaing their sizes.
unsigned int calculate_buckets_by_size(size_t exp_size, int version, int is_small_tab, int index_gap);
#endif
