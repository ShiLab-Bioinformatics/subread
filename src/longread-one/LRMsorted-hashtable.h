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
#include "LRMconfig.h"

#define LRMGEHASH_DEFAULT_SIZE	2000000000
#define LRMGEHASH_BUCKET_LENGTH	2291

#define LRMinit_gene_vote(a) {memset((a)->items, 0, LRMGENE_VOTE_TABLE_SIZE*sizeof( *((a)->items)));  }

size_t LRMgehash_go_q(LRMgehash_t * the_table, LRMgehash_key_t key, int offset, int read_len, int is_reversed, LRMgene_vote_t * vote, int indel_tolerance, int subread_number);
size_t LRMgehash_go_tolerance(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, LRMgehash_t * the_table, LRMgehash_key_t key, int offset, int read_len, int is_reversed, LRMgene_vote_t * vote, int indel_tolerance, int subread_number, int max_MM);

void LRMgehash_destory(LRMgehash_t * the_table);
void LRMfinalise_vote(LRMgene_vote_t * vote);
int LRMgehash_load(LRMgehash_t * the_table, const char fname []);
void LRMassign_best_vote(LRMgene_vote_t * vote, int i, int j);

void LRMselect_best_vote(LRMgene_vote_t * vote);
void LRMgehash_sort(LRMgehash_t * the_table);
int LRMgehash_load_option(const char fname [], int option_no, int * result);
void LRMprint_v(LRMcontext_t * context, LRMread_iteration_context_t * iteration_context, int min_votes);

size_t LRMgehash_go_QQ(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, LRMgehash_t * the_table, LRMgehash_key_t key, int offset, int read_len, int is_reversed, LRMgene_vote_t * vote, int indel_tolerance, int subread_number);
#endif
