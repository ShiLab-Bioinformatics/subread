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
  
  
#include "sorted-hashtable.h"
#include "gene-algorithms.h"
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

#define _gehash_hash(k) ((unsigned int)(k))


int gehash_create(gehash_t * the_table, size_t expected_size, char is_small_table) {
	return gehash_create_ex(the_table, expected_size, is_small_table, SUBINDEX_VER0, 3, 0);
}

unsigned int calculate_buckets_by_size( size_t expected_size, int version_number, int is_small_table, int index_gap ){
	if(expected_size ==0)
		expected_size = GEHASH_DEFAULT_SIZE;
	int expected_bucket_number = expected_size / GEHASH_BUCKET_LENGTH; 

	if(index_gap >= 3) expected_bucket_number /=3;

	if(SUBINDEX_VER1 > version_number) {
		if(expected_bucket_number < 10111 && !is_small_table)expected_bucket_number = 10111;
		else if(is_small_table)	expected_bucket_number = 4;
	} else 
		if(expected_bucket_number <= 0x3ffff )expected_bucket_number = 0x3ffff+4;


	for (;;expected_bucket_number++) {
		int j, valid_v;
		
		valid_v = 1;
		int test_prime = 13;
		if(SUBINDEX_VER1 > version_number && is_small_table) test_prime = 3;
		for(j=2; j<=test_prime;j++)
		{
			if (expected_bucket_number % j == 0)
				valid_v = 0;
		}

		if (valid_v)
			break;
	}

	//SUBREADprintf("EXP_BUCKS=%u\n",expected_bucket_number);

	return expected_bucket_number;
}

int gehash_create_ex(gehash_t * the_table, size_t expected_size, char is_small_table, int version_number, int index_gap, int padding) {
	int expected_bucket_number;
	int i;

	memset(the_table, 0, sizeof(gehash_t));
	if(expected_size ==0)
		expected_size = GEHASH_DEFAULT_SIZE;

	expected_bucket_number = calculate_buckets_by_size(expected_size, version_number, is_small_table, index_gap);

	// calculate the number of buckets for creating the data structure


	the_table -> version_number = version_number;
	the_table -> current_items = 0;
	the_table -> is_small_table = is_small_table;
	the_table -> buckets_number = expected_bucket_number;
	the_table -> buckets = (struct gehash_bucket *) malloc(
			  expected_bucket_number *
  			  sizeof(struct gehash_bucket )
			);
	the_table -> padding = padding;

	if(!the_table -> buckets)
	{
		SUBREADputs(MESSAGE_OUT_OF_MEMORY);
		return 1;
	}

	for(i=0; i<expected_bucket_number; i++){
		the_table -> buckets [i].item_keys = NULL;
		the_table -> buckets [i].current_items = 0;
		the_table -> buckets [i].space_size = 0;
	}

	the_table -> index_gap = index_gap;

	return 0;
}

int _gehash_resize_bucket(gehash_t * the_table , int bucket_no, char is_small_table)
{
	int new_bucket_length;
	struct gehash_bucket * current_bucket = &(the_table -> buckets [bucket_no]);
	gehash_key_t * new_item_keys = NULL;
	short * new_new_item_keys = NULL;
	gehash_data_t * new_item_values;

	if (current_bucket->space_size<1)
	{
		if(is_small_table)
			new_bucket_length = 5 ;
		else
			new_bucket_length = 2 ;

		if(the_table->version_number == SUBINDEX_VER0)
			new_item_keys = (gehash_key_t *) malloc(sizeof(gehash_key_t) * new_bucket_length);
		else
			new_new_item_keys = (short *) malloc(sizeof(short)*new_bucket_length);

		new_item_values = (gehash_data_t *) malloc(sizeof(gehash_data_t) * new_bucket_length);

		if(!((new_item_keys|| new_new_item_keys) && new_item_values))
		{
			SUBREADputs(MESSAGE_OUT_OF_MEMORY);
			return 1;
		}

		/*if(the_table->version_number == SUBINDEX_VER0)
			bzero(new_item_keys, sizeof(gehash_key_t) * new_bucket_length);
		else
			bzero(new_new_item_keys, sizeof(short) * new_bucket_length);

		bzero(new_item_values, sizeof(gehash_data_t) * new_bucket_length);
		*/


		if(the_table->version_number == SUBINDEX_VER0)
			current_bucket->item_keys = new_item_keys;
		else
			current_bucket->new_item_keys = new_new_item_keys;

		current_bucket->item_values = new_item_values;
		current_bucket->space_size = new_bucket_length;


	}
	else
	{
		if(is_small_table)
			//new_bucket_length = (int)max(15,current_bucket->space_size*1.5+1);
			new_bucket_length = (int)(current_bucket->space_size*5);
		else
			new_bucket_length = max((int)(current_bucket->space_size*1.5), (int)(current_bucket->space_size+4));	

		if(the_table->version_number == SUBINDEX_VER0)
			current_bucket->item_keys = (gehash_key_t *) realloc(current_bucket->item_keys ,  sizeof(gehash_key_t) * new_bucket_length);
		else
			current_bucket->new_item_keys = (short *) realloc(current_bucket->new_item_keys, sizeof(short)*new_bucket_length);

		 current_bucket->item_values = (gehash_data_t *) realloc(current_bucket->item_values, sizeof(gehash_data_t) * new_bucket_length);


		if(!((current_bucket->new_item_keys||current_bucket->new_item_keys) && current_bucket->item_values))
		{
			SUBREADputs(MESSAGE_OUT_OF_MEMORY);
			return 1;
		}

		current_bucket->space_size = new_bucket_length;
	}
	return 0;
}

#define _gehash_get_bucket(tab, key)  ( (tab) -> buckets + _gehash_hash( key ) % (tab) -> buckets_number )
#define _gehash_get_bucketNO(tab, key)  ( _gehash_hash( key ) % (tab) -> buckets_number )

void gehash_try_insert_measure(unsigned int * bucket_sizes, int bucket_no, gehash_key_t key){
	unsigned int buck_this = _gehash_hash(key)%bucket_no;
	//if(buck_this == 16305862)SUBREADprintf("Insert bucket : %u/%u -> %u at PTR=%p\n", buck_this, bucket_no, bucket_sizes[buck_this] , bucket_sizes);
	bucket_sizes[buck_this]++;
}

int gehash_insert(gehash_t * the_table, gehash_key_t key, gehash_data_t data, unsigned int * bucket_sizes)
{
	struct gehash_bucket * current_bucket;
	int is_fault = 0;

	//SUBREADprintf("Insert subread : %u/%u\n",   _gehash_hash(key)% the_table->buckets_number,  the_table->buckets_number);
	current_bucket = _gehash_get_bucket (the_table, key);

	if(the_table->version_number == SUBINDEX_VER0)
	{
		if (current_bucket->current_items >= current_bucket->space_size) {
			int bucket_number;
			if(bucket_sizes){
				SUBREADprintf("Bucket size was wrongly calculated.\n");
				return 1;
			}

			bucket_number = _gehash_hash(key) % the_table -> buckets_number;
			is_fault = _gehash_resize_bucket(the_table, bucket_number, the_table->is_small_table);
			if(is_fault)
				return 1;
		}
		current_bucket->item_keys[current_bucket->current_items] = key;
	}
	else
	{
		short step_number = key/the_table -> buckets_number;
		if(current_bucket-> new_item_keys == NULL && bucket_sizes){
			//this will go only once when all the bucktes are NULL.

			unsigned int * memory_sizes = malloc(sizeof(int) * the_table -> buckets_number);
			memset(memory_sizes, 0xff, sizeof(int) * the_table -> buckets_number);

			unsigned int grp_bytes = 0, grps = 0, this_grp_bucks = 0, xk1, bucks_per_grp = the_table -> buckets_number / GEHASH_MEM_PTR_NO +2;
			for(xk1=0; xk1<  the_table -> buckets_number ; xk1++){
				unsigned int buck_size =(sizeof(short) + sizeof(gehash_data_t)) * bucket_sizes[xk1];
				grp_bytes += buck_size;
				this_grp_bucks ++;
				if(  this_grp_bucks >= bucks_per_grp ){
					memory_sizes[grps++] = grp_bytes;
					grp_bytes=0;
					this_grp_bucks=0;
				}
			}
			if(this_grp_bucks)
				memory_sizes[grps] = grp_bytes;

			for(xk1=0; xk1 < GEHASH_MEM_PTR_NO; xk1++){
				unsigned int current_bytes = memory_sizes[xk1];
				if(current_bytes<0xff000000u){
					if(the_table -> malloc_ptr[xk1]!=NULL){
						SUBREADprintf("ERROR : no-NULL ptr : %p\n",the_table -> malloc_ptr[xk1]);
					}
					the_table -> malloc_ptr[xk1] = malloc(current_bytes);
					if(!the_table -> malloc_ptr[xk1]){
						SUBREADputs(MESSAGE_OUT_OF_MEMORY);
						return 1;
					}
				}
			}

			grp_bytes = 0; grps = 0; this_grp_bucks = 0;
			for(xk1=0; xk1<  the_table -> buckets_number ; xk1++){
				struct gehash_bucket * mem_bucket = the_table -> buckets+xk1;
				memset(mem_bucket, 0, sizeof(struct gehash_bucket));
				mem_bucket -> new_item_keys = (short*)((char*)the_table -> malloc_ptr[grps] + grp_bytes);
				mem_bucket -> item_values = (gehash_data_t*)((char *)mem_bucket -> new_item_keys + sizeof(short) * bucket_sizes[xk1]);
				mem_bucket -> space_size = bucket_sizes[xk1];

				unsigned int buck_size =(sizeof(short) + sizeof(gehash_data_t)) * bucket_sizes[xk1];
				grp_bytes += buck_size;
				this_grp_bucks ++;
				if(  this_grp_bucks >= bucks_per_grp ){
					grps++;
					grp_bytes=0;
					this_grp_bucks=0;
				}
			}

			the_table -> free_item_only = 2;
			free(memory_sizes);
		}

		if (current_bucket->current_items >= current_bucket->space_size) {
			if(bucket_sizes){
				SUBREADprintf("Bucket [%d] size was wrongly calculated : %d >= %u.\n",  _gehash_get_bucketNO (the_table, key), current_bucket->current_items,  current_bucket->space_size);
				return 1;
			}

			int bucket_number;
			bucket_number = _gehash_hash(key) % the_table -> buckets_number;

			is_fault = _gehash_resize_bucket(the_table, bucket_number, the_table->is_small_table);
			if(is_fault)
				return 1;
		}
		current_bucket->new_item_keys[current_bucket->current_items] = step_number;
	}

	current_bucket->item_values[current_bucket->current_items] = data;
	current_bucket->current_items ++;
	the_table ->current_items ++;
	return 0;
}

int gehash_insert_limited(gehash_t * the_table, gehash_key_t key, gehash_data_t data, int limit, int replace_prob)
{
	struct gehash_bucket * current_bucket;
	int occurance = 0, xk1 ;

	assert(the_table->version_number == SUBINDEX_VER0);
	current_bucket = _gehash_get_bucket (the_table, key);
	

	for(xk1=0; xk1<current_bucket->current_items ; xk1++)
		if(current_bucket->item_keys[xk1]==key) occurance++;

	if(occurance >= limit) 
	{

		if(myrand_rand()%32768 < replace_prob)
			return 1;
		int replacement_id = myrand_rand()%occurance;
		occurance = 0;
		for(xk1=0; xk1<current_bucket->current_items ; xk1++)
			if(current_bucket->item_keys[xk1]==key)
			{
				if(occurance == replacement_id)
				{
					current_bucket->item_values[xk1]=data;
					return 0;
				}
				occurance++;
			}

	}

	gehash_insert(the_table, key, data, NULL);
	return 0;
}
	



#define INDEL_SEGMENT_SIZE 5

#define _index_vote(key) (((unsigned int)(key))%GENE_VOTE_TABLE_SIZE)
#define _index_vote_tol(key) (((unsigned int)(key)/INDEL_SEGMENT_SIZE)%GENE_VOTE_TABLE_SIZE)


#define is_quality_subread(scr)	((scr)>15?1:0)

size_t gehash_go_q_tolerable(gehash_t * the_table, gehash_key_t key, int offset, int read_len, int is_reversed, gene_vote_t * vote,gene_vote_number_t weight, gene_quality_score_t quality ,int max_match_number, int indel_tolerance, int subread_number, int max_error_bases, int subread_len, unsigned int low_border, unsigned int high_border)
{
	if(max_error_bases >= 10) return 0;
	char error_pos_stack[10];	// max error bases = 10;

	gehash_key_t mutation_key;
	int ret = 0;

	ret+=gehash_go_q(the_table, key, offset, read_len, is_reversed, vote, indel_tolerance, subread_number, low_border, high_border);
	int error_bases ;
	for (error_bases=1; error_bases <= max_error_bases; error_bases++)
	{
		int i, j;
		for(i=0; i<error_bases; i++)
			error_pos_stack [i] = i;
		while (1)
		{

			char mutation_stack[10];
			memset(mutation_stack, 0 , error_bases);
			while(1)
			{
				int is_end2 = 1;
				for (i = error_bases-1; i>=0; i--)
					if (mutation_stack[i] <= 3)
					{
						int base_to_change=-1;
						mutation_key = key;
		
						for (j = 0; j < error_bases; j++)
						{
							base_to_change = error_pos_stack[j];
							int new_value = mutation_stack[j];
							mutation_key = mutation_key & ~(0x3 << (2*base_to_change));
							mutation_key = mutation_key | (new_value << (2*base_to_change));
						}
						for (j = i+1; j<error_bases; j++)
							mutation_stack[j] = 0;
						is_end2 = 0;
						if(key != mutation_key)
							ret+=gehash_go_q(the_table, mutation_key, offset, read_len, is_reversed, vote, indel_tolerance, subread_number, low_border, high_border);
						mutation_stack[i] ++;
						break;
					}
				if (is_end2)break;
			}


			int is_end = 1;
			for (i = error_bases-1; i>=0; i--)
				if (error_pos_stack[i] < subread_len - (error_bases - i))
				{
					error_pos_stack[i] ++;
					for (j = i+1; j<error_bases; j++)
						error_pos_stack[j] = error_pos_stack[i] + (j-i);
					is_end = 0;
					break;
				}

			if(is_end) break;
		}
	}
	return ret;
}


size_t gehash_go_q_CtoT(gehash_t * the_table, gehash_key_t key, int offset, int read_len, int is_reversed, gene_vote_t * vote, gene_vote_number_t weight,int max_match_number, int indel_tolerance, int subread_number, int max_error_bases, unsigned int low_border, unsigned int high_border)
{
	assert(max_error_bases < 5);
	int error_pos_stack[10];	// max error bases = 10;
	int C_poses[16];
	int  availableC=0;

	gehash_key_t mutation_key;
	int ret = 0,i;

	ret+=gehash_go_q(the_table, key, offset, read_len, is_reversed, vote, indel_tolerance, subread_number, low_border, high_border);
	int error_bases ;

	for(i=0; i<16; i++)
	{
		if(((key >> ((15-i)*2) )&3) == 3)	// if it is 'T' at i-th base
		{
			 C_poses[availableC++] = i;
		}
	}



	for (error_bases= min(availableC,max_error_bases); error_bases <= min(availableC,max_error_bases); error_bases++)
	{
		assert(error_bases < 5);
		int j;
		for(i=0; i<error_bases; i++)
		{
			error_pos_stack [i] = i;
		}

		while (1)
		{
			char mutation_stack[10];
			memset(mutation_stack, 0 , error_bases);
			while(1)
			{
				int is_end2 = 1;
				for (i = error_bases-1; i>=0; i--)
					if (mutation_stack[i] <1)
					{
						int base_to_change=-1;
						mutation_key = key;
		
						for (j = 0; j < error_bases; j++)
						{
							base_to_change = C_poses[error_pos_stack[j]];
							int new_value = mutation_stack[j]?2:3; //(C or T)
							mutation_key = mutation_key & ~(0x3 << (2*(15-base_to_change)));
							mutation_key = mutation_key | (new_value << (2*(15-base_to_change)));
						}
						for (j = i+1; j<error_bases; j++)
							mutation_stack[j] = 0;
						is_end2 = 0;

						if(key != mutation_key){
							int retone;
							retone = gehash_go_q(the_table, mutation_key, offset, read_len, is_reversed, vote, indel_tolerance, subread_number, low_border, high_border);
							ret+=retone;
						}
						mutation_stack[i] ++;
						break;
					}
				if (is_end2)break;
			}


			int is_end = 1;
			for (i = error_bases-1; i>=0; i--)
				if (error_pos_stack[i] < availableC - (error_bases - i))
				{
					error_pos_stack[i] ++;
					for (j = i+1; j<error_bases; j++)
						error_pos_stack[j] = error_pos_stack[i] + (j-i);
					is_end = 0;
					break;
				}

			if(is_end) break;
		}
	}
	return ret;
}


void assign_best_vote(gene_vote_t * vote, int i, int j)
{
	vote->max_mask = vote->masks[i][j];
	vote->max_vote = vote->votes[i][j];
	vote->max_position =vote->pos[i][j];
	vote->max_coverage_start = vote->coverage_start [i][j];
	vote->max_coverage_end = vote->coverage_end [i][j];
	memcpy(vote->max_indel_recorder, vote->indel_recorder[i][j], 3*MAX_INDEL_TOLERANCE * sizeof(*vote->max_indel_recorder));
}


size_t gehash_go_q(gehash_t * the_table, gehash_key_t raw_key, int offset, int read_len, int is_reversed, gene_vote_t * vote, int indel_tolerance, int subread_number, unsigned int low_border, unsigned int high_border)
{
	//SUBREADprintf("Q=%u, OFFSET=%d, B=%u ~ %u\n", raw_key, offset, low_border, high_border);

	if(the_table->version_number == SUBINDEX_VER0)
	{
		gehash_key_t key = raw_key;
		struct gehash_bucket * current_bucket;
		int i=0, items;

		gehash_key_t  *current_keys;//, *endp12;

		current_bucket = _gehash_get_bucket (the_table, key);
		items = current_bucket -> current_items;
		current_keys = current_bucket -> item_keys;
		
		if(!items) return 0;

		#define SPEED_UP_DENOMINAOR 3
		int jump_step = items / SPEED_UP_DENOMINAOR;
		int last_accepted_index = 0;

		if(jump_step<1) jump_step=1;

		if(key > current_keys[0])
		{
			while(1)
			{
				while(1)
				{
					int next_p = last_accepted_index + jump_step;
					if(next_p >= items) break;
					if(current_keys[next_p]>=key) break; 
					last_accepted_index = next_p;
				}
				if(jump_step>SPEED_UP_DENOMINAOR)
					jump_step /= SPEED_UP_DENOMINAOR;
				else if(jump_step>1)
					jump_step = 1;
				else
					break;
			}
			last_accepted_index++;
		}

		if(current_keys[last_accepted_index]!=key) return 0;

		short offset_from_5 = offset;
		//short offset_from_5 = is_reversed?(read_len - offset - 16):offset ; 

		if(*(current_bucket -> item_values+last_accepted_index) > 0xffff0000)	// no position should be greater than this.
		{
			// assumed to be non-informative subread.
			vote -> noninformative_subreads++;
			return 0;
		}

		if (indel_tolerance <1)
			for (; last_accepted_index<items && current_keys[last_accepted_index] == key ; last_accepted_index++)
			{
				unsigned int kv = current_bucket->item_values[last_accepted_index] - offset;
				int offsetX = _index_vote(kv);
				int datalen = vote -> items[offsetX];
				unsigned int * dat = vote -> pos[offsetX];

				if (kv > 0xffff0000)
					continue;
				for (i=0;i<datalen;i++)
				{
					if (dat[i] == kv && (subread_number + 1 > vote -> last_subread_cluster[offsetX][i]))
					{
						gene_vote_number_t test_max = (vote->votes[offsetX][i]);
						test_max += 1;
						vote->votes[offsetX][i] = test_max;
						if (offset_from_5 <  vote->coverage_start [offsetX][i])
							vote->coverage_start [offsetX][i] = offset_from_5;
						if (offset_from_5 +16 > vote->coverage_end [offsetX][i])
							vote->coverage_end [offsetX][i] = offset_from_5+16;

						vote -> last_subread_cluster[offsetX][i] = subread_number + 1;

						vote->max_vote = max(vote->max_vote , test_max);
						i = 9999999;
					}
				}

				if (i < 9999999 && datalen<GENE_VOTE_SPACE)
				{

					if (kv < low_border || kv > high_border)
						continue;


					vote -> items[offsetX] ++;
					dat[i] = kv;
					vote->votes[offsetX][i]=1;
					vote->masks[offsetX][i]= (is_reversed?IS_NEGATIVE_STRAND:0);
					vote->coverage_start [offsetX][i] = offset_from_5;
					vote->coverage_end [offsetX][i] = offset_from_5+16;
					vote -> last_subread_cluster[offsetX][i] = subread_number + 1;

					if(vote->max_vote==0)
						vote->max_vote = 1;
				}
			}
		else
		{
			// We duplicated all codes for indel_tolerance >= 1 for the minimal impact to performance.
			//int ii_end = (indel_tolerance % INDEL_SEGMENT_SIZE)?(indel_tolerance - indel_tolerance%INDEL_SEGMENT_SIZE+INDEL_SEGMENT_SIZE):indel_tolerance;

			for (; last_accepted_index<items && current_keys[last_accepted_index] == key ; last_accepted_index++)
			{
				unsigned int kv = current_bucket->item_values[last_accepted_index] - offset;
				int ii_end = (indel_tolerance % INDEL_SEGMENT_SIZE)?(indel_tolerance - indel_tolerance%INDEL_SEGMENT_SIZE+INDEL_SEGMENT_SIZE):indel_tolerance;
				int iix;
				i=0;

				for(iix = 0; iix<=ii_end; iix = iix>0?-iix:(-iix+INDEL_SEGMENT_SIZE))
				{
					int offsetX = _index_vote_tol(kv+iix);
					int datalen = vote -> items[offsetX];
					if(!datalen)continue;

					unsigned int * dat = vote -> pos[offsetX];

					for (i=0;i<datalen;i++)
					{
						int di = dat[i];
						int dist0 = kv-di;
						if( dist0 >= -indel_tolerance && dist0 <= indel_tolerance && is_reversed  == (0!=(vote -> masks[offsetX][i]&IS_NEGATIVE_STRAND))) {
							int toli =  vote -> toli[offsetX][i];
	
							if(di >= 46494104 && di <= 46496104 ){
								SUBREADprintf("VOTES: at %u, subread_no = %d , last_cluster = %d , toli = %d\n", di, subread_number , vote -> last_subread_cluster[offsetX][i] , toli);
							}

							if( toli > 0 && subread_number + 1 == vote -> last_subread_cluster[offsetX][i] ){
								int move_dist = 0;
								if( toli >=3 ) move_dist = vote -> indel_recorder[offsetX][i][toli-3+2];
								int new_dist = move_dist;
								move_dist -= vote -> indel_recorder[offsetX][i][toli+2];
								new_dist -= dist0;
								if(abs(move_dist) > abs(new_dist)){
									toli -= 3;
									vote -> toli[offsetX][i] = toli;
									vote -> last_subread_cluster[offsetX][i]--;
									vote -> votes[offsetX][i] --;
								}
							}
							if(subread_number + 1 <= vote -> last_subread_cluster[offsetX][i])
								continue;


							gene_vote_number_t test_max = (vote->votes[offsetX][i]);
							test_max += 1;
							vote -> votes[offsetX][i] = test_max;

							/*
							if (offset_from_5 <  vote->coverage_start [offsetX][i])
							{
								vote->coverage_start [offsetX][i] = offset_from_5;
							}*/
							if (offset_from_5 +16 > vote->coverage_end [offsetX][i])
							{
								vote->coverage_end [offsetX][i] = offset_from_5+16;
							}


							if (dist0 !=  vote->current_indel_cursor[offsetX][i])
							{
								toli +=3;
								if (toli < indel_tolerance*3)
								{
									vote -> toli[offsetX][i] = toli;
									vote -> indel_recorder[offsetX][i][toli] = subread_number+1; 
									vote -> indel_recorder[offsetX][i][toli+1] = subread_number+1;
									vote -> indel_recorder[offsetX][i][toli+2] = dist0; 
										
									if(toli < indel_tolerance*3-3) vote -> indel_recorder[offsetX][i][toli+3]=0;
								}
								vote->current_indel_cursor [offsetX][i] = (char)dist0;
							}
							else
								vote -> indel_recorder[offsetX][i][toli+1] = subread_number+1;

							vote -> last_subread_cluster[offsetX][i] = subread_number + 1;

							vote->max_vote = max(vote->max_vote , test_max);
							i = 9999999;
						}
					}
					if (i==9999999){
						break;
					}
				}

				if (i < 9999999)
				{
					if (kv < low_border || kv > high_border)
						continue;

					int offsetX2 = _index_vote_tol(kv);
					int datalen2 = vote -> items[offsetX2];
					unsigned int * dat2 = vote -> pos[offsetX2];

					if (datalen2<GENE_VOTE_SPACE)
					{
						vote -> items[offsetX2] ++;
						dat2[datalen2] = kv;
						vote -> masks[offsetX2][datalen2]=(is_reversed?IS_NEGATIVE_STRAND:0);
						vote -> votes[offsetX2][datalen2]=1;
						vote -> toli[offsetX2][datalen2]=0;
						vote -> last_subread_cluster[offsetX2][datalen2] = subread_number + 1;

						// data structure of recorder:
						// {unsigned char subread_start; unsigned char subread_end, char indel_offset_from_start}
						// All subread numbers are added with 1 for not being 0.

						vote -> indel_recorder[offsetX2][datalen2][0] = vote -> indel_recorder[offsetX2][datalen2][1] = subread_number+1;
						vote -> indel_recorder[offsetX2][datalen2][2] = 0;
						vote -> indel_recorder[offsetX2][datalen2][3] = 0;
						vote->current_indel_cursor [offsetX2][datalen2] = 0;
						vote->coverage_start [offsetX2][datalen2] = offset_from_5;
						vote->coverage_end [offsetX2][datalen2] = offset_from_5+16;

						if (vote->max_vote==0)
							vote->max_vote = 1;
					}
				}
			}
		}	
			
		return 1;
		//return match_end-match_start;
	}
	else
	{

		// VER_1
		// VER_2

		struct gehash_bucket * current_bucket;
		int i = 0, items;

		short *current_keys;//, *endp12;
		short key = raw_key / the_table->buckets_number;

		current_bucket = _gehash_get_bucket (the_table, raw_key);
		items = current_bucket -> current_items;
		current_keys = current_bucket -> new_item_keys;
		
		if(!items) return 0;

//#warning "======== MAKE SURE THAT '-1' IS CORRECT ============"
		int imin=0, imax=items - 1;
		int last_accepted_index;

		while(1)
		{
			last_accepted_index=(imin+imax)/2;
			short current_key = current_keys[last_accepted_index];
			if(current_key>key)
			{
				imax = last_accepted_index - 1;
			}
			else if(current_key<key)
			{
				imin = last_accepted_index + 1;
			}
			else
				break;

			if(imax<imin)
				return 0;
			
		}

		while(last_accepted_index){
			if(current_keys[last_accepted_index-1] == key) last_accepted_index-=1;
			else break;
		}


		/*if(*(current_bucket -> item_values+last_accepted_index) > 0xffff0000)	// no position should be greater than this.
		{
			// assumed to be non-informative subread.
			vote -> noninformative_subreads++;
			return 0;
		}*/

		int subread_number_P1 =  subread_number + 1;
		int of_p_16 = offset + 16;

		{
			int ii_end = INDEL_SEGMENT_SIZE;
			if(indel_tolerance>5) ii_end=(indel_tolerance % INDEL_SEGMENT_SIZE)?(indel_tolerance - indel_tolerance%INDEL_SEGMENT_SIZE+INDEL_SEGMENT_SIZE):indel_tolerance;

			for (;  last_accepted_index<items && current_keys[last_accepted_index] == key ; last_accepted_index++)
			{
				int iix, offsetX2, offsetX, datalen, datalen2;
				unsigned int kv = current_bucket->item_values[last_accepted_index] - offset;
				offsetX2 = _index_vote_tol(kv);
				datalen = vote -> items[offsetX2];

				datalen2 = datalen;
				offsetX = offsetX2;
				unsigned int * dat2, *dat;
				dat = dat2 = vote -> pos[offsetX2];

				//SUBREADprintf("You can find KV at %u\n", kv);

				for(iix = 0; iix<=ii_end; iix = iix>0?-iix:(-iix+INDEL_SEGMENT_SIZE))
				{
					if(iix)
					{
						offsetX = _index_vote_tol(kv+iix);
						datalen = vote -> items[offsetX];
						if(!datalen)continue;

						dat = vote -> pos[offsetX];
					}else if(!datalen)continue;

					for (i=0;i<datalen;i++)
					{
						int dist0 = kv-dat[i];

	
						if( dist0 >= -indel_tolerance && dist0 <= indel_tolerance && is_reversed  == (0!= vote->masks[offsetX][i]))
						{

//							if(di >= 46494104 && di <= 46496104 ){
//								SUBREADprintf("VOTES: at %u, subread_no = %d , last_cluster = %d , toli = %d\n", di, subread_number , vote -> last_subread_cluster[offsetX][i] , toli);
//							}

							int toli =  vote -> toli[offsetX][i];
							if( subread_number_P1 == vote -> last_subread_cluster[offsetX][i]  && toli >0){
								int move_dist = 0;
								if( toli >=3 ) move_dist = vote -> indel_recorder[offsetX][i][toli-3+2];
								int new_dist = move_dist;
								move_dist -= vote -> indel_recorder[offsetX][i][toli+2];
								new_dist -= dist0;
								if(abs(move_dist) > abs(new_dist)){
									toli -= 3;
									vote -> toli[offsetX][i] = toli;
									vote -> last_subread_cluster[offsetX][i]--;
									vote -> votes[offsetX][i] --;
								}
							}

							if(subread_number_P1 <= vote -> last_subread_cluster[offsetX][i]) continue;
							gene_vote_number_t test_max = (vote->votes[offsetX][i]);
							test_max += 1;
							vote -> votes[offsetX][i] = test_max;

							if (offset +16 > vote->coverage_end [offsetX][i])
								vote->coverage_end [offsetX][i] = of_p_16;


							if (dist0 ==  vote->current_indel_cursor[offsetX][i])
								vote -> indel_recorder[offsetX][i][toli+1] = subread_number_P1;
							else {
								toli +=3;
								if (toli < MAX_INDEL_SECTIONS*3)
								{
									vote -> toli[offsetX][i] = toli;
									vote -> indel_recorder[offsetX][i][toli] = subread_number_P1; 
									vote -> indel_recorder[offsetX][i][toli+1] = subread_number_P1;
									vote -> indel_recorder[offsetX][i][toli+2] = dist0; 
										
									if(toli < MAX_INDEL_SECTIONS*3-3) vote -> indel_recorder[offsetX][i][toli+3]=0;
								}
								vote->current_indel_cursor [offsetX][i] = (char)dist0;
							}

							vote -> last_subread_cluster[offsetX][i] = subread_number_P1;
							if(vote->max_vote < test_max)vote->max_vote = test_max;
							i = 9999999;
							break;
						}
					}
					if (i==9999999)break;

				}

				if (i < 9999999)
				{
					if (kv < low_border || kv > high_border)
						continue;

					if (datalen2<GENE_VOTE_SPACE)
					{
						vote -> items[offsetX2] ++;
						dat2[datalen2] = kv;
						vote -> masks[offsetX2][datalen2]=(is_reversed?IS_NEGATIVE_STRAND:0);
						vote -> votes[offsetX2][datalen2]=1;
						vote -> toli[offsetX2][datalen2]=0;

						// data structure of recorder:
						// {unsigned char subread_start; unsigned char subread_end, char indel_offset_from_start}
						// All subread numbers are added with 1 for not being 0.

						vote -> indel_recorder[offsetX2][datalen2][0] = vote -> indel_recorder[offsetX2][datalen2][1] = subread_number_P1;
						vote -> indel_recorder[offsetX2][datalen2][2] = 0;
						vote -> indel_recorder[offsetX2][datalen2][3] = 0;
						vote->current_indel_cursor [offsetX2][datalen2] = 0;
						vote->coverage_start [offsetX2][datalen2] = offset;
						vote->coverage_end [offsetX2][datalen2] = of_p_16;
						vote -> last_subread_cluster[offsetX2][datalen2] = subread_number_P1;

						if (vote->max_vote==0)
							vote->max_vote = 1;
					}
				}
				else i=0;
			}
		}	
		return 1;
	}
}



size_t gehash_go_X(gehash_t * the_table, gehash_key_t raw_key, int offset, int read_len, int is_reversed, gene_vote_t * vote, int indel_tolerance, int subread_number, unsigned int low_border, unsigned int high_border, int run_round, unsigned int * shift_indel_locs, unsigned int * shift_indel_NO){
	if(0)if(the_table->version_number == SUBINDEX_VER0){
		SUBREADprintf("ERROR: the version of the index is too old.\n");
		assert(the_table->version_number != SUBINDEX_VER0);
		return -1;
	} 

	// VER_1
	// VER_2

	struct gehash_bucket * current_bucket;
	int i = 0, items;

	short *current_keys;//, *endp12;
	short key = raw_key / the_table->buckets_number;

	current_bucket = _gehash_get_bucket (the_table, raw_key);
	items = current_bucket -> current_items;
	current_keys = current_bucket -> new_item_keys;
	
	if(!items) return 0;

//#warning "======== MAKE SURE THAT '-1' IS CORRECT ============"
	int imin=0, imax=items - 1;
	int last_accepted_index;

	while(1)
	{
		last_accepted_index=(imin+imax)/2;
		short current_key = current_keys[last_accepted_index];
		if(current_key>key)
		{
			imax = last_accepted_index - 1;
		}
		else if(current_key<key)
		{
			imin = last_accepted_index + 1;
		}
		else
			break;

		if(imax<imin)
			return 0;
		
	}


	int subread_number_P1 =  subread_number + 1;
	int of_p_16 = offset + 16;
	is_reversed = is_reversed?IS_NEGATIVE_STRAND:0;
	int start_scan_idx = last_accepted_index, scan_step = 0;

	{
		int ii_end = INDEL_SEGMENT_SIZE;
		if(indel_tolerance>5) ii_end=(indel_tolerance % INDEL_SEGMENT_SIZE)?(indel_tolerance - indel_tolerance%INDEL_SEGMENT_SIZE+INDEL_SEGMENT_SIZE):indel_tolerance;

		//for (; last_accepted_index<items && current_keys[last_accepted_index] == key ; last_accepted_index++){
		while(1){
			int iix, offsetX2, offsetX, datalen, datalen2;
			unsigned int kv = current_bucket->item_values[last_accepted_index] - offset;

			offsetX = offsetX2 = _index_vote_tol(kv);
			datalen = datalen2 = vote -> items[offsetX];

			unsigned int * dat2, *dat;
			dat = dat2 = vote -> pos[offsetX2];
			int found_some = 0;

			//SUBREADprintf("You can find KV at %u\n", kv);

			for(iix = 0; iix<=ii_end; iix = iix>0?-iix:(-iix+INDEL_SEGMENT_SIZE)) {
				if(iix) {
					offsetX = _index_vote_tol(kv+iix);
					datalen = vote -> items[offsetX];
					if(!datalen)continue;
					dat = vote -> pos[offsetX];
				} else if(!datalen) continue;

				for (i=0;i<datalen;i++)
				{
					int dist0 = kv-dat[i];
					int applied_indel_tol = ( run_round>0 && vote -> marked_shift_indel[offsetX][i])? 0: indel_tolerance ;
					if( dist0 >= -applied_indel_tol && dist0 <= applied_indel_tol && is_reversed == vote->masks[offsetX][i]){
						int toli =  vote -> toli[offsetX][i];

						if(run_round == 0 && toli >0 && dist0 ==0 && ! vote -> marked_shift_indel[offsetX][i]){
							vote -> marked_shift_indel[offsetX][i] = 1;
							shift_indel_locs[(* shift_indel_NO)++] = dat[i];
						}

						if( subread_number_P1 == vote -> last_subread_cluster[offsetX][i]  && toli >0){
							int move_dist = 0;
							if( toli >=3 ) move_dist = vote -> indel_recorder[offsetX][i][toli-3+2];
							int new_dist = move_dist;
							move_dist -= vote -> indel_recorder[offsetX][i][toli+2];
							new_dist -= dist0;
							if(abs(move_dist) > abs(new_dist)){
								toli -= 3;
								vote -> toli[offsetX][i] = toli;
								vote -> last_subread_cluster[offsetX][i]--;
								vote -> votes[offsetX][i] --;
							}
						}
	
						if(subread_number_P1 <= vote -> last_subread_cluster[offsetX][i]) continue;
						gene_vote_number_t test_max = (vote->votes[offsetX][i]);
						test_max += 1;
						vote -> votes[offsetX][i] = test_max;

						if (offset +16 > vote->coverage_end [offsetX][i])
							vote->coverage_end [offsetX][i] = of_p_16;
	
						if (dist0 ==  vote->current_indel_cursor[offsetX][i]){
							vote -> indel_recorder[offsetX][i][toli+1] = subread_number_P1;
						} else {
							toli +=3;
							if (toli < MAX_INDEL_SECTIONS*3)
							{
								vote -> toli[offsetX][i] = toli;
								vote -> indel_recorder[offsetX][i][toli] = subread_number_P1; 
								vote -> indel_recorder[offsetX][i][toli+1] = subread_number_P1;
								vote -> indel_recorder[offsetX][i][toli+2] = dist0; 
									
								if(toli < MAX_INDEL_SECTIONS*3-3) vote -> indel_recorder[offsetX][i][toli+3]=0;
							}
							vote->current_indel_cursor [offsetX][i] = (char)dist0;
						}

						vote -> last_subread_cluster[offsetX][i] = subread_number_P1;
						if(vote->max_vote < test_max)vote->max_vote = test_max;
						found_some = 1;
						break;
					}
				}
				if(found_some) break;
			}

			if (!found_some) {
				if (kv >=low_border && kv <= high_border && datalen2<GENE_VOTE_SPACE) {
					vote -> items[offsetX2] ++;
					dat2[datalen2] = kv;
					vote -> masks[offsetX2][datalen2] = is_reversed;
					vote -> votes[offsetX2][datalen2] = 1;
					vote -> toli[offsetX2][datalen2] = 0;
					vote -> marked_shift_indel[offsetX2][datalen2] = 0;
					if(run_round>0){
						int kk;
						for(kk = 0; kk < * shift_indel_NO ; kk++){
							if( kv >= shift_indel_locs[kk] - indel_tolerance && kv <= shift_indel_locs[kk] + indel_tolerance ){
								vote -> marked_shift_indel[offsetX2][datalen2] = 1;
								break;
							}
						}
					}
	
					// data structure of recorder:
					// {unsigned char subread_start; unsigned char subread_end, char indel_offset_from_start}
					// All subread numbers are added with 1 for not being 0.

					vote -> indel_recorder[offsetX2][datalen2][0] = vote -> indel_recorder[offsetX2][datalen2][1] = subread_number_P1;
					vote -> indel_recorder[offsetX2][datalen2][2] = 0;
					vote -> indel_recorder[offsetX2][datalen2][3] = 0;
					vote->current_indel_cursor [offsetX2][datalen2] = 0;
					vote->coverage_start [offsetX2][datalen2] = offset;
					vote->coverage_end [offsetX2][datalen2] = of_p_16;
					vote -> last_subread_cluster[offsetX2][datalen2] = subread_number_P1;
	
					if (vote->max_vote==0)
						vote->max_vote = 1;
				}
			}

			if(! scan_step){
				last_accepted_index++;
				if(last_accepted_index == items || current_keys[last_accepted_index]!= key){
					scan_step = 1;
					last_accepted_index = start_scan_idx;
				}
			}
			if(scan_step){
				last_accepted_index --;
				if(last_accepted_index <0 || current_keys[last_accepted_index]!= key) break;
			}
		}
	}	
	return 1;
}




void select_best_vote(gene_vote_t * vote)
{
	int i,j;
	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			if(vote -> votes[i][j] == vote->max_vote)
			{
				vote -> max_position = vote -> pos[i][j];
				vote -> max_mask = vote -> masks[i][j];
				vote -> max_coverage_start = vote->coverage_start[i][j];
				vote -> max_coverage_end = vote->coverage_end[i][j];
			}
		}
}

short indel_recorder_copy(gene_vote_number_t *dst, gene_vote_number_t * src)
{
	short all_indels = 0;
//	memcpy(dst, src, 3*MAX_INDEL_TOLERANCE);  return;


	int i=0;
	while(src[i] && (i<3*MAX_INDEL_TOLERANCE-2))
	{
		dst[i] = src[i];
		i++;
		dst[i] = src[i];
		i++;
		dst[i] = src[i];
		all_indels = dst[i]; 
		i++;
	}
	dst[i] = 0;
	return all_indels;

}

void finalise_vote(gene_vote_t * vote)
{
	if(vote->max_tmp_indel_recorder)
	{
		indel_recorder_copy(vote->max_indel_recorder, vote->max_tmp_indel_recorder);
	}
//		memcpy(vote->max_indel_recorder, vote->max_tmp_indel_recorder,  sizeof(char)*3 * MAX_INDEL_TOLERANCE);
}

int gehash_exist(gehash_t * the_table, gehash_key_t key)
{
	struct gehash_bucket * current_bucket;
	int  items;
	gehash_key_t * keyp, *endkp;

	assert(the_table->version_number == SUBINDEX_VER0);
	current_bucket = _gehash_get_bucket (the_table, key);
	items = current_bucket -> current_items;

	if(items <1)return 0;

	keyp = current_bucket -> item_keys;
	endkp = keyp + items;

	while(1)
	{
		if(*keyp == key)
			return 1;
		if(++keyp >= endkp) break;
	}
	return 0;
}


size_t gehash_update(gehash_t * the_table, gehash_key_t key, gehash_data_t data_new)
{
	struct gehash_bucket * current_bucket;
	size_t matched;
	int items;
	gehash_key_t * keyp, *endkp;

	assert(the_table->version_number == SUBINDEX_VER0);
	current_bucket = _gehash_get_bucket (the_table, key);

	matched = 0;
	items = current_bucket -> current_items;
	keyp = current_bucket -> item_keys;
	endkp = keyp + items;

	while(1)
	{
		if(*keyp == key)
		{
			current_bucket -> item_values[keyp-current_bucket -> item_keys] = data_new;
			matched++;
		}
		keyp +=1;
		if(keyp >= endkp)
			break;
	}
	return matched;
}

size_t gehash_get(gehash_t * the_table, gehash_key_t key, gehash_data_t * data_result, size_t max_result_space)
{
	struct gehash_bucket * current_bucket;
	size_t matched;
	int items;
	gehash_key_t * keyp, *endkp;

	if(max_result_space<1)
		return 0;

	assert(the_table->version_number == SUBINDEX_VER0);
	current_bucket = _gehash_get_bucket (the_table, key);

	matched = 0;
	items = current_bucket -> current_items;
	if(items<1) return 0;

	keyp = current_bucket -> item_keys;
	endkp = keyp + items;

	while(1)
	{
		if(*keyp == key)
		{
			data_result [matched] = current_bucket -> item_values[keyp-current_bucket -> item_keys];
			matched +=1;
			if(matched >= max_result_space)
				break;
		}
		keyp +=1;
		if(keyp >= endkp)
			break;
	}
	return matched;
}


size_t gehash_remove(gehash_t * the_table, gehash_key_t key)
{
	struct gehash_bucket * current_bucket;
	int i;
	size_t removed;

	assert(the_table->version_number == SUBINDEX_VER0);
	current_bucket = _gehash_get_bucket (the_table, key);	

	if(current_bucket -> current_items < 1)
		return 0;

	removed = 0;
	for(i=0; ; i++)
	{
		while(current_bucket -> item_keys [i+removed] == key && i+removed < current_bucket -> current_items)
			removed += 1;

		if(i+removed >= current_bucket -> current_items)
			break;

		if(removed)
		{
			current_bucket -> item_keys [i] = 
				current_bucket -> item_keys [i + removed];

			current_bucket -> item_values [i] = 
				current_bucket -> item_values [i + removed];
		}

	}

	current_bucket -> current_items -= removed;
	the_table-> current_items -= removed;

	return removed;
}

size_t gehash_get_hpc(gehash_t * the_table, gehash_key_t key, gehash_data_t * data_result, size_t max_result_space)
{
	return -1;
}


void gehash_insert_sorted(gehash_t * the_table, gehash_key_t key, gehash_data_t data)
{

	SUBREADprintf("UNIMPLEMENTED! gehash_insert_sorted \n");
}



// Data Struct of dumpping:
// {
//      size_t current_items;
//      size_t buckets_number;
//      struct 
//      {
//	      size_t current_items;
//	      size_t space_size;
//	      gehash_key_t item_keys [current_items];
//	      gehash_data_t item_values [current_items]
//      } [buckets_number];
// }
//

unsigned int load_int32(FILE * fp)
{
	int ret;
	int read_length;
	read_length = fread(&ret, sizeof(int), 1, fp);
	if(read_length<=0)assert(0);
	return ret;
}

srInt_64 load_int64(FILE * fp)
{
	srInt_64 ret;
	int read_length;
	read_length = fread(&ret, sizeof(srInt_64), 1, fp);
	if(read_length<=0)assert(0);
	return ret;
}


int gehash_load_option(const char fname [], int option_no, int * result){
	char tabname[MAX_FILE_NAME_LENGTH];
	char magic_chars[8];
	int found = 0, rrtv;
	sprintf(tabname, "%s.00.b.tab", fname);
	FILE * fp = f_subr_open(tabname, "rb");
	if(fp == NULL){
		sprintf(tabname, "%s.00.c.tab", fname);
		fp = f_subr_open(tabname, "rb");
	}
	if(fp){	
		rrtv = fread(magic_chars,1,8,fp);
		if(rrtv < 8) return -1;
		if(memcmp(magic_chars, "2subindx",7)==0) {
			while(1) {
				short option_key, option_length;

				rrtv = fread(&option_key, 2, 1, fp);
				if(rrtv < 1) return -1;
				if(!option_key) break;

				rrtv = fread(&option_length, 2, 1, fp);
				if(rrtv < 1) return -1;

				if(option_key == option_no){
					*result = 0;
					rrtv = fread(result ,option_length,1,fp);
					if(rrtv < 1) return -1;
					found = 1;
				}
				else
					fseek(fp, option_length, SEEK_CUR);
			}
		}
		fclose(fp);
	}
	return found;
}

int gehash_load(gehash_t * the_table, const char fname [])
{
	int i, read_length, rrtv;
	char magic_chars[8];
	magic_chars[7]=0;

	if(0)if(sizeof( size_t ) != 8|| sizeof(long int ) != 8|| sizeof(int ) != 4){
		SUBREADprintf("LINT: %zd , INT: %zd , SIZD : %zd\n", sizeof(long int ), sizeof(int ), sizeof(size_t));
		//return -1;
	}
	
	memset(the_table -> malloc_ptr ,0, sizeof(void*) * GEHASH_MEM_PTR_NO);
	the_table -> index_gap = 0;

	FILE * fp = f_subr_open(fname, "rb");
	if (!fp)
	{
		SUBREADprintf ("Table file '%s' is not found.\n", fname);
		return 1;
	}

	rrtv = fread(magic_chars,1,8,fp);
	if(rrtv !=8){
		SUBREADprintf("Error: the index magic string cannot be found. It may contain format errors or file '%s' may be truncated.\n", fname);
		assert(rrtv==8);
	}

	if(memcmp(magic_chars, "2subindx",8)!=0)
	{
		print_in_box(80,0,0,"");
		print_in_box(80,0,0,"WARNING your reference index was built under an old version of subread");
		print_in_box(80,0,0,"        package. Please rebuild the index for the reference genome.");
		print_in_box(80,0,0,"");
	}
	if(memcmp(magic_chars+1, "subindx",7)==0)
	{
		if('1'==magic_chars[0])
			the_table -> version_number = SUBINDEX_VER1;
		else if('2'==magic_chars[0])
			the_table -> version_number = SUBINDEX_VER2;
		else	assert(0);

		if(SUBINDEX_VER2 == the_table -> version_number)
		{
			while(1)
			{
				short option_key, option_length;

				rrtv = fread(&option_key, 2, 1, fp);
				if(rrtv <1){
					SUBREADprintf("Error: the index header cannot be fully loaded. It may contain format errors or file '%s' may be truncated.\n", fname);
					return -1;
				}
				//SUBREADprintf("READOPT_TAB = %04X\n", option_key);
				
				if(!option_key) break;

				rrtv = fread(&option_length, 2, 1, fp);
				if(rrtv <1){
					SUBREADprintf("Error: the index header cannot be fully loaded. It may contain format errors or file '%s' may be truncated.\n", fname);
					return -1;
				}

				rrtv = 999;
				if(option_key == SUBREAD_INDEX_OPTION_INDEX_GAP)
					rrtv = fread(&(the_table -> index_gap),2,1,fp);
				else if (option_key ==SUBREAD_INDEX_OPTION_INDEX_PADDING)
					rrtv = fread(&(the_table -> padding),2,1,fp);
				else
					fseek(fp, option_length, SEEK_CUR);
				if(rrtv <1){
					SUBREADprintf("Error: the index header cannot be fully loaded. It may contain format errors or file '%s' may be truncated.\n", fname);
					return -1;
				}
			}
			assert(the_table -> index_gap);
		}
		else if(SUBINDEX_VER1 == the_table -> version_number)
			the_table -> index_gap = 3;

		the_table -> current_items = load_int64(fp);
		if(the_table -> current_items < 1 || the_table -> current_items > 0xffffffffllu){
			SUBREADputs("ERROR: the index format is unrecognizable.");
			assert(the_table -> current_items >0 && the_table -> current_items <=0xffffffffllu);
		}

		unsigned int * bucket_bytes = malloc(sizeof(int) * GEHASH_MEM_PTR_NO);
		memset(bucket_bytes, 0xff, sizeof(int) * GEHASH_MEM_PTR_NO);
		the_table -> buckets_number = load_int32(fp);
		the_table -> buckets = (struct gehash_bucket * )malloc(sizeof(struct gehash_bucket) * the_table -> buckets_number);
		if(!the_table -> buckets)
		{
			SUBREADputs("Creating buckets.");
			SUBREADputs(MESSAGE_OUT_OF_MEMORY);
			return 1;
		}

		srInt_64 fp_curr = ftello(fp);
		unsigned int accued_bytes = 0, grp_i = 0, curr_bucks = 0, per_group_bucks = the_table -> buckets_number/ GEHASH_MEM_PTR_NO + 2;

		for(i=0; i<the_table -> buckets_number ; i++){
			unsigned int current_items = load_int32(fp);
			load_int32(fp);//useless for loading : space size
			unsigned int current_bytes = current_items*( sizeof(short) + sizeof(gehash_data_t) );
			accued_bytes += current_bytes; 
			curr_bucks ++;
			if(curr_bucks >= per_group_bucks){
				//SUBREADprintf("Allocating %d : %d buckets : %u\n", grp_i, i, accued_bytes);
				bucket_bytes[grp_i++] = accued_bytes;
				accued_bytes = 0;
				curr_bucks = 0;
			}
			#ifdef __MINGW32__
			{
//if(i%500000 == 0)SUBREADprintf("ESTMHUGE %d/%d\n", i, the_table -> buckets_number);

				char * buffkk = malloc(current_bytes );
				fread(buffkk,current_bytes ,1, fp);
				free(buffkk);
			}
			#else
			fseeko(fp, current_bytes, SEEK_CUR);
			#endif
		}
		if(curr_bucks)
			bucket_bytes[grp_i++] = accued_bytes;
		fseeko(fp, (off_t)fp_curr, SEEK_SET);

		for(i=0; i<GEHASH_MEM_PTR_NO ; i++){

//if(i%50000 == 0)SUBREADprintf("MEMHUGE %d/%d\n", i, the_table -> buckets_number);
			unsigned int current_bytes = bucket_bytes[i];
			if(current_bytes<0xff000000u){
				the_table -> malloc_ptr[i] = malloc(current_bytes);
				//SUBREADprintf("Allocating %d buckets : %u\n", i, current_bytes);
				if(!the_table -> malloc_ptr[i]){
					SUBREADputs(MESSAGE_OUT_OF_MEMORY);
					return 1;
				}
			}
		}

		grp_i = 0;
		curr_bucks = 0;
		accued_bytes = 0;
		for(i=0; i<the_table -> buckets_number ; i++){

//if(i%50000 == 0)SUBREADprintf("FILLHUGE %d/%d\n", i, the_table -> buckets_number);

			struct gehash_bucket * current_bucket = the_table -> buckets+i;
			current_bucket -> current_items = load_int32(fp);
			load_int32(fp);//useless for lo: space size
			current_bucket -> space_size =  current_bucket -> current_items;
			unsigned int current_bytes = current_bucket -> current_items*( sizeof(short) + sizeof(gehash_data_t) );

			current_bucket -> new_item_keys = (short*)(the_table -> malloc_ptr[grp_i] + accued_bytes);
			current_bucket -> item_values = (gehash_data_t*)(the_table -> malloc_ptr[grp_i] + accued_bytes + current_bucket -> current_items * sizeof(short));
			read_length = fread(current_bucket -> new_item_keys, sizeof(short), current_bucket -> current_items, fp);
			read_length += fread(current_bucket -> item_values, sizeof(gehash_data_t), current_bucket -> current_items, fp);
			if(read_length < 2* current_bucket -> current_items){
				SUBREADprintf("ERROR: the index is incomplete.\n");
				return 1;
			}
			
			accued_bytes += current_bytes;
			curr_bucks ++;
			if(curr_bucks >= per_group_bucks){
				curr_bucks=0;
				accued_bytes=0;
				grp_i++;
			}
		}

		int rval = fread(&(the_table -> is_small_table), sizeof(char), 1, fp);
		if (rval != 1)SUBREADprintf("ERROR: cannot find the table table.\n");
		free(bucket_bytes);
		fclose(fp);
		return 0;

	}
	else
	{
		fclose(fp);
		the_table -> index_gap = 3;
		fp = f_subr_open(fname, "rb");
		the_table -> version_number = SUBINDEX_VER0;
		the_table -> current_items = load_int64(fp);
		the_table -> buckets_number = load_int32(fp);
		the_table -> buckets = (struct gehash_bucket * )malloc(sizeof(struct gehash_bucket) * the_table -> buckets_number);
		if(!the_table -> buckets)
		{
			SUBREADputs(MESSAGE_OUT_OF_MEMORY);
			return 1;
		}

		for (i=0; i<the_table -> buckets_number; i++)
		{
			struct gehash_bucket * current_bucket = &(the_table -> buckets[i]);
			current_bucket -> current_items = load_int32(fp);
			current_bucket -> space_size = load_int32(fp);
			current_bucket -> space_size = current_bucket -> current_items;
			current_bucket -> item_keys = (gehash_key_t *) malloc ( sizeof(gehash_key_t) * current_bucket -> space_size);
			current_bucket -> item_values = (gehash_data_t *) malloc ( sizeof(gehash_data_t) * current_bucket -> space_size);

			if(!(current_bucket -> item_keys&&current_bucket -> item_values))
			{
				SUBREADputs(MESSAGE_OUT_OF_MEMORY);
				return 1;

			}

			if(current_bucket -> current_items > 0)
			{
				read_length = fread(current_bucket -> item_keys, sizeof(gehash_key_t), current_bucket -> current_items, fp);
				if(read_length < current_bucket -> current_items){
					SUBREADprintf("ERROR: the index is incomplete.\n");
					return 1;
				}
				read_length = fread(current_bucket -> item_values, sizeof(gehash_data_t), current_bucket -> current_items, fp);
				if(read_length < current_bucket -> current_items){
					SUBREADprintf("ERROR: the index is incomplete.\n");
					return 1;
				}
			}

		}

		read_length = fread(&(the_table -> is_small_table), sizeof(char), 1, fp);
		if(read_length!=1){
			SUBREADprintf("ERROR: the index is incomplete.\n");
			return 1;
		}
		fclose(fp);
		return 0;
	}
}


void gehash_sort(gehash_t * the_table)
{
	int i;
	for (i=0; i<the_table -> buckets_number; i++)
	{
		struct gehash_bucket * current_bucket = &(the_table -> buckets[i]);
		int ii, jj;
		gehash_key_t tmp_key;
		gehash_data_t tmp_val;

		if(current_bucket -> current_items>=1)
		{
			for(ii=0;ii<current_bucket -> current_items -1; ii++)
			{
				for (jj = ii+1; jj < current_bucket -> current_items; jj++)
				{
					if (current_bucket -> item_keys[ii] > current_bucket -> item_keys[jj])
					{
						tmp_key = current_bucket -> item_keys[ii];
						current_bucket -> item_keys[ii] = current_bucket -> item_keys[jj];
						current_bucket -> item_keys[jj] = tmp_key;

						tmp_val = current_bucket -> item_values[ii];
						current_bucket -> item_values[ii] = current_bucket -> item_values[jj];
						current_bucket -> item_values[jj] = tmp_val;
					}
				}
			}
		}
	}
}

void write_cell(short option_key, short option_length, char * option_value, FILE * fp)
{
	fwrite(&option_key, 2, 1, fp);
	fwrite(&option_length, 2, 1, fp);
	fwrite(option_value, option_length, 1, fp);
}


void write_options(FILE * fp, gehash_t * the_table)
{
	short option_key, option_length;
	short option_value;

	option_key = SUBREAD_INDEX_OPTION_INDEX_PADDING;
	option_length = 2;
	option_value = the_table -> padding;

	write_cell(option_key, option_length, (char *) &option_value,fp);

	option_key = SUBREAD_INDEX_OPTION_INDEX_GAP;
	option_length = 2;
	option_value = the_table -> index_gap;

	write_cell(option_key, option_length, (char *) &option_value,fp);

	option_key = 0;
	fwrite(&option_key, 2, 1, fp);
}

int is_1_greater_than_2(int bucket_i, int all_buckets , short k1, gehash_data_t v1, short k2, gehash_data_t v2){
	if(k1 > k2) return 1;
	if(k1 == k2){
		unsigned int real_key = 1u*k1*1u*all_buckets;
		real_key += bucket_i;
		if((real_key%791)%2==0 && v1>v2) return 1;
		if((real_key%791)%2==1 && v1<v2) return 1;
	}
	return 0;
}

int gehash_dump(gehash_t * the_table, const char fname [])
{
	int ii, jj, xx;
	int i, scroll_counter = 0;
	FILE * fp = f_subr_open(fname, "wb");
	int maximum_bucket_size = 0;
	if (!fp)
	{
		SUBREADprintf ("Table file '%s' is not able to open.\n", fname);
		return -1;
	}

	if(the_table->version_number == SUBINDEX_VER2)
	{
		fwrite("2subindx",1,8,fp);
		write_options(fp, the_table);
	}

	assert(sizeof(srInt_64 ) == 8);
	assert(sizeof(int ) == 4);
	fwrite(& (the_table -> current_items ), sizeof(srInt_64), 1, fp);
	fwrite(& (the_table -> buckets_number), sizeof(int), 1, fp);

	print_in_box(80,0,0,"Save current index block...              ");

	for (i=0; i<the_table -> buckets_number; i++)
	{
		struct gehash_bucket * current_bucket = &(the_table -> buckets[i]);
		if(current_bucket -> current_items > maximum_bucket_size)
			maximum_bucket_size = current_bucket -> current_items;
	}

	#define SORT_LANE_NUMBER 16
	short * sort_space_new_key[SORT_LANE_NUMBER];
	gehash_data_t * sort_space_data [SORT_LANE_NUMBER];
	int items_in_sort[SORT_LANE_NUMBER] ;
	int items_in_merge[SORT_LANE_NUMBER] ;
	if(the_table->version_number > SUBINDEX_VER0)
	{
		for(xx=0;xx<SORT_LANE_NUMBER;xx++)
		{
			sort_space_new_key[xx] = (short *)malloc(sizeof(short)*(maximum_bucket_size/SORT_LANE_NUMBER+2));
			sort_space_data[xx] = (gehash_data_t *)malloc(sizeof(gehash_data_t)*(maximum_bucket_size/SORT_LANE_NUMBER+2));
		}
	}

	for (i=0; i<the_table -> buckets_number; i++)
	{
		struct gehash_bucket * current_bucket = &(the_table -> buckets[i]);
		gehash_data_t tmp_val=0;

		if(i % (the_table -> buckets_number/10) == 0)
			print_window_scrolling_bar("", 1.0*i/the_table -> buckets_number, 74, &scroll_counter);

		if(current_bucket -> current_items>=1)
		{
			if(the_table->version_number == SUBINDEX_VER0)
			{
				for(ii=0;ii<current_bucket -> current_items -1; ii++)
				{
					gehash_key_t tmp_key;
					for (jj = ii+1; jj < current_bucket -> current_items; jj++)
					{
						
						if(is_1_greater_than_2(i, the_table -> buckets_number, current_bucket -> item_keys[ii], current_bucket -> item_values[ii],
											    current_bucket -> item_keys[jj], current_bucket -> item_values[jj] )) {
							tmp_key = current_bucket -> item_keys[ii];
							current_bucket -> item_keys[ii] = current_bucket -> item_keys[jj];
							current_bucket -> item_keys[jj] = tmp_key;

							tmp_val = current_bucket -> item_values[ii];
							current_bucket -> item_values[ii] = current_bucket -> item_values[jj];
							current_bucket -> item_values[jj] = tmp_val;
						}
					}
				}
			}
			else
			{
				if(1){
					memset(items_in_sort, 0, sizeof(int)*SORT_LANE_NUMBER);
					for(ii=0;ii<current_bucket -> current_items;ii++)
					{
						int sort_lane_no = ii%SORT_LANE_NUMBER;
						int xx_th_item = items_in_sort[sort_lane_no] ++;
						sort_space_new_key[sort_lane_no][xx_th_item] = current_bucket -> new_item_keys[ii];
						sort_space_data[sort_lane_no][xx_th_item] = current_bucket -> item_values[ii];
					}

					// sort the individual lanes.
					for(xx=0;xx<SORT_LANE_NUMBER;xx++)
					{
						for(ii=0;ii<items_in_sort[xx]-1; ii++)
						{
							for(jj = ii+1; jj < items_in_sort[xx]; jj++)
							{
								short tmp_key;
								if(is_1_greater_than_2(i, the_table -> buckets_number,   sort_space_new_key[xx][ii], sort_space_data[xx][ii], sort_space_new_key[xx][jj], sort_space_data[xx][jj] )) {
									tmp_key = sort_space_new_key[xx][ii];
									sort_space_new_key[xx][ii] = sort_space_new_key[xx][jj];
									sort_space_new_key[xx][jj] = tmp_key;

									tmp_val = sort_space_data[xx][ii];
									sort_space_data[xx][ii]=sort_space_data[xx][jj];
									sort_space_data[xx][jj]=tmp_val;
								}
							}
						}
					}

					// items_in_sort : array, each element for a lane
					// items_in_merge : how many items in each lane that have been merged (ie, read_ptr for each lane)
					// sort_space_new_key : N lanes each has been sorted, each has many keys (ie subreads)
					// sort_space_data : N lanes, each has many values (ie positions)
					// current_bucket -> new_item_keys : One array containing sorted keys (write_ptr is ii)
					// current_bucket -> item_values:  One array containing sorted values (write_ptr is ii)) 
					memset(items_in_merge, 0, sizeof(int)*SORT_LANE_NUMBER);

					// tmp_key : the current key being examed. A larger key is put AFTER a smaller key.
					for(ii=0;ii<current_bucket -> current_items;ii++)
					{
						int tmp_key=0x7fffffff;
						int selected_lane = 0, trying_lane;
						for(trying_lane=0;trying_lane<SORT_LANE_NUMBER;trying_lane++)
						{
							int ii_in_trying_lane = items_in_merge[trying_lane];
							if(ii_in_trying_lane >= items_in_sort[trying_lane]) continue;

							if(tmp_key >= 0x10000 ||is_1_greater_than_2(
									i, the_table -> buckets_number, 
									sort_space_new_key[selected_lane][items_in_merge[selected_lane]], 
									sort_space_data[selected_lane][items_in_merge[selected_lane]],
									sort_space_new_key[trying_lane][items_in_merge[trying_lane]],
									sort_space_data[trying_lane][items_in_merge[trying_lane]]
								)){
								selected_lane=trying_lane;
								tmp_key = sort_space_new_key[trying_lane][ii_in_trying_lane];
							}
						}

						current_bucket -> new_item_keys[ii] = sort_space_new_key[selected_lane][items_in_merge[selected_lane]];
						current_bucket -> item_values[ii] = sort_space_data[selected_lane][items_in_merge[selected_lane]];

						items_in_merge[selected_lane]++;
						//printf("V [%d] [%d] =%d\n",i , ii, tmp_key);
					}
				}
			}
		}


		if(0) // debug
		{
			int inidx=0;
			srInt_64 real_old = -1;
			for(xx = 0; xx < current_bucket -> current_items; xx++){
				unsigned int real_key = 1u * current_bucket -> new_item_keys[xx] *1u * the_table -> buckets_number;
				real_key += i;

				if(real_old != real_key) inidx=0;

				if(1||real_key >= 0xffffff00 || real_key <5000)
				SUBREADprintf("SUBREAD_%010u POS_%010u INIDX_%04d BUCKET_%04d_ITEM_%08d\n",real_key ,current_bucket -> item_values[xx], inidx, i, xx);
				real_old = real_key;
				inidx++;
			}
		}

		int is_full = 0;
		int write_len = fwrite(& (current_bucket -> current_items), sizeof(int), 1, fp);
		if(write_len<1) is_full = 1;
		write_len = fwrite(& (current_bucket -> space_size), sizeof(int), 1, fp);
		if(write_len<1) is_full = 1;

		if(the_table->version_number == SUBINDEX_VER0)
			fwrite(current_bucket -> item_keys, sizeof(gehash_key_t), current_bucket -> current_items, fp);
		else{
			write_len = fwrite(current_bucket -> new_item_keys, sizeof(short), current_bucket -> current_items, fp);
			if(write_len < current_bucket -> current_items) is_full = 1;
		}
		write_len = fwrite(current_bucket -> item_values, sizeof(gehash_data_t), current_bucket -> current_items, fp);
		if(write_len < current_bucket -> current_items) is_full = 1;
		if(is_full){
			fclose(fp);
			SUBREADprintf("ERROR: Unable to write into the output file. Please check the disk space in the output directory.\n");
			return 1;
		}
	}

	if(the_table->version_number > SUBINDEX_VER0)
	{
		for(xx=0;xx<SORT_LANE_NUMBER;xx++)
		{
			free(sort_space_new_key[xx]);
			free(sort_space_data[xx]);
		}
	}


	int write_len = fwrite(&(the_table -> is_small_table), sizeof(char), 1, fp);
	fclose(fp);

	if(write_len < 1){
		SUBREADprintf("ERROR: Unable to write into the output file. Please check the disk space in the output directory.\n");
		return 1;
	}
	print_in_box(80,0,0,"");
	return 0;
}

void gehash_destory(gehash_t * the_table)
{
	int i, is_ptr = 0;

	for(i = 0; i < GEHASH_MEM_PTR_NO; i++) if(the_table -> malloc_ptr[i]){
		free(the_table -> malloc_ptr[i]);
		is_ptr=1;
	}

	if(!is_ptr) for (i=0; i<the_table -> buckets_number; i++)
	{
		struct gehash_bucket * current_bucket = &(the_table -> buckets[i]);
		if (current_bucket -> space_size > 0)
		{
			if(0==the_table->free_item_only)free (current_bucket -> item_keys);
			free (current_bucket -> item_values);
		}
	}

	free (the_table -> buckets);

	the_table -> current_items = 0;
	the_table -> buckets_number = 0;
}
