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
  
  
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "LRMsorted-hashtable.h"

#ifndef MACOS
#ifndef FREEBSD
#include <malloc.h>
#endif
#endif

#include<math.h>
#include "LRMfile-io.h"

#define _gehash_hash(k) ((unsigned int)(k))
#define WITHOUT_CLUSTER_ORDERING 0

struct LRMgehash_bucket * LRM_gehash_get_bucket(LRMgehash_t * the_table, LRMgehash_key_t key) {
	int bucket_number;

	bucket_number = _gehash_hash(key) % the_table -> buckets_number;
	return  the_table -> buckets +bucket_number;
}



#define INDEL_SEGMENT_SIZE 5

#define _index_vote(key) (((unsigned int)(key))%LRMGENE_VOTE_TABLE_SIZE)
#define _index_vote_tol(key) (((unsigned int)(key)/INDEL_SEGMENT_SIZE)%LRMGENE_VOTE_TABLE_SIZE)


#define is_quality_subread(scr)	((scr)>15?1:0)


size_t LRMgehash_go_q(LRMgehash_t * the_table, LRMgehash_key_t raw_key, int offset, int read_len, int is_reversed, LRMgene_vote_t * vote, int indel_tolerance, int subread_number){
	//LRMprintf("Q=%u, OFFSET=%d, B=%u ~ %u\n", raw_key, offset, low_border, high_border);

	// VER_1
	// VER_2

	struct LRMgehash_bucket * current_bucket;
	int i = 0, items;

	short *current_keys;//, *endp12;
	short key = raw_key / the_table->buckets_number;

	current_bucket = LRM_gehash_get_bucket (the_table, raw_key);
	items = current_bucket -> current_items;
	current_keys = current_bucket -> new_item_keys;
	
	if(!items) return 0;

	int imin=0, imax=items;
	int last_accepted_index = 0;

	while( imin < items )
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

	int ii_end = INDEL_SEGMENT_SIZE;
	if(indel_tolerance>5) ii_end=(indel_tolerance % INDEL_SEGMENT_SIZE)?(indel_tolerance - indel_tolerance%INDEL_SEGMENT_SIZE+INDEL_SEGMENT_SIZE):indel_tolerance;

	for (; last_accepted_index<items && current_keys[last_accepted_index] == key ; last_accepted_index++)
	{
		unsigned int kv = current_bucket->item_values[last_accepted_index] - offset;
		int iix, offsetX2, offsetX, datalen, datalen2;
		offsetX2 = offsetX = _index_vote_tol(kv);
		datalen = datalen2 = vote -> items[offsetX2];
		unsigned int * dat2, *dat;
		dat = dat2 = vote -> pos[offsetX2];

		//LRMprintf("You can find KV at %u\n", kv);

		for(iix = 0; iix<=ii_end; iix = iix>0?-iix:(-iix+INDEL_SEGMENT_SIZE))
		{
			if(iix)
			{
				offsetX = _index_vote_tol(kv+iix);
				datalen = vote -> items[offsetX];
				dat = vote -> pos[offsetX];
			}


			if(!datalen)continue;

			for (i=0;i<datalen;i++)
			{
				int di = dat[i];
				int dist0 = kv-di;
				if( dist0 >= -indel_tolerance && dist0 <= indel_tolerance )
				{
					if(is_reversed  == (0!=(vote -> masks[offsetX][i]&LRMIS_NEGATIVE_STRAND)))
					{
						if(offset < vote->coverage_end [offsetX][i] + 10){
								unsigned char test_max = (vote->votes[offsetX][i]);
								test_max += 1;
								vote -> votes[offsetX][i] = test_max;

								if (offset +16 > vote->coverage_end [offsetX][i])
									vote->coverage_end [offsetX][i] = offset+16;

								/*
								int toli =  vote -> toli[offsetX][i];

								if (dist0 !=  vote->current_indel_cursor[offsetX][i])
								{
									toli +=3;
									if (toli < LRMMAX_INDEL_SECTIONS*3)
									{
										vote -> toli[offsetX][i] = toli;
										vote -> indel_recorder[offsetX][i][toli] = subread_number+1; 
										vote -> indel_recorder[offsetX][i][toli+1] = subread_number+1;
										vote -> indel_recorder[offsetX][i][toli+2] = dist0; 
										//LRMprintf("subread=#%d ,TOLI=%d, DIST0=%d, POS=%u \n", subread_number, toli, dist0, kv);
											
										if(toli < LRMMAX_INDEL_SECTIONS*3-3) vote -> indel_recorder[offsetX][i][toli+3]=0;
									}
									vote->current_indel_cursor [offsetX][i] = (char)dist0;
								} else vote -> indel_recorder[offsetX][i][toli+1] = subread_number+1;
								*/
								i = 9999999;
						}
					}
				}
			}
			if (i>=9999999){
				break;
			}

		}

		if (i < 9999999)
		{
			if (datalen2<LRMGENE_VOTE_SPACE)
			{
				vote -> items[offsetX2] ++;
				dat2[datalen2] = kv;
				vote -> masks[offsetX2][datalen2]=(is_reversed?LRMIS_NEGATIVE_STRAND:0);
				vote -> votes[offsetX2][datalen2]=1;
				vote -> toli[offsetX2][datalen2]=0;

				// data structure of recorder:
				// {unsigned char subread_start; unsigned char subread_end, char indel_offset_from_start}
				// All subread numbers are added with 1 for not being 0.

				vote -> indel_recorder[offsetX2][datalen2][0] = vote -> indel_recorder[offsetX2][datalen2][1] = subread_number+1;
				vote -> indel_recorder[offsetX2][datalen2][2] = 0;
				vote -> indel_recorder[offsetX2][datalen2][3] = 0;
				vote -> current_indel_cursor [offsetX2][datalen2] = 0;
				vote -> coverage_start [offsetX2][datalen2] = offset;
				vote -> coverage_end [offsetX2][datalen2] = offset+16;
				//LRMprintf("subread=#%d ,NEW RECORD =%u\n", subread_number, kv);
			}
		}
		else i=0;
	}
	return 1;
}


short LRMindel_recorder_copy(unsigned short *dst, unsigned short * src)
{
	short all_indels = 0;
//	memcpy(dst, src, 3*MAX_INDEL_TOLERANCE);  return;


	int i=0;
	while(src[i] && (i<3*LRMMAX_INDEL_TOLERANCE-2))
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

unsigned int LRMload_int32(FILE * fp)
{
	int ret;
	int read_length;
	read_length = fread(&ret, sizeof(int), 1, fp);
	if(read_length<=0)assert(0);
	return ret;
}

long long int LRMload_int64(FILE * fp)
{
	long long int ret;
	int read_length;
	read_length = fread(&ret, sizeof(long long int), 1, fp);
	if(read_length<=0)assert(0);
	return ret;
}


int LRMgehash_load_option(const char fname [], int option_no, int * result){
	char tabname[LRMMAX_FILENAME_LENGTH];
	char magic_chars[8];
	int found = 0;
	sprintf(tabname, "%s.00.b.tab", fname);
	FILE * fp = fopen(tabname, "rb");
	if(fp == NULL){
		sprintf(tabname, "%s.00.c.tab", fname);
		fp = fopen(tabname, "rb");
	}
	if(fp){	
		int rlen = fread(magic_chars,1,8,fp);
		if(rlen < 8){
			LRMprintf("ERROR: Unable to read-in the index.\n");
			return -1;
		}
		if(memcmp(magic_chars, "2subindx",7)==0) {
			while(1) {
				short option_key, option_length;

				rlen = fread(&option_key, 2, 1, fp);
				if(rlen<1){
					LRMprintf("ERROR: Unable to read-in the index.\n");
					return -1;
				}

				if(!option_key) break;

				rlen = fread(&option_length, 2, 1, fp);
				if(rlen<1){
					LRMprintf("ERROR: Unable to read-in the index.\n");
					return -1;
				}


				if(option_key == option_no){
					*result = 0;
					rlen = fread(result ,option_length,1,fp);
					if(rlen<1){
						LRMprintf("ERROR: Unable to read-in the index.\n");
						return -1;
					}
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

int LRMgehash_load(LRMgehash_t * the_table, const char fname [])
{
	int i, read_length;
	char magic_chars[8];
	magic_chars[7]=0;

	the_table -> index_gap = 0;

	FILE * fp = fopen(fname, "rb");
	if (!fp)
	{
		LRMprintf ("Table file '%s' is not found.\n", fname);
		return 1;
	}

	int rlen = fread(magic_chars,1,8,fp);
	if(rlen<8){
		LRMprintf("ERROR: Unable to read-in the index.\n");
		return -1;
	}
	
	if(memcmp(magic_chars+1, "subindx",7)==0)
	{
		if('2'==magic_chars[0])
			the_table -> version_number = LRMSUBINDEX_VER2;
		else	assert(0);

		while(1)
		{
			short option_key, option_length;

			rlen = fread(&option_key, 2, 1, fp);
			if(rlen<1){
				LRMprintf("ERROR: Unable to read-in the index.\n");
				return -1;
			}
			if(!option_key) break;

			rlen = fread(&option_length, 2, 1, fp);
			if(rlen<1){
				LRMprintf("ERROR: Unable to read-in the index.\n");
				return -1;
			}

			rlen = 999;
			if(option_key == LRMSUBREAD_INDEX_OPTION_INDEX_GAP)
				rlen = fread(&(the_table -> index_gap),2,1,fp);
			else if (option_key == LRMSUBREAD_INDEX_OPTION_INDEX_PADDING)
				rlen = fread(&(the_table -> padding),2,1,fp);
			else
				fseek(fp, option_length, SEEK_CUR);

			if(rlen<1){
				LRMprintf("ERROR: Unable to read-in the index.\n");
				return -1;
			}
		}
		assert(the_table -> index_gap);
	
		the_table -> current_items = LRMload_int64(fp);
		if(the_table -> current_items < 1 || the_table -> current_items > 0xffffffffllu){
			LRMputs("ERROR: the index format is unrecognizable.");
			return 1;
		}
		the_table -> buckets_number = LRMload_int32(fp);
		the_table -> buckets = (struct LRMgehash_bucket * )malloc(sizeof(struct LRMgehash_bucket) * the_table -> buckets_number);
		if(!the_table -> buckets)
		{
			LRMputs("Error: out of memory");
			return 1;
		}

		for (i=0; i<the_table -> buckets_number; i++)
		{
			struct LRMgehash_bucket * current_bucket = &(the_table -> buckets[i]);
			current_bucket -> current_items = LRMload_int32(fp);
			current_bucket -> space_size = LRMload_int32(fp);
			current_bucket -> space_size = current_bucket -> current_items;
			current_bucket -> new_item_keys = (short *) malloc ( sizeof(short) * current_bucket -> space_size);
			current_bucket -> item_values = (LRMgehash_data_t *) malloc ( sizeof(LRMgehash_data_t) * current_bucket -> space_size);

			if(!(current_bucket -> new_item_keys&&current_bucket -> item_values))
			{
				LRMputs("Error: out of memory");
				return 1;

			}

			if(current_bucket -> current_items > 0)
			{
				read_length = fread(current_bucket -> new_item_keys, sizeof(short), current_bucket -> current_items, fp);
				if(read_length < current_bucket -> current_items){
					LRMprintf("ERROR: the index is incomplete : %d < %u.\n",read_length, current_bucket -> current_items);
					return 1;
				}
				read_length = fread(current_bucket -> item_values, sizeof(LRMgehash_data_t), current_bucket -> current_items, fp);
				if(read_length < current_bucket -> current_items){
					LRMprintf("ERROR: the index value is incomplete : %d < %u.\n",read_length, current_bucket -> current_items);
					return 1;
				}
			}

		}

		read_length = fread(&(the_table -> is_small_table), sizeof(char), 1, fp);
		assert(read_length>0);
		fclose(fp);
		return 0;

	}
	else assert(0);
	return 0;
}

void LRMtest2key(unsigned int kk, char * obuf){
	int xx,oo=0;
	for(xx=0; xx<32;xx++){
		obuf[oo++] = (kk & (1<<xx))  ?'1':'0';
		if(xx%2 == 1 && xx < 31) obuf[oo++]=' ';
	}
	obuf[oo]=0;
}

int LRMtest2key_dist(unsigned int k1, unsigned int k2){
	int xx, ret = 0;
	for(xx=0; xx<16;xx++){
		int b1 = (k1 >> (xx*2)) & 3;
		int b2 = (k2 >> (xx*2)) & 3;
		if(b1!=b2) ret++;
	}
	return ret;
}

size_t LRMgehash_go_QQ(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, LRMgehash_t * the_table, LRMgehash_key_t raw_key, int offset, int read_len, int is_reversed, LRMgene_vote_t * vote, int indel_tolerance, int subread_number){

	struct LRMgehash_bucket * current_bucket;
	int i = 0, items;

	short *current_keys;//, *endp12;
	short key = raw_key / the_table->buckets_number;

	current_bucket = LRM_gehash_get_bucket (the_table, raw_key);
	items = current_bucket -> current_items;
	current_keys = current_bucket -> new_item_keys;
	
	if(!items) return 0;

	int imin=0, imax=items;
	int last_accepted_index = 0;

	while( imin < items )
	{
		last_accepted_index=(imin+imax)/2;
		short current_key = current_keys[last_accepted_index];
		if(current_key>key)
			imax = last_accepted_index - 1;
		else if(current_key<key)
			imin = last_accepted_index + 1;
		else
			break;

		if(imax<imin)
			return 0;
		
	}

	while(last_accepted_index){
		if(current_keys[last_accepted_index-1] == key) last_accepted_index-=1;
		else break;
	}

	for (; last_accepted_index<items && current_keys[last_accepted_index] == key ; last_accepted_index++)
	{
		unsigned int kv = current_bucket->item_values[last_accepted_index] - offset;
		int offsetX, datalen;
		offsetX = _index_vote(kv);
		datalen = vote -> items[offsetX];
		unsigned int *dat;
		dat = vote -> pos[offsetX];


		for (i=0;i<datalen;i++)
		{
			int di = dat[i];
			if( kv == di && is_reversed  == (0!=(vote -> masks[offsetX][i]&LRMIS_NEGATIVE_STRAND))  && offset < vote->coverage_end [offsetX][i] + 14){
				vote -> votes[offsetX][i]++;

				if (offset +16 > vote->coverage_end [offsetX][i])
					vote->coverage_end [offsetX][i] = offset+16;
				i = 9999999;
			}
		}

		if (i < 9999999)
		{
			if (datalen<LRMGENE_VOTE_SPACE)
			{
				vote -> items[offsetX] ++;
				dat[datalen]=kv;
				vote -> masks[offsetX][datalen]=(is_reversed?LRMIS_NEGATIVE_STRAND:0);
				vote -> votes[offsetX][datalen]=1;
				vote -> coverage_start [offsetX][datalen] = offset;
				vote -> coverage_end [offsetX][datalen] = offset+16;
			}
		}
		else i=0;
	}
	return 1;
}

size_t LRMgehash_go_tolerance(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, LRMgehash_t * the_table, LRMgehash_key_t key, int offset, int read_len, int is_reversed, LRMgene_vote_t * vote, int indel_tolerance, int subread_number, int max_MM){
	if(max_MM >= 10) return 0;	
	int ret = 0;

	ret+=LRMgehash_go_q(the_table, key, offset, read_len, is_reversed, vote, indel_tolerance, subread_number);
	int error_bases ;
	for (error_bases=1; error_bases <= max_MM; error_bases++)
	{
		assert(error_bases<5);
		int i, j;
		char error_pos_stack[10];	// max error bases = 10;
		LRMgehash_key_t mutation_key;

		for(i=0; i<error_bases; i++)
			error_pos_stack [i] = i;
		while (1)
		{

			char mutation_stack[10];
			memset(mutation_stack, 0 , error_bases);
			while(1)
			{
				int base_to_change=-1;
				mutation_key = key;

				for (j = 0; j < error_bases; j++)
				{
					base_to_change = error_pos_stack[j];
					int old_value = (key >> 2*base_to_change) & 3;
					int new_index = mutation_stack[j];

					int new_value;
					if(old_value <= new_index) new_value = 1+ new_index;
					else new_value = new_index;

					mutation_key = mutation_key & ~(0x3 << (2*base_to_change));
					mutation_key = mutation_key | (new_value << (2*base_to_change));
				}

				if(key != mutation_key ){
					int dret=LRMgehash_go_q(the_table, mutation_key, offset, read_len, is_reversed, vote,  indel_tolerance, subread_number);
					ret += dret;
					if(0 && LRMFIXLENstrcmp("@39076b1f-df29-4487-be51-4c30bf6c1cc4_Basecall_Alignment_template", iteration_context->read_name)==0){
						char bin_mutation_key[53], bin_key[53];
						LRMtest2key(mutation_key, bin_mutation_key);
						LRMtest2key(key, bin_key);
						LRMprintf("GO_TOLE_TEST: %s + %d DIST = %d HITS = %d :\nNEWKY:      %s\nORGKY:      %s\n",  iteration_context->read_name , offset, LRMtest2key_dist(mutation_key,key), dret, bin_mutation_key, bin_key);
					}
				}
				// increase one in the mutation_stack
				mutation_stack[error_bases-1]++;
				for(i = error_bases-1; i>=0; i--){
					if( mutation_stack[i]>2 ) {
						if(i == 0)break;
						else{
							mutation_stack[i] = 0;
							mutation_stack[i-1]++;
						}
					}
				}
				if(mutation_stack[0]>2)break;
			}


			int is_end = 1;
			for (i = error_bases-1; i>=0; i--)
				if (error_pos_stack[i] < 16 - (error_bases - i))
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

void LRMgehash_destory(LRMgehash_t * the_table)
{
	int i;

	for (i=0; i<the_table -> buckets_number; i++)
	{
		struct LRMgehash_bucket * current_bucket = &(the_table -> buckets[i]);
		if (current_bucket -> space_size > 0)
		{
			free (current_bucket -> new_item_keys);
			free (current_bucket -> item_values);
		}
	}

	free (the_table -> buckets);

	the_table -> current_items = 0;
	the_table -> buckets_number = 0;
}

void LRMprint_v(LRMcontext_t * context, LRMread_iteration_context_t * iteration_context, int min_votes){
	LRMgene_vote_t * v = &iteration_context-> vote_table;

	LRMprintf(" ==== VOTING TABLE ========================= \n");
	int iii, jjj;
	for(iii=0; iii < LRMGENE_VOTE_TABLE_SIZE; iii++){
		for(jjj = 0; jjj <  v -> items[iii]; jjj++){
			unsigned int this_pos = v -> pos[iii][jjj];
			int votes = v -> votes[iii][jjj];
			if(votes >= min_votes){
				char postxt[100];
				int conf_start = v -> coverage_start[iii][jjj] , conf_end = v -> coverage_end [iii][jjj];
				LRMpos2txt(context, this_pos, postxt);
				LRMprintf("     %d (%s) %3d - %3d : %16s   ", votes,(v -> masks[iii][jjj]&LRMIS_NEGATIVE_STRAND)?"NEG":"POS", conf_start, conf_end, postxt);
				int ix;
				for(ix=0; v->indel_recorder[iii][jjj][ix]; ix+=3){
					LRMprintf(" %d ~ %d : %d ",  v->indel_recorder[iii][jjj][ix] ,  v->indel_recorder[iii][jjj][ix+1],  v->indel_recorder[iii][jjj][ix+2]);
				}
				LRMprintf("\n");
			}
		}
	}
	LRMprintf(" =========================================== \n\n");
}
