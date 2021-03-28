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
  
  
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "gene-value-index.h"
#include "input-files.h"


int gvindex_init(gene_value_index_t * index, unsigned int start_point, unsigned int total_bases_estimated)
{
	index->start_point = start_point;
	index->memory_block_size = total_bases_estimated/4+10;
	index->values = calloc(index->memory_block_size,1);
	if(!index->values)
	{
		SUBREADputs(MESSAGE_OUT_OF_MEMORY);
		return 1;
	}
	index -> start_base_offset = index -> start_point - index -> start_point%4;
	return 0;
}

#define gvindex_baseno2offset_m(base_number, index, offset_byte, offset_bit)	{offset_byte =  (base_number - index -> start_base_offset) >>2; offset_bit = base_number % 4 * 2;}


void gvindex_baseno2offset(unsigned int base_number, gene_value_index_t * index, unsigned int * offset_byte, unsigned int * offset_bit)
{
	// the base number corrsponding to the 0-th bit in the whole value array;

	unsigned int offset = (base_number - index -> start_base_offset);

	* offset_byte = offset >>2 ;
	* offset_bit = base_number % 4 * 2;
}

int is_offset_in_chro(gene_value_index_t * offsets, gehash_data_t linear){
	int ret = 1;
	if(offsets -> appendix1 && offsets -> appendix2){
		gene_offset_t * chros = offsets -> appendix1;
	
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
				unsigned int pos;
	
				if (n==0)
					pos = linear;
				else
					pos = linear - chros->read_offsets[n-1];
	
				if( pos < chros -> padding  || pos >= ( chros->read_offsets[n] - chros->read_offsets[n-1] - chros -> padding ))
					ret = 0;
					SUBREADprintf("INCHRO:%d ; POS:%d\n", ret, pos);
				break;
			}
		}
	}
	return ret;
}

// return 'A', 'G', 'T' and 'C'
int gvindex_get(gene_value_index_t * index, gehash_data_t offset)
{
	unsigned int offset_byte, offset_bit;
	//if(!is_offset_in_chro( index, offset ))return -1;

	gvindex_baseno2offset_m(offset, index , offset_byte, offset_bit);
	if(offset_byte >= index-> values_bytes -1)return 'N';
	unsigned int one_base_value = (index->values [offset_byte]) >> (offset_bit);
	//SUBREADprintf("RECV_BASE=%d (%d -  %d)\n",one_base_value & 3, offset_byte , offset_bit);

	return int2base(one_base_value & 3);
}

int gvindex_match(gene_value_index_t * index, gehash_data_t offset, gehash_key_t base_values)
{
	unsigned int offset_byte, offset_bit;

	gvindex_baseno2offset_m(offset, index , offset_byte, offset_bit);
	int i, ret = 0;

	for (i=0; i<16; i++)
	{
		unsigned char mask = 0x3 << (offset_bit);
		unsigned char one_base_value = (index->values [offset_byte] & mask) >> (8-offset_bit);
		if ( ((base_values >> (30 - i*2)) & 0x3) == one_base_value)
			ret |= 1 << i;

		offset_bit +=2;
		if(offset_bit >=8)
		{
			offset_bit = 0;
			offset_byte ++;
		}
	}

	return ret;

}

void gvindex_set(gene_value_index_t * index, gehash_data_t offset, gehash_key_t base_values, int padding)
{
	unsigned int offset_byte, offset_bit;
	gvindex_baseno2offset(offset, index , &offset_byte, &offset_bit);
	int i;

	unsigned int offset_byte_margin = offset_byte + (padding / 8 +1);
	if(index -> memory_block_size <= offset_byte_margin + 2)
	{
		assert(index -> memory_block_size> offset_byte_margin+2);
		index -> memory_block_size *= 1.5;
		index->values = realloc(index->values, index -> memory_block_size);
	}
	for (i=0; i<16; i++)
	{
		// 11110011
		//     ^^ base
		unsigned char mask = 0xff << (offset_bit+2) | 0xff >> (8-offset_bit);
		index->values [offset_byte] &= mask;
		index->values [offset_byte] |= ((base_values >> (30 - i*2))&0x03) << (offset_bit);

		//SUBREADprintf("PUT_BASE=%d (%d - %d)\n", (base_values  >> (30 - i*2))&0x03, offset_byte, offset_bit);
		offset_bit +=2;
		if(offset_bit >=8)
		{
			offset_bit = 0;
			offset_byte ++;
		}
	}

	index -> length = offset + 16 - index -> start_point + padding;
}

int gvindex_dump(gene_value_index_t * index, const char filename [])
{
	FILE * fp = f_subr_open(filename, "wb");
	int is_full = 0;
	int write_len = fwrite(&index->start_point,4,1, fp);

	write_len += fwrite(&index->length, 4, 1, fp);
	if(write_len < 2){
		is_full = 1;
	}

	unsigned int useful_bytes, useful_bits;
	gvindex_baseno2offset (index -> length+ index -> start_point, index,&useful_bytes,&useful_bits);

	write_len = fwrite(index->values, 1, useful_bytes+1, fp);
	if(write_len <= useful_bytes) is_full = 1;
	fclose(fp);
	if(is_full) SUBREADprintf("ERROR: unable to writeinto the output file. Please check the disk space in the output directory.\n");
	return is_full;
}


int gvindex_load(gene_value_index_t * index, const char filename [])
{
	memset(index,0, sizeof(gene_value_index_t));
	FILE * fp = f_subr_open(filename, "rb");
	int read_length;
	read_length = fread(&index->start_point,4,1, fp);
	if(read_length<1){
		SUBREADprintf("ERROR: the array index is incomplete : %d\n", read_length );
		return 1;
	}
	read_length = fread(&index->length,4,1, fp);
	if(read_length<1){
		SUBREADputs("ERROR: the index is incomplete.");
		return 1;
	}
	//SUBREADprintf ("\nBINDEX %s : %u ~ +%u\n",filename, index->start_point, index->length );

	unsigned int useful_bytes, useful_bits;
	index -> start_base_offset = index -> start_point - index -> start_point%4;
	gvindex_baseno2offset (index -> length+ index -> start_point, index ,&useful_bytes,&useful_bits);
	index -> values = malloc(useful_bytes+1);
	index -> values_bytes = useful_bytes+1;
	if(!index->values)
	{
		SUBREADputs(MESSAGE_OUT_OF_MEMORY);
		return 1;
	}
	

	read_length =fread(index->values, 1, useful_bytes+1, fp);
	if(read_length < useful_bytes){
		SUBREADprintf("ERROR: the array index is incomplete : %d < %d.\n", read_length, useful_bytes+1 );
		return 1;
	}

	fclose(fp);
	return 0;

}


int match_chro_wronglen(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int space_type, int * left_match_bases, int * right_match_bases)
{
	int ret = 0;
	int i;
	char last_char='A';
	int left_correct_end = 0;
	if(left_match_bases) *left_match_bases=0;
	if(right_match_bases) *right_match_bases=0;

	if (space_type == GENE_SPACE_COLOR)
		last_char = (pos <= index -> start_point)?'A': gvindex_get(index,pos-1);

	for (i=0;i<test_len;i++)
	{
		char tt = gvindex_get (index, pos +i);
		int newv;
		if(space_type == GENE_SPACE_COLOR)
		{

			newv = read[i] == '0'+chars2color(last_char, tt); 
			last_char = tt;
		}
		else
			newv =read[i] == tt; 

		//if(left_wrong_bases)
		//	SUBREADprintf("I=%d, *LWB=%d, LWE=%d\n", i, *left_wrong_bases, left_wrong_end);

		if(left_match_bases && (newv) && (!left_correct_end ))
			(*left_match_bases)++;
		else if (!newv)left_correct_end=1;

		if(right_match_bases && (newv))
			(*right_match_bases) ++;
		else if (right_match_bases)
			(*right_match_bases) =0;

		ret += newv;
	}

	return ret;
}

#define INDEL_TEST_WINDOW 5
#define MIN_INDEL_SEARCH_MATCH_SCORE 8500

int match_indel_chro_to_front(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int * indels, int * indel_point, int max_indel_number, int max_test_length)
{
	int offset = 0;
	int i;
	int ret = 0;
	int bast_match_score_remailing=-1;

	for(i=0; i < test_len; i++)
	{
		unsigned int npos = pos + i +  offset;
		char tt;
		tt = gvindex_get (index, npos);

		//if(i-min(0,offset) >= max_test_length)break;

	//	printf("Xi=%d  POS=%u INDEL=%d MATCH=%d TT vs R = %c , %c ; IPOS=%d\n", i,  npos, offset, read[i] ==tt, tt,  read[i], *indel_point);
		if(read[i]==tt)ret++;
		else if((i + offset < test_len - INDEL_TEST_WINDOW - 3 && i >0))
		{
			// if there is a base unmatched, it is potentially an indel from here.
			int window_match = match_chro(read+i, index, pos+i+ offset, INDEL_TEST_WINDOW ,0,GENE_SPACE_BASE);

			if(window_match < INDEL_TEST_WINDOW -1)
			{
				// if the window is badly matched, it is very likely to be an indel from this base.
				int indel_test_i;
				for(indel_test_i =0; indel_test_i < 7; indel_test_i++)
				{
					int indel_test = (indel_test_i+1)/2*(indel_test_i%2?1:-1);
					if(abs(indel_test)>max_indel_number) continue;

					if(indel_test < 0)	//insertion 
					{
						int matched_tail = match_chro(read+i-indel_test, index, pos+i, test_len - i+indel_test,0,GENE_SPACE_BASE);
						int matched_score = matched_tail * 10000 / ( test_len - i + indel_test);
						if(matched_score >  bast_match_score_remailing &&  matched_score > MIN_INDEL_SEARCH_MATCH_SCORE)
						{
							offset = indel_test;
							bast_match_score_remailing = matched_score;
						}
					}else	// deletion
					{
						int matched_tail = match_chro(read+i, index, pos+i+indel_test, test_len - i ,0,GENE_SPACE_BASE);
						int matched_score = matched_tail * 10000 / (test_len - i);
						if(matched_score >  bast_match_score_remailing &&  matched_score > MIN_INDEL_SEARCH_MATCH_SCORE)
						{
							offset = indel_test;
							bast_match_score_remailing = matched_score;
						}
					}
				}
			}
			if(bast_match_score_remailing>0)
			{
				if(offset > 0)//deletion
				{
					*indel_point  = i;
					npos = pos + i + offset;
					tt = gvindex_get (index, npos);
					ret += read[i] == tt;
				}
				else if(offset<0)
				{
					*indel_point  = i;
					i-=offset+1;
				}
			}
		}
	}
	*indels = offset;
	return ret;

}


int find_all_indels(HashTable *indel_table, unsigned long pos, indel_record_t * ret, int search_to_back)
{
	int max_test_indels = 16;

	int xi, iret=0;

	
	if(!search_to_back)
	{
		int mask_pos = (pos>>3)& 0x1fffffff;
		if(!((char *)indel_table->appendix1)[mask_pos]) return 0;

		int mask_off = pos %8;
		if(!(((char *)indel_table->appendix1)[mask_pos] & (1<<mask_off)))
			return 0;
	}

	for(xi = -max_test_indels; xi<=max_test_indels; xi++)
	{
		if(!xi)continue;

		unsigned int search_pos = pos;

		if(search_to_back)
		{
			if(xi>0)search_pos-= xi;
			search_pos--;
		}


		if(search_to_back)
		{
			int mask_pos = (search_pos>>3)& 0x1fffffff;
			if(!((char *)indel_table->appendix1)[mask_pos]) continue;

			int mask_off = search_pos %8;
			if(!(((char *)indel_table->appendix1)[mask_pos] & (1<<mask_off)))
				continue;
		}



		indel_record_t key;
		key.pos = search_pos;
		key.len = xi;


		indel_result_t * res = HashTableGet(indel_table, &key);
		if(res)
		{
			memcpy(ret+iret, &key, sizeof(indel_record_t));
			iret++;
			if(iret>9)break;
		}
	}
	//printf("%u has %d indel for %d\n", pos, iret, search_to_back);
	return iret;
}


void match_indel_table_to_front_in(HashTable * indel_table , char * read, int read_pos_first_base, gene_value_index_t * index, unsigned int first_base_pos, int test_len, short * total_indels, short * indel_point, int max_indel_number, int min_test_offset, short * indels, short * indel_poses, int sofar_matched_bases, short * best_indels, short * best_indel_poses, int * best_matching_bases, int level);
// This function iterates all possible combinations of indels, and select the one that has the highest matched base number as the final result.
// The final result will be written in *sec

int match_indel_table_to_front(HashTable * indel_table , char * read, gene_value_index_t * index, unsigned int first_base_pos, int test_len, short * indels, short * indel_point, int max_indel_number, int min_test_offset, struct explorer_section_t * sec)
{
	short tmp_indel_rec[10], tmp_indel_pos_rec[10];
	memset(tmp_indel_rec, 0, sizeof(short)*10);
	memset(tmp_indel_pos_rec, 0, sizeof(short)*10);

	short *best_indel_rec=NULL, *best_indel_pos_rec=NULL;

	if(sec)
	{
		best_indel_rec = sec -> all_indels ;
		best_indel_pos_rec = sec -> all_indel_poses;
	
		memset(best_indel_rec, 0, sizeof(short)*10);
		memset(best_indel_pos_rec, 0, sizeof(short)*10);
	}
	int best_matching_bases = -1;

	match_indel_table_to_front_in(indel_table, read, 0, index, first_base_pos, test_len , indels , indel_point , max_indel_number , min_test_offset , tmp_indel_rec , tmp_indel_pos_rec , 0, best_indel_rec, best_indel_pos_rec , &best_matching_bases , 0);

	return best_matching_bases;

}

#define TEST_TARGET "AGGAAAGTGAGACCATGCTCAGGAAAACAGAAGGCAACAGAAACTATGACAAATTCCAGACATCGGATTAACAGACAAAGACAG"


void match_indel_table_to_front_in(HashTable * indel_table , char * read, int read_pos_first_base, gene_value_index_t * index, unsigned int first_base_pos, int test_len, short * total_indels, short * indel_point, int max_indel_number, int min_test_offset, short * indels, short * indel_poses, int sofar_matched_bases, short * best_indels, short * best_indel_poses, int * best_matching_bases, int level)
{

	unsigned int xi;
	int xj, xk;
	unsigned int index_low_boundary = index -> start_base_offset;
	unsigned int index_up_boundary = index -> start_base_offset + index -> length;

	if(first_base_pos <= index_low_boundary)return;
	if(first_base_pos + test_len >= index_up_boundary)return;


	if(level >9)return;

	for(xi = 1; xi < test_len; xi++)
	{
		indel_record_t found_records[20];
		int section_matched_bases = 0 ,  num_records = find_all_indels(indel_table , xi + first_base_pos, found_records, 0);
		if(num_records >0)
			section_matched_bases = match_chro(read, index, first_base_pos, xi ,0,GENE_SPACE_BASE);
//		if(num_records) printf("XI %d   FOUND %d INDELS AT %u ; SECMATCH=%d\n",xi, num_records , first_base_pos + xi, section_matched_bases);
		for(xj=0; xj<num_records; xj++)
		{
			indel_record_t * this_record = &found_records[xj];

			unsigned int next_first_base_pos = xi + first_base_pos;
			if(this_record -> len >0) next_first_base_pos += this_record -> len; // if this is a deletion, jump the gap.

			char * next_read = read + xi;
			if(this_record -> len <0) next_read -= this_record -> len ;

			int remainder_testlen = test_len - xi;
			if(this_record -> len <0) remainder_testlen += this_record -> len ;

			indels[level] = this_record -> len;
			indel_poses[level] = read_pos_first_base + xi;

			int next_read_pos_first_base = read_pos_first_base + xi;
			if(this_record -> len <0) next_read_pos_first_base -= this_record -> len;

			if(remainder_testlen>0 && (remainder_testlen<test_len -5))
				match_indel_table_to_front_in(indel_table, next_read, next_read_pos_first_base, index, next_first_base_pos , remainder_testlen , total_indels, indel_point, max_indel_number, min_test_offset, indels, indel_poses, sofar_matched_bases + section_matched_bases, best_indels, best_indel_poses, best_matching_bases, level+1) ;

			for(xk=level;xk<10;xk++)
			{
				indels[xk]=0;
				indel_poses[xk]=0;
			}
		}
	} 
	
	
	// test the number of matched bases if no indels are in this section


	int onepiece_matched_bases = match_chro(read, index, first_base_pos, test_len ,0,GENE_SPACE_BASE);

	//if(onepiece_matched_bases * 1. / test_len < 0.85) onepiece_matched_bases = 0;



	if(onepiece_matched_bases + sofar_matched_bases > *best_matching_bases)
	{
		if(best_indels)
		{
			memcpy(best_indels, indels, sizeof(short)*10);
			memcpy(best_indel_poses, indel_poses, sizeof(short)*10);
		}

		*best_matching_bases = onepiece_matched_bases + sofar_matched_bases;
		xk=0;

		int ret_indels=0;
		while(indels[xk])
			ret_indels += indels[xk++];
		*total_indels = ret_indels;

		if(indels[0])
			*indel_point = indel_poses[0];
		else
			*indel_point = 0;
	}
}



void match_indel_table_to_back_in(HashTable * indel_table , char * read, gene_value_index_t * index, unsigned int last_base_pos, int test_len, short * total_indels, short * indel_point, int max_indel_number, int min_test_offset, short * indels, short * indel_poses, int sofar_matched_bases, short * best_indels, short * best_indel_poses, int * best_matching_bases, int level);
int match_indel_table_to_back(HashTable * indel_table, char * read, gene_value_index_t * index, unsigned int last_base_pos, int test_len, short * indels, short * indel_point, int max_indel_number, int min_test_offset, struct explorer_section_t * sec)
{
	short tmp_indel_rec[10], tmp_indel_pos_rec[10];
	memset(tmp_indel_rec, 0, sizeof(short)*10);
	memset(tmp_indel_pos_rec, 0, sizeof(short)*10);

	short best_indel_rec [10] , best_indel_pos_rec [10];
	memset(best_indel_rec, 0, sizeof(short)*10);
	memset(best_indel_pos_rec, 0, sizeof(short)*10);
	int best_matching_bases = -1;

	match_indel_table_to_back_in(indel_table, read, index, last_base_pos, test_len , indels , indel_point , max_indel_number , min_test_offset , tmp_indel_rec , tmp_indel_pos_rec , 0, best_indel_rec, best_indel_pos_rec , &best_matching_bases , 0);

	int xi;
	if(sec)
	{
		int retv=0;
		while(best_indel_rec[retv])retv++;
		//printf("BACK-OUT NO.INDEL=%d\n", retv);

		for(xi=0;xi<retv;xi++)
		{
			sec->all_indels[xi] = best_indel_rec[retv-xi-1];
			sec->all_indel_poses[xi] = best_indel_pos_rec[retv-xi-1];
		}
		for(xi=retv;xi<10;xi++)
			sec->all_indels[xi] = 0;

	}
	return best_matching_bases;
}

void match_indel_table_to_back_in(HashTable * indel_table , char * read, gene_value_index_t * index, unsigned int last_base_pos, int test_len, short * total_indels, short * indel_point, int max_indel_number, int min_test_offset, short * indels, short * indel_poses, int sofar_matched_bases, short * best_indels, short * best_indel_poses, int * best_matching_bases, int level)
{

	unsigned int xi;
	int xj, xk;
	unsigned int index_low_boundary = index -> start_base_offset;
	unsigned int index_up_boundary = index -> start_base_offset + index -> length;

	if(last_base_pos - test_len <= index_low_boundary)return;
	if(last_base_pos >= index_up_boundary)return;

	if(level >9)return;
	if(0&&memcmp(read, TEST_TARGET, 60)==0)
		printf("BACK-OUT: %s, LEV=%d, LEN=%d, TAILPOS=%u\n",read, level, test_len, last_base_pos);
	for(xi = 1 ; xi < test_len-1 ; xi++)
	{
		indel_record_t found_records[20];
		int section_matched_bases = 0, num_records = find_all_indels(indel_table , last_base_pos - xi, found_records, 1);

		if(num_records>0)
			section_matched_bases = match_chro(read + test_len - xi, index, last_base_pos - xi, xi ,0,GENE_SPACE_BASE);
	//	if(num_records) printf("XI %d   FOUND %d INDELS AT %u ; SECMATCH=%d\n",xi, num_records , last_base_pos - xi, section_matched_bases);

		for(xj=0; xj<num_records; xj++)
		{
			indel_record_t * this_record = &found_records[xj];

			unsigned int next_last_base_pos = last_base_pos - xi -1;
			if(this_record -> len >0) next_last_base_pos -= this_record -> len; // if this is a deletion, jump the gap.

			int remainder_testlen = test_len - xi -1;
			if(this_record -> len <0) remainder_testlen += this_record -> len ;


			indels[level] = this_record -> len;
			indel_poses[level] = test_len - xi -1;
			if(this_record -> len <0) indel_poses[level] +=  this_record -> len ;

			if(remainder_testlen> 0  &&remainder_testlen<test_len -5)
				match_indel_table_to_back_in(indel_table, read, index, next_last_base_pos , remainder_testlen , total_indels, indel_point, max_indel_number, min_test_offset, indels, indel_poses, sofar_matched_bases + section_matched_bases, best_indels, best_indel_poses, best_matching_bases, level+1) ;

			for(xk=level;xk<10;xk++)
			{
				indels[xk]=0;
				indel_poses[xk]=0;
			}
		}
	} 
	
	
	// test the number of matched bases if no indels are in this section


	int onepiece_matched_bases = match_chro(read, index, last_base_pos - test_len, test_len ,0,GENE_SPACE_BASE);
	//if(onepiece_matched_bases *1. / test_len < 0.85)onepiece_matched_bases = 0;




	if(onepiece_matched_bases + sofar_matched_bases > *best_matching_bases)
	{
	if(0&&memcmp(read, TEST_TARGET, 60)==0)
	{
		for(xk=0; xk<level; xk++)
		{
			printf("%+d@%d ; ", indels[xk], indel_poses[xk]);
		}


		printf("BACK FINAL LEVEL %d : MATCHED=%d+%d=%d\n", level, onepiece_matched_bases , sofar_matched_bases, onepiece_matched_bases + sofar_matched_bases );
	}


		memcpy(best_indels, indels, sizeof(short)*10);
		memcpy(best_indel_poses, indel_poses, sizeof(short)*10);
		*best_matching_bases = onepiece_matched_bases + sofar_matched_bases;
		xk=0;

		int ret_indels=0;
		while(best_indels[xk])
			ret_indels += best_indels[xk++];
		*total_indels = ret_indels;

		if(best_indels[0])
			*indel_point = best_indel_poses[0];
		else
			*indel_point = 0;
	}
}





// "pos" here is the expected position of the head of the read , but it is not the final position. If indels are found in the read, the head position must be offset.
// Only certain is pos+test_len is the EXACT position of the TAIL of the read.
int match_indel_chro_to_back(char * read, gene_value_index_t * index, unsigned int pos /*_of_the_first_base (assumed)*/, int test_len, int * indels, int * indel_point, int max_indel_number, int min_test_offset)
{
	//return  match_chro(read, index, pos, test_len, 0, 1);
	//SUBREADprintf("TEST_INDEL_CHRO %s VS %u LEN=%d\n", read, pos, test_len);
	int offset = 0;
	int i;
	int ret = 0;
	unsigned int tail_pos = pos  + test_len;
	int bast_match_score_remailing=-1;


	if(pos > 0xffff0000 || pos + test_len>= index-> length + index->start_point) 
	{
		*indels = 0;
		return 0;
	}
	

	// i is the distance from the tail.
	for(i=1 ; i < test_len /*max(0,offset)*/; i++)
	{
		unsigned int npos = tail_pos - i - 1 -  offset;
		char tt;

		tt = gvindex_get (index, npos);

		//printf("POS=%u INDEL=%d MATCH=%d\n",  npos, offset, read[test_len - i - 1] ==tt);

		if(read[test_len - i - 1] ==tt)ret++;
		else if(test_len - i - 1 > INDEL_TEST_WINDOW +1 && i > 1)
		{
			// if there is a base unmatched, it is potentially an indel from here.
			int window_match = match_chro(read + test_len - i  - INDEL_TEST_WINDOW, index, tail_pos - i  - offset - INDEL_TEST_WINDOW, INDEL_TEST_WINDOW ,0,GENE_SPACE_BASE);

			if(window_match < INDEL_TEST_WINDOW - 1)
			{	
				// if the window is badly matched, it is very likely to be an indel from this base.
				int indel_test_i;
				for(indel_test_i =0; indel_test_i < 7; indel_test_i++)
				{
					int indel_test = (indel_test_i+1)/2*(indel_test_i%2?1:-1);
					if(abs(indel_test)>max_indel_number) continue;
					if(indel_test > 0)	//deletion 
					{
						int matched_tail = match_chro(read , index, tail_pos - test_len - indel_test , test_len - i ,0,GENE_SPACE_BASE);
						int matched_score = matched_tail * 10000 / ( test_len - i );

						if(matched_score >  bast_match_score_remailing &&  matched_score > MIN_INDEL_SEARCH_MATCH_SCORE)
						{
							offset = indel_test;
							bast_match_score_remailing = matched_score;
						}
					}else	//insertion 
					{
						int matched_tail = match_chro(read, index, tail_pos - test_len - indel_test , test_len - i + indel_test ,0,GENE_SPACE_BASE);
						int matched_score = matched_tail * 10000 / (test_len - i + indel_test);

						if(matched_score >  bast_match_score_remailing &&  matched_score > MIN_INDEL_SEARCH_MATCH_SCORE)
						{
							offset = indel_test;
							bast_match_score_remailing = matched_score;
						}
					}
				}
			}
			if(bast_match_score_remailing>0)
			{
				if(offset < 0)//insertion
				{
					//ret += read[i+offset] == tt;
					*indel_point  = test_len - i + offset ;
					i -= offset + 1;
				}
				else if(offset >0)	//deletion
				{
					npos = (tail_pos - i - offset) ;
					tt = gvindex_get (index, npos);
					ret += read[test_len - i - 1] == tt;
					*indel_point  = test_len - i ;
				}
			}
		}
	}
	//if(memcmp(read, "AACCCCTTGCAGAAAA", 15)==0)
	//	SUBREADprintf("\nINDEL_SPOT %d\n", offset);
	*indels = offset;
	return ret;
}



float match_chro_support(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int is_negative_strand, int space_type, char * qual_txt, int qual_format)
{
	int i;
	int all_qual = 0;
	int supported_qual = 0; 
				

	if (is_negative_strand)
	{

		for (i=test_len -1;i>=0;i--)
		{
			char tt = gvindex_get(index, pos+test_len-1-i);
			int is_correct = 0;
			switch(tt)
			{
				case 'A': is_correct = read[i] == 'T'; break;
				case 'T': is_correct = read[i] == 'A'; break;
				case 'G': is_correct = read[i] == 'C'; break;
				case 'C': is_correct = read[i] == 'G'; break;
			}

			int base_p = 0.;

			if(qual_txt[0])
			{
				if(FASTQ_PHRED64 == qual_format)
				{
					base_p = get_base_error_prob64i(qual_txt[i]);
				}
				else
				{
					base_p = get_base_error_prob33i(qual_txt[i]);
				}
			}

			if(base_p > 300000) continue;
			base_p=0;
			all_qual += (1000000-base_p);

			if(is_correct)
				supported_qual += (1000000-base_p);
		}
	}
	else
	{
		if(qual_txt[0])
		{
			if(FASTQ_PHRED33 == qual_format)
			{
				for (i=0;i<test_len;i++)
				{
					char tt = gvindex_get (index, pos +i);
					int is_correct =read[i] == tt; 
					int base_p = get_base_error_prob33i(qual_txt[i]);

					if(base_p > 300000) continue;
					base_p=0;

					all_qual += (1000000-base_p);
	
					if(is_correct)
						supported_qual += (1000000-base_p);
				}
			}
			else
			{
				for (i=0;i<test_len;i++)
				{
					char tt = gvindex_get (index, pos +i);
					int is_correct =read[i] == tt; 
					int base_p = get_base_error_prob64i(qual_txt[i]);

					if(base_p > 300000) continue;
					base_p=0;

					all_qual += (1000000-base_p);
	
					if(is_correct)
						supported_qual += (1000000-base_p);
				}
			}
		}
		else
			for (i=0;i<test_len;i++)
			{
				char tt = gvindex_get (index, pos +i);
				int is_correct =read[i] == tt; 
				

				//printf("i=%d/%d\t%c==%c  : %d\n",i, test_len, tt, read[i], is_correct);

				all_qual += (1000000);

				if(is_correct)
					supported_qual += (1000000);
			}	
	}

	//SUBREADprintf("%d\n", test_len);
	//if(all_qual < 3100000) return 0;
	return supported_qual*1. / all_qual * test_len;
}




int match_chro(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int is_negative_strand, int space_type)
{

	
	int ret = 0;
	int i;
	char last_char='A';
	
	if ((unsigned int)(pos + test_len) >= index -> length + index -> start_point) return 0;
	if (pos > 0xffff0000) return 0;

	if (is_negative_strand)
	{

		if (space_type == GENE_SPACE_COLOR)
		{
			pos++;
			last_char = (pos+test_len>= index -> length + index -> start_point)?'A': gvindex_get(index,pos+test_len);
			for (i=test_len -1;i>=0;i--)
			{
				char tt = gvindex_get (index, pos+test_len-1-i);
				ret += read[i] == '0'+chars2color(tt, last_char); 
				last_char = tt;
			}
		}
		else
		{

			for (i=test_len -1;i>=0;i--)
			{
				char tt = gvindex_get (index, pos+test_len-1-i);
				switch(tt)
				{
					case 'A': ret += read[i] == 'T'; break;
					case 'T': ret += read[i] == 'A'; break;
					case 'G': ret += read[i] == 'C'; break;
					case 'C': ret += read[i] == 'G'; break;
				}
			}
	


		}
	}	
	else
	{
		if(space_type == GENE_SPACE_BASE)
		{
			unsigned int offset_byte, offset_bit;
	//puts("P1");

			gvindex_baseno2offset_m(pos, index , offset_byte, offset_bit);

	//printf("OB=%d; V=%016llX\n", offset_byte, (long long )index->values);
			if(offset_byte >= index-> values_bytes)return 0;
			char int_value = index->values [offset_byte];

			for (i=0;i<test_len;i++)
			{
				char tt = (int_value >> offset_bit) & 3;
				char tv = read[i];
				switch(tv){
					case 'A':
						ret += tt==0;
						break;
					case 'G':
						ret += tt==1;
						break;
					case 'C':
						ret += tt==2;
						break;
					case 0:
						//SUBREADprintf("NON-ATGC-CHAR:%d\n", tv);
						//assert(0);
						break;
					default:
						ret += tt==3;

				}
				offset_bit+=2;
				if(offset_bit==8)
				{
					offset_byte++;
					if(offset_byte == index-> values_bytes)return 0;
					int_value = index->values [offset_byte];
					offset_bit = 0;
				}
			}
	//puts("P3");
		}
		else
		{
			last_char = (pos <= index -> start_point)?'A': gvindex_get(index,pos);
			for (i=0;i<test_len;i++)
			{
				char tt = gvindex_get (index, pos +i+1);
				ret += read[i] == '0'+chars2color(last_char, tt);
				last_char = tt;
			}
		}
		
	}
	return ret;
}


int match_chro_slow(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int is_negative_strand, int space_type)
{
	if(is_negative_strand || space_type == GENE_SPACE_COLOR)
		return match_chro_slow(read, index, pos, test_len, is_negative_strand, space_type);


	unsigned int i;
	unsigned int test_end = test_len+pos- index->start_base_offset;
	int ret = 0, offset, bits;
	for(i=(pos - index->start_base_offset); i< test_end ; i++)
	{
		offset = i/4;
		bits = i%4*2;
		ret += ((index -> values[offset] >> bits) & 0x3 )== base2int(*read);
		read++;
	}

	return ret;
}

unsigned int match_chro_range(char * read, gene_value_index_t * index, unsigned int pos, int read_len, int search_length, int search_to_back)
{
	short key[4];
	int i, j;
	for(i = 0; i < 4; i++)
	{
		key[i]=0;
		for(j = i+7; j >= i; j--)
		{
			key[i] = key[i] << 2 | base2int(read[j]);
		}
	}


	unsigned int offset_byte, offset_bit, search_dist;

	gvindex_baseno2offset(pos, index , &offset_byte, &offset_bit);

	/*
	if(offset_byte > 0x35000000)
		SUBREADprintf("POS=%u, OFFBYTE=%u, index->start_base_off=%u\n" , pos, offset_byte, index -> start_base_offset);
	*/
	
	search_dist = search_length /4;
	if(search_to_back == SEARCH_BACK)
	{
		if(search_dist > offset_byte - 500)search_dist = offset_byte- 500;
	}else
		if(search_dist + offset_byte >=  index -> values_bytes - 500) search_dist = index -> values_bytes - offset_byte -501;

	for (i=2; i<search_dist; i++)
	{
		unsigned long test_offset = offset_byte;
		test_offset += (search_to_back == SEARCH_BACK)?-i:i;
		short tv = *(short *)(index->values +test_offset);

		for(j=0; j<4; j++)
		{
			if(tv == key[j])
			{
				//test the whole read
				unsigned int hit_pos = test_offset*4;
				hit_pos += index -> start_base_offset - j;// + index -> start_point; 
				int retv = match_chro_maxerror(read, index, hit_pos, read_len, 0, 0, 0);
				if(retv >0)
				{
					//SUBREADprintf("POS=%u, TOFF=%u, STARTBASE=%d, j=%d , RETV=%d\n" , test_offset*4 + index -> start_base_offset - j , test_offset , index -> start_base_offset , j , retv);
					return hit_pos; 
				}
			}
		}
	}
	return 0xffffffff;
}

int match_chro_maxerror(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int is_negative_strand, int space_type, int max_error)
{
	int ret = 0;
	int i;
	char last_char='A';
	

	if (is_negative_strand)
	{

		if (space_type == GENE_SPACE_COLOR)
		{
			pos++;
			last_char = (pos+test_len>= index -> length + index -> start_point)?'A': gvindex_get(index,pos+test_len);
		}
		for (i=test_len -1;i>=0;i--)
		{
			char tt = gvindex_get (index, pos+test_len-1-i);
			if(space_type == GENE_SPACE_COLOR)
			{
				ret += read[i] != '0'+chars2color(tt, last_char); 
				last_char = tt;
			}
			else
				switch(tt)
				{
					case 'A': ret += read[i] != 'T'; break;
					case 'T': ret += read[i] != 'A'; break;
					case 'G': ret += read[i] != 'C'; break;
					case 'C': ret += read[i] != 'G'; break;
				}
			if(ret>max_error)return 0;
		}
	}
	else
	{
		if (space_type == GENE_SPACE_COLOR)
			last_char = (pos <= index -> start_point)?'A': gvindex_get(index,pos-1);
		for (i=0;i<test_len;i++)
		{
			char tt = gvindex_get (index, pos +i);
			if(space_type == GENE_SPACE_COLOR)
			{
				ret += read[i] != '0'+chars2color(last_char, tt);
				last_char = tt;
			}
			else
				ret +=read[i] != tt; 
			//SUBREADprintf("RET=%d\n",ret);
			if(ret>max_error)return 0;
		}
	}
	return test_len - ret;
}




int gvindex_match_base(gene_value_index_t * index, gehash_data_t offset, const char base_int_value)
{
	unsigned int offset_byte, offset_bit;

	gvindex_baseno2offset(offset, index, &offset_byte, &offset_bit);

	unsigned char mask = 0x3 << (offset_bit);

	if(offset_byte >= index->values_bytes)
		return 0;
//		SUBREADprintf("\nERROR: %u > %u\n", offset_byte, index->values_bytes);

	char reference_base = ((index->values [offset_byte] & mask) >> offset_bit);

	return  (reference_base == base_int_value)?1:0 ;
}

void gvindex_destory(gene_value_index_t * index)
{
	free(index -> values);
}


void gvindex_get_string(char *buf, gene_value_index_t * index, unsigned int pos, int len, int is_negative_strand)
{
	int i;
	if (is_negative_strand)
		for (i=len-1;i>=0;i--)
		{
			buf[i] = gvindex_get (index, pos + len - 1 - i);
			switch(buf[i])
			{
				case 'A': buf[i] = 'T'; break;
				case 'T': buf[i] = 'A'; break;
				case 'G': buf[i] = 'C'; break;
				case 'C': buf[i] = 'G'; break;
			}
		}
	else
		for (i=0;i<len;i++)
			buf[i] = gvindex_get (index, pos +i);
}


