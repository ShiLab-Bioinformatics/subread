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
#include "LRMconfig.h"
#include "LRMbase-index.h"
#include "LRMfile-io.h"


#define LRMgvindex_baseno2offset_m(base_number, index, offset_byte, offset_bit)	{offset_byte =  (base_number - index -> start_base_offset) >>2; offset_bit = base_number % 4 * 2;}


void LRMgvindex_baseno2offset(unsigned int base_number, LRMgene_value_index_t * index, unsigned int * offset_byte, unsigned int * offset_bit)
{
	// the base number corrsponding to the 0-th bit in the whole value array;

	unsigned int offset = (base_number - index -> start_base_offset);

	* offset_byte = offset >>2 ;
	* offset_bit = base_number % 4 * 2;
}

// return 'A', 'G', 'T' and 'C'
int LRMgvindex_get(LRMgene_value_index_t * index, LRMgehash_data_t offset)
{
	unsigned int offset_byte, offset_bit;
	LRMgvindex_baseno2offset_m(offset, index , offset_byte, offset_bit);

	if(offset_byte >= index-> values_bytes -1)return 'N';

	unsigned int one_base_value = (index->values [offset_byte]) >> (offset_bit);


	//LRMprintf("RECV_BASE=%d (%d -  %d)\n",one_base_value & 3, offset_byte , offset_bit);

	return LRMint2base(one_base_value & 3);
}

int LRMgvindex_match(LRMgene_value_index_t * index, LRMgehash_data_t offset, LRMgehash_key_t base_values)
{
	unsigned int offset_byte, offset_bit;

	LRMgvindex_baseno2offset_m(offset, index , offset_byte, offset_bit);
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

int LRMgvindex_load(LRMgene_value_index_t * index, const char filename [])
{
	FILE * fp = fopen(filename, "rb");
	int read_length;
	read_length = fread(&index->start_point,4,1, fp);
	if(read_length<1){
		LRMprintf("ERROR: the array index is incomplete : %d", read_length );
		return 1;
	}
	read_length = fread(&index->length,4,1, fp);
	if(read_length<1){
		LRMprintf("Bad index\n");
		return 1;
	}
	//LRMprintf ("\nBINDEX %s : %u ~ +%u\n",filename, index->start_point, index->length );

	unsigned int useful_bytes, useful_bits;
	index -> start_base_offset = index -> start_point - index -> start_point%4;
	LRMgvindex_baseno2offset (index -> length+ index -> start_point, index ,&useful_bytes,&useful_bits);
	index -> values = malloc(useful_bytes+1);
	index -> values_bytes = useful_bytes+1;
	if(!index->values)
	{
		LRMprintf("Out of memory\n");
		return 1;
	}
	

	read_length =fread(index->values, 1, useful_bytes+1, fp);
	if(read_length < useful_bytes){
		LRMprintf("ERROR: the array index is incomplete : %d < %d.", read_length, useful_bytes+1 );
		return 1;
	}

	fclose(fp);
	return 0;

}

void LRMgvindex_get_string(char *buf, LRMgene_value_index_t * index, unsigned int pos, int len, int is_negative_strand){
	int i;
	if (is_negative_strand)
		for (i=len-1;i>=0;i--)
		{
			buf[i] = LRMgvindex_get (index, pos + len - 1 - i);
			switch(buf[i])
			{
				case 'A': buf[i] = 'T'; break;
				case 'G': buf[i] = 'C'; break;
				case 'C': buf[i] = 'G'; break;
				default: buf[i] = 'A';
			}
		}
	else
		for (i=0;i<len;i++)
			buf[i] = LRMgvindex_get (index, pos +i);
}


int LRMvalidate_mapping(LRMcontext_t * context, char * read, char * cigar, LRMgene_value_index_t * index, unsigned int pos, int neg, int * maplen, int show_txt){
	unsigned int chro_cursor = pos;
	int read_chrsor = 0, all_matched = 0, all_mismatched = 0;
	int tmpi = 0,cigar_i, nch, tmpi_sign = 1, x1;

	if(neg) LRMreverse_read(read, strlen(read));

	if(show_txt){
		char postxt[100];
		LRMpos2txt(context, chro_cursor, postxt);
		LRMprintf("Starting Pos : Read + %d ( %s )\n", read_chrsor, postxt);
	}

	for(cigar_i = 0; (nch = cigar[cigar_i])!=0; cigar_i++){
		if(nch == '-') tmpi_sign = -1;
		else if(nch >='0' && nch <='9'){
			tmpi = tmpi * 10 + (nch - '0');
		}else{
			tmpi *= tmpi_sign;
			if(nch == 'M'){
				int this_matched = LRMmatch_chro(read + read_chrsor, index, chro_cursor, tmpi, 0);
				if(show_txt){
					unsigned txt_chro_cursor = chro_cursor;
					int txt_read_chrsor = read_chrsor;
					for(x1 = 0; x1 < tmpi ; x1++){
						int knownval = LRMgvindex_get(index, txt_chro_cursor);
						int readval = read[txt_read_chrsor];
						txt_chro_cursor ++;
						txt_read_chrsor ++;

						LRMprintf("%c[3%dm%c", 27, knownval == readval ? 7:1, readval);
					}
				}

				if(0 && abs(tmpi) > 22 && this_matched < tmpi * 0.6){
					LRMprintf("Too many mismatched (%d%c) : %d / %d : read + %d\n", tmpi, nch, this_matched, tmpi, read_chrsor);
				}

				all_matched += this_matched;
				all_mismatched += ( tmpi - this_matched );
				(*maplen) += tmpi;
			}	

			if(nch == 'M' || nch == 'I' || nch == 'S'){
				if(nch == 'I' && show_txt) for(x1 = 0; x1 < tmpi ; x1++) LRMprintf("%c[32m%c", 27, read[ read_chrsor + x1 ]);
				if(nch == 'S' && show_txt) for(x1 = 0; x1 < tmpi ; x1++) LRMprintf("%c[4m%c%c[0m", 27, read[ read_chrsor + x1 ], 27);
				read_chrsor += tmpi;
			}
			if(nch == 'M' || nch == 'D' || nch == 'N' || nch == 'S'){
				if(( nch == 'N' || nch == 'D' ) && show_txt) LRMprintf("%c[36m//%c[37m", 27, 27);
				if(nch != 'S')chro_cursor += tmpi;
			}

			tmpi = 0;
			tmpi_sign = 1;
		}
	}
	if(show_txt){
		char postxt[100];
		LRMpos2txt(context, chro_cursor, postxt);
		LRMprintf("%c[37m\n", 27);
		LRMprintf("Ending Pos : Read + %d ( %s )\n", read_chrsor, postxt);
	}
	if(neg) LRMreverse_read(read, strlen(read));
	return all_matched;
}

int LRMmatch_chro(char * read, LRMgene_value_index_t * index, unsigned int pos, int test_len, int is_negative_strand){
	int ret = 0;
	int i;
	
	if ((unsigned int)(pos + test_len) >= index -> length + index -> start_point) return 0;
	if (pos > 0xffff0000) return 0;

	if (is_negative_strand)
	{
		for (i=test_len -1;i>=0;i--)
		{
			char tt = LRMgvindex_get (index, pos+test_len-1-i);
			switch(tt)
			{
				case 'A': ret += read[i] == 'T'; break;
				case 'T': ret += read[i] == 'A'; break;
				case 'G': ret += read[i] == 'C'; break;
				case 'C': ret += read[i] == 'G'; break;
			}
		}
	}else{
		unsigned int offset_byte, offset_bit;

		LRMgvindex_baseno2offset_m(pos, index , offset_byte, offset_bit);

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
	}
	return ret;
}

void LRMgvindex_destory(LRMgene_value_index_t * index)
{
	free(index -> values);
}
