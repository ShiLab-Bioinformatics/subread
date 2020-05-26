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
  
  
#ifndef __LRMBASE_INDEX_H_
#define __LRMBASE_INDEX_H_

#include "LRMconfig.h"

int LRMgvindex_load(LRMgene_value_index_t * index, const char filename []);

void LRMgvindex_destory(LRMgene_value_index_t * index);

void LRMgvindex_baseno2offset(unsigned int base_number, LRMgene_value_index_t * index, unsigned int * offset_byte, unsigned int * offset_bit);

int LRMgvindex_get(LRMgene_value_index_t * index, LRMgehash_data_t offset);

int LRMmatch_chro(char * read, LRMgene_value_index_t * index, unsigned int pos, int test_len, int is_negative_strand);

void LRMgvindex_get_string(char *buf, LRMgene_value_index_t * index, unsigned int pos, int len, int is_negative_strand);

int LRMvalidate_mapping(LRMcontext_t * context, char * read, char * cigar, LRMgene_value_index_t * index, unsigned int pos, int rev, int * mapped_length, int show_txt);
#endif
