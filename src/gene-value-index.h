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
  
  
#ifndef __GENE_VALUE_INDEX_
#define __GENE_VALUE_INDEX_

#include "hashtable.h"
#include "subread.h"
#include "sorted-hashtable.h"

int gvindex_init(gene_value_index_t * index, unsigned int start_point, unsigned int all_bases_estm);

void gvindex_set (gene_value_index_t * index, gehash_data_t offset, gehash_key_t base_value, int padding);

int gvindex_dump(gene_value_index_t * index, const char filename []);

int gvindex_load(gene_value_index_t * index, const char filename []);

void gvindex_destory(gene_value_index_t * index);

void gvindex_baseno2offset(unsigned int base_number, gene_value_index_t * index, unsigned int * offset_byte, unsigned int * offset_bit);

// returns a 16-bit bitmap showing if each base is matched.
int gvindex_match(gene_value_index_t * index, gehash_data_t offset, gehash_key_t base_value);

int gvindex_match_base(gene_value_index_t * index, gehash_data_t offset, const char base_int_value);

int gvindex_get(gene_value_index_t * index, gehash_data_t offset);

int match_chro(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int is_negative_strand, int space_type);
float match_chro_support(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int is_negative_strand, int space_type, char * qual_str, int qual_format);


int match_chro_maxerror(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int is_negative_strand, int space_type, int maxerror);

void gvindex_get_string(char *buf, gene_value_index_t * index, unsigned int pos, int len, int is_negative_strand);

int match_chro_wronglen(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int space_type, int * left_match_bases, int * right_match_bases);

int match_indel_chro_to_front(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int * indels, int * indel_point,int max_indel_number, int max_test_len);
int match_indel_chro_to_back(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int * indels, int * indel_point,int max_indel_number, int min_test_offset);

int match_indel_table_to_front(HashTable * indel_table , char * read, gene_value_index_t * index, unsigned int pos, int test_len, short * indels, short * indel_point, int max_indel_number, int min_test_offset, struct explorer_section_t * sec);
int match_indel_table_to_back(HashTable * indel_table , char * read, gene_value_index_t * index, unsigned int pos, int test_len, short * indels, short * indel_point, int max_indel_number, int min_test_offset, struct explorer_section_t * sec);

unsigned int match_chro_range(char * read, gene_value_index_t * index, unsigned int pos, int read_len, int search_length, int search_to_back);
#endif
