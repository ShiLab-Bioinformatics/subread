#ifndef _INPUT_BLC_H
#define _INPUT_BLC_H

#include <stdio.h>
#include <stdlib.h>

#include "subread.h"
#include "seek-zlib.h"
#include "hashtable.h"

int input_BLC_init( input_BLC_t * blc_input , char * data_dir );
int input_BLC_next_read( input_BLC_t * blc_input, char * readname , char * read, char * qual );
int input_BLC_tell ( input_BLC_t * blc_input , input_BLC_pos_t * pos );
int input_BLC_seek ( input_BLC_t * blc_input , input_BLC_pos_t * pos );
void input_BLC_close (input_BLC_t * blc_input);

// "cached BCL" maintains a chunk of reads; it decompresses 
int cacheBCL_init( cache_BCL_t * cache_input, char * data_dir, int reads_in_chunk, int all_threads );
int cacheBCL_next_read(  cache_BCL_t * cache_input, char * read_name, char * seq, char * qual, srInt_64 * read_number_in_all);
int cacheBCL_go_chunk_start( cache_BCL_t * blc_input );
int cacheBCL_go_chunk_end( cache_BCL_t * blc_input );
void cacheBCL_close ( cache_BCL_t * blc_input );

// it returns a hashtable : sample name -> [ ( "ACGTAATT", 1 ), ( "CGTTATGG", 2 ), ... ]
// The hashtable can be simply destroyed and all contents are deallocated automatically
// It returns NULL if no sample sheet is found
HashTable * input_BLC_parse_SampleSheet(char * fname);

// It returns a list of barcodes (not the barcode conversation table)
// It returns NULL if no list is found
// The arraylist can be simply destroyed and all contents are deallocated automatically
ArrayList * input_BLC_parse_CellBarcodes(char * fname);
int hamming_dist_ATGC_max2(char* s1, char* s2 );

// returns -1 if error, or 0 if no error.
int cacheBCL_quality_test(char * datadir, HashTable * sample_sheet_table, ArrayList * cell_barcode_list, int testing_reads, int * tested_reads, int * valid_sample_index, int * valid_cell_barcode);
#endif
