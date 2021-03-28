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
  
  
#ifndef _SAMBAM_FILE_H_
#define _SAMBAM_FILE_H_

#include <zlib.h>
#include "HelperFunctions.h"

typedef unsigned char BS_uint_8;
typedef unsigned short BS_uint_16;
typedef unsigned int BS_uint_32;

#define BAM_MAX_CHROMOSOME_NAME_LEN 200 
#define BAM_MAX_CIGAR_LEN (30000)
#define BAM_MAX_READ_NAME_LEN 256
#define BAM_MAX_READ_LEN 3000

#define SAMBAM_FILE_SAM	10
#define SAMBAM_FILE_BAM 20

#define BAM_FILE_STAGE_HEADER 10
#define BAM_FILE_STAGE_ALIGNMENT 20


#define SAMBAM_COMPRESS_LEVEL_FASTEST Z_NO_COMPRESSION 
#define SAMBAM_COMPRESS_LEVEL_NORMAL Z_BEST_SPEED

#define SAMBAM_GZIP_WINDOW_BITS -15
#define SAMBAM_INPUT_STREAM_SIZE 140000

#define TEST_BAD_BAM_CHUNKS 9999925 

typedef struct
{
	char read_name[BAM_MAX_READ_NAME_LEN];
	char * chro_name;
	unsigned int chro_offset;
	unsigned short flags;
	char * mate_chro_name;
	unsigned int mate_chro_offset;
	int templete_length;
	unsigned char mapping_quality;
	int NH_number;

	char cigar[BAM_MAX_CIGAR_LEN];
	char sequence[BAM_MAX_READ_LEN];
	char seq_quality[BAM_MAX_READ_LEN];

	char buff_for_seq[BAM_MAX_READ_LEN*2];

} SamBam_Alignment;


#define SB_FETCH(a)  if((a) -> input_binary_stream_write_ptr - (a) -> input_binary_stream_read_ptr < 3000){int test_rlen_2 = SamBam_fetch_next_chunk(a); if(test_rlen_2 == -2){(a)->is_bam_broken = 1;}}
#define SB_EOF(a)  ((a)-> is_eof && (  (a) -> input_binary_stream_write_ptr <= (a) -> input_binary_stream_read_ptr ))
#define SB_READ(a)  ((a) -> input_binary_stream_buffer + (a) -> input_binary_stream_read_ptr - (a) -> input_binary_stream_buffer_start_ptr)
#define SB_RINC(a, len)   ((a) -> input_binary_stream_read_ptr) += len

typedef struct
{
	FILE * os_file;

	int file_type;
	int bam_file_stage;

	unsigned long long bam_file_next_section_start;
	unsigned long long input_binary_stream_read_ptr;
	unsigned long long input_binary_stream_write_ptr;
	unsigned long long input_binary_stream_buffer_start_ptr;
	unsigned long long header_length;

	SamBam_Reference_Info * bam_chro_table; 
	int bam_chro_table_size;
	SamBam_Alignment aln_buff;

	char * input_binary_stream_buffer;
	int is_eof;
	int is_paired_end;
	int is_bam_broken;
} SamBam_FILE;

struct SamBam_sorted_compressor_st{
	char plain_text[66000];
	char zipped_bin[70000];
	int text_size, bin_size;
	unsigned int CRC32_plain;
	z_stream strm;
	pthread_t thread_stub;
	srInt_64 bam_block_no;
	int last_job_done;
};

typedef struct
{
	FILE * bam_fp;
	FILE * BAI_fp;
	long long current_BAM_pos;
	char tmpf_prefix[MAX_FILE_NAME_LENGTH];
	z_stream output_stream;
	char * chunk_buffer;
	char * compressed_chunk_buffer;
	char * header_plain_text_buffer;
	int header_plain_text_buffer_used;
	int header_plain_text_buffer_max;
	long long chunk_buffer_used;
	long long chunk_buffer_max_size;
	int writer_state;
	int is_internal_error;
	int sort_reads_by_coord;
	int fastest_compression;
	int sorted_batch_id;
	unsigned int crc0;

	int threads;
	z_stream * threads_output_stream;
	char ** threads_chunk_buffer;
	char ** threads_chunk_buffer_compressed;
	long long * threads_chunk_buffer_used;
	long long * threads_chunk_buffer_max_size;
	
	HashTable * chromosome_name_table;
	HashTable * chromosome_id_table;
	HashTable * chromosome_len_table;
	subread_lock_t thread_bam_lock;

	worker_master_mutex_t sorted_notifier;
	HashTable * block_no_p1_to_vpos_tab;
	//int sorted_compress_plain_text_used;
	int sorted_compress_this_thread_no;
	srInt_64 this_bam_block_no;
	struct SamBam_sorted_compressor_st * writer_threads;
} SamBam_Writer;

// This function reads the next BAM section from the bam_fp. The buffer has a variable length but should be at least 64K bytes.
// I recommend you to allocate 80KB of memory.
// This function returns the size of the compressed data ( CDATA ). It returns < 0 if EOF.
int PBam_get_next_zchunk(FILE * bam_fp, char * buffer, int buffer_length, unsigned int * real_len);

// load the header of a BAM file (the header is important to load BAM reads)
// this function puts the File Pointer to the first read chunk in the BAM.
// It returns 0 if finished loading, or non-zero if wrong.
// If the chunk contains read data after the chromosome table, the read data is copied into remainder_read_data, and its lengtb is returned in remainder_read_data_len.
int PBum_load_header(FILE * bam_fp, SamBam_Reference_Info** chro_tab, char * remainder_read_data, int * remainder_read_data_len);


// load a new line from the BAM buffer (chunk) at chunk_ptr.
// if seq_needed==0, then no sequence nor quality str will be loaded.
// it returns the length (without "\0" after the tail) of the SAM string.
int PBam_chunk_gets(char * chunk, int *chunk_ptr, int chunk_limit, SamBam_Reference_Info * bam_chro_table, char * buff , int buff_len, SamBam_Alignment*aln, int seq_needed);

// This function opens a file, either SAM or BAM, in read-only mode.
// The "file_type" parameter specifies which type of file it is: SAMBAM_FILE_BAM or SAMBAM_FILE_SAM.
SamBam_FILE * SamBam_fopen(char * fname , int file_type);

// This function closes any opened file and releases memory footprint. It works just like "fclose()".
void SamBam_fclose(SamBam_FILE * fp);

// This function tells if a file is exhausted.
// Note that a non-exhausted file can still contain no more alignment results.
// Hence, it is recommended to check the return value of SamBam_fgets() to tell if the file has reached its end.
int SamBam_feof(SamBam_FILE * fp);

// This function works like fgets except it decode the BAM file.
// If the buffer is not long enough to store the line, the remainder of this line is omitted and the next call will read the next alignment.
// A very important difference between fgets and SamBam_fgets is that this function returns NULL when there are no more lines.
// It is recommended to use the return value as the indicator of EOF like:
/**
 * SamBam_FILE * fp = SamBam_fopen("my.bam", SAMBAM_FILE_BAM);
 * while(1)
 * {
 *   char buf[3000];
 *   char * ret = SamBam_fgets(fp, buf, 3000);
 *   if(ret) puts(buf);
 *   else break;
 * }
 * SamBam_fclose(fp);
 */
char * SamBam_fgets(SamBam_FILE * fp , char * buff , int buff_len , int seq_needed);

int SamBam_writer_create(SamBam_Writer * writer, char * BAM_fname, int threads, int sort_reads_by_coord, int is_tmp_BAM, char * tmpfname);

int SamBam_writer_close(SamBam_Writer * writer);

int SamBam_writer_add_header(SamBam_Writer * writer, char * header_text, int add_chro);

int SamBam_writer_add_chromosome(SamBam_Writer * writer, char * chro_name, unsigned int chro_length, int add_header_too);

int SamBam_writer_add_read_bin(SamBam_Writer * writer, int thread_no, char * rbin, int committable);

int SamBam_writer_calc_cigar_span(char * bin);

int SamBam_writer_add_read_fqs_scRNA(gzFile * outfp, char * bambin);

int SamBam_writer_add_read(SamBam_Writer * writer, int threadno, char * read_name, unsigned int flags, char * chro_name, unsigned int chro_position, int mapping_quality, char * cigar, char * next_chro_name, unsigned int next_chro_pos, int temp_len, int read_len, char * read_text, char * qual_text, char * additional_columns, int can_submit);

void SamBam_writer_optimize_bins(HashTable *bin_tab, ArrayList *bin_arr, HashTable ** new_tab, ArrayList ** new_arrs);

int is_badBAM(char * fn);

int SamBam_unzip(char * out, int out_max_len , char * in , int inlen, int sync_only);

int SamBam_fetch_next_chunk(SamBam_FILE *fp);

int SamBam_compress_cigar(char * cigar, int * cigar_int, int * ret_coverage, int max_secs);
char cigar_op_char(int ch);
void SamBam_read2bin(char * read_txt, char * read_bin);

int convert_BAM_binary_to_SAM(SamBam_Reference_Info * chro_table, char * bam_bin, char * sam_txt);
int is_paired_end_BAM(char * fn);
void SamBam_writer_finalise_thread(SamBam_Writer * writer, int thread_id);
void SamBam_writer_finish_header( SamBam_Writer * writer );
void SamBam_writer_finalise_one_thread(SamBam_Writer * writer);
int SamBam_writer_add_read_line(SamBam_Writer * writer, int thread_no, char * rline, int commitable);
char *duplicate_TAB_record_field(char * rline, int fld_no, int toend);
#endif
