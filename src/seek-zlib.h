#ifndef __SEEK_ZLIB_H_
#define __SEEK_ZLIB_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "subread.h"

#define PARALLEL_GZIP_TXT_BUFFER_MARGIN (2* MAX_FC_READ_LENGTH + 500)
#define PARALLEL_GZIP_TXT_BUFFER_SIZE (1024*1024)
#define PARALLEL_GZIP_ZIPPED_BUFFER_SIZE (PARALLEL_GZIP_TXT_BUFFER_SIZE *9/8 )

typedef struct{
	int thread_no;
	int in_buffer_used;
	int out_buffer_used;
	unsigned int CRC32;
	unsigned int zipped_CRC32;
	int plain_length;
	char in_buffer[PARALLEL_GZIP_TXT_BUFFER_SIZE];
	char out_buffer[PARALLEL_GZIP_ZIPPED_BUFFER_SIZE];
	z_stream zipper;
} parallel_gzip_writer_thread_t;

typedef struct{
	int threads;
	srInt_64 plain_length;
	unsigned int CRC32;
	FILE * os_file;
	parallel_gzip_writer_thread_t * thread_objs;
} parallel_gzip_writer_t;

void parallel_gzip_writer_init(parallel_gzip_writer_t * pzwtr, char * output_filename, int total_threads);
void parallel_gzip_writer_add_text(parallel_gzip_writer_t * pzwtr, char * text, int tlen, int thread_no);
// because we have to keep sync between three fastq files, the flush function has to be manually called three times at the same time point.
// otherwise R1, I2 and R2 files will have inconsistent read orders.
// the outer program has to check if any of the three in_buffers is full.
void parallel_gzip_zip_texts(parallel_gzip_writer_t * pzwtr, int thread_no, int for_eof_marker);
void parallel_gzip_writer_flush(parallel_gzip_writer_t * pzwtr, int thread_no);
void parallel_gzip_writer_close(parallel_gzip_writer_t * pzwtr);
int parallel_gzip_writer_add_read_fqs_scRNA(parallel_gzip_writer_t**outfps, char * bambin, int thread_no);

// returns 0 if OK; returns 1 if the file is not indexable; returns -1 if file doesn't exist.
int seekgz_open(const char * fname, seekable_zfile_t * fp, FILE * old_fp);

// returns length in bytes if OK (length includes the line break at the end); returns 0 if EOF
int seekgz_gets(seekable_zfile_t * fp, char * buf, int buf_size);

void seekgz_tell(seekable_zfile_t * fp, seekable_position_t * pos);

void seekgz_seek(seekable_zfile_t * fp, seekable_position_t * pos);

// Diff: seekgz_next_char returns EOF for EOF but seekgz_next_int8 returns -1 for EOF
int seekgz_next_char(seekable_zfile_t * fp);
int seekgz_next_int8(seekable_zfile_t * fp);

void seekgz_close(seekable_zfile_t * fp);

// returns length in bytes if OK (length includes the line break at the end); returns 0 if EOF
int autozip_gets(autozip_fp * fp, char * buf, int buf_size);


// return -1 for EOF
int autozip_getch(autozip_fp * fp);

void autozip_close(autozip_fp * fp);

// return -1 if error, return 0 if plain text, return 1 if gzipped 
int autozip_open(const char * fname, autozip_fp * fp);

void autozip_rewind(autozip_fp * fp);

int seekgz_preload_buffer( seekable_zfile_t * fp , subread_lock_t * read_lock);

// returns length in bytes if OK (length includes the line break at the end); returns 0 if EOF
int seekgz_gets(seekable_zfile_t * fp, char * buff, int buff_len);
#endif
