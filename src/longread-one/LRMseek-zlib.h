#ifndef __SEEK_ZLIB_H_LRM_
#define __SEEK_ZLIB_H_LRM_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define SEEKGZ_ZLIB_WINDOW_SIZE (32*1024)


typedef struct {
	FILE * gz_fp;
	char * current_chunk_txt;
	char * current_chunk_bin;
	z_stream stem;
	int current_chunk_txt_size;
	unsigned int in_pointer;
	unsigned int in_chunk_offset;
	unsigned int in_block_offset;
	//unsigned int txt_buffer_size;
	unsigned int txt_buffer_used;
	unsigned long long block_start_in_file_offset;
	unsigned int block_start_in_file_bits;

	unsigned long long next_block_file_offset;
	unsigned int next_block_file_bits;

	int is_the_last_chunk;
	int internal_error;

	unsigned int dict_window_pointer;
	unsigned int dict_window_used;
	char dict_window[SEEKGZ_ZLIB_WINDOW_SIZE];

	unsigned int block_dict_window_size;
	char block_dict_window[SEEKGZ_ZLIB_WINDOW_SIZE];
} seekable_zfile_t;

typedef struct{
	char dict_window[SEEKGZ_ZLIB_WINDOW_SIZE];
	unsigned long long block_gzfile_offset;
	unsigned int block_gzfile_bits;
	unsigned int block_dict_window_size;

	unsigned int in_block_text_offset;
} seekable_position_t;
#define SEEKZLIBmin(a,b) ( (a)<(b)?(a):(b) )

// returns 0 if OK; returns 1 if the file is not indexable; returns -1 if file doesn't exist.
int LRMseekgz_open(const char * fname, seekable_zfile_t * fp);

// returns length in bytes if OK (length includes the line break at the end); returns 0 if EOF
int LRMseekgz_gets(seekable_zfile_t * fp, char * buf, int buf_size);

void LRMseekgz_tell(seekable_zfile_t * fp, seekable_position_t * pos);

void LRMseekgz_seek(seekable_zfile_t * fp, seekable_position_t * pos);

int LRMseekgz_next_char(seekable_zfile_t * fp);

void LRMseekgz_close(seekable_zfile_t * fp);
#endif
