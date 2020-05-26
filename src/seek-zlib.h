#ifndef __SEEK_ZLIB_H_
#define __SEEK_ZLIB_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "subread.h"

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
