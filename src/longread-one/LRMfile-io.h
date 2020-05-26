#ifndef __LRMFILEIO_H_
#define __LRMFILEIO_H_

#define LRMGENE_INPUT_FASTQ 1
#define LRMGENE_INPUT_GZIP_FASTQ 51
#define LRMMAX_LINE_LENGTH (LRMMAX_READ_LENGTH + 10)

#include "LRMconfig.h"

int LRMhash_strcmp(const void * s1, const void * s2);
srUInt_64 LRMhash_strhash(const void * sv);

int LRMgenekey2int(char key []);

int LRMgeinput_open(const char * filename, LRMgene_input_t * input);

void LRMsambam_write_header(LRMcontext_t * context, LRMthread_context_t * thread_context);
void LRMbam_generate_tail_binary(LRMcontext_t * context, LRMthread_context_t * thread_context);

int LRMwrite_chunk_check_buffer_write(LRMcontext_t * context, LRMthread_context_t * thread_context, int force_write);
// return 0 if successful
int LRMfetch_next_read(LRMcontext_t * context, LRMthread_context_t * thread_context, unsigned int *read_len, char * read_name, char * read_text, char * qual_text, unsigned int * read_no_in_chunk);
// Return the length of this read or -1 if EOF.
int LRMgeinput_next_read(LRMgene_input_t * input, char * read_name, char * read_string, char * quality_string);
void LRMgeinput_close(LRMgene_input_t * input);

// returns read length
int LRMgeinput_readline(LRMgene_input_t * input, int buf_len, char * linebuffer) ;
void LRMreverse_read(char * ReadString, int Length);
void LRMreverse_quality(char * QualtyString, int Length);
int LRMload_offsets(LRMcontext_t * context);
int LRMlocate_gene_position(LRMcontext_t * context, unsigned int linear, char ** chro_name, int * pos);

void LRMpos2txt(LRMcontext_t * context, unsigned int linear, char * txt);
int LRMlocate_chro_length(LRMcontext_t * context, unsigned int linear, char ** chro_name, long long * chro_len);
#endif
