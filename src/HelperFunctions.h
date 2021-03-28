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
  
  
#ifndef __HELPER_FUNCTIONS_H_
#define __HELPER_FUNCTIONS_H_

#include <pthread.h>
#include "subread.h"
#include "hashtable.h"

#define PARSE_STATUS_TAGNAME 1
#define PARSE_STATUS_TAGTYPE 2
#define PARSE_STATUS_TAGVALUE 3

typedef struct{
	int workers;
	int * mutex_with_master;
	pthread_cond_t *conds_worker_wait;
	pthread_mutex_t *mutexs_worker_wait;
	int all_terminate;
	int * worker_is_working;
} worker_master_mutex_t;


/**** How to run master and worker
 *    Master: prepare_next_job -> wait_for_last_job_done -> collect_results -> notify_worker_do_next -> goto 0 
 *    Worker: wait_for_next_job -> do_next_job -> goto 0
 *    The result is collected from worker while it is "wait_for_next_job": it releases the lock when "cond_wait" is called.
 *    Master can prepare the next job while the worker is running on the last job, hence parallel.
 */

void worker_master_mutex_init(worker_master_mutex_t * wmt, int all_workers);
void worker_thread_start(worker_master_mutex_t * wmt, int worker_id);
void worker_master_mutex_destroy(worker_master_mutex_t * wmt);

// 0 : a job is allocated; 1 : worker should terminate
int worker_wait_for_job(worker_master_mutex_t * wmt, int worker_id);

void master_wait_for_job_done(worker_master_mutex_t * wmt, int worker_id);
// master collects results between master_wait_for_job_done and master_notify_worker
void master_notify_worker(worker_master_mutex_t * wmt, int worker_id);

// this function must be called when the master thread has all the worker locks in control.
// ie all the workers should be in the "wait_for_job" status.
void terminate_workers(worker_master_mutex_t * wmt);

typedef struct{
	HashTable * contig_table;
	HashTable * size_table;
} fasta_contigs_t;

#ifndef MAKE_STANDALONE

typedef struct{
	ArrayList * message_queue;
	int is_thread_mode;
	subread_lock_t queue_lock;
	subread_lock_t queue_notifier;
	int is_finished;
} message_queue_t;

extern message_queue_t mt_message_queue;
#endif

void msgqu_init();
void msgqu_destroy();
void msgqu_main_loop();
void msgqu_notifyFinish();
void msgqu_printf(const char * fmt, ...);

int read_contig_fasta(fasta_contigs_t * tab, char * fname);
int get_contig_fasta(fasta_contigs_t * tab, char * chro, unsigned int pos, int len, char * out_bases);
void destroy_contig_fasta(fasta_contigs_t * tab);

// This function parses CIGAR_Str and extract the relative starting points and lengths of all sections (i.e., the sections of read that are separated by 'N').
// CIGAR_Str is a CIGAR string containing 'S', 'M', 'I', 'D' and 'N' operations. Other operations are all ignored. The length of CIGAR_Str should be always less than 100 bytes or "-1" is returned.
// Staring_Points and Section_Length are empty arrays to write the sections. The minimum length of each array is 6 items.
// The length of a section is its length on the chromosome, namely 'I' is ignored but 'D' is added into the length.
// This function ignores all sections from the 7-th.

// This function returns the number of sections found in the CIGAR string. It returns -1 if the CIGAR string cannot be parsed.

int RSubread_parse_CIGAR_string(char * chro , unsigned int first_pos, const char * CIGAR_Str, int max_M, char ** Section_Chromosomes, unsigned int * Section_Start_Chro_Pos,unsigned short * Section_Start_Read_Pos, unsigned short * Section_Chro_Length, int * is_junction_read);


int RSubread_parse_CIGAR_Extra_string(int FLAG, char * MainChro, unsigned int MainPos, const char * CIGAR_Str, const char * Extra_Tags, int max_M, char ** Chros, unsigned int * Staring_Chro_Points, unsigned short * Section_Start_Read_Pos, unsigned short * Section_Length, int * is_junction_read);

// This function try to find the attribute value of a given attribute name from the extra column string in GTF/GFF.
// If the value is found, it returns the length of the value (must be > 0 by definition), or -1 if no attribute is found or the format is wrong.

int GTF_extra_column_value(const char * Extra_Col, const char * Target_Name, char * Target_Value, int TargVal_Size);


// Replacing `rep' with `with' in `orig'. 
// Rhe return value must be freed if it is not NULL.
char *str_replace(char *orig, char *rep, char *with) ;


// rule: the string is ABC123XXXXXX...
// // This is the priroity:
// // First, compare the letters part.
// // Second, compare the pure numeric part.
// // Third, compare the remainder.
int strcmp_number(char * s1, char * s2);

unsigned int reverse_cigar(unsigned int pos, char * cigar, char * new_cigar);
unsigned int find_left_end_cigar(unsigned int right_pos, char * cigar);
int mac_or_rand_str(char * char_14);

double fast_fisher_test_one_side(unsigned int a, unsigned int b, unsigned int c, unsigned int d, long double * frac_buffer, int buffer_size);
int load_features_annotation(char * file_name, int file_type, char * gene_id_column, char * transcript_id_column, char * used_feature_type,
 void * context, int do_add_feature(char * gene_name, char * transcript_id, char * chrome_name, unsigned int start, unsigned int end, int is_negative_strand, void * context)  );

HashTable * load_alias_table(char * fname) ;

char * get_short_fname(char * lname);

// Rebuild a string containing the command line.
// Return the string length (without the terminating \0)
// You need to free(*lineptr) after all.
int rebuild_command_line(char ** lineptr, int argc, char ** argv);


// Calculate a full round of MD5 or SHA256. 
void Helper_md5sum(char * plain_txt, int plain_len, unsigned char * bin_md5_buff);

typedef unsigned int HelpFuncMD5_u32plus;

typedef struct {
	HelpFuncMD5_u32plus lo, hi;
	HelpFuncMD5_u32plus a, b, c, d;
	unsigned char buffer[64];
	HelpFuncMD5_u32plus block[16];
} HelpFuncMD5_CTX;
 

void HelpFuncMD5_Init(HelpFuncMD5_CTX *ctx);
void HelpFuncMD5_Update(HelpFuncMD5_CTX *ctx, const void *data, unsigned long size);
void HelpFuncMD5_Final(unsigned char *result, HelpFuncMD5_CTX *ctx);


void Helper_sha256sum(char * plain_txt, int plain_len, unsigned char * bin_md5_buff);
unsigned long long plain_txt_to_long_rand(char * plain_txt, int plain_len);

// give me a p, I give you the value such that Pr( x < value ) == p in a 0/1 normal distribution.
double inverse_sample_normal(double p);


// big number functions
// retrived from https://github.com/kokke/tiny-TNbignum-c
// "This is free and unencumbered software released into the public domain."




#include <stdint.h>
#include <assert.h>


/* This macro defines the word size in bytes of the array that constitues the big-number data structure. */
#ifndef WORD_SIZE
  #define WORD_SIZE 4
#endif

/* Size of big-numbers in bytes */
#define BN_MAXIMUM_BITS  4096
#define BN_ARRAY_SIZE    ( BN_MAXIMUM_BITS / 8 / WORD_SIZE )


/* Here comes the compile-time specialization for how large the underlying array size should be. */
/* The choices are 1, 2 and 4 bytes in size with uint32, uint64 for WORD_SIZE==4, as temporary. */
#ifndef WORD_SIZE
  #error Must define WORD_SIZE to be 1, 2, 4
#elif (WORD_SIZE == 1)
  /* Data type of array in structure */
  #define DTYPE                    uint8_t
  /* bitmask for getting MSB */
  #define DTYPE_MSB                ((DTYPE_TMP)(0x80))
  /* Data-type larger than DTYPE, for holding intermediate results of calculations */
  #define DTYPE_TMP                uint32_t
  /* sprintf format string */
  #define SPRINTF_FORMAT_STR       "%.02x"
  #define SSCANF_FORMAT_STR        "%2hhx"
  /* Max value of integer type */
  #define MAX_VAL                  ((DTYPE_TMP)0xFF)
#elif (WORD_SIZE == 2)
  #define DTYPE                    uint16_t
  #define DTYPE_TMP                uint32_t
  #define DTYPE_MSB                ((DTYPE_TMP)(0x8000))
  #define SPRINTF_FORMAT_STR       "%.04x"
  #define SSCANF_FORMAT_STR        "%4hx"
  #define MAX_VAL                  ((DTYPE_TMP)0xFFFF)
#elif (WORD_SIZE == 4)
  #define DTYPE                    uint32_t
  #define DTYPE_TMP                uint64_t
  #define DTYPE_MSB                ((DTYPE_TMP)(0x80000000))
  #define SPRINTF_FORMAT_STR       "%.08x"
  #define SSCANF_FORMAT_STR        "%8x"
  #define MAX_VAL                  ((DTYPE_TMP)0xFFFFFFFF)
#endif
#ifndef DTYPE
  #error DTYPE must be defined to uint8_t, uint16_t uint32_t or whatever
#endif


/* Custom assert macro - easy to disable */
#define require(p, msg) assert(p && #msg)


/* Data-holding structure: array of DTYPEs */
struct bn
{
  DTYPE array[BN_ARRAY_SIZE];
};



/* Tokens returned by TNbignum_cmp() for value comparison */
enum { SMALLER = -1, EQUAL = 0, LARGER = 1 };



/* Initialization functions: */
void TNbignum_init(struct bn* n);
void TNbignum_from_int(struct bn* n, DTYPE_TMP i);
int  TNbignum_to_int(struct bn* n);
void TNbignum_from_string(struct bn* n, char* str, int nbytes);

// warning: maxsize MUST >= 1026
void TNbignum_to_string(struct bn* n, char* str, int maxsize);

/* Basic arithmetic operations: */
void TNbignum_add(struct bn* a, struct bn* b, struct bn* c); /* c = a + b */
void TNbignum_sub(struct bn* a, struct bn* b, struct bn* c); /* c = a - b */
void TNbignum_mul(struct bn* a, struct bn* b, struct bn* c); /* c = a * b */
void TNbignum_div(struct bn* a, struct bn* b, struct bn* c); /* c = a / b */
void TNbignum_mod(struct bn* a, struct bn* b, struct bn* c); /* c = a % b */
void TNbignum_divmod(struct bn* a, struct bn* b, struct bn* c, struct bn* d); /* c = a/b, d = a%b */

/* Bitwise operations: */
void TNbignum_and(struct bn* a, struct bn* b, struct bn* c); /* c = a & b */
void TNbignum_or(struct bn* a, struct bn* b, struct bn* c);  /* c = a | b */
void TNbignum_xor(struct bn* a, struct bn* b, struct bn* c); /* c = a ^ b */
void TNbignum_lshift(struct bn* a, struct bn* b, int nbits); /* b = a << nbits */
void TNbignum_rshift(struct bn* a, struct bn* b, int nbits); /* b = a >> nbits */

/* Special operators and comparison */
int  TNbignum_cmp(struct bn* a, struct bn* b);               /* Compare: returns LARGER, EQUAL or SMALLER */
int  TNbignum_is_zero(struct bn* n);                         /* For comparison with zero */
void TNbignum_inc(struct bn* n);                             /* Increment: add one to n */
void TNbignum_dec(struct bn* n);                             /* Decrement: subtract one from n */
void TNbignum_pow(struct bn* a, struct bn* b, struct bn* c); /* Calculate a^b -- e.g. 2^10 => 1024 */
void TNbignum_isqrt(struct bn* a, struct bn* b);             /* Integer square root -- e.g. isqrt(5) => 2*/
void TNbignum_assign(struct bn* dst, struct bn* src);        /* Copy src into dst -- dst := src */
size_t get_sys_mem_info(char * keyword);
int get_free_total_mem(size_t * total, size_t * free_mem);
void *windows_memmem(const void *haystack_start, size_t haystack_len, const void *needle_start, size_t needle_len);
#endif
