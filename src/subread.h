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
  
  
#ifndef _SUBREAD_H_
#define _SUBREAD_H_

#include <stdlib.h>
#include <pthread.h>
#include <stdio.h>
#include <zlib.h>

#ifndef MAKE_STANDALONE
#ifndef RUNNING_ENV
#include <R.h>
#endif
#endif

#include "hashtable.h" 

#define MAX_SCRNA_FASTQ_FILES 256 
#define SCRNA_FASTA_SPLIT1 "|Rsd:cCounts:mFQs|"
#define SCRNA_FASTA_SPLIT2 "|Rsd:cCounts:1mFQ|"


#define SAM_FLAG_PAIRED_TASK	0x01
#define SAM_FLAG_FIRST_READ_IN_PAIR 0x40
#define SAM_FLAG_SECOND_READ_IN_PAIR 0x80
#define SAM_FLAG_MATE_UNMATCHED 0x08
#define SAM_FLAG_MATCHED_IN_PAIR 0x02
#define SAM_FLAG_REVERSE_STRAND_MATCHED 0x10
#define SAM_FLAG_MATE_REVERSE_STRAND_MATCHED 0x20
#define SAM_FLAG_SECONDARY_MAPPING 0x100
#define SAM_FLAG_DUPLICATE 0x400
#define SAM_FLAG_UNMAPPED 0x04

#define SUBREAD_MAX_ULONGLONG 0xffffffffffffffffllu
#define SUBREAD_MAX_LONGLONG 0x7fffffffffffffffll

#define FUSION_BREAK_POINT	2
#define FUSION_JUNCTION		1
#define SPLICING_JUNCTION	0

#define RUN_ALIGN 		0
#define RUN_FINAL 		1



#define MAX_THREADS 40
#define FC_MAX_THREADS 64
#define MAX_EVENTS_IN_READ 8

//#warning "============== REMOVE '* 15' FROM THE NEXT LINE ================"
#define MAX_READ_LENGTH ( 1210 )
#define MAX_READ_NAME_LEN 200 
#define MAX_CHROMOSOME_NAME_LEN 200 
#define MAX_FILE_NAME_LENGTH (1000)
#define FEATURE_NAME_LENGTH 256 
#define INPUT_BLC_MAX_READS 20
#define MAX_BARCODE_LEN 32

//#warning "============== REMOVE '*1.2' FROM THE NEXT LINE ================"
#define MULTI_THREAD_OUTPUT_ITEMS  (4096 * 3/5 *3)
#define EXON_LONG_READ_LENGTH 160 
#define EXON_MAX_CIGAR_LEN 256 
#define FC_CIGAR_PARSER_ITEMS 11 
#define FC_LONG_READ_RECORD_HARDLIMIT (8*1024*1024)

#define MAX_INDEL_SECTIONS 7
//#define XBIG_MARGIN_RECORD_SIZE 24 
#define MAX_INSERTION_LENGTH 200
#define MAX_DELETION_LENGTH 1000
//#define BASE_BLOCK_LENGTH 15000000
//#define NEED_SUBREAD_STATISTIC


#define IS_MIN_POS_NEGATIVE_STRAND 4
#define IS_MAX_POS_NEGATIVE_STRAND 12 
#define IS_PAIRED_HINTED 16
#define IS_R1_CLOSE_TO_5 1
#define IS_REVERSED_HALVES 2
#define	IS_PROCESSED_READ 32
#define	IS_PROCESSED_READ_R2 64
#define IS_PAIRED_MATCH 128
#define IS_NEGATIVE_STRAND_R1 256 
#define IS_NEGATIVE_STRAND_R2 512 
#define IS_FUSION 1024 
#define IS_NEGATIVE_STRAND 2048
#define IS_RECOVERED_JUNCTION_READ 4096
#define IS_FINALISED_PROCESSING 8192
#define IS_RECOVERED_JUNCTION_READ_STEP4 (8192*2)
#define	IS_BREAKEVEN_READ (8192*4)
#define IS_R1R2_EQUAL_LEN 1024

#define USE_POSIX_MUTEX_LOCK

#if defined(MACOS) || defined(FREEBSD) || defined(USE_POSIX_MUTEX_LOCK)
typedef pthread_mutex_t subread_lock_t;
#define pthread_spinlock_t pthread_mutex_t
#define pthread_spin_lock pthread_mutex_lock
#define pthread_spin_unlock pthread_mutex_unlock
#define pthread_spin_init(a, b) pthread_mutex_init(a, NULL)
#define pthread_spin_destroy(a) pthread_mutex_destroy(a) 
#define strnlen(a,l) strlen(a)
#else
typedef pthread_spinlock_t subread_lock_t;
#endif

#ifndef SRINT_64_DEFINED
#define SRINT_64_DEFINED
typedef long long srInt_64;
typedef unsigned long long srUInt_64;
#endif

#ifdef __MINGW32__
#define ftello ftello64
#define fseeko fseeko64
#endif

#if defined(MAKE_STANDALONE) || defined(RUNNING_ENV)
#define STANDALONE_exit(i) exit(i)
#define SUBREADprintf(...) fprintf(stderr, __VA_ARGS__)
#define SUBREADputs(x) fprintf(stderr, "%s\n", x)
#define SUBREADputchar(x) fputc(x, stderr)
#define SUBREADfflush(x) fflush(x)
#define CORE_SOFT_BR_CHAR '\n'
#else
#define STANDALONE_exit(i) return i;

int safeRprintf(char *fmt, ...);
void msgqu_printf(const char * fmt, ...);

#define SUBREADprintf  msgqu_printf
#define SUBREADputs(x) msgqu_printf("%s\n",(x))
#define SUBREADputchar(X) msgqu_printf("%c",(X)) 
#define SUBREADfflush(X) 
#define CORE_SOFT_BR_CHAR '\n'

#endif

#ifndef NONONO_DONOTDEF

#define QUALITY_KILL	198
#define QUALITY_KILL_SUBREAD	160
#define MAX_QUALITY_TO_CALL_JUNCTION 2195
#define MAX_QUALITY_TO_EXPLORER_JUNCTION 209

#else

#define TEST_TARGET ""

#endif

#define SNP_CALLING_ONLY_HIGHQUAL 1

#define MESSAGE_OUT_OF_MEMORY "Out of memory. If you are using Rsubread in R, please save your working environment and restart R. \n"
#define fatal_memory_size(a) SUBREADputs(MESSAGE_OUT_OF_MEMORY);

//#define QUALITY_KILL	175
//#define QUALITY_KILL_SUBREAD	150


typedef long long subread_read_number_t;
typedef unsigned int gehash_key_t;
typedef unsigned int gehash_data_t;
//typedef float gene_quality_score_t;
typedef int gene_quality_score_t;
//typedef unsigned char gene_vote_number_t;
typedef short gene_vote_number_t;


#define XOFFSET_TABLE_SIZE 250000

#define ANCHORS_NUMBER 259
#define MAX_ALIGNMENT_PER_ANCHOR 2

#define BEXT_RESULT_LIMIT 16

#define SEARCH_BACK 0
#define SEARCH_FRONT 1

//#define LARGE_GENE_VOTE_TABLE
#ifdef LARGE_GENE_VOTE_TABLE
#warning "Using LARGE_GENE_VOTE_TABLE"
#define GENE_VOTE_SPACE 173 
#define GENE_VOTE_TABLE_SIZE 331
#else
#define GENE_VOTE_SPACE 24
#define GENE_VOTE_TABLE_SIZE 30 
#endif

#define MAX_ANNOTATION_EXONS 30000 
#define MAX_EXONS_PER_GENE 400 
#define MAX_EXON_CONNECTIONS 10

#define MAX_GENE_NAME_LEN 128
#define MAX_INDEL_TOLERANCE 7

#define SUBINDEX_VER0 100
#define SUBINDEX_VER1 200
#define SUBINDEX_VER2 201

#define SUBREAD_INDEX_OPTION_INDEX_GAP 0x0101
#define SUBREAD_INDEX_OPTION_INDEX_PADDING 0x0102


#define CHAR_ESC 27

//#define base2int(c) ((c)=='A'?0:((c)=='T'?3:((c)=='C'?2:1)))
#define base2int(c) ((c)<'G'?((c)=='A'?0:2):((c)=='G'?1:3))

/*
#define base2int(c) (("\x3\x3\x3\x3\x3\x3\x3\x3" "\x3\x3\x3\x3\x3\x3\x3\x3"  "\x3\x3\x3\x3\x3\x3\x3\x3"  "\x3\x3\x3\x3\x3\x3\x3\x3"  "\x3\x3\x3\x3\x3\x3\x3\x3"  "\x3\x3\x3\x3\x3\x3\x3\x3"  "\x3\x3\x3\x3\x3\x3\x3\x3"  "\x3\x3\x3\x3\x3\x3\x3\x3"    "\x3"\
	 "\x0\x3\x2\x3\x3\x3\x1\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3")[(int)(c)])
         // A  B  C  D  E  F  G
*/

//#define int2base(c) ((c)==1?'G':((c)==0?'A':((c)==2?'C':'T')))
//#define int2base(c) ("AGCT"[(c)]) 
#define int2base(c) (1413695297 >> (8*(c))&0xff)
#define color2int(c) ((c) - '0')
#define int2color(c) ("0123"[(c)])
#define remove_backslash(str) { int xxxa=0; while(str[xxxa]){ if(str[xxxa]=='/'){str[xxxa]='\0'; break;} xxxa++;} /* SUBREADprintf("BSRRR %s\n", str);*/ }

/*
#define get_base_error_prob64(a) ((a) < '@'-1?1:pow(10., -0.1*((a)-'@')))
#define get_base_error_prob33(a) ((a) < '!'-1?1:pow(10., -0.1*((a)-'!'))) 

*/
#define SUBREAD_malloc(a) malloc(a)

#define FASTQ_PHRED33 1
#define FASTQ_PHRED64 0

#define IS_DEBUG 0


typedef struct {
  char gene_name [MAX_GENE_NAME_LEN]; 
  // The chromosome name is not stored in this data structure
  // All coordinates are translated into the linear location in the entire referenced genome, usually 0 ~ 3.2G 
  unsigned int start_offset;
  unsigned int end_offset;

  // All exons are marked with the linear location in the entire referenced genome, usually 0 ~ 3.2G
  // This marks the end of the list: exon_ends [total_number_of_exons] = 0
  // It shouldn't be equal to 0, should it be?

  unsigned int exon_starts [MAX_EXONS_PER_GENE];
  unsigned int exon_ends [MAX_EXONS_PER_GENE];
} gene_t;



struct gehash_bucket {
	int current_items;
	int space_size;
	union
	{
		short * new_item_keys;
		gehash_key_t * item_keys;
	};
	gehash_data_t * item_values;
};


#ifdef __MINGW32__
#define GEHASH_MEM_PTR_NO (64) 
#else
#define GEHASH_MEM_PTR_NO (64*1024) 
#endif
typedef struct {
	int version_number;
	unsigned long long int current_items;
	int buckets_number;
	char is_small_table;
	struct gehash_bucket * buckets;
	int index_gap;
	int padding;
	char * malloc_ptr [GEHASH_MEM_PTR_NO];
	int free_item_only;
} gehash_t;


typedef struct{
	unsigned int memory_block_size;
	unsigned int start_base_offset;
	unsigned int start_point;
	unsigned int length;
	unsigned char * values;
	unsigned int values_bytes;
	void * appendix1;
	void * appendix2;
} gene_value_index_t;


typedef struct {
	gene_vote_number_t max_vote;
	gehash_data_t max_position;
	gene_quality_score_t max_quality;
	gene_vote_number_t max_indel_recorder[MAX_INDEL_TOLERANCE*3];
	gene_vote_number_t * max_tmp_indel_recorder;
	int max_mask;
	gene_vote_number_t noninformative_subreads;

	unsigned short items[GENE_VOTE_TABLE_SIZE];
	unsigned int pos [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	int masks [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	int marked_shift_indel[GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	gene_vote_number_t votes [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	gene_quality_score_t quality [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	gene_vote_number_t last_subread_cluster [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	gene_vote_number_t indel_recorder [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE][MAX_INDEL_TOLERANCE*3];
	char current_indel_cursor[GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	char toli[GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];

	#ifdef MAKE_FOR_EXON
	short coverage_start [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	short coverage_end [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	short max_coverage_start;
	short max_coverage_end;
	//#warning Switch "MAKE_FOR_EXON" is turned on. It may cost more time. Do not turn it on unless you want to detect junction reads.
	#endif
} gene_vote_t ;

typedef struct{
	unsigned int pos;
	char len;
} indel_record_t;

typedef struct{
	int count;
	int support;
} indel_result_t;

typedef struct{
	unsigned char best_len;
	unsigned int offsets [BEXT_RESULT_LIMIT];
	unsigned char is_reverse [BEXT_RESULT_LIMIT];
} gene_best_record_t;



typedef struct{
	int max_len;
	unsigned int * max_positions;
	unsigned char * is_counterpart;
	gene_vote_number_t * max_votes;
	gene_quality_score_t * max_quality;
	gene_quality_score_t * max_final_quality;
	short * masks;
	char * max_indel_recorder;
	char * span_coverage;
#ifdef REPORT_ALL_THE_BEST
	gene_best_record_t * best_records;
#endif
	char max_indel_tolerance;
	short indel_recorder_length;

	unsigned char *repeated_regions;

} gene_allvote_t;


typedef struct{
	int total_offsets;
        char *read_names;
        unsigned int *read_offsets;
	HashTable * read_name_to_index;
	int padding;
} gene_offset_t;


#define EXON_BUFFER_SIZE 3000 

struct thread_input_buffer {
	char read_names [EXON_BUFFER_SIZE][121];
	char read [EXON_BUFFER_SIZE][1201];
	char quality [EXON_BUFFER_SIZE][1201];
	int rl[EXON_BUFFER_SIZE];
	int write_pointer;
	int read_pointer;

	unsigned int read_id[EXON_BUFFER_SIZE];

};

#define SEEKGZ_CHAIN_BLOCKS_NO 15
#define SEEKGZ_ZLIB_WINDOW_SIZE (32*1024)
#define PAIRER_GZIP_WINDOW_BITS -15
#define PAIRER_DEFAULT_MEM_LEVEL 8

typedef struct {
	unsigned long long block_start_in_file_offset;
	unsigned int block_start_in_file_bits;

	char block_dict_window[SEEKGZ_ZLIB_WINDOW_SIZE]; // copied from the rolling window before this block is decompressed.
	unsigned int block_dict_window_size;

	char * block_txt;
	unsigned int * linebreak_positions;
	int linebreaks;
	int block_txt_size;
} seekable_decompressed_block_t;

typedef struct {
	FILE * gz_fp;
	z_stream stem;
	char * in_zipped_buffer;
	unsigned int in_zipped_buff_read_ptr;

	unsigned int current_block_txt_read_ptr;
	int blocks_in_chain;
	int has_multi_thread_accessed;
	int block_chain_current_no;
	seekable_decompressed_block_t block_rolling_chain[SEEKGZ_CHAIN_BLOCKS_NO];
	
	int internal_error;
	subread_lock_t write_lock;

	unsigned int rolling_dict_window_used;
	char rolling_dict_window[SEEKGZ_ZLIB_WINDOW_SIZE];

	unsigned long long next_block_file_offset; // for the next block after ALL blocks in the chain
	unsigned int next_block_file_bits;  // for the next block after ALL blocks in the chain
} seekable_zfile_t;

typedef struct{
	char dict_window[SEEKGZ_ZLIB_WINDOW_SIZE];
	unsigned long long block_gzfile_offset;
	unsigned int block_gzfile_bits;
	unsigned int block_dict_window_size;

	unsigned int in_block_text_offset;
} seekable_position_t;

typedef struct {
	unsigned long long read_number;
	int lane_id;

	union{
		seekable_position_t ** pos_of_bclgzs;
		unsigned long long * pos_of_bcls;
	};
	union{
		seekable_position_t * pos_of_filtergz;
		unsigned long long pos_of_filter;
	};
	int is_EOF;
} input_BLC_pos_t;

typedef struct {
	unsigned long long read_number;
	int total_bases_in_each_cluster;
	int single_read_lengths[INPUT_BLC_MAX_READS+1];
	int single_read_is_index[INPUT_BLC_MAX_READS];
	int current_lane, bcl_is_gzipped, filter_is_gzipped;
	char bcl_format_string[MAX_FILE_NAME_LENGTH];
	char filter_format_string[MAX_FILE_NAME_LENGTH];
	union{
		seekable_zfile_t ** bcl_gzip_fps; 
		FILE ** bcl_fps; 
	};
	union{
		seekable_zfile_t *  filter_gzip_fp;
		FILE *  filter_fp;
	};
	subread_lock_t read_lock;
	int is_EOF;
} input_BLC_t;

typedef struct{
	union{
		srInt_64 pos_file1, pos_file2, pos_file3;
		seekable_position_t zpos_file1;
	};
	seekable_position_t zpos_file2;
	seekable_position_t zpos_file3;
	int current_file_no;
	srInt_64 current_read_no;
} input_mFQ_pos_t;

typedef struct {
	char filename[MAX_FILE_NAME_LENGTH+1];

	int is_plain;
	FILE * plain_fp;
	seekable_zfile_t gz_fp;
	int is_first_chars;
	unsigned char first_chars[2];
} autozip_fp;



typedef struct {
	char ** files1;
	char ** files2;
	char ** files3;
	int total_files;
	int current_file_no;
	int current_guessed_lane_no;
	srInt_64 current_read_no;
	autozip_fp autofp1;
	autozip_fp autofp2;
	autozip_fp autofp3;
} input_mFQ_t;


typedef struct
{
        char chro_name[MAX_CHROMOSOME_NAME_LEN];
        unsigned int chro_length;
} SamBam_Reference_Info;


typedef struct{
	int current_BAM_file_no;
	srInt_64 section_start_pos;
	int in_section_offset;
	srInt_64 current_read_no;
} input_scBAM_pos_t;

typedef struct {
	FILE * os_file;
	char * BAM_file_names[MAX_SCRNA_FASTQ_FILES];
	char section_buff[66000];
	char align_buff[FC_LONG_READ_RECORD_HARDLIMIT];
	int current_BAM_file_no;
	int total_BAM_files;
	int in_section_offset;
	int section_bin_bytes;
	int chro_table_size;
	SamBam_Reference_Info * chro_table;
	srInt_64 section_start_pos;
	srInt_64 current_read_no;
	subread_lock_t read_lock;
} input_scBAM_t;

typedef struct {
	int read_no_in_chunk;
	int reads_available_in_chunk; // -1 : EOF of all input reads : no next chunk available.
								  // This can be set to -1 only when calling cacheBCL_netx_chunk().
	int chunk_no;
	int chunk_start_lane;
	int chunk_end_lane;
	int reads_per_chunk;
	int last_chunk_in_cache;
	int total_bases_in_each_cluster;
	int single_read_lengths[INPUT_BLC_MAX_READS+1];
	int single_read_is_index[INPUT_BLC_MAX_READS];
	int current_lane, bcl_is_gzipped, filter_is_gzipped;
	int all_threads;
	char bcl_format_string[MAX_FILE_NAME_LENGTH];
	char filter_format_string[MAX_FILE_NAME_LENGTH];
	int bcl_no_is_used[MAX_READ_LENGTH];
	autozip_fp * bcl_gzip_fps; 
	autozip_fp   filter_fp; 
	subread_lock_t read_lock;
	char ** bcl_bin_cache;
	int  flt_bin_cache_size;
	char *  flt_bin_cache;
	char * lane_no_in_chunk;
	int is_EOF;
} cache_BCL_t;

typedef struct{
	union{
		unsigned long long simple_file_position;
		seekable_position_t seekable_gzip_position;
		input_BLC_pos_t BCL_position;
		input_mFQ_pos_t mFQ_position;
		input_scBAM_pos_t scBAM_position;
	};
	char gzfa_last_name[MAX_READ_NAME_LEN];
} gene_inputfile_position_t;


typedef struct{
	unsigned int small_key;
	unsigned int big_key;
} paired_exon_key;

typedef struct{
	unsigned int supporting_reads;
	char is_fusion;
	char big_pos_neg;
	char small_pos_neg;
} fusion_record;

double miltime();

typedef struct{
	char chromosome_name[MAX_CHROMOSOME_NAME_LEN];
	unsigned long known_length;
} chromosome_t;

typedef struct{
	unsigned char record_type;
	unsigned char  mapping_quality;
	unsigned short read_pos;
	unsigned short read_len;
	unsigned short flags;
	unsigned int read_number;
	unsigned int pos;
	unsigned short mapped_segment_in_read;
	char strand;	// 0 = positive, 1 = negative
} base_block_temp_read_t;

typedef struct{
	unsigned char record_type;
	unsigned int pos;
	short type;	
} VCF_temp_read_t;

typedef struct {
	char filename [MAX_FILE_NAME_LENGTH * MAX_SCRNA_FASTQ_FILES * 3];
	int space_type ;
	int file_type ;
	void * input_fp;   // can be system (FILE * sam or fastq or fasta), (seekable_zfile_t *)
	char gzfa_last_name[MAX_READ_NAME_LEN];
	unsigned long long read_chunk_start;
	union{
		cache_BCL_t bcl_input;
		input_mFQ_t scRNA_fq_input;
		input_scBAM_t scBAM_input;
	};
} gene_input_t;



struct explorer_section_t
{
	unsigned int start_pos;

	short read_pos_start;
	short read_pos_end;

	char is_neg_strand;
	short indels;
	short indel_pos;

	short all_indel_poses[10];
	short all_indels[10];
};

struct explorer_record_t
{
	struct explorer_section_t cigar_record[6];
	short b_search_tail;
};

FILE * f_subr_open(const char * fname, const char * mode);
void myrand_srand(unsigned long long seed);
int myrand_rand();
#define abs(a) 	  ((a)>=0?(a):-(a))
#define max(a,b)  ((a)<(b)?(b):(a))
#define min(a,b)  ((a)>(b)?(b):(a))


#endif
