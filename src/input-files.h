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
  
  
#ifndef __INPUT_FILES_H_
#define __INPUT_FILES_H_

#include "subread.h"
#include "hashtable.h"

#define GENE_SPACE_BASE 1
#define GENE_SPACE_COLOR 2

#define GENE_INPUT_PLAIN 0
#define GENE_INPUT_FASTQ 1
#define GENE_INPUT_FASTA 2
#define GENE_INPUT_BCL   3
#define GENE_INPUT_SCRNA_FASTQ 4
#define GENE_INPUT_SCRNA_BAM 5
#define GENE_INPUT_GZIP_FASTQ 51
#define GENE_INPUT_GZIP_FASTA 52

#define GENE_INPUT_SAM_SINGLE   93
#define GENE_INPUT_SAM_PAIR_1   94
#define GENE_INPUT_SAM_PAIR_2   95


#define MAX_LINE_LENGTH 3000
#define MIN_FILE_POINTERS_ALLOWED 50

#define FILE_TYPE_SAM     50
#define FILE_TYPE_BAM     500
#define FILE_TYPE_FAST_   100
#define FILE_TYPE_FASTQ   105
#define FILE_TYPE_FASTA   110
#define FILE_TYPE_GZIP_FAST_   1000
#define FILE_TYPE_GZIP_FASTQ   1105
#define FILE_TYPE_GZIP_FASTA   1110
#define FILE_TYPE_UNKNOWN 999
#define FILE_TYPE_EMPTY   999990
#define FILE_TYPE_NONEXIST 999999
#define FILE_TYPE_RSUBREAD 10
#define FILE_TYPE_GTF 100

#define FEATURECOUNTS_BUFFER_SIZE (1024*1024*12)



#include <stdlib.h>
#include <stdio.h>
#include "core-indel.h"
#include "hashtable.h"


#define SAM_SORT_BLOCKS 229
#define SAM_SORT_BLOCK_SIZE 512333303LLU
//#define SAM_SORT_BLOCK_SIZE 11123333LLU
//
typedef struct {
	srInt_64 output_file_size;
	srInt_64 current_chunk_size;
	unsigned int current_chunk;
	srInt_64 written_reads;
	srInt_64 unpaired_reads;

	FILE * current_block_fp_array [SAM_SORT_BLOCKS];
	FILE * all_chunks_header_fp;

	FILE * out_fp;
	char tmp_path[MAX_FILE_NAME_LENGTH];
} SAM_sort_writer;


typedef struct {
	int thread_id;

	char * input_buff_SBAM;
	int input_buff_SBAM_used;
	int input_buff_SBAM_ptr;
	int reads_in_SBAM;
	subread_lock_t SBAM_lock;

	srInt_64 input_buff_SBAM_file_start;
	srInt_64 input_buff_SBAM_file_end;

	unsigned int chunk_number;
	unsigned int readno_in_chunk;

	unsigned char * input_buff_BIN;
	int input_buff_BIN_used;
	int input_buff_BIN_ptr;
	int input_buff_BIN_capacity;
	int orphant_block_no;
	int need_find_start;
	srInt_64 orphant_space;
	z_stream strm;

	char immediate_last_read_bin[FC_LONG_READ_RECORD_HARDLIMIT];
	char immediate_last_read_full_name[MAX_READ_NAME_LEN*2 +80 ];
	int immediate_last_read_flags;
	int immediate_last_read_bin_len;
	int immediate_last_read_name_len;

	HashTable * orphant_table;
	pthread_t thread_stab;
} SAM_pairer_thread_t;

typedef struct {
	FILE * input_fp;
	int input_is_BAM;
	int tiny_mode;
	int display_progress;
	int is_bad_format;
	int is_single_end_mode;
	int force_do_not_sort;
	int long_cigar_mode;
	int long_read_minimum_length;
	int is_finished;
	int merge_level_finished;
	int max_file_open_number;
	subread_lock_t input_fp_lock;
	subread_lock_t SAM_BAM_table_lock;
	subread_lock_t unsorted_notification_lock;

	srInt_64 total_input_reads;
	srInt_64 total_orphan_reads;

	HashTable * unsorted_notification_table;
	HashTable * sam_contig_number_table;
	HashTable * bam_margin_table;

	int total_threads;
	int input_chunk_no;
	int input_buff_SBAM_size;
	int input_buff_BIN_size;
	char tmp_file_prefix[MAX_FILE_NAME_LENGTH+1];
	char in_file_name[MAX_FILE_NAME_LENGTH+1];

	SAM_pairer_thread_t * threads;
	int BAM_header_parsed;
	unsigned int BAM_l_text;
	unsigned int BAM_n_ref;
	int is_unsorted_notified;
	int is_incomplete_BAM;
	int need_read_group_tag;
	int is_internal_error;
	int is_final_run;

	void (* reset_output_function) (void * pairer);
	int (* output_function) (void * pairer, int thread_no, char * bin1, char * bin2); 
	int (* output_header) (void * pairer, int thread_no, int is_text, unsigned int items, char * bin, unsigned int bin_len); 
	void (* unsorted_notification) (void * pairer, char * bin1, char * bin2); // it is called only once
	// reserved for the application passing its own data to the output function.
	void * appendix1;
	void * appendix2;
	void * appendix3;
	void * appendix4;
	void * appendix5;

} SAM_pairer_context_t;


#define SAM_PAIRER_WRITE_BUFFER ( 64000 )
typedef struct {
	unsigned char BIN_buffer[SAM_PAIRER_WRITE_BUFFER];
	int BIN_buffer_ptr;
	z_stream strm;

} SAM_pairer_writer_thread_t;

typedef struct {
	SAM_pairer_writer_thread_t * threads;	
	int all_threads;
	int compression_level;
	int has_dummy;
	FILE * bam_fp;
	char bam_name[MAX_FILE_NAME_LENGTH];
	subread_lock_t output_fp_lock;
} SAM_pairer_writer_main_t;


void fastq_64_to_33(char * qs);

// the caller is in charge of deallocation
char * memstrcpy(char * in);

int chars2color(char c1, char c2);

int genekey2color(char last_base,char * key);

// Convert a string key into an integer key
int genekey2int(char * key, int space_type);

// Open a read file. This function automatically decides its type.
// Return 0 if successfully opened, or 1 if error occurred
int geinput_open(const char * filename, gene_input_t * input);

// Open a sam file. Parameter half_no indicates which read of the pair is concerned.
// half_no = 1: the first read; half_no = 2: the last read; half_no = 0: single-end
// Return 0 if successfully opened, or 1 if error occurred
int geinput_open_sam(const char * filename, gene_input_t * input, int half_no);

// Read a line from the input and reward the pointer to the last position.
int geinput_readline_back(gene_input_t * input, char * linebuffer) ;

// Get the next read from the input file
// Return the length of this read or -1 if EOF. 
// The memory space for read_string must be at least 512 bytes.
int geinput_next_read(gene_input_t * input, char * read_name, char * read_string, char * quality_string);
int geinput_next_read_trim(gene_input_t * input, char * read_name, char * read_string, char * quality_string, short trim_5, short trim_3, int * is_secondary);

void geinput_jump_read(gene_input_t * input);

// Close the input file
void geinput_close(gene_input_t * input);

// return the next ATGC char from a input file.
// return 0 if this read segment reaches the end; return -1 if EOF; return -2 if error
int geinput_next_char(gene_input_t * input);

// line buffer has to be at least 300 bytes
// it returns the length of reading
int geinput_readline(gene_input_t * input, char * linebuffer, int conv_to_upper) ;



// read a line into the buff,
// the line should not be longer than 300 chars or the remaining part will be discarded.
// therefore the buff has to be at least 300 chars.
int read_line(int max_len, FILE * fp, char * buff, int conv_to_upper);

// count the number of reads in a flie
double guess_reads_density(char * fname, int is_sam) ;

// guess the size of the chromosome lib
// return the number of bases, or (-index-1) if the file at the index is not found.
srInt_64 guess_gene_bases(char ** files, int file_number);


void reverse_read(char * ReadString, int Length, int space_type);

void reverse_quality(char * QualtyString, int Length);

unsigned int read_numbers(gene_input_t * input);

//This function returns 0 if the line is a mapped read; -1 if the line is in a wrong format and 1 if the read is unmapped.
int parse_SAM_line(char * sam_line, char * read_name, int * flags, char * chro, unsigned int * pos, char * cigar, int * mapping_quality, unsigned int * pair_dist, char * sequence , char * quality_string, int * rl, int * repeated);

#define reverse_char(c)	((c)=='A'?'T':((c)=='G'?'C':((c)=='C'?'G':'A')))

int find_subread_end(int len, int  TOTAL_SUBREADS,int subread) ;

int break_SAM_file(char * in_SAM_file, int is_BAM, char * temp_file_prefix, unsigned int * real_read_count, int * block_no, chromosome_t * known_chromosomes, int is_sequence_needed, int base_ignored_head_tail, gene_value_index_t *array_index, gene_offset_t * offsets, srInt_64 * all_Mapped_bases , HashTable * event_table_ptr, char * VCF_file, srInt_64 * all_mapped_reads, int do_fragment_filtering, int push_to_read_head, int use_soft_clipped_bases);

int get_known_chromosomes(char * in_SAM_file, chromosome_t * known_chromosomes);


int load_exon_annotation(char * annotation_file_name, gene_t ** output_genes, gene_offset_t* offsets );

int is_in_exon_annotations(gene_t *output_genes, unsigned int offset, int is_start);

int does_file_exist (char * filename);

double guess_reads_density_format(char * fname, int is_sam, int * min_phred, int * max_phred, int * tested_reads);

FILE * get_temp_file_pointer(char *temp_file_name, HashTable* fp_table, int * close_immediately);

int write_read_block_file(FILE *temp_fp , unsigned int read_number, char *read_name, int flags, char * chro, unsigned int pos, char *cigar, int mapping_quality, char *sequence , char *quality_string, int rl , int is_sequence_needed, char strand, unsigned short read_pos, unsigned short read_len, unsigned short M_seg);

int get_read_block(char *chro, unsigned int pos, char *temp_file_suffix, chromosome_t *known_chromosomes, unsigned int * max_base_position);
int my_strcmp(const void * s1, const void * s2);

void destroy_cigar_event_table(HashTable * event_table);


int is_SAM_unsorted(char * SAM_line, char * tmp_read_name, short * tmp_flag, srInt_64 read_no);
int sort_SAM_add_line(SAM_sort_writer * writer, char * SAM_line, int line_len);
int sort_SAM_finalise(SAM_sort_writer * writer);
int sort_SAM_create(SAM_sort_writer * writer, char * output_file, char * tmp_path);
void colorread2base(char * read_buffer, int read_len);

int warning_file_type(char * fname, int expected_type);
char color2char(char clr, char c1);

int is_certainly_bam_file(char * fname, int * is_firstread_PE, srInt_64  * SAMBAM_header_length);

srUInt_64 sort_SAM_hash(char * str);

char * fgets_noempty(char * buf, int maxlen, FILE * fp);

char * gzgets_noempty(void * fp, char * buf, int maxlen);
int probe_file_type(char * fname, int * is_first_PE);
int probe_file_type_fast(char * fname);
void geinput_seek(gene_input_t * input, gene_inputfile_position_t * pos);
void geinput_tell(gene_input_t * input, gene_inputfile_position_t * pos);
srInt_64 geinput_file_offset( gene_input_t * input);
int probe_file_type_EX(char * fname, int * is_first_read_PE, srInt_64 * SAMBAM_header_length);


int SAM_pairer_create(SAM_pairer_context_t * pairer, int all_threads, int bin_buff_size_per_thread, int BAM_input, int is_Tiny_Mode, int is_single_end_mode, int force_do_not_sort, int need_read_group_tag, int display_progress, char * in_file, void (* reset_output_function) (void * pairer), int (* output_header_function) (void * pairer, int thread_no, int is_text, unsigned int items, char * bin, unsigned int bin_len), int (* output_function) (void * pairer, int thread_no, char * bin1, char * bin2), char * tmp_path, void * appendix1, int long_read_minimum_length) ;
int SAM_pairer_run( SAM_pairer_context_t * pairer);
void SAM_pairer_destroy(SAM_pairer_context_t * pairer);
void SAM_pairer_writer_reset(void * pairer);
void SAM_pairer_set_unsorted_notification(SAM_pairer_context_t * pairer, void (* unsorted_notification) (void * pairer, char * bin1, char * bin2));

int SAM_pairer_multi_thread_output( void * pairer, int thread_no, char * bin1, char * bin2 );
int SAM_pairer_multi_thread_header (void * pairer_vp, int thread_no, int is_text, unsigned int items, char * bin, unsigned int bin_len);
int SAP_pairer_skip_tag_body_len(char * tag_start);
int SAM_pairer_writer_create( SAM_pairer_writer_main_t * bam_main , int all_threads, int has_dummy , int BAM_output, int BAM_compression_level, char * out_file);
void SAM_pairer_writer_destroy( SAM_pairer_writer_main_t * bam_main ) ;
int SAM_pairer_iterate_int_tags(unsigned char * bin, int bin_len, char * tag_name, int * saved_value);
int SAM_pairer_iterate_tags(unsigned char * bin, int bin_len, char * tag_name, char * data_type, char ** saved_value);
int SAM_pairer_warning_file_open_limit();
int SAM_pairer_get_tag_bin_start(char * bin);
void *delay_realloc(void * old_pntr, size_t old_size, size_t new_size);
int is_comment_line(const char * l, int file_type, unsigned int lineno);
void warning_hash_hash(HashTable * t1, HashTable * t2, char * msg);
int geinput_preload_buffer(gene_input_t * input, subread_lock_t* read_lock);
int geinput_open_scRNA_fqs(char * fnames,  gene_input_t * input, int reads_per_chunk, int threads );
int geinput_open_scRNA_BAM(char * fnames,  gene_input_t * input, int reads_per_chunk, int threads );
int geinput_open_bcl( const char * dir_name,  gene_input_t * input, int reads_in_chunk, int threads );
char *strtokmm(char *str, const char *delim, char ** next);
#endif
