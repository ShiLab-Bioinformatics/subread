#ifndef __LRMconfig_h_
#define __LRMconfig_h_


#ifndef SRINT_64_DEFINED
#define SRINT_64_DEFINED
typedef long long srInt_64;
typedef unsigned long long srUInt_64;
#endif

#include <pthread.h>
#include <zlib.h>
#include "LRMseek-zlib.h"
#include "LRMhashtable.h"
#ifndef MAKE_STANDALONE
#ifndef RUNNING_ENV
#include <R.h>
#endif
#endif

#define LRMMAX_FILENAME_LENGTH 500
#define LRMMAX_READ_NAME_LEN 256
#define LRMMAX_CHROMOSOME_NAME_LEN 256
#define LRMMAX_EVENTS_PER_SITE 3

#define LRMMAX_READ_LENGTH 1200000
#define LRMMAX_THREADS 64
#define LRMMAX_SUBREAD_PER_READ_HARDLIMIT 1200000

#define LRMMAX_INDEL_TOLERANCE 16
#define LRMMAX_INDEL_SECTIONS 7
#define LRMMAX_CIGAR_OPTS_IN_READ 65500 
#define LRMSEGMENT_CIGAR_SIZE 1000
#define LRMSOFTCLIPPING_MAX_MISMATCH 2
#define LRMJUMP_MISMATCH_TOLERANCE 4
#define LRMDYNAMIC_MAXIMUM_GAP_LENGTH (15000)

// " * 250" is for RNA-seq mapping -- a much larger dynamic programming space.
#define LRMINDEL_DYNAMIC_CHANNEL_TOLERANCE (150 * 250 / 250 )


#define LRMSUBREAD_INDEX_OPTION_INDEX_GAP 0x0101
#define LRMSUBREAD_INDEX_OPTION_INDEX_PADDING 0x0102
#define LRMMULTI_THREAD_OUTPUT_ITEMS  (4096 * 3/5 *3)


#define LRMIS_REVERSED_HALVES 1
#define LRMIS_PAIRED_MATCH 2
#define LRMIS_NEGATIVE_STRAND 4
#define LRMIS_BREAKEVEN_READ 8
#define LRMIS_IN_KNOWN_EXONS 16

#define LRM_EVENT_IS_GT_AT_DONOR 1
#define LRM_EVENT_IS_ANNOTATED 2

#define LRMEVENT_TYPE_JUNCTION 30
#define LRMEVENT_TYPE_INDEL 20

#define LRMBAM_COMPRESS_LEVEL Z_BEST_SPEED
#define LRMGZIP_WINDOW_BITS -15
#define LRMSUBINDEX_VER2 201
#define LRMRUNNING_STEP_VOTING 10

#define LRMthread_lock_t pthread_mutex_t
#define LRMthread_lock  pthread_mutex_lock
#define LRMthread_lockrelease pthread_mutex_unlock
#define LRMthread_lockinit(a) pthread_mutex_init(a, NULL)

#define LRMGENE_VOTE_SPACE (131 - 80 ) // 131 gave the bast recall & accuracy & speed
#define LRMGENE_VOTE_TABLE_SIZE 64973
#define LRMMAX_LOCATIONS_PER_READ_HARDLIMIT (LRMGENE_VOTE_SPACE * LRMGENE_VOTE_TABLE_SIZE)

typedef char LRMgene_quality_score_t;
typedef unsigned int LRMgehash_key_t;
typedef unsigned int LRMgehash_data_t;

typedef struct {
	unsigned short items[LRMGENE_VOTE_TABLE_SIZE];
	unsigned int pos [LRMGENE_VOTE_TABLE_SIZE][LRMGENE_VOTE_SPACE];
	unsigned short votes [LRMGENE_VOTE_TABLE_SIZE][LRMGENE_VOTE_SPACE];
	short masks [LRMGENE_VOTE_TABLE_SIZE][LRMGENE_VOTE_SPACE];
	short indel_recorder [LRMGENE_VOTE_TABLE_SIZE][LRMGENE_VOTE_SPACE][LRMMAX_INDEL_SECTIONS*3];
	char current_indel_cursor[LRMGENE_VOTE_TABLE_SIZE][LRMGENE_VOTE_SPACE];
	unsigned short toli[LRMGENE_VOTE_TABLE_SIZE][LRMGENE_VOTE_SPACE];

	unsigned int coverage_start [LRMGENE_VOTE_TABLE_SIZE][LRMGENE_VOTE_SPACE];
	unsigned int coverage_end [LRMGENE_VOTE_TABLE_SIZE][LRMGENE_VOTE_SPACE];
} LRMgene_vote_t ;

#define LRMMAX_SORTED_WINDOW_NO_HARD_LIMIT 3

typedef struct {
	unsigned int read_no_in_chunk;
	unsigned int read_length;
	char read_name[LRMMAX_READ_NAME_LEN];
	char read_text[LRMMAX_READ_LENGTH];
	char qual_text[LRMMAX_READ_LENGTH];
	int is_reversed;
	unsigned int extracted_subreads;
	LRMgene_vote_t vote_table;

	unsigned int sorting_total_votes;
	unsigned int sorting_vote_locations[LRMMAX_LOCATIONS_PER_READ_HARDLIMIT];
	unsigned int sorting_subread_nos[LRMMAX_LOCATIONS_PER_READ_HARDLIMIT];
	unsigned short sorting_subread_votes[LRMMAX_LOCATIONS_PER_READ_HARDLIMIT];
	unsigned int sorting_is_negative_strand[LRMMAX_LOCATIONS_PER_READ_HARDLIMIT];

	unsigned int sorted_window_vote_start[LRMMAX_SORTED_WINDOW_NO_HARD_LIMIT];
	unsigned int sorted_window_vote_stop[LRMMAX_SORTED_WINDOW_NO_HARD_LIMIT];
	unsigned int sorted_window_total_votes[LRMMAX_SORTED_WINDOW_NO_HARD_LIMIT];
	unsigned int sorted_window_is_negative_strand[LRMMAX_SORTED_WINDOW_NO_HARD_LIMIT];

	unsigned int chain_tosmall_items;
	unsigned int chain_tolarge_items;
	unsigned int chain_total_items;
	unsigned int chain_cov_start[LRMMAX_SUBREAD_PER_READ_HARDLIMIT];
	unsigned int chain_cov_end[LRMMAX_SUBREAD_PER_READ_HARDLIMIT];
	unsigned int chain_chro_at_cov_start[LRMMAX_SUBREAD_PER_READ_HARDLIMIT];

	ArrayList * chain_used_gaps;
	double subread_extract_gap;
	int subread_total;
	int total_matched_bases ;
	int current_alignment_no;
	unsigned int final_mapping_location;
} LRMread_iteration_context_t;


typedef struct {
	char filename [LRMMAX_FILENAME_LENGTH];
	int file_type ;
	void * input_fp;   // can be system (FILE * sam or fastq or fasta), (seekable_zfile_t *)
	unsigned long long read_chunk_start;
} LRMgene_input_t;

typedef struct{
	unsigned int memory_block_size;
	unsigned int start_base_offset;
	unsigned int start_point;
	unsigned int length;
	unsigned char * values;
	unsigned int values_bytes;
} LRMgene_value_index_t;

struct LRMgehash_bucket {
	int current_items;
	int space_size;
	short * new_item_keys;
	LRMgehash_data_t * item_values;
};

typedef struct {
	int version_number;
	unsigned long long int current_items;
	int buckets_number;
	char is_small_table;
	struct LRMgehash_bucket * buckets;
	int index_gap;
	int padding;
} LRMgehash_t;



typedef struct{
	int thread_id;

	// SAM BAM output
	int output_buffer_item, output_buffer_pointer;
	LRMthread_lock_t output_lock;
	char * out_SAMBAM_buffer;
	int out_buff_used;
	int out_buff_capacity;

	z_stream bam_file_output_stream;
	
	// Chro events
	HashTable * events_fresh;		// only used in voting step ; event_name "I:chrX:99881122:-5" -> supporting_reads
	unsigned int * events_realignment_counts; // event_ids -> n_supp

	unsigned int * dynamic_programming_score_buffer;
	char * dynamic_programming_movement_buffer;
	char * dynamic_programming_indel_movement_buf;
	char * final_cigar_string;
	unsigned int dynamic_programming_indel_movement_start;
	unsigned int mapped_reads;
	
} LRMthread_context_t;

typedef struct{
	unsigned int small_side, large_side;
	unsigned int supporting_reads;
	short masks;
	char event_type;
	char indel_length; // >0: del; <0:ins
} LRMevent_t;


typedef struct {
	char user_command_line[10000];
	
	char input_file_name [LRMMAX_FILENAME_LENGTH];
	char output_file_name [LRMMAX_FILENAME_LENGTH];
	char index_prefix [LRMMAX_FILENAME_LENGTH];
	
	int threads;
	int is_SAM_output;
	int max_dynamic_indel_length;
	int max_read_indel_length;
	int multi_best_read_alignments;
	int is_Phred_64;
	int max_junction_distance;
	int max_mismatched_bases_in_subread;
	int min_voting_number;
	int min_matched_bases_in_alignment;
	int result_merge_tolerance;
	int show_read_validation;
	int max_cigars_in_read;
	int is_RNAseq_mode;
	
	// RUNNING STATUS
	
	double start_running_time;
	
	pthread_t running_threads[LRMMAX_THREADS];
	LRMthread_context_t thread_contexts[LRMMAX_THREADS];
	LRMthread_lock_t input_lock;
	LRMgene_input_t input_file;
	LRMgehash_t current_index;
	int current_index_padding;
	LRMgene_value_index_t current_base_index;
	seekable_position_t last_saved_zlib_pos;
	unsigned long long int last_saved_raw_pos;

	int input_exhausted;
	int processed_reads_in_chunk;
	int all_processed_reads;
	
	// Mapping results
	unsigned int mapped_reads;
	
	// Output
	char bam_file_tail_binary[200];
	int bam_file_tail_length;
	int sam_bam_file_header_written;
	FILE * sam_bam_file;

	ArrayList * chromosome_size_list;
	HashTable * chromosome_size_table;

	LRMthread_lock_t sam_bam_file_lock;
	HashTable * sam_bam_chromosome_table;
	ArrayList * sam_bam_chromosome_list;
	
	// Chro Events	
	LRMthread_lock_t event_space_lock;
	HashTable * events_realignment;	// only used in realignment step ; entry_position -> [max_N_size, N, event_id1, event_id2, ..., event_idN]
	LRMevent_t * event_space;
	unsigned int event_space_size;
	unsigned int event_number;	
	
	// Dynamic Programming
	int dynamic_programming_score_match, dynamic_programming_score_mismatch, dynamic_programming_score_create_gap, dynamic_programming_score_extend_gap;
} LRMcontext_t;


#ifdef MAKE_STANDALONE
#define LRMprintf printf
#define LRMputs puts

#else
#define LRMprintf Rprintf
#define LRMputs(x) Rprintf("%s\n",(x))
#endif

#define SEEKZLIBprintf LRMprintf

#define abs(a)    ((a)>=0?(a):-(a))
#define max(a,b)  ((a)<(b)?(b):(a))
#define min(a,b)  ((a)>(b)?(b):(a))
#define LRMint2base(c) (1413695297 >> (8*(c))&0xff)
#define LRMbase2int(c) ((c)<'G'?((c)=='A'?0:2):((c)=='G'?1:3))

#define LRMpthread_create pthread_create
#define LRMpthread_join pthread_join



int LRMshow_conf(LRMcontext_t* context);
int LRMrun_task(LRMcontext_t* context);
int LRMfinalise(LRMcontext_t* context);
int LRMdestroy_context(LRMcontext_t* context);
int LRMprint_mapping_summary(LRMcontext_t* context);

int LRMstart_thread(LRMcontext_t * context, int task );
int LRMwait_threads( LRMcontext_t * context );
int LRMmerge_threads(  LRMcontext_t * context , int step);

void LRMset_default_values_context(LRMcontext_t * context);
int LRMinput_has_finished( LRMcontext_t * context );
int LRMload_index(LRMcontext_t* context);
int LRMsave_input_pos(LRMcontext_t* context);
int LRMrewind_input_pos(LRMcontext_t* context);
int LRMmap_chunk_reads(LRMcontext_t* context);
int LRMrealign_write_chunk_reads(LRMcontext_t* context);
int LRMfinalise_chunk_reads(LRMcontext_t* context);
int LRMiterate_reads( LRMcontext_t * context, int task );
double LRMmiltime();

void * LRMmap_chunk_reads_run(void * thread_args);
void * LRMrealign_write_chunk_reads_run(void * thread_args);
int LRMchunk_read_iteration(LRMcontext_t * context, int thread_id, int task);
int LRMwrite_chunk_add_buffered_output(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, int flags,  char * chro_name, unsigned int chro_pos, int map_quality, char * cigar, int mis_matched);
//void LRMwrite_chunk_destroy_write_buffer(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iter_context, LRMwritting_buffer_t * buf);




int LRMFIXLENstrcmp(char * fixed_len, char * rname);

#endif
