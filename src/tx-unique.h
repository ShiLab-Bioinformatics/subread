#ifndef __TX_UNIQUE_H_
#define __TX_UNIQUE_H_

#include "subread.h"
#include "hashtable.h"
#include "input-files.h"
#include "HelperFunctions.h"

typedef struct{
	char chro_name[MAX_CHROMOSOME_NAME_LEN];
	unsigned int exon_start;
	unsigned int exon_stop;
	int is_negative_strand;
} txunique_exon_t;

typedef struct{
	char transcript_id[FEATURE_NAME_LENGTH];
	ArrayList * exon_list;
} txunique_transcript_t;

typedef struct{
	char gene_name[FEATURE_NAME_LENGTH];
	ArrayList * transcript_list;
} txunique_gene_t;

typedef struct{
	char input_GTF_file_name[MAX_FILE_NAME_LENGTH];
	char output_file_name[MAX_FILE_NAME_LENGTH];
	char gene_name_column_name[FEATURE_NAME_LENGTH];
	char transcript_id_column_name[FEATURE_NAME_LENGTH];
	char used_feature_type[FEATURE_NAME_LENGTH];

	HashTable * gene_table;    // gene_id => array of transcripts [ (transcript_id, list_of_exons) ]
	HashTable * result_table;	// "$gene_id\t$transcript_id\nALL|UNIQUE" => NULL + number + 1
} txunique_context_t;

#endif
