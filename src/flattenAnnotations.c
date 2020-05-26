/***************************************************************

   The Subread and Rsubread software packages are free
   software packages:

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


#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>



#ifndef MAKE_STANDALONE
  #include <R.h>
#endif

#include <zlib.h>
#include <math.h>
#include <pthread.h>
#include <getopt.h>
#include "subread.h"
#include "interval_merge.h"
#include "core.h"
#include "gene-algorithms.h"
#include "sambam-file.h"
#include "input-files.h"
#include "hashtable.h"
#include "seek-zlib.h"
#include "HelperFunctions.h"

#define MERGING_MODE_CHOPPING 100   // keep all the edges : [100,199] + [150,249] = [100,149], [150,199], [200,249] 
#define MERGING_MODE_MERGING 0    // just merge into one: [100,199] + [150,249] = [100, 249]

typedef struct{
	char GTF_gene_id_column[MAX_READ_NAME_LEN];
	char GTF_wanted_feature_type[MAX_READ_NAME_LEN];
	char GTF_file_name[MAX_FILE_NAME_LENGTH];
	char output_file_name[MAX_FILE_NAME_LENGTH];
	FILE * output_FP;
	int merging_mode;

	HashTable * gene_to_chro_strand_table;
	HashTable * gene_chro_strand_to_exons_table;
} flatAnno_context_t;

static struct option long_options[] =
{
	{0, 0, 0, 0}
};


void flatAnno_print_usage(){
	SUBREADprintf("flattenGTF Version %s\n\n", SUBREAD_VERSION);
	SUBREADputs("  Flatten features included in a GTF annotation and save the modified annotation");
	SUBREADputs("  to a SAF format file.");
	SUBREADputs("");
	SUBREADputs("Usage:");
	SUBREADputs("  ./flattenGTF [options] -a <input_file> -o <output_file>");
	SUBREADputs("");
	SUBREADputs("## Mandatory arguments: ");
	SUBREADputs("");
	SUBREADputs("  -a <file>      Name of an annotation file in GTF/GFF format.");
	SUBREADputs("");
	SUBREADputs("  -o <file>      Name of output file.");
	SUBREADputs("");
	SUBREADputs("## Optional arguments: ");
	SUBREADputs("");
	SUBREADputs("  -t <string>    Specify feature type in a GTF annotation. 'exon' by default.");
	SUBREADputs("                 Features with the specified feature type are extracted from the");
	SUBREADputs("                 annotation for processing.");
	SUBREADputs("");
	SUBREADputs("  -g <string>    Specify attribute type in GTF annotation. 'gene_id' by default.");
	SUBREADputs("                 This attribute type is used to group features into meta-");
	SUBREADputs("                 features.");
	SUBREADputs("");
	SUBREADputs("  -C             Merging overlapping exons into multiple non-overlapping exons but");
	SUBREADputs("                 all the edges are kept.");
	SUBREADputs("");
}

int flatAnno_finalise(flatAnno_context_t * context){
	HashTableDestroy(context -> gene_to_chro_strand_table);
	HashTableDestroy(context -> gene_chro_strand_to_exons_table);
	fclose(context -> output_FP );
	SUBREADputs("Finished.\n");
	return 0;
}

int flatAnno_do_anno_1R(char * gene_name, char * transcript_name, char * chro_name, unsigned int start, unsigned int end, int is_negative_strand, void * Vcontext){
	flatAnno_context_t * context = Vcontext;
	ArrayList * chro_strand_list_for_gene = HashTableGet(context -> gene_to_chro_strand_table, gene_name);
	if(NULL == chro_strand_list_for_gene){
		char * mem_gene = malloc(strlen(gene_name)+1);
		strcpy(mem_gene, gene_name);

		chro_strand_list_for_gene = ArrayListCreate(3);
		ArrayListSetDeallocationFunction(chro_strand_list_for_gene, free);
		HashTablePut(context -> gene_to_chro_strand_table, mem_gene, chro_strand_list_for_gene);
	}
	int i, found = 0;
	char chro_strand[MAX_CHROMOSOME_NAME_LEN+10+FEATURE_NAME_LENGTH];
	sprintf( chro_strand, "%s\t%s\t%c", gene_name,chro_name,is_negative_strand?'-':'+');
	for(i=0; i<chro_strand_list_for_gene->numOfElements; i++){
		char * old_ch_st = ArrayListGet(chro_strand_list_for_gene,i);
		if(strcmp(old_ch_st, chro_strand) == 0){
			found=1;
			break;
		}
	}
	if(!found){
		char * mem_ch_st = strdup(chro_strand);
		ArrayListPush(chro_strand_list_for_gene, mem_ch_st);
	}

	ArrayList * ge_ch_st_exon_list = HashTableGet(context -> gene_chro_strand_to_exons_table, chro_strand);
	if(NULL == ge_ch_st_exon_list){
		ge_ch_st_exon_list =  ArrayListCreate(3);
		ArrayListSetDeallocationFunction(ge_ch_st_exon_list,free);
		HashTablePut(context -> gene_chro_strand_to_exons_table, strdup(chro_strand), ge_ch_st_exon_list);
	}
	int * mem_start_end = malloc(sizeof(int)*2);
	mem_start_end[0] = (int)start;
	mem_start_end[1] = (int)end;
	ArrayListPush(ge_ch_st_exon_list, mem_start_end);
	return 0;
}

int flatAnno_do_anno_merge_one_array_compare(void * vL, void * vR, ArrayList *me){
	int * iL = vL, *iR = vR;
	if((*iL)>(*iR))return 1;
	if((*iL)<(*iR))return -1;
	return 0;
}

void flatAnno_do_anno_chop_one_array(void * key, void * hashed_obj, HashTable * tab){
	ArrayList * this_list = hashed_obj; // each element is a pair of integers defining start and end of an exon.
	ArrayList * edge_before_me_List = ArrayListCreate( this_list->numOfElements * 2 );
	int i;
	for(i=0; i<this_list -> numOfElements; i++){
		int * curr_2i = this_list -> elementList[ i ];
		int findi, endi;
		for(endi = 0; endi < 2; endi++){
			assert(curr_2i[endi]>0);
			int search_tag = curr_2i[endi] + endi, found=0;

			for(findi = 0; findi < edge_before_me_List->numOfElements; findi++){
				if(ArrayListGet(edge_before_me_List,findi)-NULL == search_tag){
					found = 1;
					break;
				}
			}
			if(0==found) ArrayListPush(edge_before_me_List, NULL+search_tag);
		}
	}
	ArrayListSort(edge_before_me_List, NULL); 
	char * continue_after_an_edge = malloc(edge_before_me_List -> numOfElements -1);
	memset(continue_after_an_edge,0, edge_before_me_List -> numOfElements -1);
	int missed_gaps = edge_before_me_List->numOfElements -1;
    for(i=0; i<edge_before_me_List->numOfElements -1; i++){
		int an_edge_before_me = edge_before_me_List-> elementList[i]-NULL;
		int j;
		for(j=0; j<this_list->numOfElements; j++){
			int * curr_2i = this_list -> elementList[ j ];
			if( curr_2i[0] <= an_edge_before_me && curr_2i[1] >= an_edge_before_me ) {
				continue_after_an_edge[i]=1;
				missed_gaps --;
				break;
			}
		}
	}

	long old_capacity = this_list -> capacityOfElements, old_items = this_list -> numOfElements;
	if(edge_before_me_List->numOfElements -1 - missed_gaps > old_capacity){
		this_list -> elementList = realloc( this_list -> elementList, sizeof(void*) * edge_before_me_List->numOfElements -1 -missed_gaps );
		this_list -> capacityOfElements = edge_before_me_List->numOfElements -1 -missed_gaps ;
	}

	int written_i = 0;
    for(i=0; i<edge_before_me_List->numOfElements -1; i++){
		int * curr_2i = NULL;
		//if(edge_before_me_List->elementList[i] - NULL == 108395906)SUBREADprintf("CV_I %d = %ld ~ %ld ; CTN=%d\n", i, edge_before_me_List->elementList[i] - NULL,  edge_before_me_List -> elementList[i+1] -NULL-1, continue_after_an_edge[i]);
		if(0==continue_after_an_edge[i]){
		//	SUBREADprintf("NC_I %d ~ [ 0 , %ld ]\n", i, edge_before_me_List->numOfElements -2);
			assert(i>0 && i<edge_before_me_List->numOfElements -2);
			assert(continue_after_an_edge[i+1] && continue_after_an_edge[i-1] ); // otherwise it is a singleton edge.
			continue;
		}

		if(written_i < old_items){
			curr_2i = this_list -> elementList[ written_i ];
		}else{
			curr_2i = malloc(sizeof(int)*2);
			this_list -> elementList[ written_i ] = curr_2i;
		}
		curr_2i[0] = edge_before_me_List -> elementList[i] -NULL;
		curr_2i[1] = edge_before_me_List -> elementList[i+1] -NULL-1;
		written_i++;
	}
	assert(written_i == edge_before_me_List->numOfElements -1 - missed_gaps );
	for(i=written_i; i<old_items; i++) free( this_list -> elementList[i]);
	this_list -> numOfElements = written_i;
	ArrayListDestroy(edge_before_me_List);
	free(continue_after_an_edge);
}

void flatAnno_do_anno_merge_one_array(void * key, void * hashed_obj, HashTable * tab){
	ArrayList * this_list = hashed_obj;
	ArrayListSort(this_list, flatAnno_do_anno_merge_one_array_compare);
	int i, n1_items = 0;
	for(i=1; i<this_list -> numOfElements; i++){
		int * last_2i = this_list -> elementList[ n1_items ];
		int * curr_2i = this_list -> elementList[ i ];
		if(last_2i[1] >= curr_2i[1]) continue;
		if(last_2i[1] >= curr_2i[0] -1){
			last_2i[1]=curr_2i[1];
			continue;
		}
		n1_items++;
		if(n1_items< i){
			last_2i = this_list -> elementList[ n1_items ];
			last_2i[0] = curr_2i[0];
			last_2i[1] = curr_2i[1];
		}
	}
	for(i=n1_items+1; i<this_list -> numOfElements; i++)free(this_list -> elementList[i]);
	this_list -> numOfElements = n1_items+1;
}

int flatme_strcmp(void * L, void * R, ArrayList * me){
	return strcmp(L, R);
}

int flatAnno_do_anno_merge_and_write(flatAnno_context_t * context){
	context -> gene_chro_strand_to_exons_table -> appendix1 = context;

	HashTableIteration(context -> gene_chro_strand_to_exons_table, context -> merging_mode == MERGING_MODE_CHOPPING? flatAnno_do_anno_chop_one_array :flatAnno_do_anno_merge_one_array);
	ArrayList * all_chro_st_list = HashTableKeyArray(context -> gene_chro_strand_to_exons_table);
	ArrayListSort(all_chro_st_list, flatme_strcmp);

	fprintf( context -> output_FP , "GeneID\tChr\tStart\tEnd\tStrand\n");
	int i,j;
	for(i = 0; i< all_chro_st_list -> numOfElements; i++){
		char * ge_chro_strand = ArrayListGet(all_chro_st_list,i);
		char * local_ge_chro_strand = strdup(ge_chro_strand);
		char * strand_ptr = local_ge_chro_strand;

		for(j=0; j<2; strand_ptr++)
			if(*strand_ptr=='\t')j++;
		strand_ptr[-1] = 0;

		ArrayList * exon_in_chro_strand = HashTableGet(context -> gene_chro_strand_to_exons_table, ge_chro_strand);

		for(j=0; j< exon_in_chro_strand  -> numOfElements; j++){
			int * start_end_2i = ArrayListGet(exon_in_chro_strand,j);
			fprintf( context -> output_FP ,"%s\t%d\t%d\t%s\n", local_ge_chro_strand, start_end_2i[0], start_end_2i[1], strand_ptr);
		}
		free(local_ge_chro_strand);
	}
	ArrayListDestroy(all_chro_st_list );
	return 0;
}

int flatAnno_do_anno(flatAnno_context_t * context){
	int loaded_features = load_features_annotation(context -> GTF_file_name, FILE_TYPE_GTF, context -> GTF_gene_id_column, NULL, context ->  GTF_wanted_feature_type, context,  flatAnno_do_anno_1R);
	if(loaded_features<0)SUBREADputs("ERROR: Unable to open the GTF file.");
	if(loaded_features==0)SUBREADprintf("ERROR: No '%s' feature was found in the GTF file. (the '%s' attribute is required)\n", context ->  GTF_wanted_feature_type,  context -> GTF_gene_id_column);
	if(loaded_features<=0) return -1;

	return flatAnno_do_anno_merge_and_write(context);
}


int flatAnno_start(flatAnno_context_t * context){
	SUBREADputs("");
	if(!context -> GTF_file_name[0]){
		flatAnno_print_usage();

		if(context -> output_file_name[0]) SUBREADprintf("Error: no input file is specified.\n");
		return -1;
	}

	if(!context -> output_file_name[0]){
		flatAnno_print_usage();
		SUBREADprintf("Error: no output file is specified.\n");
		return -1;
	}
	SUBREADprintf("Flattening GTF file: %s\n", context -> GTF_file_name);
	SUBREADprintf("Output SAF file: %s\n", context -> output_file_name);
	context -> output_FP = fopen(context -> output_file_name, "w");
	if(NULL == context -> output_FP ){
		SUBREADprintf("Error: unable to open the output file.\n");
		return -1;
	}

	SUBREADprintf("\nLooking for '%s' features... (grouped by '%s')\n\n", context ->  GTF_wanted_feature_type, context -> GTF_gene_id_column);
	context -> gene_to_chro_strand_table = StringTableCreate(30000);
	HashTableSetDeallocationFunctions(context -> gene_to_chro_strand_table, free, (void(*)(void *))ArrayListDestroy);
	context -> gene_chro_strand_to_exons_table =  StringTableCreate(30000);
	HashTableSetDeallocationFunctions(context -> gene_chro_strand_to_exons_table, free, (void(*)(void *))ArrayListDestroy);
	return 0;
}

#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int R_flattenAnnotations(int argc, char ** argv)
#endif
{
	flatAnno_context_t context;
	memset(&context, 0, sizeof(context));
	strcpy(context.GTF_gene_id_column, "gene_id");
	strcpy(context.GTF_wanted_feature_type, "exon");

	int option_index = 0,c;
	optind=0;
	opterr=1;
	optopt=63;
	while ((c = getopt_long (argc, argv, "Ct:g:a:o:v?", long_options, &option_index)) != -1)
		switch(c) {
			case 'C':
				context.merging_mode = MERGING_MODE_CHOPPING;
				break;
			case 'o':
				strcpy(context.output_file_name,  optarg);
				break;
			case 'a':
				strcpy(context.GTF_file_name,  optarg);
				break;
			case 'g':
				strcpy(context.GTF_gene_id_column, optarg);
				break;
			case 't':
				strcpy(context.GTF_wanted_feature_type, optarg);
				break;

			case '?':
			default :
				flatAnno_print_usage();
				return -1;
		}

	int ret = 0;
	ret = ret || flatAnno_start(&context);
	ret = ret || flatAnno_do_anno(&context);
	ret = ret || flatAnno_finalise(&context);
	return ret;
}
