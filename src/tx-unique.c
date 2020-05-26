#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "tx-unique.h"

void txunique_gene_free(void * v_gene){
	txunique_gene_t * gene = v_gene;

	ArrayListDestroy(gene -> transcript_list);
	free(gene);
}

int txunique_init_context(txunique_context_t * context){
	memset(context, 0 , sizeof(txunique_context_t));
	strcpy(context -> gene_name_column_name, "gene_id");
	strcpy(context -> transcript_id_column_name, "transcript_id");
	strcpy(context -> used_feature_type, "exon");

	context -> gene_table = HashTableCreate(62333);
	HashTableSetKeyComparisonFunction(context -> gene_table, my_strcmp);
	HashTableSetHashFunction(context -> gene_table, HashTableStringHashFunction);
	HashTableSetDeallocationFunctions(context -> gene_table, NULL, txunique_gene_free);

	context -> result_table = HashTableCreate(62333);
	HashTableSetKeyComparisonFunction(context -> result_table, my_strcmp);
	HashTableSetHashFunction(context -> result_table, HashTableStringHashFunction);
	HashTableSetDeallocationFunctions(context -> result_table, free, NULL);
	return 0;
}

void txunique_free_trans(void * v_trans){
	txunique_transcript_t * trans = v_trans;

	ArrayListDestroy(trans -> exon_list);
	free(trans);
}

int txunique_do_add_exon(char * gene_name, char * transcript_id, char * chrome_name, unsigned int start, unsigned int end, int is_negative_strand, void * v_context){
	int tx_i;

	txunique_context_t * context = v_context;
	
	txunique_gene_t * tag_gene = HashTableGet(context -> gene_table, gene_name);
	if(!tag_gene){
		tag_gene = malloc(sizeof(txunique_gene_t));
		strncpy(tag_gene -> gene_name, gene_name, FEATURE_NAME_LENGTH-1);
		tag_gene -> transcript_list = ArrayListCreate(6);
		ArrayListSetDeallocationFunction(tag_gene -> transcript_list, txunique_free_trans);
		HashTablePut(context -> gene_table, tag_gene -> gene_name, tag_gene);
	}

	txunique_transcript_t * tag_tx = NULL;
	//SUBREADprintf("NEW EXON: %s of %s : %s %u ~ %u ; txs = %lu in %lu\n", gene_name, transcript_id, chrome_name, start, end, tag_gene -> transcript_list -> numOfElements, tag_gene -> transcript_list -> capacityOfElements);
	for(tx_i = 0; tx_i < tag_gene -> transcript_list -> numOfElements; tx_i++){
		txunique_transcript_t * try_tx = ArrayListGet(tag_gene -> transcript_list , tx_i);
		if(strcmp(try_tx -> transcript_id, transcript_id) == 0){
			tag_tx = try_tx;
			break;
		}
	}

	if(!tag_tx){
		tag_tx = malloc(sizeof(txunique_transcript_t));
		strncpy(tag_tx -> transcript_id, transcript_id, FEATURE_NAME_LENGTH-1);
		tag_tx -> exon_list = ArrayListCreate(6);
		ArrayListSetDeallocationFunction(tag_tx -> exon_list, free);
		ArrayListPush(tag_gene -> transcript_list, tag_tx);
	}

	txunique_exon_t * tag_exon = malloc(sizeof(txunique_exon_t));
	strncpy(tag_exon -> chro_name, chrome_name, MAX_CHROMOSOME_NAME_LEN -1);
	tag_exon -> exon_start = start;
	tag_exon -> exon_stop = end;
	tag_exon -> is_negative_strand = is_negative_strand;
	ArrayListPush(tag_tx -> exon_list, tag_exon);
	
	return 0;
}

int txunique_load_annotation(txunique_context_t * context){
	int loaded_features = load_features_annotation(context -> input_GTF_file_name, FILE_TYPE_GTF, context -> gene_name_column_name, context -> transcript_id_column_name, context -> used_feature_type, context, txunique_do_add_exon);
	if(loaded_features<1) return -1;
	return 0;
}

int txunique_process_flat_comp( void * ex1p, void * ex2p, ArrayList *me ){
	txunique_exon_t * ex1 = ex1p;
	txunique_exon_t * ex2 = ex2p;

	if(ex1 -> exon_start < ex2 -> exon_start)return -1;
	if(ex1 -> exon_start > ex2 -> exon_start)return  1;
	return 0;
}

void debug_print_exs(ArrayList * exs){
	int ex_i;
	for(ex_i=0; ex_i<exs -> numOfElements; ex_i++){
		txunique_exon_t *ex = ArrayListGet(exs,ex_i);
		SUBREADprintf("   %s (%s) : %u ~ %u\n", ex->chro_name, ex->is_negative_strand?"NEG":"POS", ex->exon_start, ex->exon_stop);
	}
}

ArrayList * txunique_process_flat_exons(ArrayList * exs){
	ArrayList * ret = ArrayListCreate(5);
	ArrayListSetDeallocationFunction(ret, free);
	if(exs->numOfElements < 1) return ret;

	//SUBREADputs("Before Sorting");
	//debug_print_exs(exs);
	int ex_i;
	ArrayListSort(exs , txunique_process_flat_comp);

	//SUBREADputs("Before Flatten");
	//debug_print_exs(exs);

	txunique_exon_t * memex = malloc(sizeof(txunique_exon_t));
	memcpy(memex, ArrayListGet(exs,0), sizeof(txunique_exon_t));
	ArrayListPush(ret, memex);

	for(ex_i=1; ex_i<exs -> numOfElements; ex_i++){
		txunique_exon_t * lastex = ArrayListGet(ret, ret -> numOfElements - 1);
		txunique_exon_t * tryex = ArrayListGet(exs, ex_i);
		if(tryex -> exon_start > lastex -> exon_stop + 1){
			txunique_exon_t * memex = malloc(sizeof(txunique_exon_t));
			memcpy(memex, tryex, sizeof(txunique_exon_t));
			ArrayListPush(ret, memex);
		}
		else
			lastex -> exon_stop = max(tryex -> exon_stop, lastex -> exon_stop);
	}

	//SUBREADputs("After Flatten");
	//debug_print_exs(ret);
	return ret;
}

struct _txunique_tmp_edges{
	int is_exon_start;
	int nsupp;
	unsigned int base_open_end;
};

void debug_print_edges(ArrayList * exs){
	int ex_i;
	for(ex_i=0; ex_i<exs -> numOfElements; ex_i++){
		struct _txunique_tmp_edges *ex = ArrayListGet(exs,ex_i);
		SUBREADprintf("   %u : %s - nsup=%d\n", ex->base_open_end, ex->is_exon_start?"START":"END  ", ex->nsupp);
	}
}

int txunique_process_gene_edge_comp(void * e1p, void * e2p, ArrayList *me){
	struct _txunique_tmp_edges * e1 = e1p;
	struct _txunique_tmp_edges * e2 = e2p;
	if(e1 -> base_open_end < e2 -> base_open_end) return -1;
	if(e1 -> base_open_end > e2 -> base_open_end) return  1;

	if(e1 -> is_exon_start && !e2 -> is_exon_start) return -1;
	if(e2 -> is_exon_start && !e1 -> is_exon_start) return  1;
	return 0;
}

void txunique_process_gene_chro(txunique_context_t * context, char * chro, int strand_mode, txunique_gene_t * gene){
	int tx_i, ex_i;

	ArrayList **flatten_trans = malloc(gene -> transcript_list -> numOfElements*sizeof(void *));
	assert(flatten_trans);
	ArrayList * edge_list = ArrayListCreate(6);
	ArrayListSetDeallocationFunction(edge_list, free);

	for(tx_i = 0; tx_i < gene -> transcript_list -> numOfElements; tx_i ++){
		txunique_transcript_t * try_tx = ArrayListGet(gene -> transcript_list , tx_i);
		ArrayList * used_exons = ArrayListCreate(6);
		ArrayListSetDeallocationFunction(used_exons, free);

		for(ex_i = 0; ex_i < try_tx -> exon_list -> numOfElements; ex_i ++){
			txunique_exon_t * try_ex = ArrayListGet(try_tx -> exon_list, ex_i);
			if(strand_mode != try_ex -> is_negative_strand|| strcmp(try_ex -> chro_name, chro)!=0)continue;
			txunique_exon_t * memex = malloc(sizeof(txunique_exon_t));
			memcpy(memex, try_ex, sizeof(txunique_exon_t));
			ArrayListPush(used_exons, memex);
		}

		//SUBREADprintf("===== Process : %s : %s : %s : %s\n", gene -> gene_name, try_tx -> transcript_id, chro, strand_mode?"NEG":"POS");
		flatten_trans[tx_i] = txunique_process_flat_exons(used_exons);
		ArrayListDestroy(used_exons);

		for(ex_i = 0; ex_i < flatten_trans[tx_i] -> numOfElements; ex_i ++){
			txunique_exon_t * try_ex = ArrayListGet(flatten_trans[tx_i], ex_i);
			struct _txunique_tmp_edges * edge_start = malloc(sizeof(struct _txunique_tmp_edges)), *edge_end = malloc(sizeof(struct _txunique_tmp_edges));
			edge_start -> is_exon_start = 1;
			edge_start -> base_open_end = try_ex -> exon_start;
			edge_start -> nsupp = 0;
			edge_end -> is_exon_start = 0;
			edge_end -> base_open_end = try_ex -> exon_stop + 1;
			edge_end -> nsupp = 0;

			ArrayListPush(edge_list, edge_start);
			ArrayListPush(edge_list, edge_end);
		}
	}

	if(edge_list -> numOfElements >0){
			ArrayListSort(edge_list, txunique_process_gene_edge_comp);
			ArrayList * combined_edge_list = ArrayListCreate(6);
			ArrayListSetDeallocationFunction(combined_edge_list, free);

			struct _txunique_tmp_edges * merged_edge = ArrayListGet(edge_list, 0);
			merged_edge -> nsupp = 1;
			for(ex_i = 1; ex_i <= edge_list -> numOfElements; ex_i ++){
				struct _txunique_tmp_edges * tmpedge = NULL;
				if(ex_i < edge_list -> numOfElements) tmpedge = ArrayListGet(edge_list, ex_i);
				
				if(NULL == tmpedge || merged_edge -> is_exon_start != tmpedge -> is_exon_start || merged_edge -> base_open_end != tmpedge -> base_open_end){
					struct _txunique_tmp_edges * memedge = malloc(sizeof(struct _txunique_tmp_edges));
					memcpy(memedge, merged_edge, sizeof(struct _txunique_tmp_edges));
					ArrayListPush(combined_edge_list, memedge);
					if(tmpedge){
							merged_edge = tmpedge;
							merged_edge -> nsupp = 1;
					}
				}else merged_edge -> nsupp++;
			}

			//debug_print_edges(combined_edge_list);
			for(tx_i = 0; tx_i < gene -> transcript_list -> numOfElements; tx_i ++){
				txunique_transcript_t * try_tx = ArrayListGet(gene -> transcript_list , tx_i);
				ArrayList * flatten_exons = flatten_trans[tx_i];

				unsigned int total_bases = 0, unique_bases = 0, unique_start = 0, total_start = 0;
				int overlapping_count = 0, txex_i = 0, is_on = 0;
				for(ex_i = 0; ex_i < combined_edge_list -> numOfElements; ex_i ++){
					txunique_exon_t * txex = NULL;
					if(txex_i < flatten_exons -> numOfElements) txex = ArrayListGet(flatten_exons, txex_i);
					struct _txunique_tmp_edges *edge = ArrayListGet(combined_edge_list, ex_i);

					//SUBREADprintf("  TXN %s, Edge: %u (%s), Exon [%d]: %u ~ %u, %s, depth=%d, uniq_base=%u, all_base=%u\n", try_tx->transcript_id, edge -> base_open_end, edge -> is_exon_start?"START":"END  ", txex_i, txex?txex->exon_start:0, txex?txex->exon_stop:0 , is_on?"ON":"OFF", overlapping_count, unique_bases, total_bases);
					if(total_start<1) assert(!is_on);
					if(!is_on) assert(total_start < 1);
					if(edge -> is_exon_start){
						overlapping_count += edge -> nsupp;

						if(txex && edge -> base_open_end == txex -> exon_start)
							is_on = 1;

						if(unique_start >0){
							unique_bases += edge -> base_open_end - unique_start;
							unique_start = 0;
							assert(overlapping_count > 1);
						}else if(overlapping_count == 1&& is_on) unique_start = edge -> base_open_end;

						if(total_start<1 && is_on)total_start = edge -> base_open_end;
					}

					if(is_on) assert(overlapping_count>0);

					if(!edge -> is_exon_start){
						if(txex && edge -> base_open_end == txex -> exon_stop+1) is_on = 0;
						overlapping_count -= edge -> nsupp;

						if(unique_start >0){
							unique_bases += edge -> base_open_end - unique_start;
							unique_start =0;
							assert(overlapping_count ==0);
						} else if(overlapping_count == 1 && is_on) unique_start = edge -> base_open_end;

						if(total_start && !is_on){
							total_bases += edge -> base_open_end - total_start;
							total_start = 0;
						}
					}

					if(overlapping_count<1) assert(!is_on);
					if(!is_on) assert(total_start < 1);
					if(total_start<1) assert(!is_on);

					if(txex && edge-> is_exon_start==0){
							if(edge-> base_open_end > txex -> exon_stop) txex_i++;
					}

					//SUBREADprintf("      %s, Edge: %u (%s), Exon [%d]: %u ~ %u, %s, depth=%d, uniq_base=%u, all_base=%u\n\n", try_tx->transcript_id, edge -> base_open_end, edge -> is_exon_start?"START":"END  ", txex_i, txex?txex->exon_start:0, txex?txex->exon_stop:0 , is_on?"ON":"OFF", overlapping_count, unique_bases, total_bases);
				}
				assert(overlapping_count == 0);
				assert(total_start == 0);
				assert(unique_start == 0);
				assert(is_on == 0);

				char * hash_key = malloc(strlen(try_tx->transcript_id) + strlen(gene -> gene_name)+20);
				sprintf(hash_key, "%s\t%s\nALL", gene -> gene_name, try_tx->transcript_id);
				int old_all_bases = HashTableGet(context -> result_table, hash_key) - NULL;
				if(old_all_bases < 1)old_all_bases = 1;
				HashTablePut(context -> result_table, hash_key, NULL + old_all_bases + total_bases);

				hash_key = malloc(strlen(try_tx->transcript_id) + strlen(gene -> gene_name)+20);
				sprintf(hash_key, "%s\t%s\nUNIQUE", gene -> gene_name, try_tx->transcript_id);
				old_all_bases = HashTableGet(context -> result_table, hash_key) - NULL;
				if(old_all_bases < 1)old_all_bases = 1;
				HashTablePut(context -> result_table, hash_key, NULL + old_all_bases + unique_bases);
			}

			ArrayListDestroy(combined_edge_list);

		}
	ArrayListDestroy(edge_list);
	free(flatten_trans);
}

void txunique_process_write_gene(void * key, void * hashed_obj, HashTable * tab){
	txunique_context_t * context = tab -> appendix1;
	FILE * out_fp = tab -> appendix2;
	int tx_i;

	txunique_gene_t * gene = hashed_obj;
	for(tx_i = 0; tx_i < gene -> transcript_list -> numOfElements; tx_i ++){
		txunique_transcript_t * try_tx = ArrayListGet(gene -> transcript_list , tx_i);
		char hash_key [ FEATURE_NAME_LENGTH * 2 + 20];
		sprintf(hash_key, "%s\t%s\nALL", gene->gene_name, try_tx -> transcript_id);
		int all_bases = HashTableGet(context -> result_table, hash_key)-NULL-1;
		sprintf(hash_key, "%s\t%s\nUNIQUE", gene->gene_name, try_tx -> transcript_id);
		int unique_bases = HashTableGet(context -> result_table, hash_key)-NULL-1;
		fprintf(out_fp, "%s\t%s\t%d\t%d\n",  gene->gene_name, try_tx -> transcript_id, unique_bases, all_bases);
	}

}

void txunique_process_gene(void * key, void * hashed_obj, HashTable * tab){
	txunique_context_t * context = tab -> appendix1;
	txunique_gene_t * gene = hashed_obj;
	int tx_i, ex_i, ch_i;

	ArrayList * chro_list = ArrayListCreate(5);
	for(tx_i = 0; tx_i < gene -> transcript_list -> numOfElements; tx_i ++){
		txunique_transcript_t * try_tx = ArrayListGet(gene -> transcript_list , tx_i);
		for(ex_i = 0; ex_i < try_tx -> exon_list -> numOfElements; ex_i ++){
			txunique_exon_t * try_ex = ArrayListGet(try_tx -> exon_list, ex_i);
			int found_chro = 0;
			for(ch_i = 0; ch_i < chro_list->numOfElements; ch_i++){
				char * t_chro =  ArrayListGet(chro_list, ch_i);
				if(strcmp(t_chro, try_ex -> chro_name)==0){
					found_chro = 1;
						break;
				}
			}
			if(!found_chro) ArrayListPush(chro_list, try_ex -> chro_name);
		}
	}

	for(ch_i = 0; ch_i < chro_list -> numOfElements; ch_i ++){
		int is_negative_strand;
		char * chro = ArrayListGet(chro_list, ch_i);
		for(is_negative_strand = 0; is_negative_strand<2; is_negative_strand++)
				txunique_process_gene_chro(context, chro, is_negative_strand, gene);
	}
}

int txunique_find_unique_bases(txunique_context_t * context){
	context -> gene_table -> appendix1 = context;
	HashTableIteration(context -> gene_table, txunique_process_gene);
	return 0;
}

int txunique_write_output_file(txunique_context_t * context){
	FILE * out_fp = fopen(context -> output_file_name, "w" );
	fprintf(out_fp, "Gene_ID\tTranscript_ID\tUnique_Bases\tAll_Bases\n");
	if(!out_fp){
		SUBREADprintf("ERROR: unable to write output file : '%s'\n", context -> output_file_name);
		return 1;
	}

	context -> gene_table -> appendix1 = context;
	context -> gene_table -> appendix2 = out_fp;
	HashTableIteration(context -> gene_table, txunique_process_write_gene);
	fclose(out_fp);
	return 0;
}

int txunique_destroy_context(txunique_context_t * context){
	HashTableDestroy(context -> gene_table);
	HashTableDestroy(context -> result_table);
	return 0;
}

int txunique_parse_options(txunique_context_t * context, int argc, char ** argv){
	int c;

	optind = 1;
	opterr = 1;
	optopt = 63;

	while ((c = getopt (argc, argv, "a:o:g:t:f:h"))!=-1){
		switch(c){
			case 'a':
				strcpy(context -> input_GTF_file_name, optarg);
			break;

			case 'o':
				strcpy(context -> output_file_name, optarg);
			break;

			case 'g':
				strcpy(context -> gene_name_column_name, optarg);
			break;

			case 't':
				strcpy(context -> transcript_id_column_name, optarg);
			break;

			case 'f':
				strcpy(context -> used_feature_type, optarg);
			break;

			default:
				SUBREADputs("./txUnique -a <GTF_file> -o <output_text> { -g <gene_id_column> -t <transcript_id_column> -f <feature_type> }");
			break;
		}
	}

	if( context -> input_GTF_file_name[0]==0 || context -> output_file_name[0] == 0 ){
		SUBREADputs("The GTF file name and the output file name must be specified.");
		return 1;
	}
	return 0;
}

#ifdef MAKE_STANDALONE
int main(int argc, char ** argv){
#else
int TxUniqueMain(int argc, char ** argv){
#endif
		txunique_context_t context;
		int ret = 0;

		ret = ret || txunique_init_context(&context);
		ret = ret || txunique_parse_options(&context, argc, argv);
		ret = ret || txunique_load_annotation(&context);
		ret = ret || txunique_find_unique_bases(&context);
		ret = ret || txunique_write_output_file(&context);
		ret = ret || txunique_destroy_context(&context);
		if(!ret) SUBREADputs("All finished.");
		return ret;
}
