#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <math.h>
#include "subread.h"
#include "HelperFunctions.h"
#include "hashtable.h"
#include "input-files.h"
#include "gene-algorithms.h"

typedef struct{
	char gene_name[FEATURE_NAME_LENGTH+2];
	char chro_name[MAX_CHROMOSOME_NAME_LEN+2];
	unsigned int start, end;
	char is_negative_strand;
} DTCexon_t;

typedef struct {
	ArrayList * exon_table;
	ArrayList * chro_list;
	HashTable * chro_4bit_table;
	HashTable * chro_length_table;
	HashTable * gene_to_GC_TOTAL_table;
	HashTable * sam_chro_to_anno_chr_alias;
	FILE * out_FP_genes, * out_FP_bins;

	char out_file_name[MAX_FILE_NAME_LENGTH];
	char anno_file_name[MAX_FILE_NAME_LENGTH];
	char alias_file_name[MAX_FILE_NAME_LENGTH];
	int anno_file_format;
	char gene_id_column_name[MAX_FILE_NAME_LENGTH];
	char transcript_id_column_name[MAX_FILE_NAME_LENGTH];
	char used_feature_type[MAX_FILE_NAME_LENGTH];
	char fasta_file_name[MAX_FILE_NAME_LENGTH];

	int use_intron_bins;
	int filled_bins;
} DTCcontext_t;

int DTCdo_add_feature(char * gene_name, char * transcript_name, char * chro_name, unsigned int start, unsigned int end, int is_negative_strand, void * vcontext){
	DTCcontext_t * context = vcontext;
	DTCexon_t * new_exon = malloc(sizeof(DTCexon_t));

	if(context -> sam_chro_to_anno_chr_alias){
		char * sam_chro = get_sam_chro_name_from_alias(context -> sam_chro_to_anno_chr_alias, chro_name);
		if(sam_chro!=NULL) chro_name = sam_chro;
	}

	strncpy(new_exon -> gene_name, gene_name, FEATURE_NAME_LENGTH);
	strncpy(new_exon -> chro_name, chro_name, MAX_CHROMOSOME_NAME_LEN);
	new_exon -> start = start;
	new_exon -> end = end;
	new_exon -> is_negative_strand = is_negative_strand;

	ArrayListPush(context -> exon_table , new_exon);
	return 0;
}

int DTCcompare_merge_genes(void * L_elem, void * R_elem, ArrayList *me){
	unsigned int * le = L_elem, * re = R_elem;
	if(le[1] > re[1]) return 1;
	if(le[1] < re[1]) return -1;
	return 0;
}
int DTCcompare_exons(void * L_elem, void * R_elem, ArrayList *me){
	DTCexon_t * le = L_elem, * re = R_elem;
	int chname_cmp = strcmp(le -> chro_name, re -> chro_name);
	if(chname_cmp) return chname_cmp;
	long long lll = le -> start;
	lll -= re -> start;
	if(lll > 0) return 1;
	if(lll < 0) return -1;
	return 0;
}

#define BIT2_START_SIZE 1024
void DTCload_fasta(DTCcontext_t * context){
	FILE * fp = fopen(context -> fasta_file_name, "r");
	assert(fp);
	char * base_2bits = NULL;//malloc(BIT2_START_SIZE);
	unsigned int base_2bits_capasity = BIT2_START_SIZE*2;
	unsigned int base_current_chro_len = 0;
	char * curr_chro_mem = NULL;

	while(!feof(fp)){
		char fasta_line[501];
		char * fl = fgets(fasta_line, 501, fp);
		if(!fl) break;
		if(fl[0]=='>'){
			if(curr_chro_mem){
				ArrayListPush(context -> chro_list, curr_chro_mem);
				HashTablePut(context -> chro_length_table, curr_chro_mem , NULL + base_current_chro_len);
				HashTablePut(context -> chro_4bit_table, curr_chro_mem , base_2bits);
			}

			curr_chro_mem = malloc(strlen(fl));
			int ii;
			for(ii = 1; fl[ii] != '\n' && fl[ii] != '\r' && fl[ii] != 0 ; ii++){
				curr_chro_mem[ii-1] = fl[ii];
			}
			curr_chro_mem[ii-1] = 0;
			base_current_chro_len = 0;
			base_2bits = malloc(BIT2_START_SIZE);
			base_2bits_capasity = BIT2_START_SIZE*2;

		}else{
			int ii;
			for(ii=0; fl[ii] != '\n' && fl[ii] != '\r' && fl[ii] != 0 ; ii++){
				if(base_2bits_capasity <= base_current_chro_len){
					if(base_2bits_capasity < 1000000)
						base_2bits_capasity = base_2bits_capasity*2;
					else
						base_2bits_capasity = base_2bits_capasity*3/2;
					base_2bits = realloc(base_2bits, base_2bits_capasity/2 );
					HashTablePutReplaceEx(context -> chro_4bit_table, curr_chro_mem , base_2bits,0,0,0);
				}

				int ubase = toupper(fl[ii]), b2v=15;
				if(ubase == 'A')b2v = 0;
				else if(ubase == 'C')b2v = 1;
				else if(ubase == 'G')b2v = 2;
				else if(ubase == 'T')b2v = 3;
				int new_bch = base_2bits[ base_current_chro_len/2 ];
				if( base_current_chro_len%2 ) new_bch = ( b2v << 4 ) | (new_bch&15);
				else new_bch = b2v;
				base_2bits[ base_current_chro_len/2 ] = new_bch;
				base_current_chro_len ++;
			}
		}
	}
	if(curr_chro_mem){
		ArrayListPush(context -> chro_list, curr_chro_mem);
		HashTablePut(context -> chro_length_table, curr_chro_mem , NULL + base_current_chro_len);
		HashTablePut(context -> chro_4bit_table, curr_chro_mem , base_2bits);
	}
	SUBREADprintf("%lld items loaded from FASTA.\n", context -> chro_4bit_table -> numOfElements);
	fclose(fp);
}

void DTCadd_annotation(DTCcontext_t * context, char * gene_name, char * chro_name, unsigned int start, unsigned int stop, int is_negative_strand, unsigned int GCs, unsigned int totals, int is_bins){
	float fc = NAN;
	if(totals>0) fc =(float)GCs/totals;
	
	fprintf(is_bins? context -> out_FP_bins :context -> out_FP_genes , "%s\t%s\t%u\t%u\t%c\t%.5f\n", gene_name, chro_name, start, stop, is_negative_strand?'-':'+', fc);
}


int DTCwrite_annotations(char * gene_name, char * transcript_name, char * chro_name, unsigned int start, unsigned int end, int is_negative_strand, void * vcontext){
	DTCcontext_t * context = vcontext;
	unsigned int * gc_2int = HashTableGet(context ->gene_to_GC_TOTAL_table, gene_name);
	if(!gc_2int) gc_2int = (unsigned int[]){0,0,0};
	DTCadd_annotation(context, gene_name, chro_name, start, end, is_negative_strand, gc_2int[0], gc_2int[1], 0);
	return 0;
}

void DTCprint_lentab_items(void * key, void * hashed_obj, HashTable * tab){
	SUBREADprintf("%s => %ld\n", (char*)key, (long)(hashed_obj - NULL));
}

// start and end are 1-based
// chro 4bit bin are 0-based
unsigned int DTCtotal_GC_bases(DTCcontext_t * context, char * chro, unsigned int start, unsigned int end, unsigned int * total_bases){
	unsigned int ret = 0;
	char * bin2 = HashTableGet(context -> chro_4bit_table, chro);
	unsigned int chro_len = HashTableGet(context -> chro_length_table, chro) - NULL;
	
	//SUBREADprintf("TCHO: %s : %p , len = %u , range = %u ~ %u\n", chro, bin2, chro_len, start, end);
	if(NULL == bin2 && strlen(chro) <3){
		char tna [MAX_CHROMOSOME_NAME_LEN];
		sprintf(tna, "chr%s", chro);
		bin2 = HashTableGet(context -> chro_4bit_table, tna);
		chro_len = HashTableGet(context -> chro_length_table, tna) - NULL;
		//HashTableIteration(context -> chro_length_table, DTCprint_lentab_items);
	}

	if(NULL == bin2 && strlen(chro) >3 && memcmp("chr", chro, 3)==0){
		char *tna = chro + 3;
		bin2 = HashTableGet(context -> chro_4bit_table, tna);
		chro_len = HashTableGet(context -> chro_length_table, tna) - NULL;
	}

	if(chro_len < end) return 0;
	unsigned int ii;
	for(ii = start - 1; ii < end; ii++){
		int byte_no = ii / 2;
		int baseint = (ii%2)?( (bin2[byte_no]>> 4) & 0xf ):( bin2[byte_no] & 0xf );
		ret += (baseint == 1 || baseint == 2);
		(*total_bases) += (baseint<4);
	}
	return ret;
}

#define DTC_GAP_BIN_MIN_LENGTH 2000
void DTCanno_fill_gap(DTCcontext_t * context, char * chro_name, unsigned int start, unsigned int stop){
	//SUBREADprintf("FILL: %s - %u ~ %u\n", chro_name, start, stop);
	unsigned int cur;
	for(cur = start; ; cur += DTC_GAP_BIN_MIN_LENGTH){
		unsigned int bin_end = cur + DTC_GAP_BIN_MIN_LENGTH;
		if(bin_end > stop + 1 - DTC_GAP_BIN_MIN_LENGTH) bin_end = stop +1;
		unsigned int binlens = 0;
		unsigned int GCs = DTCtotal_GC_bases(context, chro_name, cur, bin_end - 1, &binlens);
		char bin_name[40];

		sprintf(bin_name, "_fill_bin_%07d", context -> filled_bins++);
		DTCadd_annotation(context, bin_name, chro_name, cur, bin_end - 1, 0, GCs, binlens, 1);
		if(bin_end >= stop) break;
	}
}

void DTCprocess_gene_tab(void * key, void * hashed_obj, HashTable * tab){
	ArrayList * merged_gene_list = tab -> appendix1;
	ArrayListPush(merged_gene_list , hashed_obj);
}

typedef struct{
	ArrayList * gene_3int_ptr_list;
	unsigned int start, end;
} MergedExon_t;

void DTC_free_ME_list(void * mel){
	MergedExon_t * me = mel;
	ArrayListDestroy(me -> gene_3int_ptr_list);
	free(mel);
}

void DTC_arrayListDestroy(void * l){
	ArrayListDestroy(( ArrayList *)l);
}

int DTCparse_GTF_and_Genome(DTCcontext_t * context){
	DTCload_fasta(context);

	load_features_annotation(context -> anno_file_name, context -> anno_file_format, context -> gene_id_column_name,
		context -> transcript_id_column_name, context -> used_feature_type, context, DTCdo_add_feature);
	SUBREADprintf("%lld items loaded from annotation file.\n", context -> exon_table -> numOfElements);

	long chri, ii;

	HashTable * chro_to_exons_table = HashTableCreate(1000);
	HashTableSetDeallocationFunctions(chro_to_exons_table, free, DTC_arrayListDestroy);
	HashTableSetKeyComparisonFunction(chro_to_exons_table, my_strcmp);
	HashTableSetHashFunction(chro_to_exons_table, fc_chro_hash);
	HashTable * gene_to_GC_TOTAL_table = HashTableCreate(10000);
	HashTableSetDeallocationFunctions(gene_to_GC_TOTAL_table, free, free);
	HashTableSetKeyComparisonFunction(gene_to_GC_TOTAL_table, my_strcmp);
	HashTableSetHashFunction(gene_to_GC_TOTAL_table, fc_chro_hash);

	for(ii = 0; ii < context -> exon_table -> numOfElements; ii++){
		DTCexon_t * tmpexon = ArrayListGet(context -> exon_table , ii);
		ArrayList * in_chro_exons = HashTableGet(chro_to_exons_table, tmpexon -> chro_name);
		if(NULL == in_chro_exons){
			in_chro_exons = ArrayListCreate(2000);
			char * memky = strdup(tmpexon -> chro_name);
			HashTablePut(chro_to_exons_table, memky, in_chro_exons);
		}
		ArrayListPush(in_chro_exons, tmpexon);
	}

	for(chri = 0; chri < context -> chro_list -> numOfElements ; chri++){
		char * current_chro = ArrayListGet(context -> chro_list, chri);
		ArrayList * in_chro_exons = HashTableGet(chro_to_exons_table , current_chro);
		if(NULL == in_chro_exons && memcmp("chr", current_chro, 3) == 0){
			char * tchro = current_chro + 3;
			in_chro_exons = HashTableGet(chro_to_exons_table , tchro);
		}
		if(NULL == in_chro_exons && strlen(current_chro) < 3){
			char tchro[MAX_CHROMOSOME_NAME_LEN];
			sprintf(tchro, "chr%s", current_chro);
			in_chro_exons = HashTableGet(chro_to_exons_table , tchro);
		}
		if(NULL == in_chro_exons) continue;

		HashTable * genename_to_3int_table = HashTableCreate(10000);
		HashTableSetDeallocationFunctions(genename_to_3int_table, free, free);
		HashTableSetKeyComparisonFunction(genename_to_3int_table, my_strcmp);
		HashTableSetHashFunction(genename_to_3int_table, fc_chro_hash);

		ArrayListSort(in_chro_exons, DTCcompare_exons);
		assert(in_chro_exons -> numOfElements > 0);
		SUBREADprintf("Annotations for %s (%lld) were sorted.\n", current_chro, in_chro_exons -> numOfElements);
	
		unsigned int gene_number = 1;
	
		char * current_chro_name_in_anno = NULL;
		for(ii = 0; ii < in_chro_exons -> numOfElements; ii++){
			DTCexon_t * tmpexon = ArrayListGet(in_chro_exons , ii);

			if(ii==0) current_chro_name_in_anno = tmpexon -> chro_name;

			char * tmpky = tmpexon -> gene_name;
			unsigned int * tmp_gene_info = HashTableGet(genename_to_3int_table, tmpky);
			unsigned int * gc_total_2int = HashTableGet(gene_to_GC_TOTAL_table, tmpky);

			if(tmp_gene_info == NULL){
				tmp_gene_info = malloc(sizeof(unsigned int) * 3);
				tmp_gene_info[0] = gene_number++;
				tmp_gene_info[1] = tmpexon -> start;
				tmp_gene_info[2] = tmpexon -> end;

				if(NULL == gc_total_2int){
					gc_total_2int = malloc(sizeof(int)*2);
					gc_total_2int[0] = gc_total_2int[1] = 0;
					HashTablePut(gene_to_GC_TOTAL_table, strdup(tmpky), gc_total_2int);
				}
				gc_total_2int[0] += DTCtotal_GC_bases(context, current_chro, tmpexon -> start, tmpexon -> end, &gc_total_2int[1]);
	
				char * memky = strdup(tmpky);
				HashTablePut(genename_to_3int_table , memky , tmp_gene_info);
			} else {
				assert(tmpexon -> start >= tmp_gene_info[1]);
				unsigned int new_start = max(tmpexon -> start, tmp_gene_info[2] +1);
				if(new_start <= tmpexon -> end){
					gc_total_2int[0] += DTCtotal_GC_bases(context, current_chro, new_start, tmpexon -> end, &gc_total_2int[1]);
					tmp_gene_info[2] = max(tmp_gene_info[2], tmpexon -> end);
				}
			}
		}

		/*for(ii = 0; ii < in_chro_exons -> numOfElements; ii++){
			DTCexon_t * tmpexon = ArrayListGet(in_chro_exons , ii);
			unsigned int * tmp_gene_info = HashTableGet(genename_to_3int_table , tmpexon -> gene_name);
			DTCadd_annotation(context, tmpexon -> gene_name, current_chro, tmpexon -> start, tmpexon -> end, tmpexon -> is_negative_strand, tmp_gene_info[3], tmp_gene_info[4]);	
		}*/

		if(context -> use_intron_bins){
			DTCexon_t * tmpexon = ArrayListGet(in_chro_exons , 0);
			unsigned int top_start = tmpexon -> start;
			unsigned int top_end = tmpexon -> end;

			ArrayList * current_ME_3int_ptrs = ArrayListCreate(5);
			ArrayList * merged_exon_list = ArrayListCreate(5000);
			ArrayListSetDeallocationFunction(merged_exon_list, DTC_free_ME_list);
			ArrayListPush(current_ME_3int_ptrs, HashTableGet(genename_to_3int_table, tmpexon -> gene_name));

			for(ii = 1; ii < in_chro_exons -> numOfElements; ii++){
				tmpexon = ArrayListGet(in_chro_exons , ii);
				if(top_end < tmpexon -> start){
					MergedExon_t * new_ME = malloc(sizeof(MergedExon_t));
					new_ME -> gene_3int_ptr_list = current_ME_3int_ptrs;
					new_ME -> start = top_start;
					new_ME -> end = top_end;
					ArrayListPush(merged_exon_list, new_ME);

					top_start = tmpexon -> start;
					top_end = tmpexon -> end;
					current_ME_3int_ptrs = ArrayListCreate(5);
					ArrayListPush(current_ME_3int_ptrs, HashTableGet(genename_to_3int_table, tmpexon -> gene_name));

				}else if(top_end < tmpexon -> end){
					top_end = tmpexon -> end;
					ArrayListPush_NoRepeatedPtr(current_ME_3int_ptrs, HashTableGet(genename_to_3int_table, tmpexon -> gene_name));
				}else{
					ArrayListPush_NoRepeatedPtr(current_ME_3int_ptrs, HashTableGet(genename_to_3int_table, tmpexon -> gene_name));
				}
			}
			if(1){
				MergedExon_t * new_ME = malloc(sizeof(MergedExon_t));
				new_ME -> gene_3int_ptr_list = current_ME_3int_ptrs;
				new_ME -> start = top_start;
				new_ME -> end = top_end;
				ArrayListPush(merged_exon_list, new_ME);
			}

			unsigned int last_exon_end = 0;
			unsigned int max_exon_end = 0;
			for(ii = 0; ii < merged_exon_list -> numOfElements; ii++){
				long ingene;

				MergedExon_t * current_ME = ArrayListGet(merged_exon_list, ii);

				if(last_exon_end>0){
					assert(last_exon_end < current_ME -> start);
					DTCanno_fill_gap(context, current_chro_name_in_anno, last_exon_end + 1, current_ME -> start -1);
				}
				for(ingene = 0; ingene < current_ME -> gene_3int_ptr_list -> numOfElements; ingene++){
					unsigned int *int3_gene = ArrayListGet(current_ME -> gene_3int_ptr_list, ingene);
					max_exon_end = max(max_exon_end, int3_gene[2]);
				}

				//SUBREADprintf("MAX_EXON_END = %u , THIS = %u - %u\n" , max_exon_end, current_ME -> start, current_ME -> end);
				if( max_exon_end > current_ME -> end ){
					last_exon_end = current_ME -> end;
				} else last_exon_end = 0;
			}
			
			ArrayListDestroy(merged_exon_list);
		}else{
			ArrayList * merged_gene_list = ArrayListCreate(10000);
			genename_to_3int_table -> appendix1 = merged_gene_list;
			HashTableIteration(genename_to_3int_table ,DTCprocess_gene_tab);
			ArrayListSort(merged_gene_list, DTCcompare_merge_genes);

			SUBREADprintf("%s has %lld merged genes from %lld item hashtable\n", current_chro, merged_gene_list -> numOfElements, genename_to_3int_table -> numOfElements);
			assert(merged_gene_list -> numOfElements>0);

			unsigned int * geneints = ArrayListGet(merged_gene_list, 0);
			unsigned int top_start = geneints[1];
			unsigned int top_end = geneints[2];

			DTCanno_fill_gap(context, current_chro_name_in_anno, 1, top_start);
			for(ii = 1; ii < merged_gene_list -> numOfElements; ii++){
				geneints = ArrayListGet(merged_gene_list, ii);
				if(top_end < geneints[1]){
					DTCanno_fill_gap(context, current_chro_name_in_anno, top_end + 1, geneints[1] - 1);
					top_start = geneints[1];
					top_end = geneints[2];
				}else if(top_end < geneints[2])
					top_end = geneints[2];
			}
			DTCanno_fill_gap(context, current_chro_name_in_anno, top_end + 1, HashTableGet(context -> chro_length_table, current_chro) - NULL - 1);
		}
	
		HashTableDestroy(genename_to_3int_table);
	}

	context -> gene_to_GC_TOTAL_table = gene_to_GC_TOTAL_table;
	load_features_annotation(context -> anno_file_name, context -> anno_file_format, context -> gene_id_column_name,
		context -> transcript_id_column_name, context -> used_feature_type, context, DTCwrite_annotations);
	HashTableDestroy(chro_to_exons_table);
	HashTableDestroy(gene_to_GC_TOTAL_table);
	return 0;
}

int DTCinit_context(DTCcontext_t ** context, int argc, char ** argv){
	DTCcontext_t * ret = malloc(sizeof(DTCcontext_t));
	memset(ret, 0, sizeof(DTCcontext_t));

	ret -> filled_bins = 0;
	ret -> anno_file_format = FILE_TYPE_GTF; 
	ret -> alias_file_name [0] = 0;

	strcpy(ret -> gene_id_column_name, "gene_id");
	strcpy(ret -> transcript_id_column_name, "transcript_id");
	strcpy(ret -> used_feature_type, "exon");

	ret -> chro_list = ArrayListCreate(100);
	ret -> exon_table = ArrayListCreate(100);
	ArrayListSetDeallocationFunction(ret -> exon_table, free);

	ret -> chro_4bit_table = HashTableCreate(100);
	HashTableSetDeallocationFunctions(ret -> chro_4bit_table, free, free);
	HashTableSetKeyComparisonFunction(ret -> chro_4bit_table, my_strcmp);
	HashTableSetHashFunction(ret -> chro_4bit_table, fc_chro_hash);

	ret -> chro_length_table = HashTableCreate(100); // chro_length_table is deallocated when chro_4bit_table is deallocated.
	HashTableSetKeyComparisonFunction(ret -> chro_length_table, my_strcmp);
	HashTableSetHashFunction(ret -> chro_length_table, fc_chro_hash);

	int c = 0;
	optind = 0;
	opterr = 1;
	optopt = 63;
	while ((c = getopt (argc, argv, "a:G:o:F:A:g:e:I"))!=-1){
		switch(c){
			case 'e':
				strcpy(ret -> used_feature_type, optarg);
				break;
			case 'g':
				strcpy(ret -> gene_id_column_name, optarg);
				break;
			case 'A':
				strcpy(ret -> alias_file_name, optarg);
				break;
			case 'F':
				if(strcmp("GTF", optarg)==0)
					ret -> anno_file_format = FILE_TYPE_GTF;
				else if(strcmp("SAF", optarg)==0)
					ret -> anno_file_format = FILE_TYPE_RSUBREAD;
				else{
					SUBREADprintf("ERROR: unknown annotation format: %s\n", optarg);
					return 1;
				}
			break;
			case 'o':
				strcpy(ret -> out_file_name, optarg);
			break;
			case 'a':
				strcpy(ret -> anno_file_name, optarg);
			break;
			case 'G':
				strcpy(ret -> fasta_file_name, optarg);
			break;
			case 'I':
				ret -> use_intron_bins = 1;
			break;
		}
	}

	if(ret -> alias_file_name[0]) ret -> sam_chro_to_anno_chr_alias = load_alias_table(ret -> alias_file_name);
	else ret -> sam_chro_to_anno_chr_alias = NULL;
	
	ret -> out_FP_genes = fopen(ret -> out_file_name,"w");
	char binfn[MAX_FILE_NAME_LENGTH+12];
	sprintf(binfn,"%s-bins", ret -> out_file_name);
	ret -> out_FP_bins = fopen(binfn,"w");

	fprintf(ret -> out_FP_genes, "GeneID\tChr\tStart\tEnd\tStrand\tGCfraction\n");
	fprintf(ret -> out_FP_bins, "GeneID\tChr\tStart\tEnd\tStrand\tGCfraction\n");
	*context = ret;
	return 0;
}

int DTCdestroy_context(DTCcontext_t *context){
	ArrayListDestroy(context -> exon_table);
	ArrayListDestroy(context -> chro_list);

	if(context -> sam_chro_to_anno_chr_alias)HashTableDestroy(context -> sam_chro_to_anno_chr_alias);
	HashTableDestroy(context -> chro_length_table);
	HashTableDestroy(context -> chro_4bit_table);
	fclose(context -> out_FP_genes);
	fclose(context -> out_FP_bins);
	free(context);
	return 0;
}

#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int detectionCall_Main(int argc, char ** argv)
#endif
{
	DTCcontext_t * context;
	int ret = DTCinit_context(&context, argc, argv);
	ret = ret || DTCparse_GTF_and_Genome(context);
	ret = ret || DTCdestroy_context(context);
	return ret;
}
