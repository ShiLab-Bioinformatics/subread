#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <dirent.h>
#include <getopt.h>
#include "subread.h" 
#include "core.h" 
#include "hashtable.h" 
#include "gene-algorithms.h" 
#include "HelperFunctions.h"

char output_file_name[MAX_FILE_NAME_LENGTH];
FILE * output_file_pointer;

static struct option propm_long_options[] =
{
	{"out",  required_argument, 0, 'o'},
	{0, 0, 0, 0}
};

typedef struct{
	char ** keys;
} VCF_sort_context;


void VCF_sort_getv(char * ky, char * chro, int * pos)
{
	int cursor, nch, tabs =0, chro_cursor=0;
	(*pos)= 0;
	for(cursor = 0;0!=(nch=ky[cursor]); cursor++){
		if(nch == '\t')
		{
			tabs ++;
			if(tabs == 3)break;
		}
		else if(tabs == 1)
			chro[chro_cursor++]=nch;
		else if(tabs == 2)
			(*pos) = (*pos)*10 + (nch - '0');
	}
	chro[chro_cursor]=0;

} 

int compare_VCF_rows(void * arri, int l, int r)
{
	VCF_sort_context * arr = (VCF_sort_context*) arri;
	char * lk = arr -> keys[l];
	char * rk = arr -> keys[r];

	char chro_1[MAX_CHROMOSOME_NAME_LEN+1];
	char chro_2[MAX_CHROMOSOME_NAME_LEN+1];
	int pos1 =0, pos2 = 0;

	VCF_sort_getv(lk,chro_1, &pos1);
	VCF_sort_getv(rk,chro_2, &pos2);
//	printf("LK=%s; CHRO=%s; POS=%d\n", lk, chro_1, pos1);
	int ret = strcmp_number(chro_1, chro_2);
	if(ret==0)
	{
		if(pos1>pos2) ret = 1;
		if(pos1<pos2) ret = -1;
	}

//	printf("%s %d %s\n\n", lk, ret, rk);
	return ret;
}

void exchange_VCF_rows(void * arri, int l, int r)
{
        VCF_sort_context * arr = (VCF_sort_context*) arri;
	char * lk = arr -> keys[l];
	char * rk = arr -> keys[r];
	arr -> keys[l] = rk;
	arr -> keys[r] = lk;
}

void merge_VCF_rows(void * arri, int start, int items, int items2)
{
        VCF_sort_context * arr = (VCF_sort_context*) arri;

	int x1, first_cursor = start, second_cursor = start+items;
	char ** merge_tmp = malloc(sizeof(char *)* (items+items2));

	
	for(x1=0;x1<items+items2;x1++)
	{
		int choose_2 = 0;
		if(first_cursor >= start+items)
			choose_2 = 1;
		else if(second_cursor < start+items+items2)
			choose_2 = (0 < compare_VCF_rows(arr, first_cursor, second_cursor));
		if(choose_2)
			merge_tmp[x1] = arr -> keys[second_cursor++];
		else	
			merge_tmp[x1] = arr -> keys[first_cursor++];
	}

	memcpy(arr -> keys + start, merge_tmp, sizeof(char *) * (items+items2));
	free(merge_tmp);
}


int endswith(const char *str, const char *suffix)
{
    if (!str || !suffix)
	return 0;
    size_t lenstr = strlen(str);
    size_t lensuffix = strlen(suffix);
    if (lensuffix >  lenstr)
	return 0;
    return strncmp(str + lenstr - lensuffix, suffix, lensuffix) == 0;
}

void set_hash_funcs(HashTable * tab, int is_key_freed, int is_value_string)
{
	HashTableSetDeallocationFunctions(tab,is_key_freed?free:NULL,is_value_string?free:NULL);
	HashTableSetHashFunction(tab , fc_chro_hash);
	HashTableSetKeyComparisonFunction(tab,fc_strcmp_chro);
}

int warning_reported_repeated;

void do_find_common(char ** file_names, int files) 
{
	int badfiles = 0;
	char * file_line = malloc(3000);
	int file_no;
	HashTable *  extra_count_table = HashTableCreate(200000);
	set_hash_funcs(extra_count_table,1,0);

	HashTable *  extra_qual_table = HashTableCreate(200000);
	set_hash_funcs(extra_qual_table,0,0);

	HashTable *  extra_info_table = HashTableCreate(200000);
	set_hash_funcs(extra_info_table,0,1);

	HashTable *  extra_id_qual_filter_table = HashTableCreate(200000);
	set_hash_funcs(extra_id_qual_filter_table,0,1);

	for(file_no=0; file_no < files; file_no++)
	{
		FILE * ifp = fopen(file_names[file_no],"r");
		if(NULL == ifp)
		{
			SUBREADprintf("Error: unable to open file '%s'\n", file_names[file_no]);
			badfiles ++;
			continue;
		}
		else
			SUBREADprintf("Process file '%s' ...\n", file_names[file_no]);

		int lines = 0;
		while(1)
		{
			char * tmp_pnt;
			char * ret = fgets(file_line, 3000, ifp);
			if(!ret)break;
			if(ret[0]=='#') continue;

			char * chro = strtok_r(ret,"\t",&tmp_pnt);
			if(chro == NULL|| tmp_pnt == NULL)
			{
				SUBREADprintf("Unknown format in the file! Only the VCF format is allowed for input files.\n");
				break;
			}

			char * pos_str = strtok_r(NULL,"\t",&tmp_pnt);
			if(pos_str == NULL|| strlen(pos_str)>10|| tmp_pnt == NULL)
			{
				SUBREADprintf("Unknown format in the file! Only the VCF format is allowed for input files.\n");
				break;
			}



			char * id_str = strtok_r(NULL,"\t",&tmp_pnt); // snp_id
			char * ref = strtok_r(NULL,"\t",&tmp_pnt); 
			char * alt_str = strtok_r(NULL,"\t",&tmp_pnt); // alt 
			char * qual_str = strtok_r(NULL,"\t",&tmp_pnt); //qual
			char * filter_str = strtok_r(NULL,"\t",&tmp_pnt); //filter
			char * info = strtok_r(NULL,"\t",&tmp_pnt); //addition columns are discarded
			if(info == NULL)
			{
				SUBREADprintf("Unknown format in the file! Only the VCF format is allowed for input files.\n");
				break;
			}



			int is_indel = 0;

			if(strstr(info, "INDEL;")) is_indel = 1;
			if(endswith(info, "INDEL")) is_indel = 1;

			char *tmp_pnt_alt = NULL;

			//SUBREADprintf("GO_INTO_STR: %s\n", alt_str);
			while(1)
			{
				char * alt_one = strtok_r(tmp_pnt_alt?NULL:alt_str,",", &tmp_pnt_alt);
				if(!alt_one) break;
				lines ++;

				char key_type = is_indel?'I':'S';

				char * ky = malloc(strlen(chro)+strlen(alt_one)+strlen(ref)+40);
				sprintf(ky, "%c\t%s\t%s\t.\t%s\t%s\t.", key_type, chro, pos_str, ref, alt_one);
				int qual = atoi(qual_str);

				if(file_no - badfiles == 0)
				{
					char * IQF_buf = malloc(strlen(qual_str)+strlen(filter_str)+strlen(id_str)+4);
					char * info_buf = malloc(strlen(info)+1);
					assert(info_buf);
					strcpy(info_buf, info);
					sprintf(IQF_buf,"%s\t%s\t%s", id_str,qual_str,filter_str);
					int repeat_times = HashTableGet(extra_count_table, ky) - NULL;
					if(0 == repeat_times)
					{
						HashTablePut(extra_count_table, ky, NULL+1);
						HashTablePut(extra_id_qual_filter_table, ky, IQF_buf);
						HashTablePut(extra_info_table, ky, info_buf);
						HashTablePut(extra_qual_table, ky, NULL + qual + 1);
					}
					else
					{
						if(!warning_reported_repeated)SUBREADprintf("Warning: repeated rows are found in the first input file.\n");
						warning_reported_repeated=1;
						free(ky);
					}
				}
				else{
					int repeat_times = HashTableGet(extra_count_table, ky) - NULL;
					if(repeat_times >0)
					{
						int old_qual = HashTableGet(extra_qual_table, ky) - NULL - 1;
						repeat_times ++;
						HashTablePutReplace(extra_count_table, ky, NULL + repeat_times, 0);
						if(old_qual > qual)
						{
							char * info_buf = malloc(strlen(info)+1);
							char * IQF_buf = malloc(strlen(qual_str)+strlen(filter_str)+strlen(id_str)+4);
							assert(info_buf);
							strcpy(info_buf, info);
							sprintf(IQF_buf,"%s\t%s\t%s", id_str,qual_str,filter_str);

							HashTablePutReplace(extra_info_table, ky, info_buf, 0);
							HashTablePutReplace(extra_id_qual_filter_table, ky, IQF_buf, 0);
							HashTablePutReplace(extra_qual_table, ky, NULL + 1 + qual, 0);
						}
					}
					free(ky);
				}
				if(NULL == tmp_pnt_alt) break;
			}
		}

		SUBREADprintf("There are %d variants found in this file.\n\n", lines);
		fclose(ifp);
	}


	KeyValuePair *cursor;
	int bucket, x1;
	char ** output_list;
	int out_list_size = 300, out_list_items = 0;
	output_list = malloc(out_list_size * sizeof(char *));
	for(bucket=0; bucket < extra_count_table -> numOfBuckets; bucket++)
	{
		cursor = extra_count_table -> bucketArray[bucket];
		while (1)
		{
			if (!cursor) break;
			char * ky = (char *)cursor -> key;
			int counts = cursor -> value - NULL;
			assert(counts <= files);

			if(counts == files - badfiles)
			{
				if(out_list_items >= out_list_size)
				{
					out_list_size *= 1.5;
					output_list = realloc(output_list, out_list_size * sizeof(char *));
				}
				output_list[out_list_items] = ky;

				out_list_items ++;
			}

			cursor = cursor->next;
		}
	}


	VCF_sort_context vcf_context;
	vcf_context.keys = output_list;
	merge_sort(&vcf_context, out_list_items,compare_VCF_rows,exchange_VCF_rows,merge_VCF_rows); 

	fprintf(output_file_pointer, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO");
	for(x1 =0 ; x1< out_list_items; x1++)
	{
		char * ky = output_list[x1];
		char * info = HashTableGet(extra_info_table, ky);
		char * IQF_buffer = HashTableGet(extra_id_qual_filter_table, ky);
		char * tok_buf_tab = NULL;

		char * id_str = strtok_r(IQF_buffer,"\t", &tok_buf_tab);
		char * qual_str = strtok_r(NULL,"\t", &tok_buf_tab);
		char * filter_str = strtok_r(NULL,"\t", &tok_buf_tab);

		char * chro_str = strtok_r(ky+2,"\t",&tok_buf_tab);
		char * pos_str = strtok_r(NULL,"\t",&tok_buf_tab);
		strtok_r(NULL,"\t",&tok_buf_tab);// fake_id
		char * ref_str = strtok_r(NULL,"\t",&tok_buf_tab);
		char * alt_str = strtok_r(NULL,"\t",&tok_buf_tab);


		char *brk = "\n";
		if(info[strlen(info)-1]=='\n') brk="";

		fprintf(output_file_pointer, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s", chro_str, pos_str, id_str, ref_str , alt_str, qual_str , filter_str , info, brk);
	}

	SUBREADprintf("Finished! There are %d common variants from the %d input file%s.\n\n", out_list_items, files - badfiles, (files - badfiles>1)?"s":"");
	free(output_list);
	
	HashTableDestroy(extra_qual_table);
	HashTableDestroy(extra_info_table);
	HashTableDestroy(extra_id_qual_filter_table);
	HashTableDestroy(extra_count_table);

	free(file_line);
}

void common_var_usage()
{
	SUBREADprintf("\nmergeVCF v%s\n", SUBREAD_VERSION);
	SUBREADputs("");
	SUBREADputs("  Calculating the intersection between multiple VCF files containing indels");
	SUBREADputs("  and/or SNPs.");
	SUBREADputs("");
	SUBREADputs("  The supported VCF files include VCF files generated by subread, subjunc,");
	SUBREADputs("  ExactSNP, subindel and other tools.");
	SUBREADputs("");
	SUBREADputs("Usage:");
	SUBREADputs("");
	SUBREADputs("  ./mergeVCF -o <output_vcf_file> input_vcf_1 {input_vcf_2} {input_vcf_3} ...");
	SUBREADputs("");
	SUBREADputs("Required arguments:");
	SUBREADputs("");
	SUBREADputs("  -o <file> : Name of the file containing the variants in all input VCF files.");
	SUBREADputs("");
	
}


#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int findCommonVariants(int argc, char ** argv)
#endif
{
	
	int c;
	int option_index = 0;

	optind = 0;
	opterr = 1;
	optopt = 63;
	output_file_name[0]=0;
	warning_reported_repeated=0;

	while((c = getopt_long (argc, argv, "o:h", propm_long_options, &option_index)) != -1)
	{
		switch(c){
			case 'o':
				strcpy(output_file_name, optarg);
				break;
			default:
				common_var_usage();
				return -1;
		}
	}

	if(!output_file_name[0])
	{
		common_var_usage();
		return -1;
	}
	else{
		output_file_pointer = f_subr_open(output_file_name,"w");
		if(output_file_pointer)
		{
			if( argc-optind )
				do_find_common(argv + optind, argc-optind);
			else
				SUBREADprintf("At least one input file is needed.\n\n");
			fclose(output_file_pointer);
		}
		else
		{
			SUBREADprintf("Unable to open the output file: '%s'\n", output_file_name);
			return -1;
		}
	}

	return 0;
}
