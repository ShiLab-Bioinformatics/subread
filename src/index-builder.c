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
  
  
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include <signal.h>
#include <unistd.h>
#include <time.h>
#include <sys/types.h>
#include "hashtable.h"
#include "gene-value-index.h"
#include "HelperFunctions.h"
#include "gene-algorithms.h"
#include "sorted-hashtable.h"
#include "seek-zlib.h"
#include "hashtable.h"
#include "input-files.h"

#define NO_GENE_DEBUG_
#define _GENE_DEBUG_SIZE_ 40000
#define MIN_READ_SPLICING 2000000


#define MAX_KEY_MATCH GENE_VOTE_SPACE 
int GENE_SLIDING_STEP = 3;
int IS_COLOR_SPACE = 0;
int MARK_NONINFORMATIVE_SUBREADS = 0;
int IS_FORCED_ONE_BLOCK = 0;
int ignore_bar_in_seqnames = 0;
double begin00_ftime = 0.0;

#define NEXT_READ 1
#define NEXT_FILE 2
#define FULL_PARTITION 4

long long estimate_memory_peak(unsigned int * bucket_sizes, unsigned int bucket_no, unsigned int tables ){
	long long peakuse = 0;
	int tn;
	for(tn=0; tn<tables; tn++){
		long long this_peak = 0;
		unsigned int bn;
		for(bn=0; bn<bucket_no; bn++)
			this_peak += bucket_sizes[bn + tn*bucket_no];
		peakuse = max(peakuse, this_peak);
	}
	//SUBREADprintf("EST MEMORY PEAK=%lld, bucks=%u, tabs=%u\n", peakuse, bucket_no, tables);
	return peakuse *(sizeof(short)+sizeof(gehash_data_t)) + sizeof(struct gehash_bucket)*bucket_no + sizeof(int) * tables*bucket_no;
}

void print_build_log(double finished_rate, double read_per_second, double expected_seconds, unsigned long long int total_reads)
{
	print_in_box( 81,0,0,"%4d%%%%, %3d mins elapsed, rate=%.1fk bps/s", (int)(finished_rate*100), (int)(miltime()-begin_ftime)/60, read_per_second/1000 );
}


#define MAX_BASES_IN_INDEX 4294900000.0

int build_gene_index(const char index_prefix [], char ** chro_files, int chro_file_number, int threshold, HashTable * huge_table, unsigned int * chro_lens, long long actual_bases, int for_measure_buckets, unsigned int ** bucket_sizes, unsigned int expected_hash_items, unsigned int bucket_no, unsigned int * total_tables){
	int file_number, table_no, status = NEXT_FILE;
	unsigned int offset, read_no, current_tab_items;
	long long int all_bases = guess_gene_bases(chro_files,chro_file_number);

	int chro_table_maxsize=100, dump_res = 0;
	unsigned int * read_offsets = malloc(sizeof(unsigned int) * chro_table_maxsize);
	char * read_names = malloc(MAX_READ_NAME_LEN * chro_table_maxsize);
	gehash_t table;
//	gehash_t huge_table;
	gene_value_index_t value_array_index;

	gene_input_t * ginp = malloc(sizeof(gene_input_t));
	begin_ftime = miltime();

//	SUBREADprintf ("Index items per partition = %u\n\n", expected_hash_items);

	if (chro_file_number > 199)
	{
		SUBREADprintf("ERROR: There are too many FASTA files. You may merge them into one FASTA file.\n");
		return -1;
	}

	if (strlen (index_prefix) > 290)
	{
		SUBREADprintf("ERROR: The path is too long. It should not be longer than 290 chars.\n");
		return -1;
	}

	if(all_bases<1)
	{
		SUBREADprintf("ERROR: File '%s' is inaccessible.\n", chro_files[-all_bases-1]);
		return -1;
	}

	current_tab_items = 0;
	table_no = 0;
	int padding_around_contigs = MAX_READ_LENGTH;
	if(0 == for_measure_buckets){
		if(gehash_create_ex(& table, expected_hash_items, 0, SUBINDEX_VER2, GENE_SLIDING_STEP, padding_around_contigs)) return 1;
		if(gvindex_init(&value_array_index, 0, (unsigned int) actual_bases))return 1;
	}
	if(1==for_measure_buckets){
		*bucket_sizes = realloc(*bucket_sizes, sizeof(int) * (table_no+1) * bucket_no);
		memset((*bucket_sizes) + table_no * bucket_no, 0, sizeof(int) *  bucket_no);
		memset(&table,0 , sizeof(gehash_t));
		memset(&value_array_index, 0, sizeof(gene_value_index_t));
		table.padding = padding_around_contigs;
		table.index_gap = GENE_SLIDING_STEP;
	}

	file_number = 0;
	offset = table.padding;
	read_no = 0;

	char * fn = malloc(3100);
	memset(read_offsets,0, chro_table_maxsize*sizeof(int));
	sprintf(fn, "%s.files", index_prefix);
	unlink(fn);

	status = NEXT_FILE;

	if(1==for_measure_buckets)
		print_in_box(80,0,0,"Estimate the index size...");
	else
		print_in_box(80,0,0,"Build the index...");

	{
		char window [16], last_color_base=-1, last_last_color_base=-1;
		int i, read_len = 0;
		unsigned int int_key = 0, array_int_key = 0;
		int skips=0, all_skips = 0;

		//Pre-fill
		while(1)
		{
			//Subread Cycle
			char next_char;
			if (status == NEXT_FILE)
			{
				if(file_number == chro_file_number)
				{
					FILE * fp;

					geinput_close(ginp);

					//SUBREADprintf ("Processing chromosome files ...\n");

					if(0==for_measure_buckets){
						sprintf (fn, "%s.%02d.%c.tab", index_prefix, table_no, IS_COLOR_SPACE?'c':'b');
						SUBREADfflush(stdout);
	
						if(!dump_res)dump_res |= gehash_dump(&table, fn);
	
						sprintf (fn, "%s.%02d.%c.array", index_prefix, table_no, IS_COLOR_SPACE?'c':'b');
						if(!dump_res)dump_res |= gvindex_dump(&value_array_index, fn);
						gvindex_destory(&value_array_index) ;
	
						gehash_destory(&table);
					}

					read_offsets[read_no-1] = offset + table.padding;

					for(i=table_no+1; i<100; i++)
					{
						sprintf(fn, "%s.%02d.%c.tab", index_prefix, i, IS_COLOR_SPACE?'c':'b');
						unlink(fn);
						sprintf(fn, "%s.%02d.%c.array", index_prefix, i, IS_COLOR_SPACE?'c':'b');
						unlink(fn);
					}

					sprintf (fn, "%s.reads", index_prefix);
					fp = f_subr_open(fn, "wb");
					for (i=0; i<read_no; i++)
						fprintf(fp, "%u\t%s\n", read_offsets[i], read_names+i*MAX_READ_NAME_LEN);

					fclose(fp);

					break;
				}
				else
				{
					if (file_number)
						geinput_close(ginp);
					geinput_open(chro_files[file_number++], ginp);
					status = NEXT_READ;
				}
			}
			if (status == NEXT_READ)
			{

				geinput_readline(ginp, fn, 0);

				if(read_no>0){
					read_offsets[read_no-1] = offset + table.padding;
					offset += 2*table.padding;
				}

				//printf("TTTXT FN=%s\n",fn);

				for(i=0;(fn[i+1] != ' ' && fn[i+1] != '\0' && fn[i+1] != '\t' && i<MAX_CHROMOSOME_NAME_LEN - 1 && ( fn[i+1] != '|'  || !ignore_bar_in_seqnames )); i++)
					*(read_names + MAX_READ_NAME_LEN*read_no + i) = fn[i+1];

				*(read_names + MAX_READ_NAME_LEN*read_no + i) = 0;

				sprintf(fn, "%s.files", index_prefix);
				FILE * fname_fp = f_subr_open(fn, "a");
				fprintf(fname_fp, "%s\t%s\t%ld\n", read_names+read_no*MAX_READ_NAME_LEN, ginp -> filename, ftell(ginp -> input_fp));
				fclose(fname_fp);
				
				for (i=0; i<16; i++)
				{
					char nch = geinput_next_char(ginp);
					if (nch == 'N') skips = 16;
					else if (skips>0) skips--;
					window[i] = nch;
				}

				read_len = 16;
				read_no ++;

				if(read_no >= chro_table_maxsize)
				{
					read_offsets = realloc(read_offsets, 2* chro_table_maxsize * sizeof(unsigned int));
					read_names = realloc(read_names, 2* chro_table_maxsize * MAX_READ_NAME_LEN);
					chro_table_maxsize *= 2;
				}

				if(IS_COLOR_SPACE)
				{
					int_key = genekey2color('A',window);
					array_int_key = genekey2int(window,GENE_SPACE_BASE);
					last_last_color_base = -1;
					last_color_base = window[15];
				}
				else
					array_int_key = int_key = genekey2int(window,GENE_SPACE_BASE);
	
			}
			if (status == FULL_PARTITION) 
			{
				int seek_back_reads ;

				if (read_len < MIN_READ_SPLICING)
				{
					seek_back_reads = read_len;
					offset -= seek_back_reads;
					offset += 16;
				}
				else
				{
					seek_back_reads = MIN_READ_SPLICING - 10;
					seek_back_reads -= seek_back_reads%3;
					seek_back_reads += 1;
					offset -= seek_back_reads;
					offset += 16;
				}

				if(0==for_measure_buckets){
					sprintf(fn, "%s.%02d.%c.tab", index_prefix, table_no, IS_COLOR_SPACE?'c':'b');
					SUBREADfflush(stdout);

					if(!dump_res)dump_res |= gehash_dump(&table, fn);
					sprintf(fn, "%s.%02d.%c.array", index_prefix, table_no, IS_COLOR_SPACE?'c':'b');
					if(!dump_res)dump_res |= gvindex_dump(&value_array_index, fn);
				}

				table_no ++;

				if(0==for_measure_buckets){
					gehash_destory(&table);
					gvindex_destory(&value_array_index);
				}

				read_len -= seek_back_reads;
				read_len += 16;

				while(seek_back_reads)
				{
					fseek(ginp -> input_fp, -1, SEEK_CUR);
					char bnch = fgetc(ginp -> input_fp);
					if ((bnch >='A' && bnch <= 'Z' ) || (bnch >='a' && bnch <= 'z' ) || bnch == '-' || bnch == 'N' || bnch=='.')seek_back_reads--;
					fseek(ginp -> input_fp, -1, SEEK_CUR);
				}

				for (i=0; i<16; i++)
				{
					char nch = geinput_next_char(ginp);
					if (nch == 'N' ) skips = 16;
					else if (skips>0) skips--;
					window[i] = nch;
				}
			
				if(IS_COLOR_SPACE)
				{
					int_key = genekey2color('A',window);
					array_int_key = genekey2int(window,GENE_SPACE_BASE);
					last_last_color_base = -1;
					last_color_base = window[15];
				}
				else
					array_int_key = int_key = genekey2int(window, GENE_SPACE_BASE);

				current_tab_items = 0;
				if(0==for_measure_buckets){
					if(gehash_create_ex(&table, expected_hash_items, 0, SUBINDEX_VER2, GENE_SLIDING_STEP, padding_around_contigs)) return 1;
					if(gvindex_init(&value_array_index, offset, (unsigned int) actual_bases)) return 1;
				}

				if(1==for_measure_buckets){
					//unsigned int * tmpp = *bucket_sizes;
					*bucket_sizes = realloc(*bucket_sizes, sizeof(int) * (table_no+1) * bucket_no);
					assert(NULL!= *bucket_sizes);
					memset((*bucket_sizes) + table_no * bucket_no, 0,   sizeof(int) * bucket_no);
				}
			}
	
			status = 0;

			if(skips || (IS_COLOR_SPACE && last_last_color_base<0))
				all_skips ++;
			else
			{
				int is_no_info = HashTableGet(huge_table, NULL+1+int_key)-NULL;
				if(is_no_info){
					//	SUBREADprintf("NOINFO KEY=%u AT %u\n", int_key, offset);
				}else {
					if(0 == for_measure_buckets)if(gehash_insert(&table, int_key, offset - (IS_COLOR_SPACE?1:0), (*bucket_sizes)+ table_no * bucket_no))return 1;
					if(1 == for_measure_buckets) gehash_try_insert_measure((*bucket_sizes)+ table_no * bucket_no , bucket_no, int_key);
					current_tab_items++;
				}
				if( 0 == for_measure_buckets )
					gvindex_set(&value_array_index, offset - (IS_COLOR_SPACE?0:0), array_int_key, padding_around_contigs);
			}

			if((!IS_FORCED_ONE_BLOCK) && current_tab_items >= expected_hash_items && (read_len > MIN_READ_SPLICING || read_len < 32))
			{
				status = FULL_PARTITION;
				continue;
			}

			for (i=0; i<GENE_SLIDING_STEP; i++)
			{
				next_char = geinput_next_char(ginp);
				if(next_char < 0) {
					if( 0 == for_measure_buckets ) gvindex_set(&value_array_index, offset - (IS_COLOR_SPACE?0:0), array_int_key, padding_around_contigs);

					if (next_char == -1) status = NEXT_READ;
					if (next_char == -2) status = NEXT_FILE;
					if (next_char == -3) return 0;
					break;
				}
				//SUBREADprintf("NEXT_CH=%c\n", next_char);

				if (next_char == 'N' )skips = 16;
				else if (skips>0){
					skips--;
					last_color_base = -1;
				}

				int_key = int_key << 2;


				if (IS_COLOR_SPACE)
				{
					if(last_color_base>0)
						int_key += chars2color(last_color_base, next_char);
					array_int_key = array_int_key << 2;
					array_int_key += base2int (next_char);

					last_last_color_base = last_color_base;
					last_color_base = next_char;
				}
				else
				{
					int_key += base2int (next_char); 
					array_int_key = int_key;
				}

				if(actual_bases > 12){
					if (offset>0 && offset % (actual_bases / 12) == 0)
					{
						double finished_rate = offset*1.0/actual_bases;
						double base_per_second = offset / (miltime()-begin_ftime);
						double ETA = (actual_bases-offset) / base_per_second;
						print_build_log(finished_rate,base_per_second,ETA, actual_bases);
						SUBREADfflush(stdout) ;
					}
				}

				offset ++;
				read_len ++;

				if(offset > 0xFFFFFFFDU)	
				{
					SUBREADprintf("ERROR: the provided reference sequences include more than 4 billion bases.\n") ;
					return -1;
				}

			}

		}
	}
	free(read_names);
	free(read_offsets);
	if(dump_res){
		SUBREADprintf("No index was built.\n");
		sprintf(fn, "%s.files", index_prefix);
		unlink(fn);
		sprintf(fn, "%s.reads", index_prefix);
		unlink(fn);
		int index_i;
		for(index_i = 0; index_i <= 99; index_i++){
			sprintf(fn, "%s.%02d.b.tab", index_prefix, index_i);
			unlink(fn);
			sprintf(fn, "%s.%02d.c.tab", index_prefix, index_i);
			unlink(fn);
			sprintf(fn, "%s.%02d.b.array", index_prefix, index_i);
			unlink(fn);
			sprintf(fn, "%s.%02d.c.array", index_prefix, index_i);
			unlink(fn);
		}
	}
	free(fn);
	free(ginp);
	*total_tables = table_no+1;
	return 0;
}

int add_repeated_subread(gehash_t * tab , unsigned int subr, unsigned char ** huge_index)
{
	int huge_byte = (subr>>1);
	int huge_offset = (subr % 2) * 4;
	int ind_id = (huge_byte >> 24) & 127;

	unsigned int byte_value = huge_index[ ind_id ][huge_byte&0xffffff] ;

	int huge_value = (byte_value>> huge_offset) & 0xf;
	if(huge_value <15){
		huge_value ++;
		byte_value = (byte_value & (~(0xf << huge_offset))) | (huge_value << huge_offset);
		huge_index[ ind_id ][huge_byte&0xffffff] = (unsigned char)(byte_value&0xff);
		return 0;
	}
	unsigned int times = 0;
	int matched = gehash_get(tab, subr, &times, 1);
	if(matched)
		gehash_update(tab, subr, times+1);
	else
		if(gehash_insert(tab, subr, 16, NULL)) return 1;
	
	return 0;
}

int scan_gene_index(const char index_prefix [], char ** chro_files, int chro_file_number, int threshold, HashTable *huge_table, long long * actual_total_bases_inc_marging)
{
	int file_number, i ,j;
	int status = NEXT_FILE;
	unsigned int offset, read_no;
	long long int guess_all_bases = guess_gene_bases(chro_files,chro_file_number);
	int padding_around_contigs = MAX_READ_LENGTH;
	*actual_total_bases_inc_marging = padding_around_contigs;

	gehash_t occurrence_table;
	unsigned char * huge_index[128];

	for(i=0;i<128;i++) {
		huge_index[i] = calloc(1024*1024*16,1);
		if(NULL == huge_index[i]){ 
			SUBREADprintf("ERROR: No memory can be allocated.\n");
			return -1;
		}
	}

	if(gehash_create_ex(&occurrence_table, 500000, 0, SUBINDEX_VER0, 1, 0)) return 1;

	gene_input_t * ginp = malloc(sizeof(gene_input_t));

	print_in_box(80,0,0,"Scan uninformative subreads in reference sequences ...");

	if (chro_file_number > 199)
	{
		SUBREADprintf("ERROR: There are too many FASTA files. You may merge them into one FASTA file.\n");
		return -1;
	}

	if (strlen (index_prefix) > 290)
	{
		SUBREADprintf("ERROR: The path is too long. It should not be longer than 290 chars.\n");
		return -1;
	}
	if(guess_all_bases<1)
	{
		SUBREADprintf("ERROR: File '%s' is inaccessible.\n", chro_files[-guess_all_bases-1]);
		return -1;
	}


	file_number = 0;
	offset = 0;
	read_no = 0;

	char * fn = malloc(3100);
	sprintf(fn, "%s.files", index_prefix);
	unlink(fn);

	status = NEXT_FILE;


	{
		char window [16], last_color_base=-1, last_last_color_base=-1;
		int i, read_len = 0;
		unsigned int int_key = 0, array_int_key = 0;
		int skips=0, all_skips = 0;

		//Pre-fill
		while(1)
		{
			//Subread Cycle
			char next_char;
			//if(status > 0)SUBREADprintf("BDSSTATUS = %d\n", status);

			if (status == NEXT_FILE)
			{
				(*actual_total_bases_inc_marging)+=padding_around_contigs;
				if(file_number == chro_file_number)
				{
					geinput_close(ginp);

					break;
				}
				else
				{
					if (file_number)
						geinput_close(ginp);
					geinput_open(chro_files[file_number++], ginp);
					status = NEXT_READ;
				}
			}
			if (status == NEXT_READ)
			{
				(*actual_total_bases_inc_marging)+=2*padding_around_contigs;

				geinput_readline(ginp, fn, 0);
				//SUBREADprintf("HEADER_SCAN '''%s'''\n", fn);

				for (i=0; i<16; i++)
				{
					char nch = geinput_next_char(ginp);
					if (nch == 'N') skips = 16;
					else if (skips>0) skips--;
					window[i] = nch;
				}
				read_len = 16;
				read_no ++;

				if(IS_COLOR_SPACE)
				{
					int_key = genekey2color('A',window);
					array_int_key = genekey2int(window,GENE_SPACE_BASE);
					last_last_color_base = -1;
					last_color_base = window[15];
				}
				else
					array_int_key = int_key = genekey2int(window,GENE_SPACE_BASE);
	
			}
	
			status = 0;

			if(skips || (IS_COLOR_SPACE && last_last_color_base<0))
				all_skips ++;
			else
			{
				int rv = add_repeated_subread(&occurrence_table, int_key, huge_index);
				if(rv<0) {
					return -1;
				}
			}


			for (i=0; i<GENE_SLIDING_STEP; i++)
			{
				next_char = geinput_next_char(ginp);
				if(next_char < 0)
				{
					if (next_char == -1) status = NEXT_READ;
					if (next_char == -2) status = NEXT_FILE;
					if (next_char == -3) return 0;
					break;
				}

				if (next_char == 'N' )skips = 16;
				else if (skips>0){
					skips--;
					last_color_base = -1;
				}

				int_key = int_key << 2;

				if (IS_COLOR_SPACE)
				{
					if(last_color_base>0)
						int_key += chars2color(last_color_base, next_char);
					array_int_key = array_int_key << 2;
					array_int_key += base2int (next_char);

					last_last_color_base = last_color_base;
					last_color_base = next_char;
				}
				else
				{
					int_key += base2int (next_char); 
					array_int_key = int_key;
				}

				offset ++;
				(*actual_total_bases_inc_marging)++;
				read_len ++;

				if(offset > 0xFFFFFFFDU)	
				{
					SUBREADprintf("ERROR: the provided reference sequences include more than 4 billion bases.\n") ;
					return -1;
				}

			}

		}
	}



	free(fn);

	for (i=0; i<occurrence_table.buckets_number; i++)
	{
		struct gehash_bucket * current_bucket = &(occurrence_table.buckets[i]);

		if(current_bucket -> current_items>=1)
		{
			for (j=0; j<current_bucket -> current_items; j++)
			{
				if(current_bucket -> item_values [j] > threshold)
				{
					HashTablePut(huge_table,NULL+1+ current_bucket -> item_keys[j], NULL+current_bucket -> item_values [j]);
				}
			}
		}
	}


	for(i=0;i<128;i++)
		if(huge_index[i])  free(huge_index[i]);

	free(ginp);
	gehash_destory(&occurrence_table);

	if(huge_table -> numOfElements)
	{
		print_in_box(80,0,0,"%llu uninformative subreads were found.", huge_table -> numOfElements);
		print_in_box(80,0,0,"These subreads were excluded from index building.");
	}

	return 0;
}


char *rtrim(char *s)
{
	char* back = s + strlen(s);
	while(isspace(*--back));
	*(back+1) = '\0';
	return s;
}

typedef struct{
	char * filename;
	unsigned int line_no;
	int byte_in_line;
} format_check_context_t;

int ERROR_FOUND_IN_FASTA = 0;
void check_and_convert_warn(format_check_context_t * fcc, char * msg, FILE * log_fp){
	ERROR_FOUND_IN_FASTA += 1;
	fprintf(log_fp, "A format issue is found in the %u-th line in file '%s': %s.\n", fcc->line_no, fcc->filename, msg);
}
void check_and_convert_warnOLD(char * FN, long long int fpos_line_head, unsigned line_no, int line_pos, char * msg, FILE * log_fp)
{
#define CHAR_ESC 27
	int x1,brs=0;
	long long int back_search_ptr;
	char * line_buf = malloc(MAX_READ_LENGTH+1);

	ERROR_FOUND_IN_FASTA += 1;

	fprintf(log_fp,"\n");

	//fprintf(log_fp,"%c[33m", CHAR_ESC);
	for(x1=0;x1<81;x1++)
		fprintf(log_fp,"=");
	fprintf(log_fp,"\n");
	//fprintf(log_fp,"%c[32m", CHAR_ESC);
	fprintf(log_fp,"Input file '%s':\n", FN);
	//fprintf(log_fp,"%c[31m", CHAR_ESC);
	fprintf(log_fp,"%s\n", msg);
	//fprintf(log_fp,"%c[33m", CHAR_ESC);
	for(x1=0;x1<81;x1++)
		fprintf(log_fp,".");
	fprintf(log_fp,"\n");
	//fprintf(log_fp,"%c[37m", CHAR_ESC);
	

	FILE * warn_fp = f_subr_open(FN, "r");

	for(back_search_ptr = fpos_line_head - 1; back_search_ptr>=0; back_search_ptr--)
	{
		int nch;
		fseeko(warn_fp, back_search_ptr, SEEK_SET);
		nch = fgetc(warn_fp);
		if(nch == '\n') brs++;
		if(brs >2) break;
		fseeko(warn_fp, back_search_ptr, SEEK_SET);
	}

	if(back_search_ptr<=0)brs++;

	int print_line_no = line_no - brs + 1;
	while(1)
	{
		char * ret = fgets(line_buf, MAX_READ_LENGTH, warn_fp);
		if(!ret)break;

		/*
		if(ftello(warn_fp) > fpos_line_head)
			fprintf(log_fp,"%c[9m%c[31m", CHAR_ESC, CHAR_ESC);
		else
			fprintf(log_fp,"%c[29m%c[37m", CHAR_ESC, CHAR_ESC);
		*/
		fprintf(log_fp," % 9d ", print_line_no++);

		rtrim(line_buf);
		
		fprintf(log_fp,"%s%s\n",line_buf,strlen(line_buf)<16?"              ":"");
		if(ftello(warn_fp) > fpos_line_head)
			break;
	}
	for(x1=0;x1<line_pos+11;x1++)
		fprintf(log_fp," ");
	fprintf(log_fp,"^\n");

	//fprintf(log_fp,"%c[29m%c[37m", CHAR_ESC, CHAR_ESC);
	for(x1=0;x1<2;x1++)
	{
		char * ret = fgets(line_buf, MAX_READ_LENGTH, warn_fp);
		if(!ret)break;
		fprintf(log_fp," % 9d ", print_line_no++);
		fprintf(log_fp,"%s",line_buf);
	}
	fclose(warn_fp);
	//fprintf(log_fp,"%c[33m", CHAR_ESC);
	for(x1=0;x1<81;x1++)
		fprintf(log_fp,"=");
	fprintf(log_fp,"\n");
	fprintf(log_fp,"\n");

	//fprintf(log_fp,"%c[0m", CHAR_ESC);
	free(line_buf);
}

int check_and_convert_FastA(char ** input_fas, int fa_number, char * out_fa, unsigned int ** chrom_lens, FILE * log_fp, char * log_fn)
{
	int is_R_warnned = 0, is_repeated_chro= 0;
	char * line_buf = malloc(MAX_READ_LENGTH);
	char * read_head_buf = malloc(MAX_READ_LENGTH * 3);
	unsigned int inp_file_no, line_no;
	int written_chrs = 0, is_disk_full = 0;
	int chrom_lens_max_len = 100;
	int chrom_lens_len = 0;
	ERROR_FOUND_IN_FASTA = 0;
	FILE * out_fp = f_subr_open(out_fa,"wb");
	format_check_context_t fcc;
	memset(&fcc,0,sizeof(fcc));

	if(!out_fp)
	{
		SUBREADprintf("ERROR: the output directory is not writable, but the index builder needs to create temporary files in the current directory. Please change the working directory and rerun the index builder.\n");
		return -1;
	}

	(*chrom_lens) = malloc(chrom_lens_max_len*sizeof(unsigned int));
	memset((*chrom_lens), 0, chrom_lens_max_len*sizeof(unsigned int));


	HashTable * rep_name_table = HashTableCreate(9999);
	HashTableSetDeallocationFunctions(rep_name_table, free, NULL);
	HashTableSetKeyComparisonFunction(rep_name_table, (int (*) (const void *, const void *)) strcmp);
	HashTableSetHashFunction(rep_name_table , fc_chro_hash);
	
	print_in_box( 80,0,0,"Check the integrity of provided reference sequences ...");
	for(inp_file_no = 0; inp_file_no < fa_number; inp_file_no++)
	{
		autozip_fp fafp;
		memset(&fafp, 0, sizeof(fafp));
		int faret = autozip_open(input_fas[inp_file_no], &fafp);
		fcc.filename = input_fas[inp_file_no];

		if(faret <0) {
			SUBREADprintf("ERROR: Input file '%s' is not found or is not accessible. No index was built.\n", input_fas[inp_file_no]);
			HashTableDestroy(rep_name_table);
			return -1;
		}

		line_no = 0;
		int is_head_written=0;
		read_head_buf[0]=0;
		while(1) {
			unsigned int read_len = 0;
			int ret = autozip_gets(&fafp, line_buf, MAX_READ_LENGTH);
			if(ret<1) break;
			if(ret >= MAX_READ_LENGTH-2){
				SUBREADprintf("ERROR: A fasta file cannot have a line longer than 1000 bytes. You need to split a very long line into many lines.\n");
				return -1;
			}

			line_no ++;
			fcc.line_no = line_no;
			int line_buf_len = strlen(line_buf);

			for(; line_buf[line_buf_len-1] == '\r' || line_buf[line_buf_len-1]=='\n' ;line_buf_len--)
			{
				if(line_buf[line_buf_len-1]=='\r')
				{
					if(!is_R_warnned)
					{
						is_R_warnned=1;
						fcc.byte_in_line = line_buf_len-1;
						check_and_convert_warn(&fcc ,"This line ends with '\\r\\n'. It is not a problem for building the index but we suggest to use Unix-styled line breaks.", log_fp);
					}
				}
				line_buf[line_buf_len-1] =0;
			}
			//if(strchr(line_buf,'>') || strchr(line_buf,'\r') || strchr(line_buf,'\n'))SUBREADprintf("TESTING: string contains '''%s'''  at %I64u\n", line_buf, ftello(fafp.plain_fp));

			if(line_buf_len<1)
			{
				fcc.byte_in_line = -1;
				check_and_convert_warn(&fcc ,"This line is empty. This is not allowed in the FASTA file.", log_fp);
				continue;
			}

			if(line_buf[0]=='>')
			{
				if(line_no>1 &&!is_head_written)
				{
					fcc.byte_in_line = 0;
					check_and_convert_warn(&fcc,"This sequence has less than 16 bases. It is ignored in the index because no subreads can be extracted.", log_fp);
				}
				is_head_written = 0;
				read_len = 0;
				read_head_buf[0]=0;

				strcat(read_head_buf, line_buf);
				strcat(read_head_buf, "\n");

				int chro_name_end=1;
				while(line_buf[chro_name_end])
				{
					if(line_buf[chro_name_end]==' '||line_buf[chro_name_end]=='|' || line_buf[chro_name_end]=='\t')break;
					chro_name_end++;
				}

				line_buf[chro_name_end] = 0;
				int is_exist = HashTableGet(rep_name_table , line_buf+1) - NULL;
				if(is_exist)
				{
					SUBREADprintf("ERROR: repeated chromosome name '%s' is observed in the FASTA file(s).\nThe index was NOT built.\n", line_buf+1);
					is_repeated_chro=1;
					break;
				}

				char * keymem = malloc(chro_name_end);
				strcpy(keymem, line_buf+1);
				HashTablePut(rep_name_table, keymem, NULL+1);

			}
			else if(line_no<1)
			{
				fcc.byte_in_line=0;
				check_and_convert_warn(&fcc ,"This file is not started with a header line. It seems not to be a FASTA file.", log_fp);
			}
			else
			{
				int xk2;
				for(xk2=0; xk2 < line_buf_len; xk2++)
				{
					int nextch = line_buf[xk2];
					int lowerch = tolower(nextch);
					if(!( lowerch == 'a' || lowerch == 't' || lowerch == 'g' || lowerch == 'c' || nextch == '.' || nextch=='-' || lowerch=='n'))
					{	
						fcc.byte_in_line=xk2;
						check_and_convert_warn(&fcc, "The non-ACGT base was converted to an 'A'.", log_fp);
						line_buf[xk2] = 'A';
					}
					else if(nextch == '.' || nextch=='-' ||  lowerch=='n')
						line_buf[xk2] = 'A';
					else
						line_buf[xk2] = toupper(nextch);
				}

				read_len += line_buf_len;

				if(read_len > 16 && !is_head_written)
				{
					fputs(read_head_buf, out_fp);
					written_chrs++;
					is_head_written = 1;

					chrom_lens_len++;
					if((chrom_lens_max_len-1) <= chrom_lens_len)
					{
						(*chrom_lens) = realloc((*chrom_lens), 2*chrom_lens_max_len*sizeof(unsigned int));
						chrom_lens_max_len*=2;
					}
				}

				if(is_head_written)
				{
					int line_buf_len = strlen(line_buf);
					int writen_len = fprintf(out_fp,"%s\n", line_buf);
					if(writen_len < line_buf_len){
						SUBREADprintf("ERROR: unable to write into the temporary file. Please check the free space of the output directory.\n");
						is_disk_full = 1;
						break;
					}
					(*chrom_lens)[chrom_lens_len-1] = read_len;
					(*chrom_lens)[chrom_lens_len] = 0;
				}
				else
				{
					strcat(read_head_buf, line_buf);
					strcat(read_head_buf, "\n");
				}
				
			}
		}

		autozip_close(&fafp);
		if(is_disk_full) break;
	}


	HashTableDestroy(rep_name_table);
	free(line_buf);
	free(read_head_buf);
	fclose(out_fp);

	if(!written_chrs){
		SUBREADprintf("ERROR: No index was built because there were no subreads extracted. A chromosome needs at least 16 bases to be indexed.");
		return 1;
	}

	if(is_repeated_chro|| is_disk_full)
		return 1;

	if(ERROR_FOUND_IN_FASTA)
	{
		print_in_box( 80,0,0,"There were %d notes for reference sequences.", ERROR_FOUND_IN_FASTA);
		print_in_box( 89,0,0,"The notes can be found in the log file, %c[36m'%s'%c[0m.", CHAR_ESC, log_fn, CHAR_ESC);
	}
	else	print_in_box( 80,0,0,"No format issues were found");

	return 0;
}

char * tmp_file_for_signal;

void SIGINT_hook(int param)
{
	#ifdef MAKE_STANDALONE
	if(tmp_file_for_signal[0])
	{
		unlink(tmp_file_for_signal);
		SUBREADprintf("\n\nReceived a terminal signal. The temporary file was removed. The index was NOT built successfully. Please DO NOT use the new index until they are rebuilt.\n\n");
	}

	exit(param);
	#endif
}

static struct option ib_long_options[]={
	{0, 0, 0, 0}
};

#ifdef MAKE_STANDALONE
int main(int argc,char ** argv)
#else
int main_buildindex(int argc,char ** argv)
#endif
{
	int threshold = 100, optindex=0;
	int memory_limit;	// 8000 MBytes
	char output_file[MAX_FILE_NAME_LENGTH], c, tmp_fa_file[MAX_FILE_NAME_LENGTH], log_file_name[MAX_FILE_NAME_LENGTH+20];
	char *ptr_tmp_fa_file[1];
	unsigned int * chromosome_lengths;
	begin00_ftime = miltime();

	if(sizeof(char *)>4) memory_limit=8000;
	else memory_limit=3000;

	ptr_tmp_fa_file[0]=tmp_fa_file;
	output_file[0] = 0;
	tmp_fa_file[0] = 0;
	tmp_file_for_signal = tmp_fa_file;

	IS_FORCED_ONE_BLOCK = 0;
	GENE_SLIDING_STEP = 3;
	IS_COLOR_SPACE = 0;

	SUBREADprintf("\n");

	optind = 0;
	
	while ((c = getopt_long (argc, argv, "kvcBFM:o:f:Db?", ib_long_options, &optindex)) != -1)
		switch(c)
		{
			case 'b':
				ignore_bar_in_seqnames= 1;
				break;
			case 'B':
				IS_FORCED_ONE_BLOCK =1;
				break;
			case 'F':
				GENE_SLIDING_STEP =1;
				break;
			case 'v':
				core_version_number("Subread-buildindex");
				return 0;
			case 'c':
				IS_COLOR_SPACE = 1;
				break;
			case 'M':
				memory_limit = atoi(optarg);
				break;
			case 'f':
				threshold = atoi(optarg);
				break;
			case 'o':
				strncpy(output_file, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'k':
				if(memcmp(SUBREAD_VERSION , "1.3.",4)==0)
				{
					SUBREADprintf("The \"-k\" option is not supported in version 1.3.x. ");
					return -1;
				}
				MARK_NONINFORMATIVE_SUBREADS = 1;
				break;	
			case '?':
				return -1 ;
		}

	if (argc == optind || !output_file[0])
	{
		SUBREADprintf("Version %s\n\n", SUBREAD_VERSION);

	/*
		SUBREADputs("Usage:");
		SUBREADputs("");
		SUBREADputs(" ./subread-buildindex [options] -o <basename> {FASTA file1} [FASTA file2] ...");
		SUBREADputs("");
		SUBREADputs("Required arguments:");
		SUBREADputs("");
		SUBREADputs("    -o <basename>   base name of the index to be created");
		SUBREADputs("");
		SUBREADputs("Optional arguments:");
		SUBREADputs("");
		SUBREADputs("    -G              build a gapped index for the reference genome. 16bp subreads");
		SUBREADputs("");
		SUBREADputs("    -M <int>        size of requested memory(RAM) in megabytes. Index is split into blocks if necessary.");
		SUBREADputs("");
		SUBREADputs("    -f <int>        specify the threshold for removing uninformative subreads");
		SUBREADputs("                    (highly repetitive 16mers in the reference). 24 by default.");
		SUBREADputs("");
		SUBREADputs("    -c              build a color-space index.");
		SUBREADputs("");
		SUBREADputs("    -v              output version of the program.");
		SUBREADputs("");
		SUBREADputs("For more information about these arguments, please refer to the User Manual.\n");
	*/
		 SUBREADputs("Usage:");
		 SUBREADputs("");
		 SUBREADputs(" ./subread-buildindex [options] -o <basename> {FASTA[.gz] file1}\\");
		 SUBREADputs("      [FASTA[.gz] file2] ...");
		 SUBREADputs("");
		 SUBREADputs("Required arguments:");
		 SUBREADputs("");
		 SUBREADputs("    -o <basename>   base name of the index to be created");
		 SUBREADputs("");
		 SUBREADputs("Optional arguments:");
		 SUBREADputs("");
		 SUBREADputs("    -F              build a full index for the reference genome. 16bp subreads");
		 SUBREADputs("                    will be extracted from every position of the reference");
		 SUBREADputs("                    genome. Size of the index is typically 3 times the size of");
		 SUBREADputs("                    index built from using the default setting.");
		 SUBREADputs("");
		 SUBREADputs("    -B              create one block of index. The built index will not be split");
		 SUBREADputs("                    into multiple pieces. This makes the largest amount of");
		 SUBREADputs("                    memory be requested when running alignments, but it enables");
		 SUBREADputs("                    the maximum mapping speed to be achieved. This option");
		 SUBREADputs("                    overrides -M when it is provided as well.");
		 SUBREADputs("");
		 SUBREADputs("    -M <int>        size of requested memory(RAM) in megabytes, 8000 by default.");
		 SUBREADputs("");
		 SUBREADputs("    -f <int>        specify the threshold for removing uninformative subreads");
		 SUBREADputs("                    (highly repetitive 16mers in the reference). 100 by default.");
		 SUBREADputs("");
		 SUBREADputs("    -c              build a color-space index.");
		 SUBREADputs("");
		 SUBREADputs("    -v              output version of the program.");
		 SUBREADputs("");
		 SUBREADputs("For more information about these arguments, please refer to the User Manual.\n");

		return -1 ;
	}

	if(threshold < 5) {
		SUBREADprintf("ERROR: The threshold of non-informative reads cannot be less than 5.\n");
		return -1;
	}

	// **********************************************
	//	print config summary
	// **********************************************

	print_subread_logo();
	size_t free_mem=-1, total_mem=-1;
	if(get_free_total_mem( &total_mem, &free_mem )){
		free_mem=-1;
		total_mem=-1;
	}

	SUBREADputs("");
	print_in_box(80, 1, 1, "setting");
	print_in_box(80, 0, 1, "");
	print_in_box(80, 0, 0, "               Index name : %s", get_short_fname(output_file));
	print_in_box(80, 0, 0, "              Index space : %s", IS_COLOR_SPACE?"color space":"base space");

	if(IS_FORCED_ONE_BLOCK)
	{
		print_in_box(80, 0, 0, "              Index split : no-split");
		memory_limit = GENE_SLIDING_STEP==1?22000:11500;
	}
	else
	{
		if(memory_limit > 12000 && GENE_SLIDING_STEP>2)
		{
			print_in_box(80, 0, 0, "                   Memory : %u -> %u Mbytes", memory_limit, 12000);
			memory_limit = 12000;
		}
		else
			print_in_box(80, 0, 0, "                   Memory : %u Mbytes", memory_limit);
	}
	print_in_box(80, 0, 0, "         Repeat threshold : %d repeats", threshold);
	print_in_box(80, 0, 0, "             Gapped index : %s", GENE_SLIDING_STEP>1?"yes":"no");
	print_in_box(80, 0, 0, "");
	if(free_mem>0)print_in_box(80, 0, 0, "      Free / total memory : %.1fGB / %.1fGB", free_mem*1./1024/1024/1024, total_mem*1./1024/1024/1024);
	print_in_box(80, 0, 0, "");	
	print_in_box(80, 0, 0, "              Input files : %d file%s in total",  argc - optind, (argc - optind>1)?"s":"");

	int x1;
	for(x1=0;x1< argc - optind; x1++)
	{
		char * fasta_fn = *(argv+optind+x1);
		int f_type = probe_file_type_fast(fasta_fn);
		char o_char = '?';
		if(f_type == FILE_TYPE_FASTA){
			o_char = 'o';
		}
		if(f_type == FILE_TYPE_GZIP_FASTA){
			o_char = 'o';
		}
		print_in_box(94, 0, 0, "                            %c[32m%c%c[36m %s%c[0m", CHAR_ESC, o_char, CHAR_ESC,  get_short_fname(fasta_fn) , CHAR_ESC);
	}
	print_in_box(80, 0, 0, "");
	
	#ifndef __MINGW32__
	if(free_mem>0 && free_mem < 3*1024ll*1024*1024){
		print_in_box(80, 0, 0, "");
		print_in_box( 80, 0, 0, "  WARNING: the free memory is lower than 3.0GB." );
		print_in_box( 80, 0, 0, "           the program may run very slow or crash." );
		print_in_box(80, 0, 0, "");
	}
	#endif

	print_in_box(80, 2, 1, "");
	SUBREADputs("");

	print_in_box(80, 1, 1, "Running");
	print_in_box(80, 0, 0, "");

	for(x1=0;x1< argc - optind; x1++)
	{
		char * fasta_fn = *(argv+optind+x1);
		int f_type = probe_file_type_fast(fasta_fn);
		if(f_type != FILE_TYPE_GZIP_FASTA && f_type != FILE_TYPE_FASTA && f_type != FILE_TYPE_NONEXIST){
			SUBREADprintf("ERROR: '%s' is not a Fasta file.\n", fasta_fn);
			STANDALONE_exit(-1);
		}
	}

	for(x1 = strlen(output_file); x1 >=0; x1--){
		if(output_file[x1]=='/'){
			memcpy(tmp_fa_file, output_file, x1);
			tmp_fa_file[x1]=0;
			break;
		}
	}
	if(tmp_fa_file[0]==0)strcpy(tmp_fa_file, "./");

	#ifdef __MINGW32__
	sprintf(tmp_fa_file+strlen(tmp_fa_file), "/subread-index-sam-%06u-%06d", getpid(),(int)(time(NULL) % 1000000));
	#else
	sprintf(tmp_fa_file+strlen(tmp_fa_file), "/subread-index-sam-%06u-XXXXXX", getpid());
	int tmpfdd = mkstemp(tmp_fa_file);
	if(tmpfdd == -1){
		SUBREADprintf("ERROR: cannot create temp file\n");
		return -1;
	}
	#endif

	sprintf(log_file_name, "%s.log", output_file);
	FILE * log_fp = f_subr_open(log_file_name,"wb");

	signal (SIGTERM, SIGINT_hook);
	signal (SIGINT, SIGINT_hook);


	int ret = check_and_convert_FastA(argv+optind , argc - optind, tmp_fa_file, &chromosome_lengths, log_fp, log_file_name);
	if(log_fp)
		fclose(log_fp);

	if(!ret)
	{
		long long actual_bases=0;
		HashTable * huge_table;
		huge_table = HashTableCreate(200000);
		unsigned int expected_hash_items = (unsigned int)(memory_limit * 1024.0 / 8.) * 1024 ;
		unsigned int * bucket_sizes = NULL, total_tables=0;
		unsigned int bucket_no = calculate_buckets_by_size(expected_hash_items, SUBINDEX_VER2, 0, GENE_SLIDING_STEP);
		ret = ret || scan_gene_index(output_file, ptr_tmp_fa_file , 1, threshold, huge_table, &actual_bases);
		ret = ret || build_gene_index(output_file, ptr_tmp_fa_file , 1, threshold, huge_table, chromosome_lengths, actual_bases, 1, &bucket_sizes, expected_hash_items, bucket_no, &total_tables);

		long long estm_need_memory = estimate_memory_peak( bucket_sizes, bucket_no, total_tables );
		long long needed_mem =  estm_need_memory + huge_table->numOfElements * sizeof(KeyValuePair) + 560ll*1024*1024 + actual_bases / ((GENE_SLIDING_STEP > 1 ) ?  7:5);

		if(free_mem>0 && free_mem < needed_mem){
			print_in_box(80, 0, 1, "");
			print_in_box(80, 0, 1, "WARNING: available memory is lower than %.1f GB." ,needed_mem*1./1024l/1024/1024); 
			print_in_box(80, 0, 1, "         The program may run very slow.");
			print_in_box(80, 0, 1, "Build a gapped index and/or split index into blocks to reduce memory use.");
			print_in_box(80, 0, 1, "");
		}else
			print_in_box(80, 0, 0, "%.1f GB of memory is needed for index building." ,needed_mem*1./1024l/1024/1024); 


		ret = ret || build_gene_index(output_file, ptr_tmp_fa_file , 1, threshold, huge_table, chromosome_lengths, actual_bases, 0, &bucket_sizes, expected_hash_items, bucket_no, &total_tables);

		if(!ret){
			print_in_box(80, 0, 1, "Total running time: %.1f minutes.", (miltime()-begin00_ftime)/60);
			print_in_box(89, 0, 1, "Index %c[36m%s%c[0m was successfully built.", CHAR_ESC, output_file, CHAR_ESC);
		}
		HashTableDestroy(huge_table);
		free(chromosome_lengths);
		free(bucket_sizes);
	}

	unlink(tmp_fa_file);
	tmp_fa_file[0]=0;

	print_in_box(80, 0, 0, "");
	print_in_box(80, 2, 1, "");
	SUBREADputs("");
	return ret;
}

