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
  
  
#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef FREEBSD
#include <sys/timeb.h>
#endif

#include "subread.h"
#include "input-files.h"
#include "gene-algorithms.h"
#include "sorted-hashtable.h"

void non_func(const char * fmt, ...)
{
}


void subread_lock_release(subread_lock_t * lock)
{
	#ifdef MACOS
	pthread_mutex_unlock(lock);
	#else
	pthread_spin_unlock(lock);
	#endif
}
void subread_lock_occupy(subread_lock_t * lock)
{

	#ifdef MACOS
	pthread_mutex_lock(lock);
	#else
	pthread_spin_lock(lock);
	#endif
}

void subread_destroy_lock(subread_lock_t * lock) {
	#ifdef MACOS
	pthread_mutex_destroy(lock);
	#else
	pthread_spin_destroy(lock);
	#endif
}

void subread_init_lock(subread_lock_t * lock)
{
	#ifdef MACOS
	pthread_mutex_init(lock, NULL);
	#else
	pthread_spin_init(lock, PTHREAD_PROCESS_PRIVATE);
	#endif

}

//int EXON_MAX_CIGAR_LEN = 26;


//#define indel_debug

#ifdef indel_debug
#define ddprintf SUBREADprintf
#define ddfflush fflush
#else
#define ddprintf non_func
#define ddfflush(a) 
#endif

#define get_base_quality_score(quality_chr, quality_scale)  ((quality_scale)==QUALITY_SCALE_NONE?-0.01: ((quality_scale) == QUALITY_SCALE_LINEAR?(	(quality_chr) - '@')*0.01-0.03: correct_rate_table[(quality_chr) - '@']))


double begin_ftime;

unsigned char get_next_char(FILE * fp)
{
	int find_br = 0;
	while (!feof(fp))
	{
		char nch;
		nch = fgetc(fp);
		if (find_br)
		{
			if (nch == '\n')
				find_br=0;
		}
		else if (nch == '>')
			find_br = 1;
		else if (nch > 32)
			return nch;
	}
	return 0;
}

int is_valid_subread(const char * read_str)
{
	return 1;
	int i;
	for (i=0; i<16; i++)
		if(read_str[i] == 'N' || read_str[i] == '.')
			return 0;
	return 1;
}

float read_quality_score(char * qualityb, int rl , int format)
{
	int i;
	int qual = 0;
	int testlen = 0;

	char base;

	if(format==FASTQ_PHRED64)
		base = 'B';
	else	base = '#';

	for(i=0; i<rl; i++)
	{
		int testv = qualityb[i] - base;
		if(testv > 1)
		{
			qual += testv;
			testlen ++;
		}
	}
	return qual*1./testlen;
}

int min_matched_bases(char * qualityb, int rl , int format, float match_score)
{
	if(!(qualityb) || !qualityb[0]) return 0;
	int i;

	char base;
	int ret = 0;

	if(format==FASTQ_PHRED64)
		base = 'B';
	else	base = '#';

	for(i=0; i<rl; i++)
		ret += (qualityb[i] - base)<=5;

	return (rl - ret*3/4)*match_score;
}
int bad_quality_base_number(char * qualityb, int rl , int format)
{
	if(!(qualityb) || !qualityb[0]) return 0;
	int ret = 0, i;
	if(FASTQ_PHRED64==format)
	{
		for(i=0; i<rl; i++)
			if(qualityb[i] <='F') ret ++;
	}
	else
		for(i=0; i<rl; i++)
			if(qualityb[i] <='#'+4) ret ++;

		
	return ret;
}


double correct_rate_table [] ={ -1.58147375341 , -0.99684304401 , -0.69552447133 , -0.50767587370 , -0.38013040807 , -0.28926818720 , -0.22255151597 , -0.17255657291 , -0.13455196029 , -0.10536051566 , -0.08276530267 , -0.06517417320 , -0.05141827416 , -0.04062484422 , -0.03213357402 , -0.02543972753 , -0.02015436476 , -0.01597586925 , -0.01266917021 , -0.01005033585 , -0.00797499828 , -0.00632956293 , -0.00502447389 , -0.00398901727 , -0.00316728823 , -0.00251504651 , -0.00199725550 , -0.00158615046 , -0.00125971852 , -0.00100050033 , -0.00079464388 , -0.00063115648 , -0.00050131287 , -0.00039818644 , -0.00031627778 , -0.00025122020 , -0.00019954614 , -0.00015850188 , -0.00012590047 , -0.00010000500 , -0.00007943598 , -0.00006309773 , -0.00005011998 , -0.00003981151 , -0.00003162328 , -0.00002511918 , -0.00001995282 , -0.00001584906 , -0.00001258933 , -0.00001000005 , -0.00000794331 , -0.00000630959 , -0.00000501188 , -0.00000398108 , -0.00000316228 , -0.00000251189 , -0.00000199526 , -0.00000158489 , -0.00000125893 , -0.00000100000 , -0.00000079433 , -0.00000063096 , -0.00000050119 , -0.00000039811 , -0.00000031623 , -0.00000025119 , -0.00000019953 , -0.00000015849 , -0.00000012589 , -0.00000010000 , -0.00000007943 , -0.00000006310 , -0.00000005012 , -0.00000003981 , -0.00000003162 , -0.00000002512 , -0.00000001995 , -0.00000001585 , -0.00000001259 , -0.00000001000 , 0., 0.,0.,0.,0.,0.,0.,0.};

int reported_version_error = 0;
gene_quality_score_t get_subread_quality(const char * quality_str, const char * read_str, int quality_scale, int phred_version)
{
	int ret =0;
	int i;
/*
	for (i=0; i<16; i++)
		if(read_str[i] == 'N' || read_str[i] == '.')
			return -22.;
*/

	//for(i=0;i<16;i++) ret += get_base_quality_score(quality_str[i] , quality_scale);
	if(FASTQ_PHRED64 == phred_version)
		for(i=0;i<16;i++) ret += (1000000 - get_base_error_prob64i(quality_str[i])); 
	else
		for(i=0;i<16;i++) ret += (1000000 - get_base_error_prob33i(quality_str[i])); 

	/*
	if (ret <0 && !reported_version_error)
	{
		SUBREADprintf("\nWARNING: negative Phred quality score! Please verify the version of the Phred scores.\n");
		reported_version_error=1;
	}*/
	
	return ret/16000000.;
}

/*int get_base_phred(char quality_chr)
{
	return quality_chr - '@';
}

void print_running_log(double finished_rate, double read_per_second, double expected_seconds, unsigned long long int total_reads, int is_pair)
{
	char outbuff[99]; int i;
	snprintf(outbuff, 98,"completed=%0.2f%%; time used=%.1fs; rate=%.1fk reads/s; time left=%.1fs; total=%lluk %s", finished_rate*100, miltime()-begin_ftime, read_per_second/1000 ,expected_seconds, total_reads/1000, is_pair?"pairs":"reads");
	fputs(outbuff, stdout);
	for(i=strlen(outbuff); i<105; i++)
		SUBREADprintf(" ");
	SUBREADprintf("%c",CORE_SOFT_BR_CHAR);
}

*/

void print_window_scrolling_bar(char * hint, float percentage, int width, int * internal_counter)
{
	/*
	char fan = '-';
	int bar_width = width - 7 - strlen(hint) , i;
	int dash_width = (int)(bar_width * percentage+0.5);
	dash_width = min(dash_width, bar_width - 1);
	int space_width = bar_width - 1 - dash_width;
	
	char out_buf[100];

	switch ((*internal_counter) % 4)
	{
		case 0:
			fan='-';
			break;
		case 1:
			fan='\\';
			break;
		case 2:
			fan='|';
			break;
		case 3:
			fan='/';
			break;
	}

	(*internal_counter) ++;

	sprintf (out_buf," %c %s [", fan, hint);
	for(i = 0; i< dash_width; i++)
		strcat(out_buf,"=");
	strcat(out_buf,">");
	for(i = 0; i< space_width; i++)
		strcat(out_buf," ");
	strcat(out_buf,"]");
	
	print_in_box(80,0,0,"%s%c",out_buf,CORE_SOFT_BR_CHAR);

	*/
	print_in_box(81,0,0," [ %.1f%%%% finished ]    %c",percentage*100, CORE_SOFT_BR_CHAR);
}




void print_text_scrolling_bar(char * hint, float percentage, int width, int * internal_counter)
{
	char fan = '-';
	int bar_width = width - 7 - strlen(hint) , i;
	int dash_width = (int)(bar_width * percentage+0.5);
	dash_width = min(dash_width, bar_width - 1);
	int space_width = bar_width - 1 - dash_width;

	/*
	for (i=0; i<width; i++);
		putchar(' ');
	SUBREADprintf("\r");
	*/
 
	switch ((*internal_counter) % 4)
	{
		case 0:
			fan='-';
			break;
		case 1:
			fan='\\';
			break;
		case 2:
			fan='|';
			break;
		case 3:
			fan='/';
			break;
	}

	(*internal_counter) ++;

	
	char lbuf[100];
	sprintf (lbuf," %c %s [", fan, hint);
	for(i = 0; i< dash_width; i++)
		strcat(lbuf, "=");
	strcat(lbuf, ">");
	for(i = 0; i< space_width; i++)
		strcat(lbuf, " ");
	strcat(lbuf, "]");
	SUBREADprintf("%s%c",lbuf,CORE_SOFT_BR_CHAR);
	
}


void add_repeated_numbers(int qid, gene_vote_t * vote, unsigned char * repeated_regions)
{

	int i,j;

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			int v_2 ;
			if(vote-> votes[i][j] >=2)
			{
				v_2 = vote-> votes[i][j] - 2;
				if(repeated_regions[v_2 + qid  * 16] < 255) repeated_regions[v_2 + qid  * 16]++;
			}

		}
}

int remove_repeated_reads(gehash_t * table, gehash_t * huge_table, int index_threshold)
{
	int i;
	int vals[200000];
	int val_len[200000];
	int scroll_count = 0;
	unsigned int all_removed = 0;


	for (i=0; i <table->buckets_number ; i++)
	{
		struct gehash_bucket * cb = table->buckets + i;
		int j;
		int val_c = 0;

		if(i % 300 == 0)
			print_text_scrolling_bar ("Removing non-informative subreads", 1.0*i/table->buckets_number, 80, &scroll_count);

		for(j=0; j<cb->current_items; j++)
		{
			int k, found = 0;
			for(k=0;k<val_c;k++)
			{
				if (vals[k]==cb->item_keys[j])
				{
					val_len[k]++;
					found = 1;
					break;
				}
			}

			if(!found)
			{
				if(val_c>199999)
				{
					SUBREADprintf("\nThere are too many items in a bucket; you may reduce the threshold of non-informative subreads to eliminate this problem.\n");
					#ifdef MAKE_STANDALONE
					exit(-1);
					#else
					return -1;
					#endif
				}
				else
				{
					vals[val_c] = cb->item_keys[j];
					val_len[val_c] = 1;
					val_c++;
				}
			}
		}

//			SUBREADprintf ("F3\n");

		for(j=0; j<val_c; j++)
		{

			if (gehash_exist(huge_table, vals[j]))
				 gehash_remove(table, vals[j]);
			else if(val_len[j]> index_threshold)
			{
				gehash_remove(table, vals[j]);
				gehash_insert(huge_table, vals[j], 1, NULL);
				all_removed += val_len[j];
			}
		}
	}

	if(IS_DEBUG)
		SUBREADprintf ("@LOG There are %u subreads removed from the index.\n", all_removed);

	return all_removed;
}


unsigned int linear_gene_position(const gene_offset_t* offsets , char *chro_name, unsigned int chro_pos)
{
	unsigned int ret = 0 ;
	int n;


	n = HashTableGet(offsets -> read_name_to_index, chro_name)-NULL;
	//printf("09NOV: GET '%s' = %d\n", chro_name, n);
	if(n<1) return 0xffffffff;
	if(n>1)
		ret = offsets->read_offsets[n-2];
	ret += offsets -> padding ;

	return ret + chro_pos;
}

unsigned int get_gene_linear(int chrono, int offset, const unsigned int offsets [])
{
	if (chrono>1)return offsets[chrono-1]+offset;
	return offset;
}

int locate_gene_position_max(unsigned int linear, const gene_offset_t* offsets , char ** chro_name, int * pos, int * head_cut_length, int * tail_cut_length, int rl)
{
	int n = 0;

	(*chro_name) = NULL;
	(*pos) = -1;
	int low_idx=0, high_idx= offsets->total_offsets;
	while(1){
		if(high_idx <= low_idx+1) {
			n = max(low_idx-2,0);
			break;
		}
		int mid_idx = (low_idx+high_idx)/2;
		unsigned int mid_val = offsets->read_offsets[mid_idx];
		if(mid_val > linear) high_idx = mid_idx;
		else low_idx = mid_idx+1;
	}

	for (; n<offsets->total_offsets; n++)
	{
		if (offsets->read_offsets[n] > linear)
		{

			//#warning "======== COMMENT THIS LINE !! ========"
			//SUBREADprintf("max=%u <= lim=%u : ACCEPTED.\n", rl + linear , offsets->read_offsets[n] + 16);

			// the end of the read should not excess the end of the chromosome

			if (n==0)
				*pos = linear;
			else
				*pos = linear - offsets->read_offsets[n-1];


			if(tail_cut_length == NULL){ 
				if(rl + linear > offsets->read_offsets[n] + 15 - offsets -> padding){
					return 1;
				}
			} else {
				unsigned int posn1 = 0;
				if(n>0) posn1 = offsets->read_offsets[n-1];
				long long int tct =  ( linear + rl - posn1 - offsets -> padding );
				if(tct < rl)tct = rl;
				long long int chro_leng = (offsets->read_offsets[n] - posn1 - 2*offsets -> padding + 16);

				//SUBREADprintf("CHRO_LEN : %lld, READ_TAIL %lld , RL=%d\n", chro_leng, tct, rl);
				tct -= chro_leng;
				if( tct >= rl ){
					return 1;
				}
				if( tct <0    )tct=0;
				(*tail_cut_length) = tct;
			}

			if( (*pos) < offsets -> padding ) {
				if(head_cut_length == NULL || (*pos) + rl <= offsets -> padding){
					return 1;
				}else{
					(*head_cut_length) = offsets -> padding - (*pos);
					(*pos) = offsets -> padding;
				}
			}

			(*pos) -= offsets -> padding;
			*chro_name = (char *)offsets->read_names+n*MAX_CHROMOSOME_NAME_LEN;

			return 0;
		}
	}
	return -1;
}


int locate_gene_position(unsigned int linear, const gene_offset_t* offsets , char ** chro_name, int * pos)
{
	return locate_gene_position_max(linear, offsets, chro_name, pos, NULL, NULL, 0);
}

#define _index_vote(key) (((unsigned int)key)%GENE_VOTE_TABLE_SIZE)

#ifdef USE_SLOW_HASHTABLE_INDEX
int vv=0;
void add_gene_vote(gene_vote_t* vote, int key , int add_new) {
	add_gene_vote_weighted(vote, key, add_new, 1);
}

void add_gene_vote_weighted(gene_vote_t* vote, int key , int add_new, int w) {
	int offset = _index_vote(key);
	int datalen = vote -> items[offset];
	unsigned int * dat = vote -> pos[offset];
	int i;

	for (i=0;i<datalen;i++)
	{
		if (dat[i] == key)
		{

			int test_max = (vote->votes[offset][i]);
			test_max += w;
			vote->votes[offset][i] = test_max;

			if(test_max > vote->max_vote){
				vote->max_vote = test_max;
				vote->max_position = key;
			}
			return;
		}
	}

	if(!add_new || datalen >=GENE_VOTE_SPACE)
		return;

//	w += 16;

	vote -> items[offset] ++;
	dat[i] = key;
	vote->votes[offset][i]=w;
}

#endif

int max_gene_vote(gene_vote_t* vote, int * position_result, int query)
{
	int n, i, max_index=0;
	int max_votes = -1;
	for (n=0; n<GENE_VOTE_TABLE_SIZE; n++)
	{
		int itms = vote->items[n];
		gene_vote_number_t * vots = vote->votes[n];
		for (i=0; i<itms; i++)
		{

			if (vots[i] > max_votes)
			{
				max_index = n<<16|i;
				max_votes = vots[i];
			}
		}
	}

	if(max_votes == -1)
	{
		*position_result = -1;
		return 0;
	}
	else
	{
		*position_result = vote->pos[max_index>>16][max_index&0xffff];
		return max_votes;
	}
}


void clear_allvote(gene_allvote_t* allvote)
{
	memset(allvote -> max_votes,0,  allvote -> max_len* sizeof(*allvote -> max_votes));
	memset(allvote -> masks,0,  allvote -> max_len* sizeof(*allvote -> masks));
	memset(allvote -> span_coverage ,0,  allvote -> max_len* sizeof(char));
}

void destory_allvote(gene_allvote_t* allvote)
{
	free(allvote -> max_positions);
	free(allvote -> max_votes);
	free(allvote -> max_quality);
	free(allvote -> max_final_quality);
	free(allvote -> masks);
	#ifdef REPORT_ALL_THE_BEST
	free(allvote -> best_records);
	#endif
	free(allvote -> is_counterpart);
	free(allvote -> span_coverage);
	if(allvote -> max_indel_recorder)
		free(allvote -> max_indel_recorder);
}

int init_allvote(gene_allvote_t* allvote, int expected_len, int allowed_indels)
{
	int is_OK = 0;
	allvote -> max_len = expected_len; 
	allvote -> max_positions = (unsigned int *) malloc(sizeof(int)*expected_len);
	allvote -> max_votes = (gene_vote_number_t *) calloc(sizeof(gene_vote_number_t), expected_len);
	allvote -> max_quality = (gene_quality_score_t *) calloc(sizeof(gene_quality_score_t), expected_len);
	allvote -> max_final_quality = (gene_quality_score_t *) calloc(sizeof(gene_quality_score_t), expected_len);
	allvote -> masks = (short *) calloc(sizeof(short), expected_len);
#ifdef REPORT_ALL_THE_BEST
	allvote -> best_records = (gene_best_record_t *) malloc(sizeof(gene_best_record_t)* expected_len);
#endif
	allvote -> is_counterpart = (unsigned char *) malloc(expected_len);

	allvote -> max_indel_tolerance = allowed_indels;
	allvote	-> indel_recorder_length = max(3*(allowed_indels+1)+1, 28);
	allvote -> span_coverage = (char *)calloc(sizeof(char), expected_len);
	allvote -> repeated_regions = (unsigned char *)calloc(sizeof(unsigned char), 16 * expected_len);

	if((allvote -> max_quality &&  allvote -> max_positions  && allvote -> max_votes  && allvote -> max_final_quality && allvote -> masks && allvote -> is_counterpart && allvote -> span_coverage))
		is_OK = 1;

	if(allowed_indels && is_OK)
	{
		allvote -> max_indel_recorder = (char *)malloc( allvote -> indel_recorder_length *expected_len);
		is_OK = (allvote -> max_indel_recorder!=NULL);
	}
	else	allvote -> max_indel_recorder =  NULL;

	if(is_OK)
		return 0;
	else
	{
		SUBREADputs(MESSAGE_OUT_OF_MEMORY);
		return 1;
	}
}


// rl_adjust >0: "D" in read; rl_adjust<0: "I" in read
void compress_cigar(char *cigar, int total_length, char * read, int * pos_offset, int *rl_adjust)
{
	char tmp[200];
	char cigar_piece [10];

	int cigar_len = strlen(cigar);
	cigar_len = min(cigar_len, EXON_MAX_CIGAR_LEN);
	cigar[cigar_len]=0;

	int i;
	srInt_64 tmpv = 0;

	char last_operation = 'X';
	srInt_64 last_tmpv = 0;
	int is_first_M = 1;
	int delta_i = 0;
	int delta_d = 0;
	int delta_rl = 0;
	tmp[0]=0;
	int cigar_length=0;
	
	for(i=0; i < cigar_len; i++)
	{
		char cc = cigar[i];
		if(isdigit(cc))
		{
			tmpv=tmpv*10+(cc-'0');
		}
		else if(cc=='-')
		{
			last_tmpv = 0;
			cigar_length = 0;
			break;
		}
		else
		{
			if(is_first_M && pos_offset && cc=='S')
				*pos_offset = tmpv;
			else if(cc=='M') is_first_M = 0;
				
			if((cc!=last_operation) && (last_operation!='X'))
			{
				if(last_operation == 'M' || last_operation == 'S' || last_operation == 'N' ||  last_operation == 'B' || last_operation == 'J' || last_operation == 'j' || last_operation == 'b')
				{
					if(delta_i)
					{
						sprintf(cigar_piece,"%dI", delta_i); 
						strcat(tmp, cigar_piece);
					}
					delta_i = 0;

					if(delta_d)
					{
						sprintf(cigar_piece,"%dD", delta_d); 
						strcat(tmp, cigar_piece);
					}
					delta_d = 0;
	
					#ifdef __MINGW32__
					sprintf(cigar_piece,"%I64d%c", last_tmpv, last_operation); 
					#else
					sprintf(cigar_piece,"%lld%c", last_tmpv, last_operation); 
					#endif
					strcat(tmp, cigar_piece);
				}

				if(last_operation == 'M' || last_operation == 'S' || last_operation == 'I')
					cigar_length += last_tmpv;

				last_tmpv = 0;
			}
			if(cc=='I' )
			{
				delta_i += tmpv;
				delta_rl -= tmpv;
			}
			if(cc=='D' )
			{
				delta_d += tmpv;
				delta_rl += tmpv;
			}
	

			last_tmpv += tmpv;
		
			tmpv = 0;
			last_operation = cc;
		
		}
	
	}

	if(last_tmpv)
	{
		if(delta_i)
		{
			sprintf(cigar_piece,"%dI", delta_i); 
			strcat(tmp, cigar_piece);
		}
		if(delta_d)
		{
			sprintf(cigar_piece,"%dD", delta_d); 
			strcat(tmp, cigar_piece);
		}



		if(last_operation =='M' || last_operation =='S')
		{
			#ifdef __MINGW32__
			sprintf(cigar_piece,"%I64d%c", tmpv+last_tmpv, last_operation); 
			#else
			sprintf(cigar_piece,"%lld%c", tmpv+last_tmpv, last_operation); 
			#endif
			strcat(tmp, cigar_piece);
		}

		if(last_operation == 'M' || last_operation == 'S' || last_operation == 'I')
			cigar_length += tmpv+last_tmpv;
	}
	if(cigar_length == total_length)
	{
		if(rl_adjust)(*rl_adjust)=delta_rl;
		strcpy(cigar, tmp);
	}
	else
	{
		sprintf(cigar, "%dM", total_length);
	}
}

// adjust_len <0: "I" in read; adjust_len>0: "D" in read
void show_cigar(char * info, int len, int is_reversed_map, char * buf, int indel_tolerance, int total_subreads, char *read, int * pos_offset, int * adjust_len)
{
	int i;
	int last_offset = 0, cursor = 0;

	if(info[0]==-1){
		sprintf(buf, "%dM", len);
		return;
	}
	else if(info[0]==-2){
		if (strchr( info+1, '-'))
			sprintf(buf, "%dM", len);
		else
		{
			strncpy(buf, info+1, 98);
			compress_cigar(buf, len, read, pos_offset, adjust_len);
		}
		return;
	}
	else if(info[0]==-3){
		info ++;
	}

	for(i=0; i<indel_tolerance*3; i+=3)
	{
		if (!info[i])break;
		int dist = info[i+2];
		//int subread_start = info[i]-1;
		int subread_end = info[i+1]-1;

		int base_end = (i < indel_tolerance*3-3 && info[i+3])?find_subread_end(len, total_subreads, subread_end):len;

//			SUBREADprintf("BE=%d ; II+3=%d ; len=%d\n", base_end, info[i+3], len);

		int offset = last_offset-dist;
		if (base_end - cursor - (offset>0?offset:0) <0) 
		{
			buf[0]=0;

			cursor = 0;
			break;
		}
		if (i >0)
		{
			sprintf(buf+strlen(buf), "%d%c%dM", abs(offset), offset>0?'I':'D', base_end - cursor - (offset>0?offset:0));
			if(adjust_len)
				(*adjust_len) -= offset;
		}
		else
			sprintf(buf+strlen(buf), "%dM", base_end);
		last_offset = dist;
		cursor = base_end;
	}
	compress_cigar(buf, len, read, pos_offset, adjust_len);
}

void add_allvote_q(gene_allvote_t* allvote,int qid , int pos, gene_vote_number_t votes, gene_quality_score_t quality, int is_counterpart, short mask, char * max_indel_recorder, gene_value_index_t * array_index, char * read_txt, int read_len, int max_indel, int total_subreads, int space_type, int report_junction, int is_head_high_quality, char * qual_txt, int phred_version, char span_coverage, short **dynamic_programming_short ,  char ** dynamic_programming_char)
{

//	SUBREADprintf("Ospan=%d Nspan=%d\n", allvote -> span_coverage[qid], span_coverage);

}


float CORE_RECOVER_MATCHING_RATE = 0.9;
float CORE_INDEL_MATCHING_RATE_TAIL = 0.8;
float CORE_INDEL_MATCHING_RATE_HEAD = 0.8;

void find_and_explain_indel(gene_allvote_t* allvote,int qid , int pos, gene_vote_number_t votes, gene_quality_score_t quality, int is_counterpart, char mask, char * max_indel_recorder, gene_value_index_t * array_index, char * read_txt, int read_len, int max_indel, int total_subreads, int space_type, int report_junction,  int is_head_high_quality, char * qual_txt, int phred_version , short ** dynamic_programming_short,  char **dynamic_programming_char)
{
	//unsigned int pos0 = pos;
	if(allvote -> max_indel_recorder)
	{
		//PNT111
		short head_indel_pos=-1 , tail_indel_pos=-1;
		int head_indel_movement=0, tail_indel_movement=0;
		if (array_index)
		{
			int k;
			if (max_indel_recorder[0])
			{
				int cover_start = find_subread_end(read_len, total_subreads, max_indel_recorder[0]-1)- 15;

				for (k=0; max_indel_recorder[k]; k+=3);
				int cover_end = find_subread_end(read_len, total_subreads, max_indel_recorder[k-2]-1) + max(0,-max_indel_recorder[k-1]);

				int head_qual = is_head_high_quality?0x7fffffff:0, tail_qual = is_head_high_quality?0:0x7fffffff;


				if(qual_txt && qual_txt[0])
				{
					int j;
					head_qual = 0;
					tail_qual = 0;
					for (j = 0; j<cover_start; j++)
					{
						if(FASTQ_PHRED64 == phred_version)
						{
							head_qual += (1000000-get_base_error_prob64i(qual_txt[j]));
						}
						else
						{
							head_qual += (1000000-get_base_error_prob33i(qual_txt[j]));
						}
					}
					for (j = 0; j<read_len - cover_end; j++)
					{
						if(FASTQ_PHRED64 == phred_version)
						{
							tail_qual += (1000000-get_base_error_prob64i(qual_txt[read_len - j -1]));
						}
						else
						{
							tail_qual += (1000000-get_base_error_prob33i(qual_txt[read_len - j -1]));
						}
					}
				}

				int head_must_correct = 4;
				if(cover_start > 0)
				{
					if (head_qual / cover_start < 850000)
					{
						CORE_INDEL_MATCHING_RATE_HEAD = 0.75;
						head_must_correct =2;
					}
					else if (head_qual / cover_start < 950000 )
					{
						CORE_INDEL_MATCHING_RATE_HEAD = 0.85;
						head_must_correct =3;
					}
					else CORE_INDEL_MATCHING_RATE_HEAD = 0.92;
				}else	CORE_INDEL_MATCHING_RATE_HEAD = 9999;
	

				int tail_must_correct = 4;
				if(read_len-cover_end > 0)
				{
					if (tail_qual / (read_len-cover_end) < 850000)
					{
						CORE_INDEL_MATCHING_RATE_TAIL = 0.75;
						tail_must_correct =2;
					}
					else if (tail_qual / (read_len-cover_end) < 950000)
					{
						CORE_INDEL_MATCHING_RATE_TAIL = 0.85;
						tail_must_correct =3;
					}
					else CORE_INDEL_MATCHING_RATE_TAIL = 0.92;
				}else	CORE_INDEL_MATCHING_RATE_TAIL = 9999;


				

				int is_full_covered = 0;
				is_full_covered = extend_covered_region(array_index, pos, read_txt, read_len, cover_start, cover_end, 4, head_must_correct, tail_must_correct,  max_indel, space_type, max_indel_recorder[k-1], &head_indel_pos, &head_indel_movement, &tail_indel_pos, &tail_indel_movement, is_head_high_quality, qual_txt, phred_version, CORE_INDEL_MATCHING_RATE_HEAD, CORE_INDEL_MATCHING_RATE_TAIL);
				if(head_indel_movement)
				{
					pos += head_indel_movement;
					allvote -> max_positions[qid] = pos; 
				}
				if(is_full_covered == 3)
				{
					allvote -> masks[qid] &= ~IS_RECOVERED_JUNCTION_READ;
				}
				else	allvote -> masks[qid] |= IS_RECOVERED_JUNCTION_READ;
			}
		}

		*(allvote -> max_indel_recorder + qid * allvote -> indel_recorder_length) = 0xff ;


		if ((max_indel_recorder[3] || head_indel_pos>=0 || tail_indel_pos>0)){
			int head_pos = 0;
			if((head_indel_pos>=0 && report_junction) || head_indel_movement)head_pos = head_indel_pos ;//- max(0,head_indel_movement);
			int tail_pos = read_len;
			if((tail_indel_pos>0 && report_junction) || tail_indel_movement)tail_pos = tail_indel_pos ;//- min(0, tail_indel_movement);
			explain_indel_in_middle(allvote, qid , pos, max_indel_recorder,  array_index, read_txt, read_len,  max_indel, total_subreads, head_pos, tail_pos , head_indel_movement, tail_indel_movement, report_junction, dynamic_programming_short, dynamic_programming_char);
			//if(*(allvote -> max_indel_recorder + qid * allvote -> indel_recorder_length) != 0xfd ) allvote -> max_positions[qid] = pos0;
		}

	}

}

void explain_indel_in_middle(gene_allvote_t* allvote, int qid , int pos, char * max_indel_recorder, gene_value_index_t * array_index, char * read_txt, int read_len, int max_indel, int total_subreads, int head_start_point, int tail_end_point, int head_indel_movement, int tail_indel_movement, int report_junction, short ** dynamic_programming_short ,  char ** dynamic_programming_char)
{

	char indel_operations[1500];
	int i,xx;
	char tmp_cigar [EXON_MAX_CIGAR_LEN+1];

	int current_pos = head_start_point; 
	int explain_cursor = head_start_point;
	char last_operation = 0;
	int del_number = 0;
	int dynamic_delta = 0;

	tmp_cigar[0]=0;

	if(current_pos >0)
	{
		if (head_indel_movement != 0)
			sprintf(tmp_cigar, "%dM%d%c", head_start_point  - max(head_indel_movement, 0), abs(head_indel_movement), head_indel_movement>0?'I':'D');
		else if (report_junction)
			sprintf(tmp_cigar, "%dS", head_start_point);
	}

	//ddprintf ("R=%s; REC[8]=%d [9]=%d\n", read_txt, max_indel_recorder[8], max_indel_recorder[9]);

	for (i=3; max_indel_recorder[i]; i+=3)
	{
		if (i >= max_indel*3) break;

		int last_dist = max_indel_recorder[i-1];
		int black_subread_start = max_indel_recorder[i-2]-1;
		int black_subread_end = max_indel_recorder[i]-1;
		if (black_subread_end < black_subread_start+1) black_subread_end = black_subread_start+1;
		while(max_indel_recorder[i+3])
		{
			if (max_indel_recorder[i+3]- black_subread_end>2 || i >= max_indel*3-3) break;
			black_subread_end = max_indel_recorder[i+3]-1;
			i+=3;
		}
		int black_base_start = find_subread_end(read_len, total_subreads,black_subread_start )- 3 + max(0,-last_dist);
		int black_base_end = find_subread_end(read_len, total_subreads, black_subread_end)-15 + 3 + max(0,-max_indel_recorder[i+2]);

		int blackref_base_start = black_base_start - max(0,-last_dist) + max(0,last_dist);
//		SUBREADprintf("baseref_offset=%d\n", max(0,last_dist));

		int gap_end_read =  black_base_end + max(0,max_indel_recorder[i+2]);

		blackref_base_start = max(0,blackref_base_start);
		black_base_start = max(0,black_base_start);
		black_base_end = min(read_len, black_base_end);

		int vpos = strlen(tmp_cigar);
		if (vpos > EXON_MAX_CIGAR_LEN -6) break;

		int exp_indel = last_dist - max_indel_recorder[i+2];

		
		int moves = 0;
		if(black_base_end > black_base_start)
			moves = dynamic_align (read_txt + black_base_start, black_base_end-black_base_start , array_index, pos + blackref_base_start, max_indel, indel_operations, exp_indel,  -10, black_base_end - black_base_start+5, dynamic_programming_short, dynamic_programming_char);

		last_operation = 0;
		if(moves)// < (black_base_end-black_base_start) + 6 + 1.5 * exp_indel)
		{
			int pre_read_len = 0;
			current_pos = black_base_start ;
			xx=0;
			while (pre_read_len < 0 && xx < moves)
			{
				if(indel_operations[xx]!=2)pre_read_len++;
				xx++;
			}
			for (; xx<moves; xx++)
			{
				int current_operation = indel_operations[xx];
				if(current_operation == 3) current_operation = 0;

				if(current_operation != last_operation && (current_pos < gap_end_read-1+1 || current_operation == 0))// && current_pos!=explain_cursor)
				{
					int vpos = strlen(tmp_cigar);
					if (vpos>EXON_MAX_CIGAR_LEN-6) break;
					sprintf(tmp_cigar + vpos, "%d%c", last_operation==1?del_number:(current_pos - explain_cursor), last_operation==0?'M':(last_operation==1?'D':'I'));
					explain_cursor = current_pos ;
					del_number = 0;
				}
				if(current_pos >= gap_end_read-1 || xx == moves - 1){
					if(current_operation != 0 && current_pos!=explain_cursor)
					{
						int vpos = strlen(tmp_cigar);
						sprintf(tmp_cigar + vpos, "%d%c", current_operation==1?del_number:(current_pos - explain_cursor), current_operation==1?'D':'I');
						explain_cursor = current_pos;
						del_number=0;
					}
					break;
				}

				if (indel_operations[xx]==0) current_pos ++;
				else if (indel_operations[xx]==1) {del_number++; dynamic_delta++;}//"D"
				else if (indel_operations[xx]==2) {current_pos ++; dynamic_delta--;} //"I"
				else if (indel_operations[xx]==3) current_pos ++;


				last_operation = current_operation; 
			} 
		}
		else
		{
			int movement = last_dist - max_indel_recorder[i+2];
			int vpos = strlen(tmp_cigar);
			if (vpos>EXON_MAX_CIGAR_LEN-6) break;

			current_pos = find_subread_end(read_len, total_subreads,max_indel_recorder[i]-1)-15 + (max_indel_recorder[i+2] >0?max_indel_recorder[i+2]:-max_indel_recorder[i+2]);
			current_pos -= (black_base_end - black_base_start)/2 -3;

			if(movement)
				sprintf(tmp_cigar + vpos, "%dM%d%c", current_pos - explain_cursor , abs(movement), movement<0?'D':'I' );
			else if (current_pos - explain_cursor>0)
				sprintf(tmp_cigar + vpos, "%dM", current_pos - explain_cursor );

			explain_cursor = current_pos + (movement >0?movement:0);
			current_pos = find_subread_end(read_len, total_subreads,max_indel_recorder[i+1]-1);
			last_operation = 0;
			del_number = 0;
			moves = 0;
		}
		last_operation = 0;
		current_pos = black_base_end;
			
		
	}

	int vpos = strlen(tmp_cigar);
	if (vpos > EXON_MAX_CIGAR_LEN-9 )
	{
		*(allvote -> max_indel_recorder + qid * allvote -> indel_recorder_length) = 0xfd; 
		memcpy(allvote -> max_indel_recorder + qid * allvote -> indel_recorder_length+1, max_indel_recorder, 3 * max_indel);
	}
	else
	{
		if (explain_cursor<read_len)	
		{
			if(tail_end_point > explain_cursor)
				sprintf(tmp_cigar+vpos,"%dM", tail_end_point - explain_cursor);
			vpos = strlen(tmp_cigar);
			if(tail_end_point< read_len)
			{
				if (tail_indel_movement!=0) 
				{
					int tail_len_m = read_len -(tail_end_point - min(tail_indel_movement, 0));
					if(tail_len_m)
						sprintf(tmp_cigar+vpos, "%d%c%dM", abs(tail_indel_movement), tail_indel_movement<0?'I':'D', tail_len_m);
					else
						sprintf(tmp_cigar+vpos, "%d%c", abs(tail_indel_movement), tail_indel_movement<0?'I':'D');
				}
				else if (report_junction)
					sprintf(tmp_cigar+vpos, "%dS",read_len - tail_end_point );
			}
		}
	//	if (head_start_point >0 || tail_end_point < read_len)
	//		strcat(tmp_cigar,"X");

		*(allvote -> max_indel_recorder + qid * allvote -> indel_recorder_length) = 0xfe ;
		strncpy(allvote -> max_indel_recorder + qid * allvote -> indel_recorder_length+1, tmp_cigar, allvote -> indel_recorder_length - 2);
	}
}

int evaluate_piece(char * piece_str, int chron, int offset, int is_counterpart, int start_pos, int end_pos)
{
	char fname[MAX_FILE_NAME_LENGTH];
	int inner_pos = 0, i;
	FILE * fp;
	char next_char=0;
	int ret = 0;

	if (chron == 0)
		sprintf(fname, "/opt/Work2001/Gene-Search/src/GENE-LIB/%02da.fa", chron);
	else
		sprintf(fname, "/opt/Work2001/Gene-Search/src/GENE-LIB/%02d.fa", chron);

	inner_pos = offset + offset / 70;

	fp = f_subr_open(fname,"r");

	while(next_char!='\n')
		next_char=fgetc(fp);
	fseek(fp, inner_pos, SEEK_CUR);

	for(i=0 ; i<end_pos; i++)
	{
		next_char = get_next_char(fp);
		if(i < start_pos)
			continue;
		if(next_char == 'N')
			ret ++;
		else{
			if(is_counterpart)
			{
				if (piece_str[99-i] == 'A' && next_char == 'T')
					ret ++;
				else if (piece_str[99-i] == 'G' && next_char == 'C')
					ret ++;
				else if (piece_str[99-i] == 'T' && next_char == 'A')
					ret ++;
				else if (piece_str[99-i] == 'C' && next_char == 'G')
					ret ++;
			}
			else if(piece_str[i] == next_char)
				ret ++;
		}
	}

	fclose(fp);

	return ret;

}

int match_chro_indel(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int is_negative_strand, int space_type, int indel_size, gene_vote_number_t * indel_recorder, int total_subreads)
{

	int ret = 0;
	int tali = 0;
	int last_section_end = 0;

	for(tali = 0; indel_recorder[tali*3] && tali<MAX_INDEL_TOLERANCE; tali++)
	{
		int section_last_subread = indel_recorder[tali*3+1] - 1;
		int section_offset = indel_recorder[tali*3+2];
		int this_section_end = find_subread_end(test_len , total_subreads, section_last_subread); 
		if(!indel_recorder[tali*3+3]) this_section_end = test_len;

		this_section_end = min(test_len, this_section_end);
		this_section_end = max(this_section_end, last_section_end);

		ret += match_chro(read + last_section_end - min(0, section_offset), index, pos + last_section_end + max(0, section_offset), this_section_end - last_section_end + min(0, section_offset), is_negative_strand * 0, space_type);

		last_section_end = this_section_end;
	}

	return ret;

}


int match_chro_indel_old(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int is_negative_strand, int space_type, int indel_size)
{

	int i;
	int ret = 0;
	for (i = -indel_size; i<=indel_size ; i++)
	{
		if (pos +i + test_len < index -> start_base_offset + index -> length && pos > -i)
			ret += match_chro(read, index, pos + i, test_len, is_negative_strand, space_type);
	}
	return ret;

}



void mark_votes_array_index(char * read_str, int read_len, gene_vote_t * dest, gene_vote_t * src, gene_value_index_t * my_array_index, int color_space, int indel_tolerance, int min_minor, const char quality_str [], int quality_scale, int is_negative_strand)
{
	int i, j;

	dest -> max_vote = 0;
	dest -> max_quality = 0;

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
	{
		dest -> items[i] = src -> items[i];
		for(j=0; j< src->items[i]; j++)
		{
			unsigned int potential_position = src-> pos[i][j];
			float matchingness_count = 0;
			if (src -> votes [i][j]>= min_minor)
				//matchingness_count = match_read(read_str, read_len, potential_position, my_array_index, color_space, indel_tolerance, quality_str, quality_scale);
				matchingness_count = 1.*match_chro_indel(read_str, my_array_index, potential_position, read_len, 0, color_space, indel_tolerance, NULL, 0);
	
			dest -> pos[i][j] = potential_position;
			dest -> quality[i][j] = matchingness_count;
			dest -> votes [i][j] =src -> votes [i][j];
			dest -> masks[i][j] = src -> masks[i][j];
			dest -> coverage_start[i][j] = src -> coverage_start[i][j];
			dest -> coverage_end[i][j] = src -> coverage_end[i][j];
			memcpy(dest -> indel_recorder[i][j], src -> indel_recorder[i][j], MAX_INDEL_TOLERANCE*3*sizeof(char));
			// find for the `best positions'; only the best positions are replaced by the base machingness scores.
			if((matchingness_count >  dest -> max_quality && src -> votes [i][j] == dest -> max_vote) || (src -> votes [i][j] > dest -> max_vote))
			{
				memcpy(dest -> max_indel_recorder, src -> indel_recorder[i][j], MAX_INDEL_TOLERANCE*3*sizeof(char));
				dest -> max_vote  = src -> votes [i][j];
				dest -> max_mask  = src -> masks[i][j];
				dest -> max_quality = matchingness_count;
				dest -> max_position = potential_position;
				dest -> max_coverage_start = src -> coverage_start[i][j];
				dest -> max_coverage_end = src -> coverage_end[i][j];
			}
		}
	}
}

int select_positions_array(char * read1_str, int read1_len, char * read2_str, int read2_len,  gene_vote_t * vote_read1, gene_vote_t * vote_read2, gene_vote_number_t * numvote_read1, gene_vote_number_t * numvote_read2,  gene_quality_score_t * sum_quality, gene_quality_score_t * qual_r1, gene_quality_score_t * qual_r2, gene_vote_number_t * r1_recorder, gene_vote_number_t * r2_recorder, gehash_data_t * pos_read1, gehash_data_t * pos_read2, unsigned int max_pair_dest, unsigned int min_pair_dest, int min_major, int min_minor,int is_negative_strand, gene_value_index_t * my_array_index, int color_space, int indel_tolerance, int number_of_anchors_quality, const char quality_str1 [], const char quality_str2 [], int quality_scale, int max_indel_len, int * is_break_even)
{
	gene_vote_t base_vote_1 , base_vote_2;

	mark_votes_array_index(read1_str, read1_len, &base_vote_1, vote_read1, my_array_index, color_space, indel_tolerance, min_minor, quality_str1, quality_scale, is_negative_strand);
	mark_votes_array_index(read2_str, read2_len, &base_vote_2, vote_read2, my_array_index, color_space, indel_tolerance, min_minor, quality_str2, quality_scale, is_negative_strand);

//	SUBREADprintf("Q=%.3f\n", *qual_r2);

	return select_positions(&base_vote_1, &base_vote_2, numvote_read1, numvote_read2, sum_quality, qual_r1, qual_r2,  pos_read1, pos_read2, r1_recorder , r2_recorder ,max_pair_dest,  min_pair_dest,  min_major,  min_minor, is_negative_strand, number_of_anchors_quality,max_indel_len, is_break_even);
}


int select_positions(gene_vote_t * vote_read1, gene_vote_t * vote_read2, gene_vote_number_t * numvote_read1, gene_vote_number_t * numvote_read2, gene_quality_score_t * sum_quality, gene_quality_score_t * qual_r1, gene_quality_score_t * qual_r2 , gehash_data_t * pos_read1, gehash_data_t * pos_read2, gene_vote_number_t * read1_indel_recorder, gene_vote_number_t * read2_indel_recorder, unsigned int max_pair_dest, unsigned int min_pair_dest, int min_major, int min_minor,int is_negative_strand, int number_of_anchors_quality, int max_indel_len, int * is_breakeven)
{

	int k, i, j, anchors = 0;
	gehash_data_t anchors_position [ANCHORS_NUMBER];
	unsigned char anchor_read [ANCHORS_NUMBER];

	gene_vote_number_t anchor_votes = max(vote_read1->max_vote, vote_read2->max_vote);
	gene_vote_number_t anchors_votes [ANCHORS_NUMBER];
	gene_quality_score_t anchors_quality [ANCHORS_NUMBER];
	int anchors_buckets [ANCHORS_NUMBER];
	int anchors_index [ANCHORS_NUMBER];

	gene_vote_number_t minor_votes [ANCHORS_NUMBER];
	gene_quality_score_t minor_quality [ANCHORS_NUMBER];
	gehash_data_t minor_position [ANCHORS_NUMBER];
	int minor_buckets [ANCHORS_NUMBER];
	int minor_index [ANCHORS_NUMBER];
	char is_minor_breakeven [ANCHORS_NUMBER];

	if(anchor_votes < min_major)
		return 0;

	for (k=0; k<2; k++)
	{
		gene_vote_t * current_vote = k?vote_read2:vote_read1;
		if (current_vote->max_vote < anchor_votes)
			continue;

		for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
			for(j=0; j< current_vote->items[i]; j++)
			{
				if(current_vote->votes[i][j] >= anchor_votes -	    0)
				{
					anchors_position[anchors] = current_vote->pos[i][j];
					anchors_votes   [anchors] = current_vote->votes[i][j];
					anchors_quality [anchors] = current_vote->quality[i][j];
					anchor_read	[anchors] = k;
					anchors_buckets [anchors] = i;
					anchors_index   [anchors] = j;

					anchors++;

					if(anchors >= ANCHORS_NUMBER)
					{
	//					SUBREADprintf("Abnormally too many anchors [%d, %d] %f\n", i, j, anchor_votes);
						anchors --;
					}
				}
			}

	}

	memset(minor_votes,0, ANCHORS_NUMBER*sizeof(gene_vote_number_t));
	memset(minor_quality,0, ANCHORS_NUMBER*sizeof(gene_quality_score_t));
	memset(is_minor_breakeven,0, ANCHORS_NUMBER);
	for (k=0; k<2; k++)
	{
		gene_vote_t * current_vote = k?vote_read2:vote_read1;
		gene_vote_t * current_vote2 = k?vote_read1:vote_read2;

		if (current_vote2->max_vote < anchor_votes)
			continue;

		for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
			for(j=0; j< current_vote->items[i]; j++)
			{
				if(current_vote->votes[i][j] >= min_minor)
				{
					int l;
					for (l=0; l<anchors; l++)
					{

						if(anchor_read[l]==k)continue;

						long long int abdist = current_vote->pos[i][j];
						abdist -= anchors_position[l];

						if((anchor_read[l] && !is_negative_strand) ||
						    ((!anchor_read[l]) && is_negative_strand))
							abdist = -abdist;

						if (abdist < 0)
						//#warning WGSIM generates paired-end data in a wrong way; abdist is given its absolute value in such a case.
							abdist = -abdist;
						//	continue;

						long long int abdist_old = minor_position[l];
						abdist_old -=  anchors_position[l];
						if(abdist_old <0) abdist_old=-abdist_old;
						if (	(minor_votes[l] < current_vote->votes[i][j] ||
							(minor_votes[l] == current_vote->votes[i][j] && current_vote->quality [i][j] > minor_quality[l]) ||
							(minor_votes[l] == current_vote->votes[i][j] && abs(current_vote->quality [i][j] - minor_quality[l])<0.0001 && abdist < abdist_old))&& 
							(abdist <= max_pair_dest) &&
							(abdist >= min_pair_dest) 
						   )
						{
							minor_votes[l] = current_vote->votes[i][j];
							minor_position[l] = current_vote->pos[i][j];
							minor_quality[l] = current_vote->quality [i][j];
							minor_buckets[l] = i;
							minor_index[l] = j;
						}
					}
				}
			}

	}

	gene_vote_number_t selection_minor_vote = 0;
	gene_vote_number_t selection_major_vote = 0;
	gehash_data_t selection_minor_pos = 0;
	gehash_data_t selection_major_pos = 0;

	int selected_major_index=0, selected_major_bucket=0, selected_minor_bucket=0, selected_minor_index=0;

	*sum_quality = 0;

	unsigned char selection_major_read_no = 0;
	//int reason = 0;

	for (k=0; k<anchors; k++)
	{
		if(minor_votes[k]==0)
			continue;

		int need_replace = 0;


		if ((minor_votes[k]+anchors_votes[k]) > (selection_minor_vote+selection_major_vote) ||
	   	       ((minor_votes[k]+anchors_votes[k]) == (selection_minor_vote+selection_major_vote) &&
		 	 minor_quality[k] + anchors_quality[k] > *sum_quality))
		{
			need_replace = 1;
		}
		else if (minor_votes[k]+anchors_votes[k] == selection_minor_vote+selection_major_vote &&
		 	 minor_quality[k] + anchors_quality[k] == *sum_quality)
		{
			if (selection_minor_pos != minor_position[k] && minor_position[k] != selection_major_pos)
			{
				* is_breakeven = 1;
				if (selection_major_pos >  anchors_position[k]) need_replace = 1;
			}
		}

		if(need_replace){
			selection_minor_pos = minor_position[k];
			selection_major_pos = anchors_position[k];
			selection_minor_vote = minor_votes[k];
			selection_major_vote = anchors_votes[k];
			selection_major_read_no = anchor_read[k];
			selected_major_bucket = anchors_buckets[k];
			selected_minor_bucket = minor_buckets[k];
			selected_major_index = anchors_index[k];
			selected_minor_index = minor_index[k];
			* is_breakeven = is_minor_breakeven[k];
			* qual_r1 = anchor_read[k] ? minor_quality[k]:anchors_quality[k];
			* qual_r2 = anchor_read[k] ? anchors_quality[k]: minor_quality[k];
			* sum_quality = minor_quality[k] + anchors_quality[k];
	
		}
	}

	if(selection_minor_vote > 0)
	{
		if(selection_major_read_no)	// major is on read 2
		{
			*numvote_read1 = selection_minor_vote;
			*numvote_read2 = selection_major_vote;
			*pos_read1 = selection_minor_pos;
			*pos_read2 = selection_major_pos ;
			indel_recorder_copy(read1_indel_recorder, vote_read1 -> indel_recorder[selected_minor_bucket][selected_minor_index]);
			indel_recorder_copy(read2_indel_recorder, vote_read2 -> indel_recorder[selected_major_bucket][selected_major_index]);
		}
		else
		{
			*numvote_read1 = selection_major_vote;
			*numvote_read2 = selection_minor_vote;
			*pos_read1 = selection_major_pos;
			*pos_read2 = selection_minor_pos;
			indel_recorder_copy(read1_indel_recorder, vote_read1 -> indel_recorder[selected_major_bucket][selected_major_index]);
			indel_recorder_copy(read2_indel_recorder, vote_read2 -> indel_recorder[selected_minor_bucket][selected_minor_index]);
		}
		return 1;
	}

	return 0;
}

void destroy_offsets(gene_offset_t* offsets)
{
	free(offsets->read_names);
	free(offsets->read_offsets);
	HashTableDestroy(offsets->read_name_to_index);
}

int load_offsets(gene_offset_t* offsets , const char index_prefix [])
{
	char fn[MAX_FILE_NAME_LENGTH];
	FILE * fp;
	int n=0;
	int padding = 0;

	int is_V3_index = gehash_load_option(index_prefix, SUBREAD_INDEX_OPTION_INDEX_PADDING , &padding);
	if(!is_V3_index) return 1;
	
	memset(offsets, 0, sizeof( gene_offset_t));
	sprintf(fn, "%s.reads", index_prefix);

	fp = f_subr_open(fn, "r");

	if(!fp)
	{
		SUBREADprintf("file not found :%s\n", fn);
		return 1;
	}

	int current_max_n = 100;
	offsets->read_names = malloc(current_max_n * MAX_READ_NAME_LEN);
	offsets->read_offsets= malloc(current_max_n * sizeof(int));
	offsets->read_name_to_index = HashTableCreate(5000);
	offsets->padding = padding;

	HashTableSetKeyComparisonFunction(offsets->read_name_to_index, my_strcmp);
	HashTableSetHashFunction(offsets->read_name_to_index ,HashTableStringHashFunction);
	HashTableSetDeallocationFunctions(offsets->read_name_to_index ,free, NULL);



	while (!feof(fp))
	{
		int i=0, step = 0, j=0;

		read_line(MAX_FILE_NAME_LENGTH-1,fp, fn, 0);
		if (strlen(fn)<2)continue;
		while (fn[i])
		{
			if (fn[i] == '\t')
			{
				fn[i]=0;
				offsets->read_offsets[n] = (unsigned int)atoll(fn);
				step = 1;
			}
			else if (step)
			{
				if(j<MAX_CHROMOSOME_NAME_LEN-1)
				{
					*(offsets->read_names+MAX_CHROMOSOME_NAME_LEN*n + (j++)) = fn[i];
					*(offsets->read_names+MAX_CHROMOSOME_NAME_LEN*n + j) = 0;
				}
			}
			i++;
		}

		char * read_name_mem = malloc(MAX_CHROMOSOME_NAME_LEN);
		strcpy(read_name_mem, offsets->read_names + n*MAX_CHROMOSOME_NAME_LEN);

		//printf("09NOV: add_annotation_to_junctions: '%s',%u\n", read_name_mem , n);
		HashTablePut(offsets->read_name_to_index, read_name_mem , NULL + 1 + n);

		n++;
		if(n >= current_max_n)
		{
			offsets->read_names = realloc(offsets->read_names, 2*current_max_n * MAX_CHROMOSOME_NAME_LEN);
			offsets->read_offsets = realloc(offsets->read_offsets, 2*current_max_n * sizeof(int));
			current_max_n*=2;
		}


		offsets->read_offsets[n] = 0;

	}

	offsets->total_offsets=n;
	fclose(fp);
	return 0;
}

#ifndef MAKE_STANDALONE
#define CLOCK_USE_GETTIME
#endif

double miltime(){
	double ret;
	#ifdef FREEBSD
		struct timeval tp;
		struct timezone tz;
		tz.tz_minuteswest=0;
		tz.tz_dsttime=0;
		gettimeofday(&tp,&tz);
		ret = tp.tv_sec+ 0.001*0.001* tp.tv_usec; 
	#else
    	#ifdef CLOCK_USE_GETTIME
	 	struct timespec tsc;
		clock_gettime(CLOCK_REALTIME, &tsc);
	 	ret = tsc.tv_sec*1. + tsc.tv_nsec*1./1000/1000/1000;
    	#else
		struct timeb trp;
		ftime(&trp);
		ret = trp.time*1.0+(trp.millitm*1.0/1000.0);
	#endif
	#endif

	return ret;
}

#define front2(str, bias)	(*((str)+(bias))+*((str)+1+(bias)))
#define front4(str, bias)	(front2(str, bias)+front2(str, bias+2))
#define front8(str, bias)	(front4(str, bias)+front4(str, bias+4))
#define front16(str, bias)	(front8(str, bias)+front8(str, bias+8))

float match_read(const char read_str[], int read_len, unsigned int potential_position,  gene_value_index_t * my_array_index, int space_type, int indel_tolerance, const char quality_str [], int quality_scale)
{
	float ret = 0;
	int i, bias;
	char read_matchingness [7][1250];

	if(indel_tolerance>3) indel_tolerance = 3;

	for(bias=-indel_tolerance; bias<=indel_tolerance;bias++)
	{
		for(i=0;i<read_len; i++)
		{
			char base_int = base2int(read_str[i]);
			int is_matched_base =  gvindex_match_base(my_array_index, potential_position+i+bias, base_int);
			read_matchingness[bias+indel_tolerance][i] = is_matched_base; 

/*			if(i % 3 == 2)
			{
				int b2_match = (read_matchingness[bias+indel_tolerance][i -2]  && read_matchingness[bias+indel_tolerance][i -1] && read_matchingness[bias+indel_tolerance][i]);
				if(b2_match)
				{
					read_matchingness[bias+indel_tolerance][i]  = 2;
					read_matchingness[bias+indel_tolerance][i-1]  = 2;
					read_matchingness[bias+indel_tolerance][i-2]  = 2;
				}
			}*/
		}
	}

	for(i=0; i<read_len-4; i+=4)
	{
	
		int j;
		int best_movement = 0; 
		float max_matchness = -1;
		for (j=-indel_tolerance; j<=indel_tolerance; j++)
		{
			int m = front4(read_matchingness[j],i);
			if(m > max_matchness)
			{
				max_matchness = m*1.;
				best_movement = j;
			}
		}

//		SUBREADprintf("best-match=%f\n", max_matchness);
//

		if (quality_str[0])
		{
			max_matchness = 0;
		//	float qs0 = max_matchness;
			for (j=0; j<4; j++)
			{
				if(read_matchingness[best_movement][j+i])
				{
					// penalty
					max_matchness += 1.03+get_base_quality_score(quality_str[j+i], quality_scale);
				}
				else
				{
					// relief
					//max_matchness += 0.2-get_base_phred(quality_str[j+i])*0.01;
				}
			}
		//	SUBREADprintf("Delta matchingness = %f\n", qs0 - max_matchness);
		}
		
		ret += max_matchness; //front8(read_matchingness[movement], i);  // (front4(read_matchingness[movement], i)==4?4:0);
	}


	return  ret;
}


void final_matchingness_scoring(const char read_str[], const char quality_str[], int read_len, gene_vote_t * vote, gehash_data_t * max_position, gene_vote_number_t * max_vote, short *max_mask, gene_quality_score_t * max_quality, gene_value_index_t * my_array_index, int space_type, int indel_tolerance, int quality_scale, gene_vote_number_t * max_indel_recorder, int * max_coverage_start, int * max_coverage_end)
{
	int i, j;
	gene_quality_score_t max_matching = -1.0;

	*max_vote = vote -> max_vote; 

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			if( vote->votes[i][j]  < vote -> max_vote -1 ) continue;

			unsigned int potential_position = vote -> pos[i][j];
			gene_quality_score_t matchingness_count = 1.*match_chro_indel((char *)read_str, my_array_index, potential_position, read_len, 0, space_type, indel_tolerance, NULL, 0);
			if(matchingness_count > max_matching)
			{
				max_matching = matchingness_count;
				*max_position = potential_position;
				*max_mask = vote -> masks[i][j];
				*max_coverage_start = vote-> coverage_start[i][j];
				*max_coverage_end = vote -> coverage_end[i][j];
				indel_recorder_copy(max_indel_recorder, vote-> indel_recorder[i][j]);

				* max_quality = max_matching;
			}
			else if (matchingness_count == max_matching)
			{
//				SUBREADprintf("\nBREAK EVEN DETECTED AT SORTED TABLE: %u (kept) and %u\n", vote->max_position, vote -> pos[i][j]);
				(*max_mask) |= IS_BREAKEVEN_READ ;
			}
		}
}

int DPALIGN_CREATEGAP_PENALTY = -2 ;
int DPALIGN_EXTENDGAP_PENALTY = 0 ;
int DPALIGN_MISMATCH_PENALTY = 0;
int DPALIGN_MATCH_SCORE = 2;



int search_DP_branch(char * read, int read_len, gene_value_index_t * index, unsigned int begin_position, int path_i, int path_j, short ** table,char ** table_mask, int max_indel, char * movement_buffer, int expected_offset, int current_score, int out_pos, int current_offset, int init_read_offset, int shutdown_read_offset, int * all_steps)
{

	if (1499 - out_pos > (read_len <<2) || (*all_steps) > 3000 + (read_len << 3) || (*all_steps) > 33000){
		/*ddprintf("\nTOO MANY STEPS: rl = %d    len=%d    steps=%d\n", read_len, 1499 - out_pos, *all_steps);
		ddfflush(stdout) ;*/
		return 0;
	}
	if(path_j < 0 || path_i < 0)
	{
		int i;
		//ddprintf("\nFINAL TEST exp %d  real-offset %d\n", expected_offset, current_offset);
		if (expected_offset !=current_offset)return 0;
		for (i=0; i<=path_j; i++) movement_buffer[out_pos--] = 1;
		for (i=0; i<=path_i; i++) movement_buffer[out_pos--] = 2;
		return  out_pos ;
	}
	
	short upper_score;
	if (path_i>0) upper_score = table[path_i-1][path_j];
	else	      upper_score = 0;

	short left_score;
	if (path_j>0)      left_score = table[path_i][path_j-1];
	else	      left_score = 0;

	short upperleft_score;
	if (path_j>0 && path_i>0) upperleft_score = table[path_i-1][path_j-1];
	else	      upperleft_score = 0;

	char is_matched_ij =  gvindex_get(index, begin_position + path_j) == read[path_i]?DPALIGN_MATCH_SCORE :DPALIGN_MISMATCH_PENALTY;

	int found = 0;
	int table_mask_j_1 = 0, table_mask_i_1 = 0;
	int table_mask2_i_1 = 0;
	if(path_j >0 && path_i >=0)
	{
		table_mask_j_1 = left_score + (table_mask[path_i][path_j-1]? DPALIGN_EXTENDGAP_PENALTY: DPALIGN_CREATEGAP_PENALTY) == current_score;
	}
	if(path_i >0 && path_j >=0)
	{
		table_mask_i_1 = upper_score + (table_mask[path_i-1][path_j]? DPALIGN_EXTENDGAP_PENALTY: DPALIGN_CREATEGAP_PENALTY) == current_score;
		table_mask2_i_1 = upper_score + (table_mask[path_i-1][path_j]? DPALIGN_EXTENDGAP_PENALTY: DPALIGN_CREATEGAP_PENALTY) < 0 ;
	}


	if ((!found ) && ( table_mask_j_1 || (current_score == 0 && table_mask_i_1)))
	{
		movement_buffer[out_pos] = 1;
		(*all_steps) ++;
		//ddprintf ("Path_i = %d ; Path_j = %d ; Offset = %d => 1\n", path_i, path_j, current_offset);
		//ddfflush(stdout) ;
		found =search_DP_branch (read, read_len, index, begin_position, path_i , path_j -1, table  , table_mask, max_indel, movement_buffer, expected_offset, left_score, out_pos -1, current_offset - ((path_i >= init_read_offset && path_i <= shutdown_read_offset)?1:0),  init_read_offset, shutdown_read_offset, all_steps); 
	}
	if ((!found )&& ( table_mask_i_1 || (current_score == 0 && table_mask2_i_1)))
	{
		movement_buffer[out_pos] = 2;
		(*all_steps) ++;
		//ddprintf ("Path_i = %d ; Path_j = %d ; Offset = %d => 2\n", path_i, path_j, current_offset);
		//ddfflush(stdout) ;
		found =search_DP_branch (read, read_len, index, begin_position, path_i -1 , path_j , table , table_mask , max_indel, movement_buffer, expected_offset, upper_score, out_pos -1, current_offset + ((path_i >= init_read_offset && path_i <= shutdown_read_offset)?1:0),  init_read_offset, shutdown_read_offset, all_steps); 
	}
	if((!found )&& (is_matched_ij + (upperleft_score - current_score) == 0))
	{
		movement_buffer[out_pos] = is_matched_ij==2?0:3;
		//ddprintf ("Path_i = %d ; Path_j = %d ; Offset = %d => 0\n", path_i, path_j, current_offset);
		(*all_steps) ++;
		found = search_DP_branch (read, read_len, index, begin_position, path_i -1 , path_j -1, table , table_mask, max_indel, movement_buffer, expected_offset, upperleft_score, out_pos -1, current_offset, init_read_offset, shutdown_read_offset, all_steps); 
	}
	

	return found;
}

#define INDEL_WINDOW_WIDTH 4

int window_indel_align(char * read, int read_len, gene_value_index_t * index, unsigned int begin_position, int max_indel, char * movement_buffer, int expected_offset, int init_read_offset, int shutdown_read_offset)
{
	short indel_windows[32];
	int scores [32][MAX_READ_LENGTH];
	int i,j,windows_number;
	char chro_str[200];

	memset(indel_windows,0, 32*sizeof(short));

	windows_number = abs(expected_offset)+1;
	//ddprintf ("\nWindow size = %d ; ExpOffset = %d\n", windows_number , expected_offset);
	ddfflush(stdout) ;
	for (i=0; i< read_len; i++)
		for (j=0; j< windows_number; j++)
		{
			int matchingness ;
			if (j==0) chro_str[i] = gvindex_get(index, begin_position+i);
			if (expected_offset<0)
				matchingness = *(read + i) == gvindex_get(index, begin_position + j + i) ;
			else
				matchingness = *(read + i) == gvindex_get(index, begin_position - windows_number + j + 1 + i) ;
			indel_windows [j] += matchingness;

			if (i >= INDEL_WINDOW_WIDTH)// The result for i - INDEL_WINDOW_WIDTH is OK
			{
				scores[j][i - INDEL_WINDOW_WIDTH] = indel_windows [j] ;
			}
			if (i >= INDEL_WINDOW_WIDTH)// The result for i - INDEL_WINDOW_WIDTH is OK
			{

				if (expected_offset<0)
					matchingness = *(read - INDEL_WINDOW_WIDTH + i) == gvindex_get(index, begin_position + j + i - INDEL_WINDOW_WIDTH) ;
				else
					matchingness = *(read - INDEL_WINDOW_WIDTH + i) == gvindex_get(index, begin_position -windows_number + j +1 + i - INDEL_WINDOW_WIDTH) ;
				indel_windows [j] -= matchingness;
			}
		}
	chro_str[i]=0;
	j = read[read_len];
	read[read_len]=0;

	ddprintf ("CHRO=%s\nREAD=%s\n", chro_str, read);

	/*
	float contingency_list [MAX_READ_LENGTH];
	for (i=0; i< read_len - windows_number-INDEL_WINDOW_WIDTH; i++)
	{
		float contingency;
		if (expected_offset >0)
			contingency  = ((scores[0][i]*1. + 1)-(scores[windows_number-1][i]*1.+1));
		else
			contingency  = ((scores[windows_number-1][i]*1. + 1)-(scores[0][i]*1.+1));
		contingency_list[i]=contingency;
	}*/

	int max_score=-1;
	int max_pos = -1;
	if (expected_offset>0)
	{
		//insertion
		for (i=read_len-INDEL_WINDOW_WIDTH-1; i>=0; i--)
		{
			if (scores[windows_number - 1][i - expected_offset]>=2 && scores[0][i] >= max_score)
			{
				max_pos = i - windows_number+1;
				max_score = scores[0][i];
			}
		}
	}
	else
	{
		//deletion
		for (i=read_len-INDEL_WINDOW_WIDTH-1; i>=0; i--)
		{
			if (scores[windows_number -1][i] >= max_score && scores[0][i + expected_offset] >=2)
			{
				max_score = scores[windows_number -1][i];
				max_pos = i ;
			}
		}
	
	}

	max_pos = min(read_len, max(0, max_pos));

	int move_number = 0;

	for (i=0; i< read_len -INDEL_WINDOW_WIDTH; i++)
	{
		if (i==max_pos)
			for (j = 0; j < windows_number-1; j++)
				movement_buffer[move_number++] = 1+(expected_offset>0);	// 1=del; 2=ins
		if (i!=max_pos || expected_offset <0)
			movement_buffer [move_number++] = 0; // matched

		#ifdef indel_debug
		for (j=0; j< windows_number; j++)
			SUBREADprintf("%d ",scores[j][i]);
		SUBREADprintf ("     %c\n", (i==max_pos)?'*':' ') ;
		#endif
	}
	for(; i< read_len; i++)
		movement_buffer [move_number++] = 0; // matched



	read[read_len] = j;
	return move_number;
}

int dynamic_align(char * read, int read_len, gene_value_index_t * index, unsigned int begin_position, int max_indel, char * movement_buffer, int expected_offset, int init_read_offset, int shutdown_read_offset, short **table  ,  char ** table_mask)
// read must be converted to the positive strand.
// movement buffer: 0:match, 1: read-insert, 2: gene-insert, 3:mismatch
// the size of the movement buffer must be equal to the length of the read plus max_indel * 3.
{


	int i,j;

	#ifdef indel_debug
	int ii, jj;
	for(ii = 0; ii<read_len - expected_offset; ii++)
		putchar(gvindex_get(index, begin_position + ii));

	SUBREADprintf ("\n%s\n", read);
	#endif

	for (i=0; i<read_len ; i++)
	{
		for(j=0; j<read_len - expected_offset; j++)
		{
			table_mask[i][j]=0;
			if (j < i - max_indel || j > max_indel + i)
			{
				table[i][j]=-9999;
			#ifdef indel_debug
				putchar('\t');
			#endif
				continue;
			}

			short from_upper;

			if (i>0) from_upper = table[i-1][j] + (table_mask[i-1][j]?DPALIGN_EXTENDGAP_PENALTY:DPALIGN_CREATEGAP_PENALTY);
			else     from_upper = DPALIGN_CREATEGAP_PENALTY;

			short from_left;

			if (j>0) from_left = table[i][j-1] + (table_mask[i][j-1]?DPALIGN_EXTENDGAP_PENALTY:DPALIGN_CREATEGAP_PENALTY);
			else     from_left = DPALIGN_CREATEGAP_PENALTY;

			char chromo_ch = gvindex_get(index, begin_position + j);
			char is_matched_ij = (chromo_ch == read[i])?DPALIGN_MATCH_SCORE:DPALIGN_MISMATCH_PENALTY;

			short from_upperleft;
			if (i>0 && j>0) from_upperleft = table[i-1][j-1] + is_matched_ij;
			else	    from_upperleft = is_matched_ij; 
			if ((i ==0 || j ==0) && (i+j>0)) from_upperleft += DPALIGN_CREATEGAP_PENALTY;

			if (from_upperleft <=from_left || from_upperleft <= from_upper)
				table_mask[i][j]=1;

			table[i][j]=max(0,max(from_upper, max(from_left, from_upperleft)));
			#ifdef indel_debug
			SUBREADprintf("%c%c\t", chromo_ch, read[i]);
			#endif
		}
		#ifdef indel_debug
		SUBREADputs("");
		#endif

	}
	#ifdef indel_debug
	//SUBREADputs("");
	//SUBREADputs("");

	#endif

	short current_score = table[read_len-1] [read_len - expected_offset-1];
	short path_i = read_len-1 ;
	int out_pos = 1499;
	char out_tmp [1500];
	int all_steps = 0;
	j = read_len-expected_offset -1;

	#ifdef indel_debug

	for(ii=0;ii< 20; ii++)
	{
		for(jj=0; jj<20; jj++)
			SUBREADprintf("%d\t",table[ii][jj]);
		SUBREADputs("");
	}
		SUBREADputs("");
		SUBREADputs("");

	for(ii=0;ii< 20; ii++)
	{
		for(jj=0; jj<20; jj++)
			SUBREADprintf("%d\t",table_mask[ii][jj]);
		SUBREADputs("");
	}
	#endif

	out_pos = search_DP_branch(read, read_len, index, begin_position, path_i, j, table, table_mask, max_indel,  out_tmp , expected_offset,  current_score , out_pos, 0, init_read_offset,shutdown_read_offset, &all_steps); 

	if (out_pos)
	{
		memcpy(movement_buffer, out_tmp  + out_pos +1, 1499 - out_pos);
		return 1499 - out_pos;
	}
	else return 0;

}

int extend_covered_region(gene_value_index_t *array_index, unsigned int read_start_pos, char * read, int read_len, int cover_start, int cover_end, int window_size, int req_match_5end , int req_match_3end, int indel_tolerance, int space_type, int tail_indel, short * head_indel_pos, int * head_indel_movement, short * tail_indel_pos, int * tail_indel_movement, int is_head_high_quality, char * qual_txt, int qual_format, float head_matching_rate, float tail_matching_rate)
{
	int ret = 0;
	*head_indel_pos = -1;
	*tail_indel_pos = -1;
	if (cover_start >= window_size && head_matching_rate < 1.0001) 
	{
		int head_test_len =  cover_start;
		int roughly_mapped = match_chro(read, array_index, read_start_pos, head_test_len , 0, space_type);
		if (roughly_mapped >= head_test_len  * CORE_RECOVER_MATCHING_RATE - 0.0001)
		{
			ret |= 1;
		}
		else
		{
			int window_end_pos = cover_start + window_size-1;
			int not_too_bad = 1;
			int right_match_number=0;

			while (1)
			{
				int matched_bases_in_window = match_chro_wronglen(read+ window_end_pos - window_size, array_index, read_start_pos + window_end_pos - window_size, window_size, space_type, NULL, &right_match_number);
				int best_indel_movement = -1000;
				int best_indel_pos = -1;

				if (matched_bases_in_window >= req_match_5end) // this window is not bad enough so that an indel is considered
					window_end_pos--;
				else
				{
					roughly_mapped = match_chro(read, array_index, read_start_pos, window_end_pos - right_match_number , 0, space_type);
					if (roughly_mapped < (int)(0.5+ (window_end_pos - right_match_number)  * CORE_RECOVER_MATCHING_RATE))
					{
						// the quality of this window is too low (at most 1 base is matched). I must consider if it is an indel.
						int indel_movement_i ;

						int best_matched_bases_after_indel = -1;

						not_too_bad = 0;
						for (indel_movement_i = 0; indel_movement_i < 2* indel_tolerance-1 ; indel_movement_i ++)
						{
							int indel_movement = (indel_movement_i+1)/2 * (indel_movement_i %2?1:-1);
							int test_length = window_end_pos /*- 1*/ - max(0, indel_movement) -  right_match_number;

							if (test_length < window_size) continue;
							//if (test_length <= 1+ abs(indel_movement)/4) continue;
							//test_length = min(10, test_length);
							if (abs(indel_movement) > indel_tolerance) continue;

							int test_start =0;// window_end_pos - max(0, indel_movement) -  right_match_number - test_length;

							int matched_bases_after_indel = match_chro_support(read +test_start, array_index, read_start_pos + indel_movement +test_start, test_length,0, space_type, qual_txt, qual_format);
							float test_rate = head_matching_rate;
							
							if(test_length < 3) test_rate = 1;

							if(best_matched_bases_after_indel <matched_bases_after_indel  && matched_bases_after_indel >= (int)(0.5+test_length * test_rate))
							{
								not_too_bad = 1;

								best_matched_bases_after_indel = matched_bases_after_indel;
								best_indel_movement = indel_movement;
								best_indel_pos = window_end_pos - right_match_number;


								
								*head_indel_pos = best_indel_pos;
								*head_indel_movement = best_indel_movement;
							}
						}
						if(best_indel_pos<0) *head_indel_pos =  window_end_pos - right_match_number;
						break;
					}else window_end_pos--;
				}
				if (window_end_pos - window_size <= 0) break;
			}
			if(not_too_bad)
				ret |=1;
			//else
				//*head_indel_pos = -1;
		}
	}
	else ret |= 1;



	if (cover_end <= read_len - window_size && tail_matching_rate < 1.0001) 
	{
		int tail_test_len = read_len - cover_end;
		int roughly_mapped = match_chro(read + cover_end, array_index, read_start_pos + cover_end + tail_indel, tail_test_len , 0, space_type);
		if (roughly_mapped >= tail_test_len  * CORE_RECOVER_MATCHING_RATE - 0.0001)
		{
			ret |= 2;
		}
		else
		{
			int window_start_pos = cover_end - window_size +1;
			int not_too_bad = 1;

			while (1)
			{
				int left_match_number = 0;
				int matched_bases_in_window = match_chro_wronglen(read+ window_start_pos, array_index, read_start_pos + window_start_pos + tail_indel, window_size, space_type, &left_match_number, NULL);
				int best_indel_movement = -1000;
				int best_indel_pos = -1;

				if (matched_bases_in_window >= req_match_3end) // this window is not bad enough so that an indel is considered
					window_start_pos++;
				else
				{
					roughly_mapped = match_chro(read+window_start_pos + left_match_number, array_index, read_start_pos + window_start_pos + tail_indel + left_match_number, read_len - window_start_pos - left_match_number , 0, space_type);

					if (roughly_mapped < (int)(0.5+(read_len - window_start_pos -left_match_number )  * CORE_RECOVER_MATCHING_RATE))
					{
						// the quality of this window is too low (at most 1 base is matched). I must consider if it is an indel.
						int indel_movement_i ;
						int best_matched_bases_after_indel = -1;
						not_too_bad = 0;
						for (indel_movement_i = 0; indel_movement_i < 2* indel_tolerance; indel_movement_i ++)
						{
							
							int indel_adjustment = (indel_movement_i+1)/2 * (indel_movement_i %2?1:-1);
							int indel_movement = tail_indel + indel_adjustment;
							int test_length = read_len - window_start_pos  - left_match_number + min(0, indel_adjustment);
							//test_length = min(10, test_length);


							if (test_length < window_size) continue;
							//if (test_length <= 1 + abs(indel_movement)/4) continue;
							if (abs(indel_movement) > indel_tolerance) continue;

							int matched_bases_after_indel = match_chro_support(read + window_start_pos - min(0, indel_adjustment) + left_match_number, array_index, read_start_pos + window_start_pos  + max(0,indel_movement) +left_match_number , test_length,0, space_type, qual_txt +  window_start_pos - min(0, indel_adjustment) + left_match_number , qual_format);

							float test_rate = tail_matching_rate;
							
							if(test_length < 3) test_rate = 1;

							if(best_matched_bases_after_indel <matched_bases_after_indel  && matched_bases_after_indel >= (int)(0.5+test_length * test_rate))
							{
								not_too_bad = 1;
								best_matched_bases_after_indel = matched_bases_after_indel;
								best_indel_movement = indel_movement;
								best_indel_pos = window_start_pos + left_match_number ;//-1;
								*tail_indel_movement = best_indel_movement ;
							}
						}
					
						if(best_indel_pos<0)
							*tail_indel_pos =  window_start_pos + left_match_number ;
						else
							*tail_indel_pos = best_indel_pos;
						break;
					}else window_start_pos++;
				}
				if (window_start_pos + window_size >= read_len) break;
			}
			if(not_too_bad)
				ret |=2;
			//else
			//	*tail_indel_pos =-1;
		}
	}
	else ret |= 2;


	return ret;
}


float match_base_quality_cs(gene_value_index_t *array_index, char * read_txt,  unsigned int pos, char * qual_txt, int read_len, int phred_version, int * high_qual_unmatch, int * all_mismatched, int ql_kill, int head_clipped, int tail_clipped)
{
	int i;
	int ret =0;
	char lastch;
	if(pos < array_index -> start_base_offset || pos + read_len >= array_index -> start_base_offset + array_index -> length){
		//SUBREADprintf("WARNING: BASE INDEX OUT OF LIMIT: %u < %u < %u\n%s\n", array_index -> start_base_offset , pos, array_index -> start_base_offset + array_index -> length, read_txt);
	//	exit(-1);
		return (read_len - tail_clipped - head_clipped);
	}
	lastch = gvindex_get(array_index, pos);
	for(i=head_clipped; i<read_len - tail_clipped; i++)
	{
		char nch = gvindex_get(array_index, pos+i+1);
		int is_matched = read_txt[i] == '0'+chars2color(lastch, nch);
		if(is_matched){
			ret ++;
		}else
		{
			(*all_mismatched)++;
			(*high_qual_unmatch)++;
			ret --;
		}

		lastch = nch;
	}
	return ret*1.;
}

float match_base_quality(gene_value_index_t *array_index, char * read_txt,  unsigned int pos, char * qual_txt, int read_len, int is_negative, int phred_version, int * high_qual_unmatch, int * all_mismatched, int ql_kill, int head_clipped, int tail_clipped)
{
	int i;
	int ret =0;
	if(pos < array_index -> start_base_offset || pos + read_len >= array_index -> start_base_offset + array_index -> length){
		//SUBREADprintf("WARNING: BASE INDEX OUT OF LIMIT: %u < %u < %u\n%s\n", array_index -> start_base_offset , pos, array_index -> start_base_offset + array_index -> length, read_txt);
	//	exit(-1);
		return (read_len - tail_clipped - head_clipped);
	}
	for(i=head_clipped; i<read_len - tail_clipped; i++)
	{
		char true_chr;
		if(is_negative)
		{
			true_chr = gvindex_get(array_index, pos + read_len - i - 1);
			if(true_chr == 'A') true_chr = 'T';
			else if(true_chr == 'G') true_chr = 'C';
			else if(true_chr == 'C') true_chr = 'G';
			else  true_chr = 'A';
		}
		else
			true_chr = gvindex_get(array_index, pos + i);

		//SUBREADprintf("%c vs %c\n", true_chr , read_txt[i]);

		if (true_chr == read_txt[i])
		{
			if(!qual_txt)
				ret += 1000000;
			else if(FASTQ_PHRED64 == phred_version)
				ret += (1000000-get_base_error_prob64i(qual_txt[i]));
			else
				ret += (1000000-get_base_error_prob33i(qual_txt[i]));
		}
		else
		{
			(*all_mismatched)++;
			if(!qual_txt)
			{
				ret -= 1000000;
				(*high_qual_unmatch)++;
			}
			else
			{
				int ql ;
				if(FASTQ_PHRED64 == phred_version)
					ql = get_base_error_prob64i(qual_txt[i]);
				else
					ql = get_base_error_prob33i(qual_txt[i]);
				/*
				#ifdef QUALITY_KILL
					#if QUALITY_KILL > 196
						#define QL_MIN 999000
					#endif
				#endif
	
				#ifndef QL_MIN
					#define QL_MIN 200000
				#endif*/

				if( ql < ql_kill) (*high_qual_unmatch)++;

				ret += ql-1000000;
			}
		}
	}
	//printf ("SECTION QUAL = %d LEN = %d\n", ret, read_len);
	return ret/1000000.;
}

float final_mapping_quality(gene_value_index_t *array_index, unsigned int pos, char * read_txt, char * qual_txt, char * cigar_txt, int phred_version, int * mismatch, int rl, char * refined_cigar, unsigned int * new_pos)
{
	int cigar_cursor = 0;
	int read_cursor = 0;
	long long int chromosome_cursor = pos, x;
	int cigar_length = strlen(cigar_txt);
	int i;
	int ret = 0;

	x= 0;
	while(cigar_cursor < cigar_length)
	{
		if(cigar_txt[cigar_cursor] =='X')
		{
			cigar_cursor++;
			continue;
		}
		if(cigar_txt[cigar_cursor]>='0' && cigar_txt[cigar_cursor]<='9')
			x = x*10+ (cigar_txt[cigar_cursor]-'0');
		else
		{
			if(cigar_txt[cigar_cursor] == 'M' || cigar_txt[cigar_cursor] == 'S') 
			{
				int all_MM=0;
				float nret = match_base_quality(array_index, read_txt + read_cursor, chromosome_cursor , (qual_txt && qual_txt[0])?qual_txt + read_cursor:NULL, x, 0, phred_version, mismatch, &all_MM, 200000,0,0);
				//printf ("%s: Q=%.6f; L=%d ; POS=%u\n",  read_txt + read_cursor, nret, x, chromosome_cursor);

				ret += (int)(nret*1000000);
				chromosome_cursor +=x;
				read_cursor +=x;
			}
			else if(cigar_txt[cigar_cursor] == 'I')
			{
				if(!qual_txt)
				{
					ret += x*1000000;
					read_cursor +=x;
				}
				else if(FASTQ_PHRED64 == phred_version)
				{
					for(i = 0; i<x; i++)
					{
						char qchar = qual_txt[read_cursor++];
						ret += (1000000-get_base_error_prob64i(qchar));
					}
				}
				else
				{
					for(i = 0; i<x; i++)
					{
						char qchar = qual_txt[read_cursor++];
						ret += (1000000-get_base_error_prob33i(qchar));
					}
				}
			}
			else if(cigar_txt[cigar_cursor] == 'D' || cigar_txt[cigar_cursor] == 'N'|| cigar_txt[cigar_cursor] == 'j' || cigar_txt[cigar_cursor] == 'J')
				chromosome_cursor +=x;
			else if(cigar_txt[cigar_cursor] == 'B' || cigar_txt[cigar_cursor] == 'b')
				chromosome_cursor -=x;

			x= 0;
		}
		cigar_cursor++;
	}
	if(read_cursor != rl)
	{
		*mismatch = 9999;
		return 0.;
	}
	//printf("S=%.5f, LEN=%d\n", ret , read_cursor);

	float match_rate = (ret/10000.) / read_cursor+100.;

	if(refined_cigar && match_rate < 190)
	{
		char read_data[1250];
		chromosome_cursor = pos;
		read_cursor = 0;
		cigar_cursor = 0;
		x=0;

		while(cigar_cursor < cigar_length)
		{
			if(cigar_txt[cigar_cursor] =='X')
			{
				cigar_cursor++;
				continue;
			}
			if(cigar_txt[cigar_cursor]>='0' && cigar_txt[cigar_cursor]<='9')
				x = x*10+ (cigar_txt[cigar_cursor]-'0');
			else
			{
				if(cigar_txt[cigar_cursor] == 'M' || cigar_txt[cigar_cursor] == 'S') 
				{
					gvindex_get_string(read_data + read_cursor  ,  array_index ,  chromosome_cursor , x , 0);
					chromosome_cursor +=x;
					read_cursor +=x;
				}
				else if(cigar_txt[cigar_cursor] == 'I')
				{
					int xk;
					for(xk=read_cursor;xk<read_cursor+x; xk++)
						read_data[cigar_cursor] = 'I';
					read_cursor +=x;
				}
				else if(cigar_txt[cigar_cursor] == 'D' || cigar_txt[cigar_cursor] == 'N'|| cigar_txt[cigar_cursor] == 'j' || cigar_txt[cigar_cursor] == 'J')
					chromosome_cursor +=x;
				else if(cigar_txt[cigar_cursor] == 'B' || cigar_txt[cigar_cursor] == 'b')
					chromosome_cursor -=x;

				x= 0;
			}
			cigar_cursor++;
		}
		
		int window_matched = 0;
		int last_confirmed_read_pos = 99999;
		int REFINE_WINDOW_SIZE = 10;
		for(read_cursor = - REFINE_WINDOW_SIZE; read_cursor < rl - REFINE_WINDOW_SIZE; read_cursor++)
		{

//			printf("MMM %c =?= %c\n", read_txt[read_cursor+REFINE_WINDOW_SIZE] , read_data[read_cursor+REFINE_WINDOW_SIZE]);
			window_matched += (read_txt[read_cursor+REFINE_WINDOW_SIZE] == read_data[read_cursor+REFINE_WINDOW_SIZE] || read_txt[read_cursor+REFINE_WINDOW_SIZE] =='I');
			if(read_cursor>=0)
			{
				if(window_matched >= REFINE_WINDOW_SIZE -1)
					if(read_txt[read_cursor+REFINE_WINDOW_SIZE] == read_data[read_cursor+REFINE_WINDOW_SIZE])
						last_confirmed_read_pos = read_cursor+REFINE_WINDOW_SIZE+1;
				window_matched -= (read_txt[read_cursor] == read_data[read_cursor] || read_txt[read_cursor] == 'I');
			}
		} 


		window_matched = 0;
		int first_confirmed_read_pos = 99999;
		for(read_cursor = rl -1 ; read_cursor >=0 ; read_cursor--)
		{

			window_matched += (read_txt[read_cursor] == read_data[read_cursor]|| read_txt[read_cursor] == 'I');
			if(read_cursor <= rl- REFINE_WINDOW_SIZE - 1)
			{
				if(window_matched >= REFINE_WINDOW_SIZE -1)
					if(read_txt[read_cursor] == read_data[read_cursor])
						first_confirmed_read_pos = read_cursor;
				window_matched -= (read_txt[read_cursor+REFINE_WINDOW_SIZE] == read_data[read_cursor+REFINE_WINDOW_SIZE] || read_txt[read_cursor + REFINE_WINDOW_SIZE] == 'I');
			}
		} 
	

		cigar_cursor = 0;
		read_cursor = 0;
		refined_cigar[0]=0;
		x=0;
		int begin_copy  =0;
		chromosome_cursor = pos;

		if(first_confirmed_read_pos<99990 && last_confirmed_read_pos<99990)
		{
			if(first_confirmed_read_pos>0)
				sprintf(refined_cigar,"%dS",first_confirmed_read_pos);
			while(cigar_cursor < cigar_length)
			{
				if(cigar_txt[cigar_cursor] =='X')
				{
					cigar_cursor++;
					continue;
				}
				if(cigar_txt[cigar_cursor]>='0' && cigar_txt[cigar_cursor]<='9')
					x = x*10+ (cigar_txt[cigar_cursor]-'0');
				else
				{
					if(cigar_txt[cigar_cursor] == 'M' || cigar_txt[cigar_cursor] == 'S' || cigar_txt[cigar_cursor] == 'I') 
					{
						int out_len = 0;
						if(read_cursor < first_confirmed_read_pos && read_cursor + x >= first_confirmed_read_pos)
							out_len = read_cursor + x - first_confirmed_read_pos;
						else if(first_confirmed_read_pos <= read_cursor)
							out_len = x;

						if(out_len)
						{
							if(!begin_copy)
								*new_pos = chromosome_cursor - read_cursor; 
							if(read_cursor >= last_confirmed_read_pos) out_len=0;
							else if(read_cursor < last_confirmed_read_pos && read_cursor +x>last_confirmed_read_pos)
								out_len -= (read_cursor +x - last_confirmed_read_pos);
						}

						if(out_len)
						{
							sprintf(refined_cigar+strlen(refined_cigar), "%d%c", out_len, cigar_txt[cigar_cursor]);
							begin_copy=1;
						}
						if(cigar_txt[cigar_cursor] == 'M' || cigar_txt[cigar_cursor] == 'S')
							chromosome_cursor +=x;

						read_cursor +=x;
						if(read_cursor >= last_confirmed_read_pos) break;
					}
					else
					{
						if(begin_copy)
							#ifdef __MINGW32__
							sprintf(refined_cigar+strlen(refined_cigar), "%I64d%c", x, cigar_txt[cigar_cursor]);
							#else
							sprintf(refined_cigar+strlen(refined_cigar), "%lld%c", x, cigar_txt[cigar_cursor]);
							#endif

						if(cigar_txt[cigar_cursor] == 'D' || cigar_txt[cigar_cursor] == 'N'|| cigar_txt[cigar_cursor] == 'j' || cigar_txt[cigar_cursor] == 'J')
							chromosome_cursor +=x;
						else if(cigar_txt[cigar_cursor] == 'B' || cigar_txt[cigar_cursor] == 'b')
							chromosome_cursor -=x;


					} 

					x= 0;
				}
				cigar_cursor++;
			}
			if(last_confirmed_read_pos<=rl-1)
				sprintf(refined_cigar+strlen(refined_cigar),"%dS",rl-last_confirmed_read_pos);
		}
		else if(new_pos)*new_pos=pos;
	} 
	else if(new_pos)*new_pos=pos;

	return match_rate;
}


void print_votes(gene_vote_t * vote, char *index_prefix)
{

	gene_offset_t offsets;
	int i,j;
	char * chrname = NULL;
	int chrpos = 0;
	load_offsets (&offsets, index_prefix);
	//locate_gene_position(vote -> max_position, &offsets, &chrname, &chrpos);

	SUBREADprintf(" ==========   Max votes = %d   ==========\n", vote->max_vote);// , Position is %s,%u\n", vote->max_vote, chrname, chrpos );
	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			locate_gene_position(vote -> pos[i][j], &offsets, &chrname, &chrpos);
			int toli = vote->toli [i][j];
			int last_offset = vote -> indel_recorder[i][j][toli+2];
			SUBREADprintf("  %s\tVote = %d , Position is %s,%d (+%u) Coverage is (%d, %d) Indel %d %s (%d)\n", vote->votes[i][j] == vote->max_vote?"***":"   ", vote->votes[i][j] , chrname, chrpos, vote -> pos[i][j], vote -> coverage_start[i][j], vote -> coverage_end[i][j],last_offset, (vote -> masks[i][j])?"NEG":"POS", vote -> masks[i][j]);

			if(1){
				int k;
				for(k=0; k<=toli; k+=3){
					SUBREADprintf("    %d - %d : D=%d    ",  vote -> indel_recorder[i][j][k], vote -> indel_recorder[i][j][k+1], vote -> indel_recorder[i][j][k+2]);
				}
				SUBREADputs("");
			}
		}
	

}


static float PROB_QUAL_TABLE[] = {1.000000000 , 0.794328235 , 0.630957344 , 0.501187234 , 0.398107171 , 0.316227766 , 0.251188643 , 0.199526231 , 0.158489319 , 0.125892541 , 0.100000000 , 0.079432823 , 0.063095734 , 0.050118723 , 0.039810717 , 0.031622777 , 0.025118864 , 0.019952623 , 0.015848932 , 0.012589254 , 0.010000000 , 0.007943282 , 0.006309573 , 0.005011872 , 0.003981072 , 0.003162278 , 0.002511886 , 0.001995262 , 0.001584893 , 0.001258925 , 0.001000000 , 0.000794328 , 0.000630957 , 0.000501187 , 0.000398107 , 0.000316228 , 0.000251189 , 0.000199526 , 0.000158489 , 0.000125893 , 0.000100000 , 0.000079433 , 0.000063096 , 0.000050119 , 0.000039811 , 0.000031623 , 0.000025119 , 0.000019953 , 0.000015849 , 0.000012589 , 0.000010000 , 0.000007943 , 0.000006310 , 0.000005012 , 0.000003981 , 0.000003162 , 0.000002512 , 0.000001995 , 0.000001585 , 0.000001259 , 0.000001000 , 0.000000794 , 0.000000631 , 0.000000501};




float get_base_error_prob33(char v)
{
	return PROB_QUAL_TABLE[v-'!'];
}
float get_base_error_prob64(char v)
{
	return PROB_QUAL_TABLE[v-'@'];
}

static int PROB_QUAL_INT_TABLE[] = { 1000000 , 794328 , 630957 , 501187 , 398107 , 316227 , 251188 , 199526 , 158489 , 125892 , 100000 , 79432 , 63095 , 50118 , 39810 , 31622 , 25118 , 19952 , 15848 , 12589 , 10000 , 7943 , 6309 , 5011 , 3981 , 3162 , 2511 , 1995 , 1584 , 1258 , 1000 , 794 , 630 , 501 , 398 , 316 , 251 , 199 , 158 , 125 , 100 , 79 , 63 , 50 , 39 , 31 , 25 , 19 , 15 , 12 , 10 , 7 , 6 , 5 , 3 , 3 , 2 , 1 , 1 , 1 , 1 , 0 , 0 , 0 , };


int get_base_error_prob33i(char v)
{
	return PROB_QUAL_INT_TABLE[v-'!'];
}
int get_base_error_prob64i(char v)
{
	return PROB_QUAL_INT_TABLE[v-'@'];
}

#ifdef SKIP_THIS_PART
void bad_reverse_cigar(char * cigar)
{
	int cigar_cursor = 0;
	srInt_64 tmpv=0;
	char ncg[100];
	ncg[0]=0;
	while(cigar[cigar_cursor])
	{
		char cc = cigar[cigar_cursor];
		if(isdigit(cc))
		{
			tmpv=tmpv*10+(cc-'0');
		}
		else if((cc>='A'&&cc<='Z')||(cc>='a'&&cc<='z')) 
		{
			char ncg2[103];
			#ifdef __MINGW32__
			sprintf(ncg2, "%I64d%c", tmpv, cc);
			#else
			sprintf(ncg2, "%lld%c", tmpv, cc);
			#endif
			strncat(ncg2, ncg, 99);
			strncpy(ncg, ncg2, 99);
			tmpv=0;
		}
		else
		{
			char ncg2[103];
			sprintf(ncg2, "%c%s",cc,ncg);
			strncpy(ncg, ncg2, 99);
			tmpv=0;
		}
		cigar_cursor++;
	}
	strcpy(cigar, ncg);
}

#ifdef TEST_WHATCIGARREVERSE
int main()
#else
int debug_main()
#endif
{
	char cg[100];
	sprintf(cg,"10M8H20M9I100N30M");
	bad_reverse_cigar(cg);
	SUBREADprintf("%s\n",cg);
	return 0;
}
#endif

void remove_indel_neighbours(HashTable * indel_table)
{

	unsigned int * to_delete;
	char * to_delete_len;
	int xi;
	KeyValuePair * cursor;
	return;

	to_delete=(unsigned int* )malloc(sizeof(int)*400000);
	to_delete_len=(char* )malloc(sizeof(char)*400000);
	int num_delete=0, bucket;

	for(bucket=0; bucket<indel_table -> numOfBuckets; bucket++)
	{
		cursor = indel_table -> bucketArray[bucket];
		while (1)
		{
			if(!cursor) break;
			indel_record_t * p = (indel_record_t *) cursor -> key;
			indel_result_t * pv = (indel_result_t *)cursor -> value;

			indel_record_t p_nb;
			p_nb.len = p->len;

			for(xi=-5;xi<=5;xi++)
			{
				if(!xi)continue;

				p_nb.pos = p->pos+xi;
				indel_result_t * p_nbv = (indel_result_t *)HashTableGet(indel_table, &p_nb);
				if(p_nbv && (p_nbv->support > pv->support || (p_nbv->support == pv->support && xi>0)))
				{
					if(num_delete<399999)
					{
						to_delete_len[num_delete]=p->len;
						to_delete[num_delete++]= p->pos+xi;
					}
					break;
				}
			}

			cursor = cursor->next;
		}
	}
	
	
	for(xi=0;xi<num_delete;xi++)
	{
		indel_record_t p;
		p.len = to_delete_len[xi];
		p.pos = to_delete[xi];
		HashTableRemove(indel_table, &p);
	}


	SUBREADprintf("\n %d low-confidence indels have been removed.", num_delete);
	free(to_delete);
	free(to_delete_len);

}

void print_version_info()
{
	SUBREADprintf("\nSubread %s\n", SUBREAD_VERSION);
	SUBREADprintf("\n");
}

int fc_strcmp_chro(const void * s1, const void * s2)
{

//	// we can compare the 3-th and 4-th bytes because we know that the buffers have enough lengths.
//	if(((char *)s1)[4] != ((char *)s2)[4] || ((char *)s1)[3] != ((char *)s2)[3] )
//		return 1;
	return strcmp((char*)s1, (char*)s2);
}



srUInt_64 fc_chro_hash(const void *key) {
	const unsigned char *str = (const unsigned char *) key;

	int xk1;
	unsigned long hashValue = 0;

	for(xk1=0; str[xk1]; xk1++)
	{
		unsigned long ch = str[xk1];
		hashValue += (ch + xk1) << (ch & 0xf);
	}
		

	return hashValue;
}


