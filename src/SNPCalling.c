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
  
  

#include <dirent.h> 
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <pthread.h>
#include <assert.h>
#include <getopt.h>
#include "subread.h"
#include "hashtable.h"
#include "HelperFunctions.h"
#include "gene-algorithms.h"
#include "SNPCalling.h"
#include "input-files.h"

#define HUGE_PILE_HEIGHT 300

#define SNP_INTERVAL_SUP 0

extern unsigned int BASE_BLOCK_LENGTH;

float LOG_10_2 ;

int (*snp_progress_report_callback) (int,int,int);

struct SNP_Calling_Parameters{

	float supporting_read_rate;
	int empty_blocks;
	int max_supporting_read_number;
	int min_supporting_read_number;
	int min_alternative_read_number;

	int neighbour_filter_testlen;
	float neighbour_filter_rate;
	int bases_ignored_head_tail;
	int all_blocks;
	int is_phred_64;
	int is_BAM_file_input;
	int is_paired_end_data;
	int is_coverage_calculation;
	int use_soft_clipping_bases;

	int fisher_exact_testlen;

	int min_phred_score;
	int excellent_phred_score;
	int final_phred_score;
	int test_two_strands;
	int neighbour_abundance_test;
	int snp_interval;
	float cutoff_multiplex ; 
	float cutoff_upper_bound ; 

	pthread_spinlock_t * output_fp_lock;

	char pile_file_name[MAX_FILE_NAME_LENGTH];
	int delete_piles;
	int disk_is_full;

	char background_input_file[MAX_FILE_NAME_LENGTH];
	char subread_index[MAX_FILE_NAME_LENGTH];
	char known_SNP_vcf[MAX_FILE_NAME_LENGTH];
	unsigned int known_SNPs_number;
	gene_offset_t * subread_index_offsets;
	gene_value_index_t * subread_index_array;
	HashTable * cigar_event_table;
	char * rebuilt_command_line;

	srInt_64 all_mapped_bases;
	unsigned int fisher_normalisation_target; 

	double start_time;
	unsigned int real_read_count;
	unsigned int reported_SNPs;
	unsigned int reported_indels;
	int supporting_read_excessed_reported;
};

#define PRECALCULATE_FACTORIAL 2000000

long double * precalculated_factorial;// [PRECALCULATE_FACTORIAL];

double factorial_float_real(int a)
{

	double ret = 0;
	while(a)
		ret += log(a--);
	return ret;
}


double factorial_float(int a)
{
	if(a<PRECALCULATE_FACTORIAL && (precalculated_factorial[a]!=0.))
		return precalculated_factorial[a]; 
	else
	{
		double ret = factorial_float_real(a);
		if(a<PRECALCULATE_FACTORIAL) precalculated_factorial[a]=ret;
		return ret;
	}
}

double fisherSub(int a, int b, int c, int d)
{
	double ret = factorial_float(a+b) + factorial_float(c+d) + factorial_float(a+c) + factorial_float(b+d) ;
	ret -= factorial_float(a) + factorial_float(b) + factorial_float(c) + factorial_float(d) + factorial_float(a+b+c+d);
	return pow(2.71828183, ret);
}




/**
 * See HELP string or run with no arguments for usage.
 * <p>
 * The code used to calculate a Fisher p-value comes originally from a
 * <a href="http://infofarm.affrc.go.jp/~kadasowa/fishertest.htm">JavaScript program</a>
 * by T. Kadosawa (kadosawa@niaes.affrc.go.jp).
 * Retrieved from http://www.users.zetnet.co.uk/hopwood/tools/StatTests.java on 3/Jul/2012
 *
 * @author David Hopwood
 * @date   2000/04/23
 */

double fisher_exact_test(int a, int b, int c, int d)
{

	if(a*1./c < b*1./d)	return 1.1;
		// the abnormal level at this base should be at least as the noise level.

	//printf("FET: %d %d %d %d\n", a, b, c, d);


	if(1){
		double ret = fast_fisher_test_one_side(a,b,c,d, precalculated_factorial, PRECALCULATE_FACTORIAL);
		//#warning ">>>>>> FOR ACCURACY MEASUREMENT ONLY <<<<<<<<"
		//SUBREADprintf("FISHER_RES %d %d %d %d %.19G %.19G\n", a,b,c,d, ret, log(ret));
		return ret;
	}else{
		int AA=a, BB=b, CC=c, DD=d;
		if (a * d > b * c) {
		    a = a + b; b = a - b; a = a - b; 
		    c = c + d; d = c - d; c = c - d;
		}

		if (a > d) { a = a + d; d = a - d; a = a - d; }
		if (b > c) { b = b + c; c = b - c; b = b - c; }

		double p_sum = 0.0;

		double p = fisherSub(a, b, c, d);
		while (a >= 0) {
		    p_sum += p;
		    if (a == 0) break;
		    --a; ++b; ++c; --d;
		    p = fisherSub(a, b, c, d);
		}

		if(1){
			// DON'T PUT IT HERE!
			double r1 = fast_fisher_test_one_side(AA,BB,CC,DD, NULL, PRECALCULATE_FACTORIAL);
			if(abs(r1-p_sum) / max(r1,p_sum) >= 0.01){
				printf("BADR: FAST = %.7f != JAVA %.7f,  %d %d %d %d\n", log(r1), log(p_sum),AA,BB,CC,DD);
			}else{
				printf("GOODR: FAST = %.7f ~= JAVA %.7f,  %d %d %d %d\n", log(r1), log(p_sum),AA,BB,CC,DD);
			}
		}

		return p_sum;
	}
}

unsigned int fisher_test_size;

void mask_snp_bitmap(char * bm, unsigned int pos)
{
	int bytes = (pos >> 3);
	int bits = pos & 0x7;

	bm [bytes] |= (1 << bits);
}

int is_snp_bitmap(char * bm, unsigned int pos)
{
	int bytes = (pos >> 3);
	int bits = pos & 0x7;

	return bm [bytes] & (1 << bits);

}
void put_hash_to_pile(HashTable * merge_table, unsigned int* snp_voting_piles, struct SNP_Calling_Parameters * parameters, char * chro, int block_start )
{

	int bucket;
	KeyValuePair * cursor;

	for(bucket=0; bucket< merge_table -> numOfBuckets; bucket++) {
		cursor = merge_table -> bucketArray[bucket];

		while (1) {
			if (!cursor) break;
			int base_N_qual = cursor->value - NULL;
			int POI_pos = cursor->key - NULL - 100;
			int base1 = (base_N_qual>>8) & 0xff;
			int qual1 = (base_N_qual&0xff)-1;

			int j;
			unsigned int supporting_bases = 0;

			for(j=0; j<4; j++) {
				supporting_bases += snp_voting_piles[POI_pos*4+j];
			}

			if(supporting_bases < parameters->max_supporting_read_number) {
				if(qual1 <  (parameters -> is_phred_64?'@':'!')+parameters->min_phred_score){	// low quality bases
					// the low-seq-quality bases are not inclued in voting.
				} else snp_voting_piles[POI_pos*4+base1]+=1;
			}else{
				if(parameters -> supporting_read_excessed_reported < 100){
					parameters -> supporting_read_excessed_reported++;
					SUBREADprintf("Warning: the depth exceeded the limit of %d at %s : %d\n", parameters->max_supporting_read_number, chro ,  block_start + POI_pos);
					if(parameters -> supporting_read_excessed_reported == 100)
						SUBREADprintf("Too many warnings.\nNo more warning messages are reported.\n");
				}
			}

			cursor = cursor->next;
		}
	}
}

int read_tmp_block(struct SNP_Calling_Parameters * parameters, FILE * tmp_fp, char ** SNP_bitmap_recorder, unsigned int * snp_voting_piles, int block_no, unsigned int reference_len, char * referenced_genome, char * chro_name, unsigned int offset)
{
	int last_read_id=-1,i;
	HashTable * merge_table = HashTableCreate(1000);
	unsigned long long int OVERLAPPED_BASES=0;
	unsigned long long int OVER_MISMA_BASES=0;
	unsigned long long int ALL_BASES=0;

	while(!feof(tmp_fp))
	{
		int type_char = fgetc(tmp_fp), rlen=-1;
		if(type_char == EOF ) break;

		fseek(tmp_fp, -1 , SEEK_CUR);
		
		if(type_char == 200)	// SNP
		{
			VCF_temp_read_t SNP_rec;

			rlen = fread(&SNP_rec, sizeof(SNP_rec),1 , tmp_fp);
			if(rlen < 1){
				SUBREADputs("ERROR: the temporary file is broken.");
				return -1;
			}
			if(!(*SNP_bitmap_recorder))
			{
				(*SNP_bitmap_recorder)=malloc((reference_len/8)+200);
				memset((*SNP_bitmap_recorder), 0 , (reference_len/8)+200);
			}

			if(SNP_rec.pos >  block_no * BASE_BLOCK_LENGTH && SNP_rec.pos <= block_no * BASE_BLOCK_LENGTH + reference_len )
   				mask_snp_bitmap((*SNP_bitmap_recorder), SNP_rec.pos - block_no * BASE_BLOCK_LENGTH - 1);
			parameters -> known_SNPs_number ++;
			//printf("SNPat: %u\n", SNP_rec.pos);
		} 
		else if(type_char == 100)	// read
		{
			
			base_block_temp_read_t read_rec;
			unsigned short read_len;
			unsigned int first_base_pos;
			char read[MAX_READ_LENGTH];
			char qual[MAX_READ_LENGTH];

			int rlen = fread(&read_rec, sizeof(read_rec), 1, tmp_fp);
			if(rlen < 1){
				SUBREADputs("ERROR: the temporary file is broken.");
				return -1;
			}

			rlen = fread(&read_len, sizeof(short), 1, tmp_fp);
			if(rlen < 1){
				SUBREADputs("ERROR: the temporary file is broken.");
				return -1;
			}

			rlen = fread(read, sizeof(char), read_len, tmp_fp);
			if(rlen < read_len){
				SUBREADputs("ERROR: the temporary file is broken.");
				return -1;
			}
			
			rlen = fread(qual, sizeof(char), read_len, tmp_fp);
			if(rlen < read_len){
				SUBREADputs("ERROR: the temporary file is broken.");
				return -1;
			} 
			first_base_pos = read_rec.pos - block_no * BASE_BLOCK_LENGTH;
			parameters->is_paired_end_data = read_rec.flags & 1;

			//SUBREADprintf("Loading bases at %u ; len=%d\n", first_base_pos, read_len);

			if(first_base_pos + read_len -1> reference_len || first_base_pos<1)
			{
				SUBREADprintf("WARNING: read length %u+%d out of boundary: %u at the %d-th block.\n", first_base_pos, read_len, reference_len, block_no);
				continue;
			}

			if(parameters->is_paired_end_data){
				if( (last_read_id >>1) != (read_rec.read_number>>1) && last_read_id>=0 && merge_table -> numOfElements > 0)
				{
					put_hash_to_pile(merge_table, snp_voting_piles, parameters,chro_name, offset);
					HashTableDestroy(merge_table);
					merge_table = HashTableCreate(1000);
				}
			}

			last_read_id = read_rec.read_number;

			for(i=0;i<read_len;i++)
			{
				char base_int = -1;

				char true_value = referenced_genome[i + first_base_pos -1];
				//int  true_value_int =  (true_value=='A'?0:(true_value=='C'?1:(true_value=='G'?2:3)));

				switch(read[i])
				{
					case 'A':
						base_int = 0;
						break;
					case 'C':
						base_int = 1;
						break;
					case 'G':
						base_int = 2;
						break;
					case 'T':
						base_int = 3;
				}
				if(base_int < 0) continue;

				//if(qual[i] < '!'+parameters->min_phred_score)continue;
				if(true_value=='N' || true_value=='.') continue;

				if(i+first_base_pos > reference_len  || i+first_base_pos<1)
				{
					SUBREADprintf("Warning: read out of boundary: %u >= %u\n", i+first_base_pos, reference_len);
					break;
				}

				if(parameters->is_paired_end_data){
					// to remove overlapped bases between the two ends of the reads.
					unsigned int POI_pos = 100 + first_base_pos+i-1;
					unsigned int original_value = HashTableGet(merge_table,NULL+POI_pos)-NULL;
					if(0 == original_value)// no value has been put to this location
					{
						int base_N_qual = base_int;
						base_N_qual = base_N_qual << 8 | (1+qual[i]);
						HashTablePut(merge_table, NULL+POI_pos, NULL+base_N_qual);
						ALL_BASES++;
					}
					else
					{
						// the value has been in the merge_table.
						ALL_BASES++;
						OVERLAPPED_BASES++;
						int old_qual = (original_value & 0xff) - 1;
						int old_base = (original_value >> 8) & 0xff;
						if(old_base != base_int)OVER_MISMA_BASES++;
						if(old_qual < qual[i])
						{
							int base_N_qual = base_int;
							base_N_qual = base_N_qual << 8 | (1+qual[i]);
							HashTablePut(merge_table, NULL+POI_pos, NULL+base_N_qual);
						}	

					}
				}
				else
				{
					// if it is single-end, just put the read into the piles.
					int j;
					unsigned int supporting_bases = 0;

					for(j=0; j<4; j++)
					{
						supporting_bases += snp_voting_piles[(first_base_pos+i-1)*4+j];
					}

					if(supporting_bases < parameters->max_supporting_read_number)
					{
						if(qual[i] <  (parameters -> is_phred_64?'@':'!')+parameters->min_phred_score) 	// low quality bases
						{
							// the low-seq-quality bases are not inclued in voting.
						}
						else
						{
							snp_voting_piles[(first_base_pos+i-1)*4+base_int]+=1;
						}
					}else{
						if(parameters -> supporting_read_excessed_reported < 100){
							parameters -> supporting_read_excessed_reported++;
							SUBREADprintf("Warning: the depth exceeded the limit of %d at %s : %d\n", parameters->max_supporting_read_number, chro_name ,  offset + first_base_pos+i);
							if(parameters -> supporting_read_excessed_reported == 100)
								SUBREADprintf("Too many warnings.\nNo more warning messages are reported.\n");
						}
					}
				}
			}
		}
		else assert(0);
	}

	if(parameters->is_paired_end_data && merge_table -> numOfElements > 0)
		put_hash_to_pile(merge_table, snp_voting_piles, parameters, chro_name, offset);
	
	HashTableDestroy(merge_table);

	return 0;
}

void set_sample_MM(unsigned int * tumor_M, unsigned int *tumor_MM , unsigned int *snp_voting_piles, int i, int true_i){
	int j;

	for(j=0;j<4;j++)
	{
		if(j==true_i)(*tumor_M)+=snp_voting_piles[4*i + j];
		else	     (*tumor_MM)+=snp_voting_piles[4*i + j];
	}
}

void fishers_test_on_POI(struct SNP_Calling_Parameters * parameters, float * snp_fisher_VS, unsigned int *  snp_voting_piles,unsigned int * snp_BGC_piles, char * referenced_genome, unsigned int  reference_len, float * tumor_P)
{
	int i;
	for(i=0; i<reference_len; i++)
	{

		if(tumor_P[i]<0) continue;
		char true_value = referenced_genome[i];
		int  true_value_int =  (true_value=='A'?0:(true_value=='C'?1:(true_value=='G'?2:3)));
		unsigned int  tumor_M = 0, tumor_MM = 0, normal_M = 0, normal_MM = 0;

		set_sample_MM(&tumor_M, &tumor_MM , snp_voting_piles, i, true_value_int);
		if(tumor_MM)
		{
			set_sample_MM(&normal_M, &normal_MM , snp_BGC_piles, i, true_value_int);
			//assert(tumor_M == normal_M);
			//assert(tumor_MM == normal_MM);

			snp_fisher_VS[i]=fisher_exact_test(tumor_MM, normal_MM, tumor_M, normal_M);
		}
	}
}

void fishers_test_on_block(struct SNP_Calling_Parameters * parameters, float * snp_fisher_raw, unsigned int * snp_voting_piles, char * referenced_genome, unsigned int  reference_len, double multiplex_base, char * SNP_bitmap_recorder, unsigned short * snp_filter_background_matched, unsigned short * snp_filter_background_unmatched, int all_result_needed, char * chro_name, unsigned int block_start) {

		int i, j; // all_SNPs_at_pois=0;
		int POIbase_MM=0, ALLbase_MM=0, POIbase_MAT=0, ALLbase_MAT=0;
		long long int reference_len_long = reference_len;
		/**          | POI  | All_Window (inc. POI)
		 * ----------+------+-------
		 * #mismatch |  a   |  b
		 * #matched  |  c   |  d
		 **/
		for(i= - parameters -> fisher_exact_testlen;i<reference_len_long; i++)
		{
			POIbase_MM=0; POIbase_MAT=0;

			if(i>=0) {
				char true_value = referenced_genome[i];
				int  true_value_int =  (true_value=='A'?0:(true_value=='C'?1:(true_value=='G'?2:3)));

				for(j=0;j<4;j++)
				{
					if(j == true_value_int)
						POIbase_MAT  =  snp_voting_piles[ i * 4 + j ] ;
					else
						POIbase_MM +=  snp_voting_piles[ i * 4 + j ] ;
				}

			}

			if(i+parameters -> fisher_exact_testlen < reference_len_long) {
				int mm=0, pm=0;
				char true_value = referenced_genome[i + parameters -> fisher_exact_testlen ];
				int  true_value_int = (true_value=='A'?0:(true_value=='C'?1:(true_value=='G'?2:3)));

				for(j=0;j<4;j++)
				{
					if(j == true_value_int)
						pm = snp_voting_piles[ (i+parameters -> fisher_exact_testlen ) *4 + j ] ;
					else
						mm += snp_voting_piles[ (i+parameters -> fisher_exact_testlen ) *4 + j ] ;
				}


				unsigned int chro_pos = block_start + i + parameters -> fisher_exact_testlen;
				if(0 && chro_pos < 16084550+12 && chro_pos > 16084550-12){
					SUBREADprintf("ADD-BASE-AT-RIGHT : %s : %u : KNOWN-SNP : %d  a b c d : %d %d %d %d  A 4 > C : %d\n" , chro_name, chro_pos,  is_snp_bitmap(SNP_bitmap_recorder,i + parameters -> fisher_exact_testlen), POIbase_MM , ALLbase_MM, POIbase_MAT, ALLbase_MAT, (POIbase_MM *4 >= POIbase_MAT));
				}

				if((!SNP_bitmap_recorder)||(!is_snp_bitmap(SNP_bitmap_recorder, i + parameters -> fisher_exact_testlen))) {
					ALLbase_MAT += pm;
					ALLbase_MM += mm;
				}
			}

			// test the middle base
			// a : mismatched at POI ; c : matched at POI
			if(i>=0 && POIbase_MM > 0){
				int Known_SNP_here = (SNP_bitmap_recorder!=NULL) && is_snp_bitmap(SNP_bitmap_recorder,i);
				//if(Known_SNP_here)
				//		SUBREADprintf("SNP known at %s:%d\n" , chro_name,  block_start + i);
				//all_SNPs_at_pois += (Known_SNP_here)?1:0;

				float observed_coverage = (ALLbase_MM + ALLbase_MAT) *1./(1. + 2 * parameters -> fisher_exact_testlen) ;
				double p_cutoff = pow(10, -(observed_coverage/multiplex_base));
				p_cutoff = min(parameters -> cutoff_upper_bound, p_cutoff);
				p_cutoff = max(1E-323, p_cutoff);

				int flanking_unmatched = ALLbase_MM;
				int flanking_matched = ALLbase_MAT;
				if(SNP_bitmap_recorder == NULL || ! Known_SNP_here){
					flanking_unmatched -= POIbase_MM;
					flanking_matched -= POIbase_MAT;
				}


				unsigned int chro_pos = block_start + i;
				if(0 && strcmp("8", chro_name)==0 && chro_pos < 144762490 + 12 && chro_pos > 144762490-12){
					SUBREADprintf("KNOWN-SNP-AT-PoI at %s : %u : KNOWN-SNP : %d  a b c d : %d %d %d %d  A 4 > C : %d\n" , chro_name, chro_pos,  is_snp_bitmap(SNP_bitmap_recorder,i), POIbase_MM, ALLbase_MM , POIbase_MAT, ALLbase_MAT, (POIbase_MM *4 >= POIbase_MAT));
				}

				float p_middle = fisher_exact_test(POIbase_MM, flanking_unmatched, POIbase_MAT, flanking_matched);

				if(0 && chro_pos < 16084550+12 && chro_pos > 16084550-12){
					SUBREADprintf("TEST: %s : %u  a,b,c,d=%d %d %d %d; FU=%d FM=%d; Goahead=%d; Tailleft=%d; p=%G; p-cut=%G\n    KNOWN-SNP : %d    A 4 > C : %d\n", chro_name, chro_pos ,POIbase_MM, ALLbase_MM,POIbase_MAT, ALLbase_MAT , flanking_unmatched, flanking_matched, 0, 0, p_middle, p_cutoff,  is_snp_bitmap(SNP_bitmap_recorder,i), (POIbase_MM *4 >= POIbase_MAT));
				}

				//#warning " ===================== If a known SNP is at POI, the Fisher's p-value is only compared with the 'upper bound' of the p-value. ====================="
				//#warning " ===================== '1 &&' is the switch ========================"
				if(all_result_needed || ( 0 &&  Known_SNP_here && p_middle < parameters -> cutoff_upper_bound ) ||  ( p_middle < p_cutoff && flanking_matched*20>(flanking_matched+ flanking_unmatched )*16)) 
					snp_fisher_raw [i] = p_middle;
				else	snp_fisher_raw [i] = -999;


				if(flanking_unmatched <0) 
					SUBREADprintf("ERROR_AB: A=%d, B=%d, C=%d, D=%d, flanking_unmatched=%d\n", POIbase_MM, ALLbase_MM , POIbase_MAT,  ALLbase_MAT,flanking_unmatched);
				assert(flanking_unmatched>=0);
				assert(flanking_matched>=0);
				//if(strcmp(chro_name, "chr20")==0 && block_no  == 3 && i == 57566366-1-3*15000000)

				if(snp_filter_background_unmatched){
					snp_filter_background_unmatched[i] = flanking_unmatched;
					snp_filter_background_matched[i] = flanking_matched;
				}
				fisher_test_size ++;
			}else if(i>=0 && all_result_needed) snp_fisher_raw[i]=1.1;

			if(i >= parameters -> fisher_exact_testlen) {
				int mm=0, pm=0;
				char true_value = referenced_genome[i - parameters -> fisher_exact_testlen];
				int  true_value_int =  (true_value=='A'?0:(true_value=='C'?1:(true_value=='G'?2:3)));
				for(j=0;j<4;j++)
				{
					if(j == true_value_int)
						pm = snp_voting_piles[ (i-parameters -> fisher_exact_testlen) *4 + j ] ;
					else
						mm += snp_voting_piles[ (i-parameters -> fisher_exact_testlen) *4 + j ] ;
				}
				if((!SNP_bitmap_recorder)||!is_snp_bitmap(SNP_bitmap_recorder, i - parameters -> fisher_exact_testlen)) {
					ALLbase_MAT-=pm;
					ALLbase_MM-=mm;
				}
			}
		}
//		SUBREADprintf("BLOCK_FINISH; SNP found = %d\n", all_SNPs_at_pois);
}

int process_snp_votes(FILE *out_fp, unsigned int offset , unsigned int reference_len, char * referenced_genome, char * chro_name , char * temp_prefix, struct SNP_Calling_Parameters * parameters)
{
	int block_no = (offset -1) / BASE_BLOCK_LENGTH, i, disk_is_full = 0;
	char temp_file_name[MAX_FILE_NAME_LENGTH];
	FILE *tmp_fp;
	unsigned int * snp_voting_piles, *snp_BGC_piles = NULL;	// offset * 4 + "A/C/G/T"[0,1,2,3]
	char * SNP_bitmap_recorder = NULL;
	double multiplex_base = parameters -> fisher_normalisation_target / parameters -> cutoff_multiplex;

	float * snp_fisher_raw, *pcutoff_list, *snp_fisher_BGC=NULL, *snp_fisher_VS=NULL;
	unsigned short * snp_filter_background_unmatched;
	unsigned short * snp_filter_background_matched;

	long long int reference_len_long = reference_len;
	char * sprint_line;

	sprintf(temp_file_name , "%s%s-%04u.bin", temp_prefix, chro_name, block_no);
	tmp_fp = f_subr_open(temp_file_name, "rb");

	// if no temp file is here, do nothing 
	if(!tmp_fp){
		parameters->empty_blocks ++;
		return 0;
	}


	sprint_line = malloc(1000);
	sprint_line[0]=0;
	LOG_10_2 = log(2) / log(10);


	//printf("I %s %d\n", chro_name, offset);

	snp_voting_piles = (unsigned int *)SUBREAD_malloc(sizeof(unsigned int) * reference_len*4); 
	snp_filter_background_matched = (unsigned short *)SUBREAD_malloc(sizeof(unsigned short) * reference_len*4);
	snp_filter_background_unmatched = (unsigned short *)SUBREAD_malloc(sizeof(unsigned short) * reference_len*4);

	pcutoff_list = (float *)SUBREAD_malloc(sizeof(float)  * reference_len);
	snp_fisher_raw = (float *)SUBREAD_malloc(sizeof(float)  * reference_len);
	if((!snp_fisher_raw))
	{
		fatal_memory_size();
		return -1;
	}

	memset(snp_voting_piles,0 ,sizeof(unsigned int) * reference_len*4);

	for(i=0; i<reference_len; i++)
	{
		snp_fisher_raw[i]=-1.;
		pcutoff_list[i]=-1.;
	}

	int read_is_error = read_tmp_block(parameters, tmp_fp,&SNP_bitmap_recorder,snp_voting_piles,block_no, reference_len, referenced_genome, chro_name, offset -  reference_len);

	fclose(tmp_fp);
	if (parameters -> delete_piles)
		unlink(temp_file_name);


	if(parameters -> background_input_file[0])
	{
		snp_BGC_piles = (unsigned int *)SUBREAD_malloc(sizeof(unsigned int) * reference_len*4);
		snp_fisher_BGC = (float *)SUBREAD_malloc(sizeof(float)  * reference_len);
		snp_fisher_VS = (float *)SUBREAD_malloc(sizeof(float)  * reference_len);
		memset(snp_BGC_piles, 0, sizeof(int) *4 * reference_len);
		memset(snp_fisher_BGC, 0, sizeof(float)  * reference_len);
		memset(snp_fisher_VS, 0, sizeof(float)  * reference_len);

		sprintf(temp_file_name , "%sBGC-%s-%04u.bin", temp_prefix, chro_name, block_no);
		tmp_fp = f_subr_open(temp_file_name, "rb");
		if(tmp_fp)
		{
			int min_phred_score_raw = parameters -> min_phred_score;
			parameters -> min_phred_score -= 3;
			read_tmp_block(parameters, tmp_fp,NULL , snp_BGC_piles,block_no, reference_len, referenced_genome, chro_name, offset);
			parameters -> min_phred_score  = min_phred_score_raw;
		}
		fclose(tmp_fp);
		unlink(temp_file_name);

	}

		// Finding SNPs from the finished voting table

	char * base_is_reliable = NULL;


	if(parameters -> fisher_exact_testlen || parameters -> test_two_strands || parameters->neighbour_abundance_test || 1)
	{
		base_is_reliable = (char *)malloc(sizeof(char) * reference_len);
		memset(base_is_reliable,1 , sizeof(char) * reference_len);
	}


	if(parameters->neighbour_abundance_test)
	{
		int window_size = 7;
		int window_no = 0;
		for(i= -window_size;i<reference_len_long; i++)
		{
			int j;

			if(i<reference_len_long-window_size)
				for(j=0;j<4;j++)
				{
					window_no+=snp_voting_piles[(window_size+i)*4+j];
				}

			if(i>=window_size)
				for(j=0;j<4;j++)
				{
					window_no-=snp_voting_piles[(i-window_size)*4+j];
				}


			
		}
		
	}

	if(parameters -> test_two_strands)
	{
		for(i= 0;i<reference_len_long; i++)
		{
			int j;
			int pos_match=0, pos_unmatch=0, neg_match=0, neg_unmatch=0;

			char true_value = referenced_genome[i + parameters -> fisher_exact_testlen];
			int  true_value_int =  (true_value=='A'?0:(true_value=='C'?1:(true_value=='G'?2:3)));

			for(j=0;j<4;j++)
			{
				if(j==true_value_int)
				{
					pos_match += snp_voting_piles[i*4+j];
				}
				else
				{
					pos_unmatch += snp_voting_piles[i*4+j];
				}

			}

			float p_diff = fisher_exact_test(pos_match, neg_match, pos_unmatch, neg_unmatch);
			if(p_diff<0.15)
				base_is_reliable[i]=0;
		}
	}



	if(parameters -> neighbour_filter_testlen )
	{
		int j;
		int a=0, b=0, c=0, d=0;
		long long int reference_len_long = reference_len;
		/**    | This | All_Window 
		 * ----+------+-------
		 * #mm |  a   |  b
		 * #pm |  c   |  d
		 **/


		for(i= - parameters -> neighbour_filter_testlen;i<reference_len_long; i++)
		{
			a=0; c=0;
			for(j=0;j<4;j++)
			{
				if(i+parameters -> neighbour_filter_testlen < reference_len_long)
				{
					char true_value = referenced_genome[i + parameters -> neighbour_filter_testlen];
					int  true_value_int =  (true_value=='A'?0:(true_value=='C'?1:(true_value=='G'?2:3)));
					if(j == true_value_int)
						d += snp_voting_piles[ (i+parameters -> neighbour_filter_testlen) *4 + j ] ;
					else
						b += snp_voting_piles[ (i+parameters -> neighbour_filter_testlen) *4 + j ] ;
				}
				
				if(i>=0)
				{
					char true_value = referenced_genome[i];
					int  true_value_int =  (true_value=='A'?0:(true_value=='C'?1:(true_value=='G'?2:3)));

					if(j == true_value_int)
						c  =  snp_voting_piles[ i * 4 + j ] ;
					else
						a +=  snp_voting_piles[ i * 4 + j ] ;
				}

			}
		
			// test the middle base
			if(i>=0 && a > 0){

				float window_correct_rate = (d-c)*1./(d+b-a-c);
				if(window_correct_rate < parameters->neighbour_filter_rate)
				{
					base_is_reliable[i] = 0;
				}
				
			}

			for(j=0;j<4;j++)
				if(i >= parameters -> neighbour_filter_testlen)
				{
					char true_value = referenced_genome[i - parameters -> neighbour_filter_testlen];
					int  true_value_int =  (true_value=='A'?0:(true_value=='C'?1:(true_value=='G'?2:3)));
					if(j == true_value_int)
						d -= snp_voting_piles[ (i-parameters -> neighbour_filter_testlen) *4 + j ] ;
					else
						b -= snp_voting_piles[ (i-parameters -> neighbour_filter_testlen) *4 + j ];
				}
	
		}
	}


	

	if(parameters -> fisher_exact_testlen)
	{
		fishers_test_on_block(parameters, snp_fisher_raw, snp_voting_piles, referenced_genome, reference_len, multiplex_base, SNP_bitmap_recorder, snp_filter_background_matched, snp_filter_background_unmatched, 0, chro_name, BASE_BLOCK_LENGTH*block_no +1);
		if(parameters->background_input_file[0])
		{
			fishers_test_on_block(parameters, snp_fisher_BGC, snp_BGC_piles, referenced_genome, reference_len, multiplex_base, SNP_bitmap_recorder, NULL, NULL, 1, chro_name, BASE_BLOCK_LENGTH*block_no +1);
			fishers_test_on_POI(parameters, snp_fisher_VS, snp_voting_piles, snp_BGC_piles, referenced_genome, reference_len, snp_fisher_raw);
		}
	}


	for(i=0;i<reference_len; i++)
	{
		char true_value = referenced_genome[i];
		int true_value_int = (true_value=='A'?0:(true_value=='C'?1:(true_value=='G'?2:3)));
		int tested_int;
		int all_reads = 0;
		char base_list[10], supporting_list[55], snps=0;

		if(true_value == 'N') continue;

		for(tested_int=0; tested_int <4; tested_int ++)
			all_reads += snp_voting_piles[i*4+tested_int] ;

		//int all_reads_index = min(287-1, all_reads/5);
		if(all_reads<parameters->min_supporting_read_number || !all_reads) continue;
		//if(base_is_reliable && !(base_is_reliable[i])) continue;

		if(1)
		{
			unsigned int POI_MM = 0;
			base_list[0]=0;
			supporting_list[0]=0;

			for(tested_int=0; tested_int <4; tested_int ++)
			{
				if(tested_int != true_value_int)
				{
					int midNexcellent_sup = snp_voting_piles[i*4+tested_int];
					POI_MM += midNexcellent_sup;
					if(midNexcellent_sup >= parameters->min_alternative_read_number)
					{
						char int_buf[12];
						sprintf(int_buf, "%u", midNexcellent_sup);
						if(snps>0)
						{
							base_list[snps*2-1] = ',';
							base_list[snps*2] = tested_int==0?'A':(tested_int == 1?'C':(tested_int == 2?'G':'T'));
							base_list[snps*2+1] = 0;
							strcat(supporting_list,",");
							strcat(supporting_list,int_buf);
						}
						else
						{
							base_list[0] = tested_int==0?'A':(tested_int == 1?'C':(tested_int == 2?'G':'T'));
							base_list[1] = 0;
							strcpy(supporting_list, int_buf);
						}
						snps++;
							
					}
				}
			}

		//	SUBREADprintf("OUTTEXT: i=%d, all_reads=%d, Fisher-p=%f, SNPs=%d, POI-MM=%d\n", i, all_reads, snp_fisher_raw[i], snps, POI_MM);
			if(snps && POI_MM *1. / all_reads >= parameters->supporting_read_rate )
			{
				if(snp_fisher_raw[i] >= 0. || parameters -> fisher_exact_testlen < 1)
				{
					float Qvalue =  -1.0*log(max(1E-40,snp_fisher_raw[i]))/log(10);
					char BGC_Qvalue_str [120];
					BGC_Qvalue_str[0]=0;

					if(snp_fisher_BGC)
					{
						float BGC_Qvalue =  -1.0*log(max(1E-40,snp_fisher_BGC[i]))/log(10);
						float VS_Qvalue =  -1.0*log(max(1E-40,snp_fisher_VS[i]))/log(10);
						int BGC_all_reads = 0;
						int BGC_alt_reads = 0;

						for(tested_int=0; tested_int <4; tested_int ++)
						{
							BGC_all_reads += snp_BGC_piles[i*4+tested_int] ;
							if(tested_int != true_value_int)
								BGC_alt_reads += snp_BGC_piles[i*4+tested_int] ;
						}
						sprintf( BGC_Qvalue_str, ";CTRL_DP=%d;CTRL_MM=%d;CTRL_QV=%.4f;VS_QV=%.4f",BGC_all_reads, BGC_alt_reads, BGC_Qvalue,max(0,VS_Qvalue));
					}

					snprintf(sprint_line,999, "%s\t%u\t.\t%c\t%s\t%.4f\t.\tDP=%d;MMsum=%d;MM=%s;BGTOTAL=%d;BGMM=%d%s\n", chro_name, BASE_BLOCK_LENGTH*block_no +1 + i, true_value,base_list, Qvalue, all_reads, POI_MM, supporting_list , snp_filter_background_matched[i]+snp_filter_background_unmatched[i], snp_filter_background_unmatched[i], BGC_Qvalue_str);
					if(parameters->output_fp_lock)
						subread_lock_occupy(parameters->output_fp_lock);
					int sprint_line_len = strlen(sprint_line);
					int wlen = fwrite(sprint_line, 1, sprint_line_len,out_fp);
					if(wlen < sprint_line_len){
						disk_is_full=1;
						break;
					}
					parameters->reported_SNPs++;
					if(parameters->output_fp_lock)
						subread_lock_release(parameters->output_fp_lock);

				}
			}


			// check the indels right after the POI location
			char event_token[100];
			char event_token2[100];
			snprintf(event_token, 99, "%s\t%u", chro_name, BASE_BLOCK_LENGTH*block_no + 1 + i);
			long long int indel_no = HashTableGet(parameters -> cigar_event_table, event_token) - NULL;
			for(tested_int = 0; tested_int < indel_no; tested_int ++)
			{
				snprintf(event_token2, 99, "%s\t%u\t%d", chro_name, BASE_BLOCK_LENGTH*block_no + 1 + i, tested_int);
				long long int tmpv = HashTableGet(parameters -> cigar_event_table, event_token2) - NULL;
				long long int event_id = 0xffffff&(tmpv >> 8);

				int indels = (tmpv & 0xff) - 0x80;

				if(parameters->output_fp_lock)
					subread_lock_occupy(parameters->output_fp_lock);

				fprintf(out_fp, "%s\t.\t", event_token);

				fwrite(referenced_genome + i, max(0, indels) + 2, 1, out_fp);

				fwrite("\t", 1, 1, out_fp);
				if(indels<0)
				{
					char ** seq_tab = (char **)parameters -> cigar_event_table-> appendix1;
					char * ins_seq = seq_tab[event_id];
					fwrite(referenced_genome + i, 1, 1, out_fp);
					fwrite(ins_seq, 1, -indels, out_fp);
				}
				else
					fwrite(referenced_genome + i, 1, 1, out_fp);
				fwrite(referenced_genome + 1 + i + max(0,indels), 1, 1, out_fp);
				unsigned short * indel_sups = parameters -> cigar_event_table-> appendix2;
				int wlen = fprintf(out_fp, "\t1.0\t.\tINDEL;DP=%d;SR=%d\n",all_reads,indel_sups[event_id]);

				if(wlen < 10){
					disk_is_full=1;
					break;
				}
				parameters->reported_indels++;
				if(parameters->output_fp_lock)
					subread_lock_release(parameters->output_fp_lock);
			}
		}
	}


	if(base_is_reliable) free(base_is_reliable);
	if(SNP_bitmap_recorder) free(SNP_bitmap_recorder);
	free(snp_filter_background_matched);
	free(snp_filter_background_unmatched);
	free(snp_voting_piles);
	if(snp_BGC_piles){
		free(snp_BGC_piles);
		free(snp_fisher_BGC);
		free(snp_fisher_VS);
	}
	free(snp_fisher_raw);
	free(pcutoff_list);
	free(sprint_line);
	//SUBREADprintf("OVERLAPPED=%llu; MISMA=%llu; ALL_BASES=%llu\n",OVERLAPPED_BASES, OVER_MISMA_BASES, ALL_BASES);
	return read_is_error || disk_is_full;
}


int run_chromosome_search(FILE *in_fp, FILE * out_fp, char * chro_name , char * temp_prefix, chromosome_t * chromosomes , struct SNP_Calling_Parameters* parameters, int * task_no, int thread_no, int all_threads)
{
	int chro_no;
	unsigned int offset=0, all_offset = 0;
	
	unsigned int chro_len=0;
	char * referenced_base;

	referenced_base = (char *) SUBREAD_malloc(sizeof(char)* BASE_BLOCK_LENGTH + 10000);
	if(!referenced_base)
	{
		fatal_memory_size();
		return -1;
	}


	for(chro_no=0;chro_no < XOFFSET_TABLE_SIZE; chro_no++)
	{
		if(!chromosomes[chro_no].chromosome_name[0]) break;
		if(strcmp(chro_name , chromosomes[chro_no].chromosome_name)==0)
		{
			chro_len = chromosomes[chro_no].known_length;
			break;
		}
	}


	//if (all_threads>1)
	//	SUBREADprintf("Thread #%d is processing chromosome %s in FASTA file; expected length is %u.\n", thread_no, chro_name, chro_len);
	if(chro_len)
	{
		//SUBREADprintf("Processing chromosome %s in FASTA file; expected length is %u.\n", chro_name, chro_len);
		while( 1 || all_offset <= chro_len - 1)
		{
			char nc = fgetc(in_fp);

			if(nc == ' ' || nc == '\r' || nc == '\n') continue;
			if(nc == '>')
				fseek(in_fp, -1, SEEK_CUR);

			if(nc != EOF && nc != '>')
			{
				nc = toupper(nc);
				referenced_base[offset++] = (nc=='A' || nc == 'G' || nc == 'C')?nc:'T';
				all_offset ++;
			}

			if((nc == '>'||nc == EOF) && all_offset < chro_len)
				SUBREADprintf("WARNING: Chromosome is shorter than expected: %s\n", chro_name);

			if(offset == BASE_BLOCK_LENGTH || nc == EOF || nc == '>')
			{
				if(nc!=EOF && nc!='>')
				{
					long long int back_seek = ftello(in_fp);
					int xk2;
					for(xk2=0; xk2<10000; xk2++)
					{
						char ncx = fgetc(in_fp);
						if(ncx == EOF || ncx == '>') break;
						if(ncx == '\r' || ncx=='\n') continue;
						referenced_base[offset+xk2] = (ncx=='A' || ncx == 'G' || ncx == 'C')?ncx:'T';
					}

					fseeko(in_fp, back_seek, SEEK_SET);
				}


				if( (*task_no) % all_threads == thread_no && all_offset <= chro_len)
				{
					//#warning "=== ONLY TEST ONE BLOCK   , USE 'if(1)' IN RELEASE          ==="
					//if(strcmp(chro_name,"chr7")==0 && all_offset == 60000000){
					if(1){
						parameters -> disk_is_full |= process_snp_votes(out_fp, all_offset, offset, referenced_base, chro_name , temp_prefix, parameters);
						print_in_box(89,0,0,"processed block %c[36m%s@%d%c[0m by thread %d/%d [block number=%d/%d]", CHAR_ESC, chro_name, all_offset, CHAR_ESC , thread_no+1, all_threads, 1+(*task_no)-parameters->empty_blocks, parameters->all_blocks);
					}
					if(parameters -> disk_is_full)break;
				}
				else if((*task_no) % all_threads == thread_no)
				{
					print_in_box(80,0,0,"Ignored in: %s@%d by thr %d/%d [tid=%d]\n", chro_name, all_offset, thread_no, all_threads, *task_no);
	//				SUBREADprintf("LEN %u > %u\n", all_offset, chro_len);
				}
				offset = 0;
				(*task_no)++;
			}

			if(nc == EOF || nc == '>'){
				if(nc == '>') fseek(in_fp, -1, SEEK_CUR);
				break;
			}

		}
		//printf("Current indels=%ld\n", parameters->cigar_event_table->numOfElements);
	}
	else
	{
		if(thread_no==0)
			print_in_box(89,0,0,"WARNING chromosome %c[36m%s%c[0m is not in the SAM header.", CHAR_ESC, chro_name, CHAR_ESC);
		while(!feof(in_fp))
		{
			long long int fp0 = ftello(in_fp);
			
			char * ret = fgets(referenced_base, 2000, in_fp);
			if(!ret) break;

			if(referenced_base[0]=='>')
			{
				fseeko(in_fp, fp0, SEEK_SET);
				break;
			}
		}
	}

	free(referenced_base);

	return 0;
}



// This scan is driven by reading the FASTA file.
// While reading the FASTA file, if you see a new sequence, then find the base blocks for that sequence.
int parse_read_lists(char * in_FASTA_file, FILE * out_fp, char * temp_prefix, chromosome_t * chromosomes, struct SNP_Calling_Parameters * parameters, int all_threads, int thread_no)
{
	char line_buffer [3000];
	int task_no = 0;

	FILE *fp = f_subr_open(in_FASTA_file,"r");

	long long int FASTA_size = ftello(fp);

	if(!fp)
	{
		SUBREADprintf("Referenced Genome not found or inaccessible: '%s'.\n", in_FASTA_file);
		return -1;
	}
	
	while(!feof(fp))
	{
		int linelen = read_line(2999, fp, line_buffer, 0);
		
		if(line_buffer[0] == '>')
		{
			char chro_name [MAX_CHROMOSOME_NAME_LEN];
			int i;
			for(i=0; i< linelen -1; i++)
			{
				if(line_buffer[i+1] == ' ' || line_buffer[i+1] == '|' || line_buffer[i+1] == '\t')
					break;
				else chro_name[i] = line_buffer[min(MAX_CHROMOSOME_NAME_LEN-1,i+1)];
			}
			chro_name[min(MAX_CHROMOSOME_NAME_LEN-1,i)]=0;
	//		if(strcmp(chro_name,"chr20")==0)
			if(run_chromosome_search(fp, out_fp, chro_name , temp_prefix, chromosomes, parameters, &task_no, thread_no, all_threads)) return -1;

			float finished_rate = ftello(fp)*1./FASTA_size;
			if(snp_progress_report_callback)
				snp_progress_report_callback(40,10,3000+finished_rate*7000);
		}
	}
	fclose(fp);

	return 0;
}

struct parse_read_thread_context
{
	int thread_id;
	int all_threads;

	struct SNP_Calling_Parameters* parameters;
	char * temp_prefix;
	char * in_FASTA_file;
	chromosome_t * chromosomes;

	pthread_spinlock_t init_lock;
	FILE * out_fp;	
};

void *parse_read_lists_wrapper(void * vparam)
{
	struct parse_read_thread_context *param = (struct parse_read_thread_context *)vparam;
	int my_thread_id = param -> thread_id;
	pthread_spin_unlock(&param -> init_lock);


	parse_read_lists(param -> in_FASTA_file, param -> out_fp, param -> temp_prefix, param -> chromosomes, param -> parameters , param -> all_threads, my_thread_id);

	return NULL;
}

int parse_read_lists_maybe_threads(char * in_FASTA_file, char * out_BED_file, char * temp_prefix, chromosome_t * chromosomes, struct SNP_Calling_Parameters* parameters, int all_threads)
{

	FILE * out_fp = f_subr_open(out_BED_file,"w");
	int ret=0;
	if(!out_fp){
		SUBREADprintf("Cannot open the output file: '%s'\n", out_BED_file);
	}
	fputs("##fileformat=VCFv4.0\n",out_fp);
	fprintf(out_fp, "##exactSNP_Commandline=%s\n", parameters -> rebuilt_command_line);
	fputs("##comment=The QUAL values for the SNPs in this VCF file are calculated as min(40, - log_10 (p_value)), where p_value is from the Fisher's Exact Test. The QUAL values for the Indels in this VCF file are always 1.0.\n", out_fp);
	fputs("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n", out_fp);
	fputs("##INFO=<ID=BGMM,Number=1,Type=Integer,Description=\"Number of mismatched bases in the background (for SNP only)\">\n", out_fp);
	fputs("##INFO=<ID=BGTOTAL,Number=1,Type=Integer,Description=\"Total number of bases in the background (for SNP only)\">\n", out_fp);
	fputs("##INFO=<ID=MM,Number=1,Type=String,Description=\"Number of supporting reads for each alternative allele (for SNP only)\">\n", out_fp);
	fputs("##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n", out_fp);
	fputs("##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Number of supporting reads (for INDEL only)\">\n", out_fp);
	fputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", out_fp);
	if(all_threads < 2)
	{
		ret =  parse_read_lists(in_FASTA_file, out_fp, temp_prefix, chromosomes, parameters , all_threads, 0);
	}
	else
	{
		pthread_t runners [200];
		struct parse_read_thread_context param;
		param.parameters = parameters;		
		pthread_spin_init(&param.init_lock, PTHREAD_PROCESS_PRIVATE);
		pthread_spin_lock(&param.init_lock);
		param.out_fp = out_fp;
		param.in_FASTA_file = in_FASTA_file;
		param.chromosomes = chromosomes;
		param.all_threads = all_threads;
		param.temp_prefix = temp_prefix;
		parameters -> output_fp_lock = malloc(sizeof(*parameters -> output_fp_lock));
		pthread_spin_init(parameters -> output_fp_lock ,  PTHREAD_PROCESS_PRIVATE);

		int i;
		for(i=0; i<all_threads; i++)
		{
			param.thread_id = i;
			pthread_create(runners+i, NULL,  parse_read_lists_wrapper, &param);
			pthread_spin_lock(&param.init_lock);
		}

		
		for (i=0; i< all_threads; i++)
		{
			pthread_join(*(runners+i), NULL);
		}
		free((void *)parameters -> output_fp_lock);
	}
	//fprintf(out_fp, "## Fisher_Test_Size=%u\n",fisher_test_size);
	fclose(out_fp);
	if(parameters -> disk_is_full){
		unlink(out_BED_file);
		SUBREADputs("ERROR: cannot write into the output VCF file. Please check the disk space in the output directory.");
		ret = 1;
	}
	return ret;

}

int guess_quality_threshold(char * qual_file, float percentage)
{
	char buff[100];
	unsigned long long readns [100], all_readn = 0;
	FILE * qfp = f_subr_open(qual_file, "r");
	int lineno = 0;

	if(!qfp)
	{
		SUBREADprintf("Unable to open file: %s\n", qual_file);
		return 0;
	}

	memset(readns, 0, 100*sizeof(long long));

	while(!feof(qfp))
	{
		int pos=0;
		unsigned long long readn = 0;
		read_line(100, qfp, buff, 0);	

		while(buff[pos]!='\t')pos++;
		pos++;

		for(;buff[pos] && buff[pos]!='\n'; pos++)
		{
			readn *=10;
			readn += buff[pos]-'0';
		}
		readns[lineno]=readn;
		all_readn += readn;

		lineno ++;
		if(lineno>64)break;
	}

	fclose(qfp);
	long long read_threshold = all_readn * percentage;

	for(lineno=99; lineno>1; lineno--)
	{
		read_threshold -= readns[lineno]; 
		if(read_threshold<=0){
			SUBREADprintf("Phred score threshold has been shifted to %d\n", lineno);
			return lineno;
		}
	}

	return 0;
}

char * _EXSNP_SNP_delete_temp_prefix = NULL;
void EXSNP_SIGINT_hook(int param)
{
	#ifndef MAKE_STANDALONE
	if(param==0)
	{
	#endif
	
		int xk1, last_slash = -1;
		if(_EXSNP_SNP_delete_temp_prefix != NULL)
		{
			char del2[MAX_FILE_NAME_LENGTH], del_suffix[MAX_FILE_NAME_LENGTH], del_name[MAX_FILE_NAME_LENGTH];
			#ifdef MAKE_STANDALONE
			if(param)
				SUBREADprintf("\n\nReceived a terminal signal. The temporary files were removed.\n");
			#endif
			for(xk1=0; _EXSNP_SNP_delete_temp_prefix[xk1]; xk1++)
			{
				if(_EXSNP_SNP_delete_temp_prefix[xk1]=='/') last_slash = xk1;
				else if(_EXSNP_SNP_delete_temp_prefix[xk1]=='\\')
				{
					SUBREADprintf("The file name is unknown.\n");
					return;
				}
			}
			if(last_slash>=0)
			{
				memcpy(del2, _EXSNP_SNP_delete_temp_prefix, last_slash);
				del2[last_slash] = 0;
				strcpy(del_suffix , _EXSNP_SNP_delete_temp_prefix + last_slash + 1);
			}
			else
			{
				strcpy(del2,".");
				strcpy(del_suffix , _EXSNP_SNP_delete_temp_prefix);
			}
		
			if(strlen(del_suffix)>8)
			{
				DIR           *d;
				struct dirent *dir;

				d = opendir(del2);
				if (d)
				{
					while ((dir = readdir(d)) != NULL)
					{
						if(strstr(dir->d_name, del_suffix))
						{
							//printf("%s\n", dir->d_name);
							strcpy(del_name, del2);
							strcat(del_name, "/");
							strcat(del_name, dir->d_name);
							unlink(del_name);
						}
					}
					closedir(d);
				}
			}
				
		}

	#ifdef MAKE_STANDALONE
	if(param)
	exit(param);
	#else
	}
	#endif
}


int SNP_calling(char * in_SAM_file, char * out_BED_file, char * in_FASTA_file, char * temp_location, unsigned int known_read_count, int threads, struct SNP_Calling_Parameters* parameters)
{
	char temp_file_prefix[MAX_FILE_NAME_LENGTH];
	chromosome_t * known_chromosomes;
	unsigned int real_read_count=0;

	int i, fpos;

	signal (SIGTERM, EXSNP_SIGINT_hook);
	signal (SIGINT, EXSNP_SIGINT_hook);

	parameters -> output_fp_lock = NULL;

	fisher_test_size = 0;

	precalculated_factorial = (long double*)malloc(sizeof(long double)*PRECALCULATE_FACTORIAL * 2);
	for(i=0; i<PRECALCULATE_FACTORIAL * 2; i++)
		precalculated_factorial[i] = 0.; 
		

	known_chromosomes = (chromosome_t *) SUBREAD_malloc(sizeof(chromosome_t) * XOFFSET_TABLE_SIZE);
	if(!known_chromosomes)
	{
		fatal_memory_size();
		return -1;
	}

	known_chromosomes[0].chromosome_name[0]=0;

	if(snp_progress_report_callback)
		snp_progress_report_callback(40,10,3000);

	if(parameters->subread_index[0])
	{
		char table_fn[MAX_FILE_NAME_LENGTH+80];
		parameters -> subread_index_offsets = (gene_offset_t*)malloc(sizeof(gene_offset_t));
		load_offsets (parameters -> subread_index_offsets, parameters->subread_index);

		sprintf(table_fn,"%s.00.b.array", parameters->subread_index);
		parameters -> subread_index_array = (gene_value_index_t*) malloc(sizeof(gene_value_index_t));
		if(gvindex_load(parameters -> subread_index_array,table_fn)) return -1;
	}

	//Step 1:The SAM file is scanned to create a number of temp files "temp-snps-chrX-21-00000000-XXXXXX" and the related read positions/sequences are written into them (in 2-bit base coding). A read can contribute to two blocks if it crosses the border of the blocks. If there are indels in a read, each continuously mapped section is individually written into the temporary file. The quality scores are written companying the bases. The hierarchy of the data is: block -> read -> sections {bases, phred scores, CIGAR string, reported mapping quality}

	if(parameters->pile_file_name[0])
	{
		fpos=0;
		while(1)
		{
			int fpos0 = fpos;
			char one_fn [MAX_FILE_NAME_LENGTH];
			while(in_SAM_file[fpos]!=',' && in_SAM_file[fpos]!=0)
				fpos++;
			strncpy(one_fn, in_SAM_file+fpos0, fpos-fpos0);
			one_fn[fpos-fpos0]=0;

			if(get_known_chromosomes(one_fn, known_chromosomes)) return -1;

			if(!in_SAM_file[fpos]) break;
			fpos++;
		}
		strcpy(temp_file_prefix, parameters->pile_file_name);
	}	
	else
	{	
		parameters-> cigar_event_table = HashTableCreate(11011);
		HashTableSetHashFunction(parameters-> cigar_event_table ,HashTableStringHashFunction);
		HashTableSetDeallocationFunctions(parameters-> cigar_event_table, NULL, NULL);
		HashTableSetKeyComparisonFunction(parameters-> cigar_event_table, my_strcmp);

		char mac_rand[13];
		mac_or_rand_str(mac_rand);

		sprintf(temp_file_prefix, "%s/temp-snps-%06u-%s-", temp_location, getpid(), mac_rand);

		_EXSNP_SNP_delete_temp_prefix = temp_file_prefix;

		print_in_box(89,0,0,"Split %s file into %c[36m%s*%c[0m ..." , parameters -> is_BAM_file_input?"BAM":"SAM" , CHAR_ESC, temp_file_prefix, CHAR_ESC);

		fpos=0;
		while(1)
		{
			int fpos0 = fpos;
			char one_fn [MAX_FILE_NAME_LENGTH];
			while(in_SAM_file[fpos]!=',' && in_SAM_file[fpos]!=0)
				fpos++;
			strncpy(one_fn, in_SAM_file+fpos0, fpos-fpos0);
			one_fn[fpos-fpos0]=0;

			if(break_SAM_file(one_fn, parameters -> is_BAM_file_input, temp_file_prefix, &real_read_count, &parameters->all_blocks, known_chromosomes, 1, parameters -> bases_ignored_head_tail, parameters->subread_index_array, parameters->subread_index_offsets, &parameters -> all_mapped_bases, parameters-> cigar_event_table, parameters->known_SNP_vcf, NULL, 0,0, parameters -> use_soft_clipping_bases)) return -1;
			if(!in_SAM_file[fpos]) break;
			fpos++;
		}
		if(parameters -> background_input_file[0])
		{
			char temp_file_prefix2[MAX_FILE_NAME_LENGTH+80];
			sprintf(temp_file_prefix2, "%sBGC-", temp_file_prefix);
			if(break_SAM_file(parameters -> background_input_file, parameters -> is_BAM_file_input, temp_file_prefix2, NULL, NULL, known_chromosomes, 1, parameters -> bases_ignored_head_tail, parameters->subread_index_array, parameters->subread_index_offsets, NULL, NULL, NULL, NULL, 0,0, parameters -> use_soft_clipping_bases)) return -1;
		}


	}

	parameters -> real_read_count = real_read_count;

	char qfname[MAX_FILE_NAME_LENGTH+40];
	sprintf(qfname, "%s.qStatic", temp_file_prefix);
	parameters -> final_phred_score = 0;
	if (parameters -> delete_piles)
		unlink(qfname);

	//Step 2:Each base blocks are load from the FASTA file into memory, then each temp file is scanned to create the SNP voting table. The sections in the temp files are scanned to create the SNP voting table. The voting table is compared with the sequences in the FASTA file to determine if each base is a SNP. The result is written into the bed file immediatly.
	if(parse_read_lists_maybe_threads(in_FASTA_file, out_BED_file, temp_file_prefix, known_chromosomes, parameters, threads)) return -1;

	free(known_chromosomes);

	free(precalculated_factorial);

	if(parameters->subread_index_array)
	{
		gvindex_destory(parameters->subread_index_array);
		free(parameters->subread_index_array);
		free(parameters->subread_index_offsets);
	}

	destroy_cigar_event_table(parameters-> cigar_event_table);

	return 0;
}


void print_usage_snp(char * myname)
{

	SUBREADprintf("\nVersion %s\n\n", SUBREAD_VERSION);
	SUBREADputs("Usage:");
	SUBREADputs("");
	SUBREADputs("  ./exactSNP [options] -i input -g reference_genome -o output ");
	SUBREADputs("");
	SUBREADputs("Required arguments:");
	SUBREADputs("");
	SUBREADputs("  -i <file>  Specify name of an input file including read mapping results. The");
	SUBREADputs(" [-b if BAM] format of input file can be SAM or BAM (-b needs to be specified");
	SUBREADputs("             if a BAM file is provided).");
	SUBREADputs("");
	SUBREADputs("  -g <file>  Specify name of the file including all reference sequences. Only");
	SUBREADputs("             one single FASTA format file should be provided.");
	SUBREADputs("");
	SUBREADputs("  -o <file>  Specify name of the output file. This program outputs a VCF format");
	SUBREADputs("             file that includes discovered SNPs.");
	SUBREADputs("");
	SUBREADputs("Optional arguments:");
	SUBREADputs("");
	SUBREADputs("  -a <file>  Provide a set of annotated SNPs (e.g. SNPs included in the dbSNP");
	SUBREADputs("             database). The supplied file should be in VCF format (gzipped file");
	SUBREADputs("             is accepted). Providing known SNPs to the program should improve");
	SUBREADputs("             its SNP calling performance. Note that the provided SNPs may or");
	SUBREADputs("             may not be called.");
	SUBREADputs("");
	SUBREADputs("  -b         Indicate the input file provided via -i is in BAM format.");
	SUBREADputs("");
	SUBREADputs("  -Q <int>   Specify the q-value cutoff for SNP calling at sequencing depth of");
	SUBREADputs("             50X. 12 by default. The corresponding p-value cutoff is 10^(-1*Q).");
	SUBREADputs("             Note that this program automatically adjusts the q-value cutoff");
	SUBREADputs("             according to the sequencing depth at each chromosomal location.");
	SUBREADputs("");
	SUBREADputs("  -f <float> Specify the minimum fraction of mis-matched bases a SNP-containing");
	SUBREADputs("             location must have. Its value must between 0 and 1. 0 by default.");
	SUBREADputs("");
	SUBREADputs("  -n <int>   Specify the minimum number of mis-matched bases a SNP-containing");
	SUBREADputs("             location must have. 1 by default.");
	SUBREADputs("");
	SUBREADputs("  -r <int>   Specify the minimum number of mapped reads a SNP-containing");
	SUBREADputs("             location must have (ie. the minimum coverage). 1 by default.");
	SUBREADputs("");
	SUBREADputs("  -x <int>   Specify the maximum depth a SNP location is allowed to have.");
	SUBREADputs("             1,000,000 reads by default. This option is useful for removing PCR");
	SUBREADputs("             artefacts.");
	SUBREADputs("");
	SUBREADputs("  -s <int>   Specify the minimum base quality scores (Phred scores) read bases");
	SUBREADputs("             must satisfy to be used for SNP calling. 13 by default. Read bases");
	SUBREADputs("             with quality scores less than 13 will be excluded from the");
	SUBREADputs("             analysis.");
	SUBREADputs("");
	SUBREADputs("  -t <int>   Specify the number of bases trimmed off from each end of the read.");
	SUBREADputs("             3 by default.");
	SUBREADputs("");
	SUBREADputs("  -T <int>   Specify the number of threads. 1 by default.");
	SUBREADputs("");
	SUBREADputs("  -v         output version of the program.");
	SUBREADputs("");
	SUBREADputs("  -C <path>  Specify the path to save the temporary files.");
	SUBREADputs("");
	SUBREADputs("Example:");
	SUBREADputs("");
	SUBREADputs("  ./exactSNP -i my-alignment.sam -g mm10.fa -o my-SNPs.txt");
	SUBREADputs("");
}


static struct option snp_long_options[]={
	{"coverage-calc", no_argument ,0,'4'},
	{0,0,0,0}
};

#ifdef MAKE_STANDALONE
int main(int argc,char ** argv)
{
	snp_progress_report_callback = NULL;
#else
int main_snp_calling_test(int argc,char ** argv)
{
#endif

	int xk1 = 0;
	int ret = 0;
	for(xk1=0;xk1<1;xk1++){


	int c;
	char in_SAM_file[5000];
	char out_BED_file[MAX_FILE_NAME_LENGTH];
	char temp_path[MAX_FILE_NAME_LENGTH];
	char in_FASTA_file[MAX_FILE_NAME_LENGTH];
	int threads, optindex=0;
	int t=0, k;
	struct SNP_Calling_Parameters parameters;
	unsigned int read_count;

	in_SAM_file [0] = 0;
	in_FASTA_file [0] = 0;
	out_BED_file[0] = 0;
	temp_path[0] = 0;

	read_count = 0;
	threads = 1;

	optind = 0;
	opterr = 1;
	optopt = 63;


	memset(&parameters, 0, sizeof(struct SNP_Calling_Parameters));
	rebuild_command_line(&parameters . rebuilt_command_line, argc, argv);
	parameters.start_time = miltime();
	parameters.empty_blocks = 0;
	parameters.reported_SNPs = 0;
	parameters.reported_indels = 0;

	parameters.supporting_read_rate = 0.;
	parameters.min_supporting_read_number = 1;
	parameters.min_alternative_read_number = 1;
	parameters.max_supporting_read_number = 1000000;
	parameters.neighbour_filter_testlen = -1; 
	parameters.neighbour_filter_rate = 0.000000001;
	parameters.min_phred_score = 13;
	parameters.fisher_exact_testlen = 5;
	parameters.test_two_strands = 0;
	parameters.snp_interval = 0;
	parameters.excellent_phred_score = 39;

	parameters.bases_ignored_head_tail = 3;
	parameters.neighbour_abundance_test = 0;

	parameters.pile_file_name[0]=0;
	parameters.delete_piles = 1;

	parameters.subread_index[0]=0;
	parameters.subread_index_offsets=NULL;
	parameters.subread_index_array=NULL;
	parameters.is_phred_64 = 0;
	parameters.is_coverage_calculation = 0;
	parameters.is_BAM_file_input = 0;

	parameters.cutoff_multiplex = 12.f;
	parameters.cutoff_upper_bound = 5E-3f;
	parameters.all_mapped_bases = 0;
	parameters.fisher_normalisation_target = 50;
	parameters.known_SNP_vcf[0]=0;
	parameters.known_SNPs_number=0;
	parameters.background_input_file[0]=0;

	if(argc<2)
	{
		print_usage_snp(argv[0]);
		return 0;
	}
	while ((c = getopt_long (argc, argv, "S7:N:C:a:i:g:o:bQ:p:f:n:r:x:w:s:t:T:v4",snp_long_options, &optindex))!=-1)
	{
		switch (c)
		{
			case 'S':
				parameters.use_soft_clipping_bases = 1;
				break;
			case 'N':
				strcpy(parameters.background_input_file, optarg);
				break;
			case 'a':
				strcpy(parameters.known_SNP_vcf, optarg);
				break;
			case 'v':
				core_version_number("exactSNP");
				return 0;
			case 'b':
				parameters.is_BAM_file_input = 1;
				break;
			case 'Q':
				parameters.cutoff_multiplex = 1.0*atoi(optarg);
				break;
			case 'p':
				parameters.cutoff_upper_bound = atof(optarg);
				break;
			case 'n':
				parameters.min_alternative_read_number = atoi(optarg);
				break;
			case '2':	// UNUSED
				parameters.neighbour_abundance_test=0;
				break;
			case 't':
				parameters.bases_ignored_head_tail = atoi(optarg);
				break;
			case 'w':
				parameters.fisher_exact_testlen = atoi(optarg);
				break;
			case 's':
				parameters.min_phred_score = atoi(optarg);
				break;
			case '3':	// UNUSED
				k=strlen(optarg);
				for(t=0;t<k;t++)
					if(optarg[t]==',')
					{
						optarg[t]=0;
						break;
					}
				if(t==k)
					SUBREADprintf("Warning: the neighbour filtering parameter is unparseable. It should be like \"-N 5,0.5\".\n");
				else
				{
					parameters.neighbour_filter_testlen = atoi(optarg);
					parameters.neighbour_filter_rate = atof(optarg+t+1);
				}
				break;

			case 'g':
				strncpy(in_FASTA_file, optarg,MAX_FILE_NAME_LENGTH-1);
				break;

			case 'i':
				strncpy(in_SAM_file, optarg,MAX_FILE_NAME_LENGTH-1);
				break;

			case 'o':
				strncpy(out_BED_file, optarg,MAX_FILE_NAME_LENGTH-1);
				break;

			case 'T':
				threads = atoi(optarg);
				if(!threads)threads=1;
				break;

			case 'f':
				parameters.supporting_read_rate = atof(optarg);
				break;

			case 'x':
				parameters.max_supporting_read_number = atof(optarg);
				break;

			case 'C':
				strncpy(temp_path, optarg,MAX_FILE_NAME_LENGTH-1);	
				break;

			case 'r':
				parameters.min_supporting_read_number = atof(optarg);
				break;

			case '8':	// UNUSED
				read_count = atoi(optarg);
				break;

			case '7':	// UNUSED
				strcpy(parameters.pile_file_name, optarg);
				break;

			case '6':	// UNUSED
				strcpy(parameters.subread_index, optarg);
				break;

			case '5':	// UNUSED
				parameters.delete_piles = 0;
				break;				

			case '4':
				parameters.is_coverage_calculation = 1;
				break;				

			case '?':
			default:
				print_usage_snp(argv[0]);
				return -1;
		}
	}

	if(argc > optind){
		SUBREADprintf("Invalid parameter '%s'\n", argv[optind]);
		return -1;
	}

	if(out_BED_file[0]==0 || in_FASTA_file[0]==0 || (parameters.pile_file_name [0] == 0 && in_SAM_file[0]==0))
	{
		SUBREADprintf("The names of the input file, the output file and the reference sequence file must be specified using -i, -o and -g options.\n");
		return -1;
	}

	FILE * tfp = f_subr_open(out_BED_file,"w");
	if(!tfp)
	{
		SUBREADprintf("Cannot write the output file. Please check if the disk has space and you can create file in the current directory.\n");
		return -1;
	}
	fclose(tfp);

	tfp = f_subr_open(in_FASTA_file,"r");
	if(!tfp)
	{
		SUBREADprintf("Cannot open the reference sequence file: %s\n", in_FASTA_file); 
		return -1;
	}
	fclose(tfp);

	if(in_SAM_file[0])
	{
		tfp = f_subr_open(in_SAM_file,"r");
		if(!tfp)
		{
			SUBREADprintf("Cannot open the input file: %s\n", in_SAM_file); 
			return -1;
		}
		fclose(tfp);
	}



	SUBREADputs("");
	print_subread_logo();
	SUBREADputs("");


	print_in_box(80,1,1,"exactSNP setting");
	print_in_box(80,0,1,"");
	print_in_box(80,0,0,"                  Input file : %s (%s)", get_short_fname(in_SAM_file), parameters.is_BAM_file_input?"BAM":"SAM");
	print_in_box(80,0,0,"                 Output file : %s", get_short_fname(out_BED_file));
	print_in_box(80,0,0,"            Reference genome : %s", get_short_fname(in_FASTA_file));
	print_in_box(80,0,0,"                   Temp path : %s", temp_path);
	print_in_box(80,0,1,"");
	print_in_box(80,0,0,"                     Threads : %d", threads);
	print_in_box(80,0,0,"        Min supporting reads : %d", parameters.min_supporting_read_number);
	print_in_box(81,0,0,"Min pct. of supporting reads : %.1f%%%%", 100* parameters.supporting_read_rate);
	print_in_box(80,0,0,"      Min base quality score : %d", parameters.min_phred_score);
	print_in_box(80,0,0,"     Number of trimmed bases : %d", parameters.bases_ignored_head_tail);
	print_in_box(80,0,1,"");
	print_in_box(80,0,0,"              Q value cutoff : 10E-%.1f", parameters.cutoff_multiplex);
	print_in_box(80,0,0,"         P value upper bound : %.5f", parameters.cutoff_upper_bound);
	print_in_box(80,0,0,"       Flanking windows size : %d", parameters.fisher_exact_testlen);
	if(parameters.known_SNP_vcf[0])
		print_in_box(80,0,0,"       Known SNP annotations : %s", get_short_fname(parameters.known_SNP_vcf));
	
	print_in_box(80,0,1,"");
	print_in_box(80,2,1,"");
	SUBREADputs("");


	char tbuf[90];
	char_strftime(tbuf);



	print_in_box(80,1,1,"Running (%s)", tbuf);
	print_in_box(80,0,1,"");
	warning_file_type(in_SAM_file, parameters.is_BAM_file_input?FILE_TYPE_BAM:FILE_TYPE_SAM);
	warning_file_type(in_FASTA_file, FILE_TYPE_FASTA);
	warning_file_limit();
	int x1;

	if(temp_path[0]==0){
		// temp path = out_dir
		for(x1 = strlen(out_BED_file); x1 >= 0; x1--){
			if(out_BED_file[x1]=='/'){
				memcpy(temp_path, out_BED_file, x1);
				temp_path[x1]=0;
				break;
			}
		}
		if(temp_path[0]==0)strcpy(temp_path, "./");
	}
	ret = SNP_calling(in_SAM_file, out_BED_file, in_FASTA_file, temp_path, read_count, threads, &parameters);
	if(ret != -1)
	{
		print_in_box(80,0,1,"");
		print_in_box(89,0,1,"%c[36mCompleted successfully.%c[0m",CHAR_ESC, CHAR_ESC);
		print_in_box(80,0,1,"");
		print_in_box(80,2,1,"");
		SUBREADputs("");

		print_in_box(80,1,1,"Summary");
		print_in_box(80,0,1,"");
		if(parameters.known_SNPs_number)
			print_in_box(80,0,0,"                  Known SNPs : %u", parameters.known_SNPs_number);
		print_in_box(80,0,0,"             Processed reads : %u", parameters.real_read_count);
		print_in_box(80,0,0,"               Reported SNPs : %u", parameters.reported_SNPs);
		print_in_box(80,0,0,"             Reported indels : %u", parameters.reported_indels);
		print_in_box(80,0,1,"");
		print_in_box(80,0,0,"                Running time : %.1f minutes", (miltime() - parameters.start_time)/60);
		print_in_box(80,0,1,"");
		print_in_box(80,2,1,"");
		SUBREADputs("");
	}


	EXSNP_SIGINT_hook(0);
	}	// XK1 END
	return ret;
	
}
