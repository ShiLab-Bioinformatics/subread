#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "subread.h" 
#include "core.h" 
#include "HelperFunctions.h"
#include "sambam-file.h" 
#include "gene-algorithms.h" 
#include "input-files.h" 
#include "hashtable.h"

#define COVERAGE_MAX_INT 0x7ffffff0 
#define MAX_FRAGMENT_LENGTH 3000
unsigned long long all_counted;
typedef unsigned int coverage_bin_entry_t;
int is_BAM_input = 0;
int max_M = 10;
int paired_end = 0;
int ignore_bad_pair = 0;
int base_values = 0;
char input_file_name[MAX_FILE_NAME_LENGTH];
char output_file_name[MAX_FILE_NAME_LENGTH];
HashTable * cov_bin_table;


static struct option cov_calc_long_options[] =
{
	{"maxMOp",required_argument, 0, 'M'},
	{"primary",no_argument, 0, 0},
	{0, 0, 0, 0}
};


void calcCount_usage()
{
	SUBREADprintf("\ncoverageCount Version %s\n\n", SUBREAD_VERSION);
	SUBREADputs("  This program calculates the coverage of mapped reads at each location on");
	SUBREADputs("the reference genome. It generates a binary file for each chromosome by concate-");
	SUBREADputs("nating the coverage levels as 4-bytes integer numbers.");
	SUBREADputs("");
	SUBREADputs("Usage");
	SUBREADputs("");
	SUBREADputs("  ./coverageCount [options] -i <input_file> -o <output_prefix>");
	SUBREADputs("");
	SUBREADputs("Required arguments:");
	SUBREADputs("");
	SUBREADputs("  -i <string>  Name of input file in SAM or BAM format.");
	SUBREADputs("");
	SUBREADputs("  -o <string>  Prefix of the output files. Each output file contains Four-byte");
	SUBREADputs("               integer numbers");
	SUBREADputs("");
	SUBREADputs("Optional arguments:");
	SUBREADputs("");
	SUBREADputs("  -p           The input file contains paired-end reads.");
	SUBREADputs("");
	SUBREADputs("  --maxMOp <int> Maximum number of 'M' operations allowed in a CIGAR string.");
	SUBREADputs("               10 by default. Both 'X' and '=' are treated as 'M' and adjacent");
	SUBREADputs("               'M' operations are merged in the CIGAR string.");
	SUBREADputs("");
}

void add_chro(char *sam_h)
{
	char *chro_name = malloc(200);
	unsigned int chro_len = 0;

	char nch;
	int cur = 0, tabs=0, state = 0, txtcur = 0;

	chro_name[0]=0;
	while(1){
		nch = sam_h[cur++];
		if(!nch || nch == '\n')break;

		if(nch == '\t')
		{
			txtcur = 0;
			tabs ++;
			state = 0;
		}
		else if(nch == ':')
		{
			state ++;
		}
		else
		{
			if(state == 1 && tabs == 1)
			{
				chro_name[txtcur++]=nch;
				chro_name[txtcur]=0;
			}
			else if(state == 1 && tabs ==2)
				chro_len = chro_len*10 + (nch-'0');
		}
	}

	if(chro_name[0]==0 || chro_len<1)
	{
		SUBREADprintf("ERROR: incorrect SAM format: %s\n", sam_h);
	}

	void ** bin_entry = malloc(sizeof(void *)*2);
	coverage_bin_entry_t * chro_bin = calloc(sizeof(coverage_bin_entry_t) , chro_len * ( base_values?5:1 ));
	if(!chro_bin)
	{
		SUBREADprintf("ERROR: cannot allocate the memory block. You need at least 4GB of memory,\n");
	}
	bin_entry[0] = (void *)chro_bin;
	bin_entry[1] = (void *)(NULL + chro_len);
	HashTablePut(cov_bin_table, chro_name, bin_entry);
	SUBREADprintf("Added a new chromosome : %s [%u]\n", chro_name, chro_len);
}

void get_read_info(char * fl, char * chro, unsigned int * pos_1base , char * cigar, int *flags, int * tlen, char ** rtext){
	char * tmp_tok = NULL;

	char * tmp_res = strtok_r(fl, "\t", &tmp_tok);

	tmp_res = strtok_r(NULL, "\t", &tmp_tok);
	(*flags) = atoi(tmp_res);

	tmp_res = strtok_r(NULL, "\t", &tmp_tok); // chro
	strcpy(chro, tmp_res);

	tmp_res = strtok_r(NULL, "\t", &tmp_tok); // pos_0base
	(*pos_1base) = atoi(tmp_res);

	tmp_res = strtok_r(NULL, "\t", &tmp_tok); // qual
	tmp_res = strtok_r(NULL, "\t", &tmp_tok); // cigar
	strcpy(cigar, tmp_res);

	tmp_res = strtok_r(NULL, "\t", &tmp_tok); // mate chro
	tmp_res = strtok_r(NULL, "\t", &tmp_tok); // mate pos
	tmp_res = strtok_r(NULL, "\t", &tmp_tok); // TLEN 
	(*tlen) = atoi(tmp_res);

	if(rtext)
		*rtext = strtok_r(NULL, "\t", &tmp_tok);
}

int base_value_to_int(char b){
	if(b=='A') return 1;
	if(b=='C') return 2;
	if(b=='G') return 3;
	if(b=='T') return 4;
	return 0;
}

int covCalc()
{

	cov_bin_table = HashTableCreate(200);
	HashTableSetHashFunction(cov_bin_table, fc_chro_hash);
	HashTableSetKeyComparisonFunction(cov_bin_table , fc_strcmp_chro);
	HashTableSetDeallocationFunctions(cov_bin_table , free, free);

	unsigned int this_hit_locs[MAX_FRAGMENT_LENGTH], reads = 0;
	char this_hit_bases[MAX_FRAGMENT_LENGTH];
	coverage_bin_entry_t * chrbin12[2];
	SamBam_FILE * in_fp = SamBam_fopen(input_file_name, is_BAM_input?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM);
	char * line_buffer = malloc(3000);
	int this_hits = 0, r1_items = 0;
	
	char last_chro[200],chro[200], *rtext = NULL;
	last_chro[0]=0;

	while(1)
	{
		char * is_ret = SamBam_fgets(in_fp, line_buffer, 2999, 1);
		if(!is_ret) break;
		if(line_buffer[0]=='@'){
			if(strstr(line_buffer,"@SQ\t"))
				add_chro(line_buffer);
		} else {
			char * Chros[ max_M ];
			unsigned int Staring_Points[max_M];
			unsigned short Staring_Read_Points[max_M];
			unsigned short Section_Lengths[max_M];
	
			int flags=0, x1, is_junc = 0, tlen = 0;
			char cigar_str[200];
			unsigned int pos_0base = 0;
			cigar_str[0]=0;

			if(!paired_end) this_hits = 0;
			if(paired_end && (reads%2) == 0){
				this_hits = 0;
				r1_items = 0;
				chrbin12[0] = chrbin12[1] = NULL;
			}
			if(paired_end && (reads%2) == 1)
				strcpy(last_chro, chro);

			chro[0]=0;

			get_read_info(line_buffer, chro, &pos_0base, cigar_str, &flags, &tlen, &rtext);
			pos_0base --;

			if(0 && FIXLENstrcmp("NS500643:166:HNCCHBGXY:1:23105:1714:15599" , line_buffer)==0)
				SUBREADprintf("INPUT READ_%d: FLAG=%d, chro=%s, pos=%u\n", 1+(reads % 2 ), flags, chro, pos_0base);
			//if((flags & 4) || ( paired_end && (tlen == 0 || abs(tlen) > 500000) )){
			//	reads++;
			//	continue;
			//}

			coverage_bin_entry_t * chrbin = NULL;
			unsigned int chrlen = 0;
			int cigar_sections = 0;
			if((flags & 4)==0){
				void ** bin_entry = HashTableGet(cov_bin_table, chro);
				if(NULL == bin_entry)
				{
					SUBREADprintf("ERROR: The chromosome name is not in header:%s\n", chro);
				}

				chrbin = (coverage_bin_entry_t*) bin_entry[0];
				chrlen = (void *)( bin_entry[1]) - NULL;
				cigar_sections = RSubread_parse_CIGAR_string(chro, pos_0base, cigar_str, max_M, Chros, Staring_Points, Staring_Read_Points, Section_Lengths, &is_junc);
			}
			if(paired_end) {
				for(x1 = 0; x1 < cigar_sections; x1++){
					unsigned int x2,x3;

					//if(strcmp( cigar_str, "8S2M1D14M4D5M1D9M1I14M1D9M4D10M1D23M1I4M1D9M1I9M4D5M4D2M" )==0){
					//	SUBREADprintf("Cigar [%d] = %u ~ + %u ; this_hits=%d\n", x1, Staring_Points[x1], Section_Lengths[x1], *this_hits);
					//}

					int read_cursor = Staring_Read_Points[x1];
					for(x2 = Staring_Points[x1]; x2<Staring_Points[x1]+Section_Lengths[x1]; x2++){
						int found = 0;

						if(reads%2 == 1 && strcmp(chro,last_chro)==0)
							for(x3=0; x3<this_hits; x3++){
								if(this_hit_locs[x3] == x2){
									found =1;
									break;
								}
							}
						if(!found){
							this_hit_bases[this_hits] = base_value_to_int(rtext[read_cursor]);
							this_hit_locs[this_hits++] = x2;
							if(reads % 2 == 0) r1_items ++;
						}
						if(this_hits >= MAX_FRAGMENT_LENGTH){
							SUBREADprintf("ERROR: read is too long : %s!\n", cigar_str);
							return -1;
						}
						read_cursor ++;
					}
				}

				chrbin12[reads%2] = chrbin;
			} else for(x1 = 0; x1 < cigar_sections; x1++) {
				unsigned int x2;
				int read_cursor = Staring_Read_Points[x1];
				for(x2 = Staring_Points[x1]; x2<Staring_Points[x1]+Section_Lengths[x1]; x2++) {
					if(x2 < chrlen) {
						if(base_values){
							chrbin[x2*5 + base_value_to_int(rtext[read_cursor])] ++;
						}
						else if(chrbin[x2] <= COVERAGE_MAX_INT)chrbin[x2] ++;
						all_counted ++;
						if(all_counted % 10000000 == 0)
							SUBREADprintf("Processed %llu bases.\n", all_counted);
					} else {
						SUBREADprintf("%s:%s %u [%s] :: %u-%u\n", line_buffer , chro, pos_0base, cigar_str, Staring_Points[x1], Staring_Points[x1]+Section_Lengths[x1]);
						SUBREADprintf("Read %s overhangs the boundary of chromosome %s (%u >= %u)\n", line_buffer, chro, x2, chrlen);
						//exit(-1);
					}
					read_cursor ++;
				}
			}
			//SUBREADprintf("SREAD_%d , R1ITEMS=%d, R1+2ITEMS=%d, BIN1=%p, BIN2=%p\n", reads, r1_items, this_hits, chrbin12[0], chrbin12[1]);
			if((reads % 2) == 1 && paired_end){
				int x1;
				for(x1 = 0; x1<this_hits;x1++){
					coverage_bin_entry_t * this_ch_binn = x1 < r1_items?chrbin12[0]:chrbin12[1];
					if(base_values) this_ch_binn[this_hit_locs[x1]*5 + this_hit_bases[x1]]++;
					else this_ch_binn[this_hit_locs[x1]]++;

					all_counted ++;
					if(all_counted % 10000000 == 0)
						SUBREADprintf("Processed %llu bases.\n", all_counted);
				}
			}
			reads ++;
		}
	}

	free(line_buffer);

	SamBam_fclose(in_fp);


	SUBREADprintf("Processed totally %llu bases.\nNow write results.\n", all_counted);

	int bucket;
	KeyValuePair *cursor;
	for(bucket=0; bucket < cov_bin_table  -> numOfBuckets; bucket++)
	{
		cursor = cov_bin_table -> bucketArray[bucket];
		while (1)
		{
			if (!cursor) break;
			coverage_bin_entry_t * chrbin = (coverage_bin_entry_t*)(((void **) cursor -> value)[0]);
			unsigned int chrlen = (((void **) cursor -> value)[1]) - NULL; 
			char * chro = (char *)(cursor -> key);
			char out_name[340];
			sprintf(out_name,"%s-%s.bin", output_file_name, chro);

			FILE * fpo = fopen(out_name,"w");
			//SUBREADprintf("PTR %s = %p\n", chro, chrbin);
			if(0 && strcmp(chro, "chr10")==0){
				SUBREADprintf("WRITE chr10 : %u ptr %p\n", chrbin [75631235-1], chrbin);
			}
			fwrite(chrbin, sizeof(coverage_bin_entry_t), chrlen * (base_values?5:1), fpo);
			fclose(fpo);
			free(chrbin);

			SUBREADprintf("Wrote bin for %s\n", chro);
			cursor = cursor->next;
		}	
	}

	HashTableDestroy(cov_bin_table);

	SUBREADprintf("Calculation finished.\n");
	return 0;
}



#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int cov_calc_main(int argc, char ** argv)
#endif
{
	int ret = 0;
	int c=0;
	int option_index=0;
	input_file_name[0]=0;
	output_file_name[0]=0;
	max_M = 10;
	is_BAM_input=0;
	all_counted = 0;

	optind=0;
	opterr=1;
	optopt=63;

	while ((c = getopt_long (argc, argv, "bCBpM:i:o:?", cov_calc_long_options, &option_index)) != -1)
		switch(c)
		{
			case 'C':
				ignore_bad_pair = 1;
				break;
			case 'M':
				if(!is_valid_digit_range(optarg, "maxMOp", 1 , 64))
					exit(-1);
				max_M = atoi(optarg);
			break;
			case 'b':
				base_values = 1;
			break;
			case 'p':
				paired_end = 1;
			break;
			case 'i':
				strcpy(input_file_name, optarg);
			break;
			case 'o':
				strcpy(output_file_name, optarg);
			break;

			case '?':
			default :
				calcCount_usage();
				return -1;
	
		}
	

	if((!output_file_name[0])||(!input_file_name[0]))
	{
		calcCount_usage();
		return 0;
	}

	int is_bam = is_certainly_bam_file(input_file_name, NULL, NULL);

	if(1==is_bam) is_BAM_input = 1;
	else if(is_bam < 0)
	{
		ret = -1;
		SUBREADprintf("Unable to open input file '%s' or the input file is empty!\n", input_file_name);
	}

	ret = ret || covCalc();

	return ret;
}
