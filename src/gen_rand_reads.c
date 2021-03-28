#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <sys/time.h>
#include <getopt.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#ifndef __MINGW32__
#include <sys/resource.h>
#endif
#include <locale.h>
#include <ctype.h>
#include <math.h>
#include <zlib.h>

#ifndef MAKE_STANDALONE
  #include <R.h>
#endif



#include "subread.h"
#include "core.h"
#include "HelperFunctions.h"
#include "seek-zlib.h"
#include "input-files.h"
#include "gene-algorithms.h"

#define MAX_SIMULATION_READ_LEN 250
#define STRATEGY_LESS_THAN_N 10
#define STRATEGY_RANDOM_ASSIGN 20
#define STRATEGY_ITERATIVE_M 30

int print_usage_gen_reads(char * pgname) {
	SUBREADputs("");
	SUBREADputs("Usage:");
	SUBREADputs("");
	SUBREADputs(" For scanning a FASTA/gz file:");
	SUBREADprintf("    %s --summarizeFasta \\\n", pgname);
	SUBREADputs("       --transcriptFasta <file> --outputPrefix <string> [--simpleTranscriptId]");
	SUBREADputs("");
	SUBREADputs(" For generating read/pairs:");
	SUBREADprintf("    %s --transcriptFasta <file>\\\n", pgname);
	SUBREADputs("       --outputPrefix <string> --expressionLevels <file> [other options]");
	SUBREADputs("");
	SUBREADputs(" --summarizeFasta           Only output the transcript names and lengths.");
	SUBREADputs("");
	SUBREADputs(" --transcriptFasta <file>   The transcript database in FASTA/gz format.");
	SUBREADputs("");
	SUBREADputs(" --outputPrefix <string>    The prefix of the output files.");
	SUBREADputs("");
	SUBREADputs(" --totalReads  <int>        Total read/pairs in output.");
	SUBREADputs("");
	SUBREADputs(" --expressionLevels <file>  Two column table delimited by <TAB>, giving the");
	SUBREADputs("                            wanted TPM values. Columns: TranscriptID and TPM");
	SUBREADputs("");
	SUBREADputs(" --readLen <int>            The length of the output reads. 100 by default.");
	SUBREADputs("");
	SUBREADputs(" --totalReads <int>         Total read/pairs in the output.");
	SUBREADputs("");
	SUBREADputs(" --randSeed <int64>         Seed to generate random numbers. UNIXTIME is used");
	SUBREADputs("                            as the random seed by default.");
	SUBREADputs("");
	SUBREADputs(" --qualityRefFile <file>    A textual file containing Phred+33 quanlity strings");
	SUBREADputs("                            for simulating sequencing errors. The quality");
	SUBREADputs("                            strings have to have the same length as the output");
	SUBREADputs("                            reads. No sequencing errors are simulated when this");
	SUBREADputs("                            option is omitted.");
	SUBREADputs("");
	SUBREADputs(" --floorStrategy            How to deal with round-up errors. 'FLOOR': generate");
	SUBREADputs("                            less than wanted reads; 'RANDOM': randomly assign");
	SUBREADputs("                            margin reads to transcripts; 'ITERATIVE': find the");
	SUBREADputs("                            best M value to have ~N reads.");
	SUBREADputs("");
	SUBREADputs(" --pairedEnd                Generate paired-end reads.");
	SUBREADputs("");
	SUBREADputs(" --insertionLenMean <float>,--insertionLenSigma <float>,--insertionLenMin <int>,");
	SUBREADputs(" --insertionLenMax <int>    Parameters of a truncated normal distribution for");
	SUBREADputs("                            deciding insertion lengths of paired-end reads.");
	SUBREADputs("                            Default values: mean=160, sigma=30, min=110, max=400");
	SUBREADputs("");
	SUBREADputs(" --simpleTranscriptId       Truncate transcript names to the first '|' or space.");
	SUBREADputs("");
	SUBREADputs(" --truthInReadNames         Encode the true locations of reads in read names.");
	SUBREADputs("");
	SUBREADputs(" --noActualReads            Do not actually generate reads in fastq.gz files.");
	SUBREADputs("");
	return 1;
}

static struct option long_options[] =
{
	{"quiet", no_argument, 0, 'Q'},
	{"truthInReadNames", no_argument, 0, 'T'},
	{"simpleTranscriptId", no_argument, 0, 'C'},
	{"summarizeFasta",  no_argument, 0, 'M'},
	{"transcriptFasta",  required_argument, 0, 't'},
	{"totalReads",  required_argument, 0, 'r'},
	{"pairedEnd",  no_argument, 0, 'p'},
	{"expressionLevels",  required_argument, 0, 'e'},
	{"qualityRefFile",  required_argument, 0, 'q'},
	{"outputPrefix",  required_argument, 0, 'o'},
	{"randSeed",  required_argument, 0, 'S'},
	{"readLen",  required_argument, 0, 'L'},
	{"floorStrategy",  required_argument, 0, 'O'},
	{"noActualReads",  no_argument, 0, 'x'},
	{"insertionLenMean",  required_argument, 0, 'F'},
	{"insertionLenMin",  required_argument, 0, 'N'},
	{"insertionLenMax",  required_argument, 0, 'X'},
	{"insertionLenSigma",  required_argument, 0, 'V'},
	{0, 0, 0, 0}
};

typedef struct {
	char random_seeds[16];
	char transcript_fasta_file[MAX_FILE_NAME_LENGTH];
	char output_prefix[MAX_FILE_NAME_LENGTH];
	char expression_level_file[MAX_FILE_NAME_LENGTH];
	char quality_string_file[MAX_FILE_NAME_LENGTH];

	unsigned long long output_sample_size;
	unsigned long long applied_M;
	int floor_strategy;
	int is_paired_end;
	int simple_transcript_names;
	int truth_in_read_names;
	float insertion_length_mean;
	int insertion_length_max;
	int insertion_length_min;
	float insertion_length_sigma;
	int quiet;
	int read_length;
	int no_actual_reads;

	ArrayList * quality_strings;
	ArrayList * transcript_hitting_space;
	ArrayList * transcript_names;
	HashTable * transcript_sequences;
	HashTable * transcript_lengths;
	HashTable * expression_levels;

	char fake_quality_string[MAX_SIMULATION_READ_LEN+3];
	char * cmd_line;
	gzFile out_fps[2];
	FILE * counts_out_fp;
} genRand_context_t;

void grc_incrand(genRand_context_t * grc){
	unsigned long long round_rand =0;
	memcpy(&round_rand, grc->random_seeds+8, sizeof(round_rand));
	round_rand++;
	memcpy(grc->random_seeds+8, &round_rand, sizeof(round_rand));
}

void grc_sequencing_error_read(char * seq, int qlen, char * qua){
	int b;
	for(b=0; b<qlen; b++){
		if(seq[b]=='N') continue;

		int qub = qua[b];
		float randv = myrand_rand()*1./RAND_MAX;
		float errorp = pow(10,3.3 - qub*.1); // Phred 33
		if(randv < errorp * 1.3333333333){// ATGC random can be the original
			seq[b]="ACGT"[myrand_rand()%4];
		}
	}
}

void gen_one_read_here(genRand_context_t * grc, char * seq, int is_PE_second, int trans_negative, unsigned long long rno, char * seq_name, int my_pos, int mate_pos){
	char read_seq [grc -> read_length+1];
	memcpy(read_seq, seq, grc -> read_length);
	read_seq[grc -> read_length]=0;
	if(trans_negative) reverse_read(read_seq , grc -> read_length , GENE_SPACE_BASE);

	char * qual_str = NULL;
	if(grc->quality_strings->numOfElements>0){
		qual_str = ArrayListRandom(grc->quality_strings);
		//SUBREADprintf("TESTQUAL\t%p\n", qual_str);
		grc_sequencing_error_read(read_seq, grc -> read_length, qual_str);
	}else{
		if(!grc->fake_quality_string[0]){
			int xx;
			for(xx = 0; xx < grc -> read_length; xx++) grc->fake_quality_string[xx]='H';
			grc->fake_quality_string[xx]=0;
		}
		qual_str = grc->fake_quality_string;
	}

	gzFile thisfp = (is_PE_second==1) ? grc -> out_fps[1] : grc -> out_fps[0];
	int R1_pos = is_PE_second?mate_pos:my_pos;
	int R2_pos = is_PE_second?my_pos:mate_pos;
	if(grc->truth_in_read_names){
		if(is_PE_second<0) gzprintf(thisfp, "@R%09llu:%s:%d\n%s\n+\n%s\n", rno, seq_name, 1+my_pos, read_seq, qual_str);
		else gzprintf(thisfp, "@R%09llu:%s:%d:%d\n%s\n+\n%s\n", rno, seq_name, 1+R1_pos, 1+R2_pos, read_seq, qual_str);
	}
	else gzprintf(thisfp, "@R%09llu\n%s\n+\n%s\n", rno, read_seq, qual_str);
}

void gen_a_read_from_one_transcript(genRand_context_t * grc, long this_transcript_no, unsigned  long long rno){
	char * trans_name = ArrayListGet(grc->transcript_names, this_transcript_no);
	char * trans_seq = HashTableGet(grc->transcript_sequences, trans_name);
	int actual_transcript_len = HashTableGet(grc->transcript_lengths, trans_name) - NULL;
	int applied_insertion_maxlen = min(grc -> insertion_length_max, actual_transcript_len);
	double rand_01 = plain_txt_to_long_rand(grc->random_seeds, 16)*1./0xffffffffffffffffllu;
	int rand_01_int = (int)(rand_01*901267351);
	myrand_srand(rand_01_int); // for generating sequencing errors. NOTE: the argument is unused in R. MYRAND module uses R's RNG for seeding.
	grc_incrand(grc);

	if(grc -> is_paired_end){
		float insertion_len = inverse_sample_normal(rand_01) * grc-> insertion_length_sigma + grc -> insertion_length_mean;
		int fraglen = (int)(min(max(insertion_len, grc -> insertion_length_min), applied_insertion_maxlen));
		rand_01 = plain_txt_to_long_rand(grc->random_seeds, 16)*1./0xffffffffffffffffllu;
		grc_incrand(grc);
		int start_pos = (actual_transcript_len - fraglen) * rand_01;
		int is_first_end_negative = rand_01_int % 2;
		if(is_first_end_negative){
			gen_one_read_here(grc, trans_seq + start_pos + fraglen - grc -> read_length, 0, 1, rno, trans_name, start_pos + fraglen - grc -> read_length, start_pos);
			gen_one_read_here(grc, trans_seq + start_pos, 1, 0, rno, trans_name, start_pos, start_pos + fraglen - grc -> read_length);
		}else{
			gen_one_read_here(grc, trans_seq + start_pos, 0, 0, rno, trans_name, start_pos, start_pos + fraglen - grc -> read_length);
			gen_one_read_here(grc, trans_seq + start_pos + fraglen - grc -> read_length, 1, 1, rno, trans_name, start_pos + fraglen - grc -> read_length, start_pos);
		}
	}else{
		int start_pos = (actual_transcript_len - grc -> read_length)*rand_01;
		int is_negative = rand_01_int % 2;
		gen_one_read_here(grc, trans_seq + start_pos, -1, is_negative, rno, trans_name, start_pos, -1);
	}
}

int grc_check_parameters(genRand_context_t * grc){
	int ret = 0;

	if(grc->read_length > MAX_SIMULATION_READ_LEN){
		SUBREADprintf("Error: the read length cannot be higher than  %d.\n", MAX_SIMULATION_READ_LEN);
		ret=1;
	}
	if(grc->is_paired_end){
		if(grc->insertion_length_min>grc->insertion_length_max){
			SUBREADprintf("Error: the minimum insertion length must be equal or higher than the maximum insertion length.\n");
			ret=1;
		}
	
		if(grc->insertion_length_min<grc->read_length){
			SUBREADprintf("Error: the minimum insertion length must be equal or higher than read length.\n");
			ret=1;
		}
	
		if(grc->insertion_length_max<1){
			SUBREADprintf("Error: the maximum insertion length must be a positive number.\n");
			ret=1;
		}
	}

	if(grc->read_length<1){
		SUBREADprintf("Error: the read length must be a positive number.\n");
		ret=1;
	}

	if(!grc->transcript_fasta_file[0]){
		SUBREADprintf("Error: a transcript file must be provide.\n");
		ret=1;
	}

	if(!grc->output_prefix[0]){
		SUBREADprintf("Error: the output prefix must be provide.\n");
		ret=1;
	}else{
		char outname[MAX_FILE_NAME_LENGTH+30];
		sprintf(outname, "%s.for_test.log",grc->output_prefix);
		FILE * test_out = fopen(outname, "w");
		if(test_out){
			fclose(test_out);
			unlink(outname);
		}else{
			SUBREADprintf("Error: cannot create the output file.\n");
			ret=1;
		}
	}

	if(!grc->expression_level_file[0]){
		SUBREADprintf("Error: the wanted expression levels must be provided.\n");
		ret=1;
	}

	if(grc->output_sample_size < 1){
		SUBREADprintf("Warning: no read number is specified. Generating one million read%s.\n", grc->is_paired_end?"-pairs":"s");
		grc->output_sample_size = 1000000;
	}

	return ret;
}

unsigned long long calc_N_from_M(genRand_context_t *grc , unsigned long long Mval){
	long ti;
	unsigned long long space_end = ArrayListGet(grc->transcript_hitting_space, grc->transcript_hitting_space->numOfElements -1)-NULL;
	unsigned long long ret =0, lastv=0;

	for(ti = 0; ti < grc->transcript_hitting_space->numOfElements; ti++){
		unsigned long long thisv = ArrayListGet(grc->transcript_hitting_space, ti) - NULL;
		unsigned long long this_space_span = thisv - lastv;
		unsigned long long expected_reads =(unsigned long long )(this_space_span *1.0/space_end * Mval);  // expected_reads = exp_level * seq_length / sum(Q*L) * M
		ret += expected_reads;
		lastv = thisv;
	}
	return ret;
}

unsigned long long itr_find_M(genRand_context_t *grc){
	long transcripts=  grc->transcript_hitting_space->numOfElements;
	unsigned long long wanted_reads = grc->output_sample_size;

	unsigned long long LL = wanted_reads, HH = wanted_reads + transcripts;
	unsigned long long MM = 0;
	
	while(1){
		MM = (LL+HH)/2;
		unsigned long long thisN_at_M = calc_N_from_M(grc, MM);
		if(thisN_at_M < wanted_reads) LL = MM+1;
		else if(thisN_at_M > wanted_reads) HH = MM-1;
		else break;
		if(LL >=HH){
			MM = max(LL,HH);
			break;
		}
	}
	return MM;
}

unsigned long long convert_hitting_space_to_num_of_reads(genRand_context_t *grc , ArrayList * num_of_frags_per_transcript, int min_seq_len){
	unsigned long long lastv = 0, current_total =0;

	ArrayList * rescure_hitting_space = ArrayListCreate(100000);
	unsigned long long to_rescure_read_top=0;
	unsigned long long read_i = 0;
	unsigned long long space_end = ArrayListGet(grc->transcript_hitting_space, grc->transcript_hitting_space->numOfElements -1)-NULL;

	if(grc -> floor_strategy == STRATEGY_ITERATIVE_M){
		grc->applied_M = itr_find_M(grc);
	}else grc->applied_M = grc->output_sample_size;
	// after M is here, the reminder is just like "LESS_THAN_N.

	{
		for(read_i = 0; read_i < grc->transcript_hitting_space->numOfElements ; read_i++){
			char *seq_name = ArrayListGet(grc->transcript_names, read_i);
			int seq_len = HashTableGet(grc-> transcript_lengths, seq_name)-NULL;
			unsigned long long thisv = ArrayListGet(grc->transcript_hitting_space, read_i) - NULL;
			unsigned long long this_space_span = thisv - lastv;
	
			// the many 9's is to aovid the sum of read numbers larger than wanted.
			double floor_chopping = grc -> floor_strategy ==   STRATEGY_RANDOM_ASSIGN ?0.999999999:1;
			unsigned long long expected_reads =(unsigned long long )((this_space_span *1.0/space_end) * grc->applied_M*floor_chopping);
			unsigned long long to_rescure_reads = 0;
			if(grc -> floor_strategy == STRATEGY_RANDOM_ASSIGN) to_rescure_reads = (unsigned long long)((this_space_span *1.0/space_end * grc->applied_M- 1.*expected_reads)*100000.);
	
			if(this_space_span < 1) to_rescure_reads=0;
			if( seq_len < min_seq_len ){
				to_rescure_reads = 0;
				expected_reads = 0;
			}
			to_rescure_read_top+= to_rescure_reads;
			assert(to_rescure_read_top < 0x5fffffffffffffffllu);
			ArrayListPush(rescure_hitting_space, NULL+to_rescure_read_top);
			ArrayListPush(num_of_frags_per_transcript, NULL+expected_reads);
			current_total += expected_reads;
	
			lastv = thisv;
		}
	
		// I'm not sure why this is so important but let it be.
		assert(current_total<=grc->applied_M);
	
		if(grc -> floor_strategy == STRATEGY_RANDOM_ASSIGN)for(read_i = current_total; read_i < grc->applied_M; read_i++){
			unsigned long long longrand = plain_txt_to_long_rand(grc->random_seeds, 16);
			grc_incrand(grc);
	
			longrand = longrand % to_rescure_read_top;
			long this_transcript_no = ArrayListFindNextDent(rescure_hitting_space, longrand);
			unsigned long long expected_reads = ArrayListGet(num_of_frags_per_transcript, this_transcript_no)-NULL;
			expected_reads++;
			num_of_frags_per_transcript->elementList[this_transcript_no] = NULL+expected_reads;
			current_total ++;
		}
	}

	ArrayListDestroy(rescure_hitting_space);
	return current_total;
}

// The 823,532,653,200th prime is 24,537,224,085,139.
//     -- https://primes.utm.edu/nthprime/index.php
#define A_LARGE_PRIME_FOR_MOD 24537224085139llu


int grc_gen( genRand_context_t *grc ){
	int ret = 0;
	unsigned long long read_i = 0;

	ArrayList * num_of_frags_per_transcript = ArrayListCreate(100000);
	int min_seq_len = grc->is_paired_end?grc->insertion_length_min:grc->read_length;

	convert_hitting_space_to_num_of_reads(grc, num_of_frags_per_transcript, min_seq_len);

	ArrayList * per_transcript_reads_hitting_space = ArrayListCreate(100000);
	unsigned long long total_read_top =0;

	for(read_i =0; read_i < num_of_frags_per_transcript -> numOfElements; read_i++) {
		char *seq_name = ArrayListGet(grc->transcript_names, read_i);
		int seq_len = HashTableGet(grc-> transcript_lengths, seq_name)-NULL;
		unsigned long long expected_reads = ArrayListGet(num_of_frags_per_transcript, read_i)-NULL;
		if(seq_len >= min_seq_len)
			#ifdef __MINGW32__
			fprintf(grc->counts_out_fp, "%s\t%d\t%I64u\n", seq_name, seq_len, expected_reads);
			#else
			fprintf(grc->counts_out_fp, "%s\t%d\t%llu\n", seq_name, seq_len, expected_reads);
			#endif
		else
			fprintf(grc->counts_out_fp, "%s\t%d\tNA\n", seq_name, seq_len);
		total_read_top+=expected_reads;
		ArrayListPush(per_transcript_reads_hitting_space, NULL+total_read_top);
	}


	if(grc -> floor_strategy == STRATEGY_RANDOM_ASSIGN) assert(total_read_top == grc->applied_M);
	else grc->applied_M = total_read_top;

	unsigned long long longrand = plain_txt_to_long_rand(grc->random_seeds, 16);
	grc_incrand(grc);

	if(grc->no_actual_reads == 0){
		unsigned long long mod_class = longrand % grc->applied_M; // an arbitratry starting point.
		for(read_i =0; read_i < grc->applied_M; read_i++) {
			mod_class += A_LARGE_PRIME_FOR_MOD;
			mod_class = mod_class % grc->applied_M;
			long this_transcript_no = ArrayListFindNextDent(per_transcript_reads_hitting_space, mod_class);
			//char * trans_name = ArrayListGet(grc->transcript_names, this_transcript_no);
			//SUBREADprintf("TESTGEN\t%s\n", trans_name);
			gen_a_read_from_one_transcript(grc, this_transcript_no, read_i);
		}
	}

	ArrayListDestroy(num_of_frags_per_transcript);
	ArrayListDestroy(per_transcript_reads_hitting_space);
	return ret;
}

int grc_finalize(genRand_context_t *grc){
	HashTableDestroy(grc->expression_levels);
	HashTableDestroy(grc->transcript_sequences);
	HashTableDestroy(grc->transcript_lengths);
	ArrayListDestroy(grc->quality_strings);
	ArrayListDestroy(grc->transcript_hitting_space);
	ArrayListDestroy(grc->transcript_names);
	gzclose(grc->out_fps[0]);
	if(grc->out_fps[1]) gzclose(grc->out_fps[1]);
	fclose(grc->counts_out_fp);
	free(grc->cmd_line);
	#ifdef __MINGW32__
	SUBREADprintf("Finished. Actual sample size : I64u\n", grc->applied_M);
	#else
	SUBREADprintf("Finished. Actual sample size : %llu\n", grc->applied_M);
	#endif
	return 0;
}

#define TRANSCRIPT_FASTA_LINE_INIT 800
#define TRANSCRIPT_MAX_EXPRESSION_LEVEL 1000001.0
#define TRANSCRIPT_FASTA_LINE_WIDTH 1000



int grc_summary_fasta(genRand_context_t * grc){
	char outname[MAX_FILE_NAME_LENGTH+30];
	autozip_fp auto_FP;

	if(!grc->output_prefix[0]){
		SUBREADprintf("Error: the output prefix must be provide.\n");
		return -1;
	}

	sprintf(outname,"%s.faSummary", grc->output_prefix);
	int ret = autozip_open(grc->transcript_fasta_file, &auto_FP);
	if(ret<0){
		SUBREADprintf("Error: cannot open the fasta file as input\n");
		return -1;
	}else ret = 0;

	FILE * sumfp = fopen(outname, "w");
	if(sumfp == NULL){
		SUBREADprintf("Error: cannot open the putput file\n");
		return -1;
	}
	fprintf(sumfp, "TranscriptID\tLength\tMD5\tUnique\tOccurrence\tDuplicated\n");

	char * seq_name = NULL;
	unsigned char md5res[16];
	int seq_len = 0, total_rep_seq = 0;
	HelpFuncMD5_CTX md5ctx;
	HelpFuncMD5_Init(&md5ctx);
	HashTable * seq_duplicate_tab = StringTableCreate(100000);
	HashTable * name_duplicate_tab = StringTableCreate(100000);
	HashTable * seq_length_tab = StringTableCreate(100000);
	HashTable * seq_md5_tab = StringTableCreate(100000);
	ArrayList * seq_name_list = ArrayListCreate(100000);
	HashTableSetDeallocationFunctions(seq_duplicate_tab, free,NULL);
	ArrayListSetDeallocationFunction(seq_name_list, free);

	while(1){
		char clinebuf[TRANSCRIPT_FASTA_LINE_WIDTH];
		int rlength = autozip_gets(&auto_FP, clinebuf, TRANSCRIPT_FASTA_LINE_WIDTH -1);
		if(rlength < 1)break;
		if(rlength >= TRANSCRIPT_FASTA_LINE_WIDTH -1 || clinebuf[rlength]!='\0' || clinebuf[rlength-1]!='\n'){
			SUBREADprintf("Error: The line width of the fasta file excessed %d bytes.\n", TRANSCRIPT_FASTA_LINE_WIDTH);
			ret = 1;
			break;
		}
		if(clinebuf[0]=='>'){
			if(seq_name){
				if(seq_len<1){
					SUBREADprintf("Error: a transcript has no sequence: '%s'\n", seq_name);
					return -1;
				}
				HelpFuncMD5_Final(md5res, &md5ctx);
				//SUBREADprintf("%s\t",seq_name);int md5i;for(md5i=0;md5i<16;md5i++)SUBREADprintf("%02X",0xff&(int)md5res[md5i]);SUBREADputs("");
				char * md5mem = malloc(33);
				int md5i;for(md5i=0;md5i<16;md5i++)sprintf(md5mem+2*md5i, "%02X", 0xff&(int)md5res[md5i]);

				long dupval = HashTableGet(seq_duplicate_tab,md5mem)-NULL;
				dupval++;
				HashTablePutReplace( seq_duplicate_tab, md5mem, NULL+dupval , 0);
				if(dupval>1){
					char * newmd5mem =HashTableGetKey(seq_duplicate_tab, md5mem);
					free(md5mem);
					md5mem = newmd5mem;
				}
				ArrayListPush( seq_name_list, seq_name);
				HashTablePut( seq_length_tab, seq_name, NULL+seq_len);
				HashTablePut( seq_md5_tab, seq_name, md5mem);
				//if(strcmp("ENST00000410691.1",seq_name)==0)SUBREADprintf("MD5CHECK: %s has %s (%p)\n", seq_name, md5mem, md5mem);
				//if(strcmp("ENST00000634833.2",seq_name)==0)SUBREADprintf("MD5CHECK: %s has %s (%p), not %s (%p)\n", seq_name, md5mem, md5mem, HashTableGet(seq_md5_tab, "ENST00000410691.1"), HashTableGet(seq_md5_tab, "ENST00000410691.1"));

				seq_len = 0;
			}
			clinebuf[rlength-1]=0;
			if(grc->simple_transcript_names){
				int xx;
				for(xx=1; xx<rlength-1; xx++) if(clinebuf[xx]=='|' || clinebuf[xx]==' ') clinebuf[xx]=0;
			}

			void * had_tab = HashTableGet(name_duplicate_tab, clinebuf+1);
			//SUBREADprintf("CHECK PTR: %s => %p\n", clinebuf+1, had_tab);
			if(had_tab){
				SUBREADprintf("Error: duplicate sequence name was found : '%s'.\n", clinebuf+1);
				return -1;
			}
			seq_name=malloc(rlength);
			strcpy(seq_name, clinebuf+1);

			HashTablePut(name_duplicate_tab, seq_name, NULL+1);
			HelpFuncMD5_Init(&md5ctx);
		}else{
			int xx; for(xx=0; xx<rlength-1; xx++) clinebuf[xx] = toupper(clinebuf[xx]);
			HelpFuncMD5_Update(&md5ctx, clinebuf, rlength-1);
			seq_len += rlength-1; // no \n
		}
	}

	if(seq_name){
		if(seq_len<1){
			SUBREADprintf("Error: a transcript has no sequence: '%s'\n", seq_name);
			return -1;
		}

		HelpFuncMD5_Final(md5res, &md5ctx);
		//SUBREADprintf("%s\t",seq_name);int md5i;for(md5i=0;md5i<16;md5i++)SUBREADprintf("%02X",0xff&(int)md5res[md5i]);SUBREADputs("");
		char * md5mem = malloc(33);
		int md5i;for(md5i=0;md5i<16;md5i++)sprintf(md5mem+2*md5i, "%02X", 0xff&(int)md5res[md5i]);
		long dupval = HashTableGet(seq_duplicate_tab,md5mem)-NULL;
		dupval++;
		HashTablePutReplace( seq_duplicate_tab, md5mem, NULL+dupval ,0);

		if(dupval>1){
			char * newmd5mem =HashTableGetKey(seq_duplicate_tab, md5mem);
			free(md5mem);
			md5mem = newmd5mem;
		}

		ArrayListPush( seq_name_list, seq_name);
		HashTablePut( seq_length_tab, seq_name, NULL+seq_len);
		HashTablePut( seq_md5_tab, seq_name, md5mem);
	}

	int seqi;
	HashTable * reprs_tab = StringTableCreate(100000);
	for(seqi = 0 ; seqi < seq_name_list->numOfElements; seqi++){
		char * seqnam = ArrayListGet(seq_name_list,seqi);
		char * md5str = HashTableGet(seq_md5_tab,seqnam);
		long md5repeated = HashTableGet(seq_duplicate_tab, md5str) - NULL;
		if(md5repeated>1) total_rep_seq++;
		long seqlen = HashTableGet(seq_length_tab,seqnam)-NULL;
		int is_reprs = HashTableGet(reprs_tab,md5str)==NULL;
		HashTablePut(reprs_tab,md5str,NULL+1);
		fprintf(sumfp, "%s\t%ld\t%s\t%s\t%ld\t%s\n", seqnam, seqlen, md5str, md5repeated>1?"FALSE":"TRUE" /*"Unique"*/, md5repeated, is_reprs?"FALSE":"TRUE" /*"Duplicated"*/);
	}
	HashTableDestroy(reprs_tab);

	if(total_rep_seq && grc->quiet==0)SUBREADprintf("Warning: %d duplicate sequences were found in the input. Please check the summary table.\n",total_rep_seq);

	ArrayListDestroy(seq_name_list);
	HashTableDestroy(seq_length_tab);
	HashTableDestroy(seq_md5_tab);
	HashTableDestroy(seq_duplicate_tab);
	HashTableDestroy(name_duplicate_tab);
	autozip_close(&auto_FP);
	fclose(sumfp);
	return ret;
}

void grc_put_new_trans(genRand_context_t *grc, char * seq_name, char * seq_str, unsigned int seq_len, unsigned long long * linear_space_top){
	if(seq_len<1){
		SUBREADprintf("Warning: a transcript, '%s', has a zero length. No read is generated from it.\n", seq_name);
	}
	HashTablePut(grc-> transcript_sequences,seq_name, seq_str);
	HashTablePut(grc-> transcript_lengths, seq_name, NULL+ seq_len);
	unsigned long long this_seq_exp_10000 = HashTableGet(grc->expression_levels, seq_name)-NULL;
	if(this_seq_exp_10000<1){
		SUBREADprintf("Warning: a transcript, '%s', has no wanted expression level. No read is generated from it.\n", seq_name);
		this_seq_exp_10000=0;
	}else this_seq_exp_10000-=1;
	//SUBREADprintf("TESTLEN\t%s\t%d\n", seq_name, seq_len);
	(*linear_space_top) += this_seq_exp_10000 * seq_len;
	ArrayListPush(grc->transcript_names, seq_name);
	ArrayListPush(grc->transcript_hitting_space, NULL+*linear_space_top);
}

int grc_load_env(genRand_context_t *grc){

	grc->expression_levels = HashTableCreate(100000);
	HashTableSetDeallocationFunctions(grc->expression_levels, free, NULL);
	HashTableSetKeyComparisonFunction(grc->expression_levels, fc_strcmp_chro);
	HashTableSetHashFunction(grc->expression_levels, fc_chro_hash);

	grc->transcript_sequences = HashTableCreate(100000);
	HashTableSetDeallocationFunctions(grc->transcript_sequences, free, free);
	HashTableSetKeyComparisonFunction(grc->transcript_sequences, fc_strcmp_chro);
	HashTableSetHashFunction(grc->transcript_sequences, fc_chro_hash);

	grc->transcript_lengths = HashTableCreate(100000);
	HashTableSetKeyComparisonFunction(grc->transcript_lengths, fc_strcmp_chro);
	HashTableSetHashFunction(grc->transcript_lengths, fc_chro_hash);

	grc -> quality_strings = ArrayListCreate(100000);
	ArrayListSetDeallocationFunction(grc -> quality_strings, free);
	grc -> transcript_hitting_space = ArrayListCreate(100000);
	grc -> transcript_names = ArrayListCreate(100000); // the names are destroyed by destroying grc->transcript_sequences

	autozip_fp auto_FP;
	int xk1;
	int ret = autozip_open(grc->expression_level_file, &auto_FP);
	if(ret<0){
		ret = 1;
		SUBREADprintf("Error: unable to open the expression level file.\n");
	}else ret = 0;
	if(ret) return ret;

	unsigned long long total_tpm = 0;
	while(1){
		char linebuf[400], * tokbuf=NULL;
		int rline = autozip_gets(&auto_FP, linebuf, 399);
		if(rline<1) break;
		if(strstr(linebuf, "ID\tTPM"))continue;
		char * seqname = strtok_r(linebuf, "\t", &tokbuf);
		char * seqexp_str = tokbuf;
		if(NULL == seqexp_str){
			SUBREADprintf("Error: expression level file format error.\n");
			ret = 1;
		}
		double seqexp = atof(seqexp_str);
		if(seqexp > TRANSCRIPT_MAX_EXPRESSION_LEVEL){
			SUBREADprintf("Error: The transcript expression level shouldn't excess %.0f\n", TRANSCRIPT_MAX_EXPRESSION_LEVEL);
		}
		
		unsigned long long seqexp_int = (unsigned long long )(seqexp*10000.);
		total_tpm += seqexp_int;
		char * seqname_buf = malloc(strlen(seqname)+1);
		strcpy(seqname_buf, seqname);

		void * had_tab = HashTableGet(grc->expression_levels, seqname_buf);
		if(had_tab){
			SUBREADprintf("Error: duplicate transcript name was found in the TMP table: '%s'.\n", seqname_buf);
			return -1;
		}
		HashTablePut(grc->expression_levels, seqname_buf, NULL+seqexp_int+1);
	}
	autozip_close(&auto_FP);
	if(ret) return ret;

	#define ROUNDUP_TOLERANCE ( 1000llu * 10000llu )
	if(total_tpm > 1000000llu * 10000llu + ROUNDUP_TOLERANCE || total_tpm < 1000000llu * 10000llu - ROUNDUP_TOLERANCE) {
		SUBREADprintf("Error: total TPM is not 1,000,000\n");
		return 1;
	}




	if(grc->quality_string_file[0]){
		ret = autozip_open(grc->quality_string_file, &auto_FP);
		if(ret<0){
			ret = 1;
			SUBREADprintf("Error: unable to open the quality string file.\n");
		}else ret = 0;
		if(ret) return ret;
		while(1){
			char linebuf[400];
			int rline = autozip_gets(&auto_FP, linebuf, 399);
			if(rline<1) break;
			if(rline==grc->read_length +1 || ( rline==grc->read_length && linebuf[rline-1]!='\n' ) ) {
				// is OK.
			}else{
				SUBREADprintf("Error: all your quality strings must be %d-byte long.\n", grc->read_length);
				ret = 1;
				break;
			}
			char * qstr = malloc(rline);
			memcpy(qstr, linebuf, rline);
			if(qstr[rline-1]=='\n') qstr[rline-1]=0;

			ArrayListPush(grc -> quality_strings, qstr);
		}
		autozip_close(&auto_FP);
	}

	if(ret) return ret;

	ret = autozip_open(grc->transcript_fasta_file, &auto_FP);
	if(ret<0){
		ret = 1;
		SUBREADprintf("Error: unable to open the transcript file.\n");
	} else ret = 0;
	if(ret) return ret;
	

	HelpFuncMD5_CTX md5ctx;
	HelpFuncMD5_Init(&md5ctx);
	HashTable * seq_duplicate_tab = StringTableCreate(100000);
	HashTableSetDeallocationFunctions(seq_duplicate_tab, free, NULL);

	unsigned long long linear_space_top = 0;
	char * lbuf = NULL, * seq_name = NULL;
	unsigned int lbuf_cap = 0, lbuf_used = 0, this_seq_len = 0, total_dup=0;
	while(1){
		
		char clinebuf[TRANSCRIPT_FASTA_LINE_WIDTH];
		int rlength = autozip_gets(&auto_FP, clinebuf, TRANSCRIPT_FASTA_LINE_WIDTH -1);
		if(rlength < 1)break;
		if(rlength >= TRANSCRIPT_FASTA_LINE_WIDTH -1 || clinebuf[rlength]!='\0' || clinebuf[rlength-1]!='\n'){
			SUBREADprintf("Error: The line width of the fasta file excessed %d bytes.\n", TRANSCRIPT_FASTA_LINE_WIDTH);
			ret = 1;
			break;
		}
		if(clinebuf[0]=='>'){
			if(NULL != seq_name){
				char * md5mem = malloc(33); unsigned char md5res[16];
				HelpFuncMD5_Final(md5res, &md5ctx);
				int md5i;for(md5i=0;md5i<16;md5i++)sprintf(md5mem+2*md5i, "%02X", 0xff&(int)md5res[md5i]);
				char * had_tab = HashTableGet(seq_duplicate_tab, md5mem);
				long seq_exp = HashTableGet(grc->expression_levels, seq_name)-NULL;
				if(had_tab && seq_exp>1) total_dup++;//SUBREADprintf("Warning: duplicate sequence was found in '%s' and '%s'.\n", seq_name, had_tab);
				if(seq_exp>1)HashTablePut(seq_duplicate_tab, md5mem, 1+NULL);	
				else free(md5mem);

				had_tab = HashTableGet(grc-> transcript_sequences, seq_name);
				if(had_tab){
					SUBREADprintf("Error: duplicate sequence names were found in the input: '%s'.\n", seq_name);
					return -1;
				}

				grc_put_new_trans(grc, seq_name, lbuf, this_seq_len, &linear_space_top);
			}

			clinebuf[rlength-1]=0;
			if(grc->simple_transcript_names)
				for(xk1=1; xk1<rlength-1; xk1++) if(clinebuf[xk1]=='|' || clinebuf[xk1]==' ') clinebuf[xk1]=0;

			seq_name = malloc(strlen(clinebuf));
			if( clinebuf[1]==0 ){
				SUBREADprintf("Error: Every transcript needs a name.\n");
				ret = 1;
				break;
			}
			strcpy(seq_name, clinebuf+1);
			lbuf_used = 0;
			lbuf = malloc(TRANSCRIPT_FASTA_LINE_INIT);
			lbuf_cap = TRANSCRIPT_FASTA_LINE_INIT;
			HelpFuncMD5_Init(&md5ctx);
		}else{
			if(NULL == seq_name){
				SUBREADprintf("Error: The fasta file did not start correctly! \n");
				ret = 1;
				break;
			}
			if(lbuf_cap - lbuf_used < rlength + 1){
				lbuf_cap = max(lbuf_cap *8/5, lbuf_cap + rlength);
				lbuf = realloc(lbuf, lbuf_cap);
			}

			int xx; for(xx=0; xx<rlength-1; xx++) clinebuf[xx] = toupper(clinebuf[xx]);
			//SUBREADprintf("STCP1 : %d used, %d len, %d cap\n", lbuf_used, strlen(clinebuf), lbuf_cap);
			HelpFuncMD5_Update(&md5ctx, clinebuf, rlength-1);
			strcpy(lbuf + lbuf_used, clinebuf );
			*(lbuf+lbuf_used+rlength-1)=0; // '\n' => 0

			lbuf_used += rlength -1;
			this_seq_len = lbuf_used;
		}
	}
	if(lbuf_used<1){
		SUBREADprintf("Error: The fasta file did not end correctly! \n");
		ret = 1;
	}
	if(NULL != seq_name && lbuf_used >0){
		char * md5mem = malloc(33); unsigned char md5res[16];
		HelpFuncMD5_Final(md5res, &md5ctx);
		int md5i;for(md5i=0;md5i<16;md5i++)sprintf(md5mem+2*md5i, "%02X", 0xff&(int)md5res[md5i]);
		long seq_exp = HashTableGet(grc->expression_levels, seq_name)-NULL;
		char * had_tab = HashTableGet(seq_duplicate_tab, md5mem);
		if(had_tab && seq_exp>1)total_dup++;// SUBREADprintf("Warning: duplicate sequence was found in '%s' and '%s'.\n", seq_name, had_tab);
		free(md5mem);

		had_tab = HashTableGet(grc-> transcript_sequences, seq_name);
		if(had_tab){
			SUBREADprintf("Error: duplicate sequence names were found in the input: '%s'.\n", seq_name);
			return -1;
		}

		grc_put_new_trans(grc, seq_name, lbuf, this_seq_len, &linear_space_top);
	}
	
	if(total_dup)SUBREADprintf("Warning: there are %d transcripts that have replicate sequences and the wanted expression levels are non-zero. You may use scanFasta() to find their names.\n", total_dup);
	autozip_close(&auto_FP);
	HashTableDestroy(seq_duplicate_tab);

	if(linear_space_top<1){
		SUBREADprintf("Error: no valid transcript found in the input. No reads can be generated.\n");
		return -1;
	}

	char outname[MAX_FILE_NAME_LENGTH+30];

	sprintf(outname,"%s.truthCounts", grc->output_prefix);
	grc->counts_out_fp = fopen(outname,"w");
	fprintf(grc->counts_out_fp, "## CMD :%s\nTranscriptID\tLength\tCount\n", grc->cmd_line);

	sprintf(outname,"%s_R1.fastq.gz", grc->output_prefix);
	grc->out_fps[0] = gzopen(outname, "wb");

	if(grc->is_paired_end){
		sprintf(outname,"%s_R2.fastq.gz", grc->output_prefix);
		grc->out_fps[1] = gzopen(outname, "wb");
	}else grc->out_fps[1]=NULL;

	return ret;
}

#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int gen_rnaseq_reads_main(int argc, char ** argv)
#endif
{
	int do_fasta_summary = 0;
	int c;
	int option_index = 0;
	genRand_context_t grc;
	memset(&grc,0,sizeof(grc));

	optind = 0;
	opterr = 1;
	optopt = 63;

	rebuild_command_line(&grc.cmd_line, argc, argv);
	// default settings of the read/insertion length: 100bp, a general illumina feeling 
	grc.insertion_length_sigma = 30.;
	grc.insertion_length_min = 110;
	grc.insertion_length_max = 400;
	grc.insertion_length_mean = 160.;
	grc.read_length = 100;

	grc.floor_strategy = STRATEGY_RANDOM_ASSIGN;

	long long seed = -1;

	while ((c = getopt_long (argc, argv, "QO:TCxS:V:N:X:F:L:q:r:t:e:o:pM?", long_options, &option_index)) != -1) {
		switch(c){
			case 'Q':
				grc.quiet=1;
				break;
			case 'x':
				grc.no_actual_reads = 1;
				break;
			case 'O':
				if(strcmp(optarg,"FLOOR")==0)
					grc.floor_strategy = STRATEGY_LESS_THAN_N;
				if(strcmp(optarg,"ITERATIVE")==0)
					grc.floor_strategy = STRATEGY_ITERATIVE_M;
				if(strcmp(optarg,"RANDOM")==0)
					grc.floor_strategy = STRATEGY_RANDOM_ASSIGN;
				break;
			case 'M':
				do_fasta_summary = 1;
				break;
			case 'V':
				grc.insertion_length_sigma = atof(optarg);
				break;
			case 'N':
				grc.insertion_length_min = atoi(optarg);
				break;
			case 'X':
				grc.insertion_length_max = atoi(optarg);
				break;
			case 'F':
				grc.insertion_length_mean = atof(optarg);
				break;
			case 'L':
				grc.read_length = atoi(optarg);
				break;
			case 'p':
				grc.is_paired_end = 1;
				break;
			case 'o':
				strcpy(grc.output_prefix, optarg);
				break;
			case 'e':
				strcpy(grc.expression_level_file, optarg);
				break;
			case 't':
				strcpy(grc.transcript_fasta_file, optarg);
				break;
			case 'T':
				grc.truth_in_read_names=1;
				break;
			case 'C':
				grc.simple_transcript_names = 1;
				break;
			case 'r':
				grc.output_sample_size = atoll(optarg);
				break;
			case 'S':
				seed = atoll(optarg);
				assert(seed>=0);
				break;
			case 'q':
				strcpy(grc.quality_string_file, optarg);
				break;
			default:
			case '?':
				print_usage_gen_reads(argv[0]);
				return 0;
		}
	} 

	#ifdef MAKE_STANDALONE
	if(seed<0){
		double timemil = miltime();
		memcpy(&seed, &timemil, sizeof(seed));
	}
	#else
	{	// in R, we use the random seed generated by R's RNG to initialize the myrand module, then use four random numbers from it to initialize the random_seed inside this program.
		// This guarantees that when R's RNG was seeded by the same seed, the behaviour of the read generator is fixed.
		myrand_srand(0);// the seed is ignored in R so give it zero.
		seed = 0;
		seed = (seed<<16) | (myrand_rand()&0xffff);
		seed = (seed<<16) | (myrand_rand()&0xffff);
		seed = (seed<<16) | (myrand_rand()&0xffff);
		seed = (seed<<16) | (myrand_rand()&0xffff);
	}
	#endif
	memcpy(grc.random_seeds, &seed, sizeof(seed));


	int ret = 0;
	if(do_fasta_summary){
		ret = grc_summary_fasta(&grc);

		if(ret && strlen(grc.output_prefix)>0){
			char delfn[30+MAX_FILE_NAME_LENGTH];
			sprintf(delfn, "%s.faSummary", grc.output_prefix);
			//SUBREADprintf("UNLINK: %s\n", delfn);
			unlink(delfn);
		}
		free(grc.cmd_line);
	}else{
		int ret = grc_check_parameters(&grc) && print_usage_gen_reads(argv[0]);
		ret =  ret || grc_load_env(&grc);
		ret =  ret || grc_gen(&grc);
		ret =  ret || grc_finalize(&grc);
	}


	return ret;
}
