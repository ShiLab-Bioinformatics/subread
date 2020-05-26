#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <zlib.h>
#include "subread.h"
#include "core.h"
#include "input-files.h"
#include "HelperFunctions.h"

void print_usage_pairer(char * cmd){
	SUBREADprintf("\nrepair Version %s\n\n", SUBREAD_VERSION);
	SUBREADputs("  Find reads that are from the same pair in the input and then place them next");
	SUBREADputs("to each other in the output. A dummy read is added for each singleton read");
	SUBREADputs("that does not have a pair. The output file is compatible with featureCounts");
	SUBREADputs("program.");
	SUBREADputs("");
	SUBREADputs("Usage:");
	SUBREADputs("");
	SUBREADputs("  ./repair [options] -i <input_file> -o <output_file>\n");
	SUBREADputs("");
	SUBREADputs("Required arguments:");
	SUBREADputs("");
	SUBREADputs("  -i <string>  Name of input file. BAM format by default.");
	SUBREADputs("");
	SUBREADputs("  -o <string>  Name of output file. The output file is in BAM format.");
	SUBREADputs("");
	SUBREADputs("Optional arguments:");
	SUBREADputs("");
	SUBREADputs("  -S           The input file is in SAM format.");
	SUBREADputs("");
	SUBREADputs("  -c           Compress the output BAM file. This will reduce the size of BAM");
	SUBREADputs("               file, but will increase the time of retrieving reads from BAM");
	SUBREADputs("               file.");
	SUBREADputs("");
	SUBREADputs("  -T <int>     Number of CPU threads. 8 by default.");
	SUBREADputs("");
	SUBREADputs("  -d           Do not add dummy reads for singleton reads.");
	SUBREADputs("");
	SUBREADputs("  -t           Do not include sequences and quality scores of reads in the");
	SUBREADputs("               output file.");
	SUBREADputs("");
}

#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int main_read_repair(int argc, char ** argv)
#endif
{
	double t0 = miltime();
	int threads = 8, is_BAM = 1;
	char c;
	char in_BAM_file[MAX_FILE_NAME_LENGTH+1];
	char out_BAM_file[MAX_FILE_NAME_LENGTH+1];
	char rand_prefix[40];
	int no_compression = 1;
	int has_dummy = 1;
	int tiny_mode = 0;
	optind = 1;
	opterr = 1;
	optopt = 63;
	int memory = 64;
	in_BAM_file[0] = out_BAM_file[0] = 0;

	while ((c = getopt(argc, argv, "i:T:M:o:vtdcS?")) != -1)
	{
		switch(c)
		{
			case '?':
			case 'v':
				print_usage_pairer(argv[0]);
				STANDALONE_exit(0);
				break;
			case 'S':
				is_BAM = 0;
				break;
			case 't':
				tiny_mode = 1;
				break;
			case 'd':
				has_dummy = 0;
				break;
			case 'o':
				strcpy(out_BAM_file, optarg);
				break;
			case 'M':
				memory = atoi(optarg);
				if(memory < 1) memory = 1;
				break;
			case 'T':
				threads = atoi(optarg);
				if(threads < 1) threads = 1;
				if(threads > MAX_THREADS) threads = MAX_THREADS;
				break;
			case 'c':
				no_compression = 0;
				break;
			case 'i':
				strcpy(in_BAM_file, optarg);
				break;
		}
	}

	if(in_BAM_file[0]==0 || out_BAM_file[0]==0){
		print_usage_pairer(argv[0]);
		//SUBREADprintf("\nNo input or output files are specified.\n");
		STANDALONE_exit(-1);
	}

	if(!is_paired_end_BAM(in_BAM_file)){
		SUBREADprintf("Error: the input file contains single-end reads. Repair can only process paired-end reads.\n");
		STANDALONE_exit(-1);
	}

	char mac_rand[13];
	mac_or_rand_str(mac_rand);

	sprintf(rand_prefix, "fsbm-p%06d-%s", getpid(), mac_rand);

	SAM_pairer_context_t pairer;
	SAM_pairer_writer_main_t writer_main;
	int ret = SAM_pairer_writer_create(&writer_main, threads, has_dummy,1, no_compression?Z_NO_COMPRESSION:Z_DEFAULT_COMPRESSION, out_BAM_file);
	if(ret){
		SUBREADprintf("Unable to open the output file. Program terminated.\n");
		return -1;
	}else{
		ret = SAM_pairer_create(&pairer, threads, memory, is_BAM, tiny_mode,0,0 , 0, 1, in_BAM_file, SAM_pairer_writer_reset, SAM_pairer_multi_thread_header, SAM_pairer_multi_thread_output, rand_prefix, &writer_main, 99999999);
		if(ret){
			SUBREADprintf("Unable to open the input file. Program terminated.\n");
			return -1;
		}else{
			SAM_pairer_run(&pairer);
			int has_error = pairer.is_bad_format;
			SAM_pairer_destroy(&pairer);
			SAM_pairer_writer_destroy(&writer_main);
			SUBREADprintf("\n%s %.2f minutes\nTotal input reads: %llu ; Unpaired reads: %llu\n\n", has_error?"Program terminated WITH ERRORS!!! Used":"All finished in", (miltime()-t0)/60, pairer.total_input_reads, pairer.total_orphan_reads);
			if(has_error){
				SUBREADprintf("No output file was generated!!!\n");
				unlink(out_BAM_file);
				return -1;
			}
			return 0;
		}
	}
}
