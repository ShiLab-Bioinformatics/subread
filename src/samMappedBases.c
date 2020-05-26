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
#include <string.h>
#include <zlib.h>
#include <getopt.h>
#include <ctype.h>
#include "subread.h"
#include "core.h"
#include "sambam-file.h"
#include "HelperFunctions.h"
#include "input-files.h"

typedef struct{
	int is_BAM;
	
}countbases_context_t;

static struct option sumb_long_options[] =
{
	{"BAMinput",  no_argument, 0, '9'},
	{"SAMinput",  no_argument, 0, '8'},
	{0, 0, 0, 0}
};

void countBases(char * fn, countbases_context_t * context)
{
	char fline[2999];
	unsigned long long int allbases = 0;
	SamBam_FILE * in_fp = SamBam_fopen(fn,context->is_BAM?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM); 
	
	while(1)
	{
		char * tok_tmp = NULL;
		char * is_ret = SamBam_fgets(in_fp, fline, 2999, 0);
		if(!is_ret) break;
		if('@' == fline[0]) continue;

		strtok_r(fline, "\t", &tok_tmp);
		char * flags_str = strtok_r(NULL, "\t", &tok_tmp);
		strtok_r(NULL, "\t", &tok_tmp);
		strtok_r(NULL, "\t", &tok_tmp);
		strtok_r(NULL, "\t", &tok_tmp);
		char * cigar = strtok_r(NULL, "\t", &tok_tmp);

		int flags = atoi(flags_str);

		if(4 & flags) continue;

		
		unsigned int Staring_Points[6];
		unsigned short Section_Length[6];

		int i, retv = RSubread_parse_CIGAR_string(cigar, Staring_Points, Section_Length);

		for(i=0;i<retv;i++)
			allbases += Section_Length[i];
		
	}
	SamBam_fclose(in_fp);

	SUBREADprintf("%s\t%llu\n", fn, allbases);
}

#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int main_mappedBases(int argc, char ** argv)
#endif
{

	int ret = 0;
	int c;
	int option_index = 0 ;
	countbases_context_t * context = calloc(sizeof(countbases_context_t),1);

	optind = 0;
	opterr = 1;
	optopt = 63;


	while((c = getopt_long (argc, argv, "98", sumb_long_options, &option_index)) != -1)
	{
		switch(c){
			case '9':
				context -> is_BAM = 1;
				break;
			case '8':
				context -> is_BAM = 0;
				break;
			default:
				return 0;
		}
	}

	int input_file;
	for(input_file = optind; input_file < argc; input_file++)
	{
		countBases(argv[input_file], context);
	}

	free(context);
	return ret;
}

