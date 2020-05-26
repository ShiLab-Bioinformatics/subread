#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "subread.h" 
#include "core.h" 
#include "sambam-file.h" 
#include "input-files.h" 


static struct option long_options[] =
{
	{"in",  required_argument, 0, 'i'},
	{"filter",  required_argument, 0, 'F'},
	{0, 0, 0, 0}
};

int main(int argc, char ** argv)
{
	int c;
	int option_index = 0;	

	optind = 1;
	opterr = 1;
	optopt = 63;

	int sort_needed = 0;
	char filter_mode[10];
	char in_name[MAX_FILE_NAME_LENGTH];
	int in_SAM = 1;
	int out_SAM = 1;

	in_name[0] = filter_mode[0]=0;


	while ((c = getopt_long (argc, argv, "i:F:", long_options, &option_index)) != -1)
	{
		switch(c){

			case 'F':
				out_SAM = (strcmp(optarg, "SAM")==0);
				break;
			case 'i':
				strcpy(in_name, optarg);
				break;
		}	
	}
}
