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
#include <signal.h>
#include <assert.h>
#include <dirent.h>
#include <getopt.h>
#include <sys/types.h>
#ifndef __MINGW32__
#include <sys/resource.h>
#endif
#include <unistd.h>
#include <sys/stat.h>
#include "subread.h" 
#include "core.h" 
#include "hashtable.h" 
#include "sambam-file.h" 
#include "input-files.h" 
#include "gene-algorithms.h"
#include "HelperFunctions.h"


typedef struct {
	char input_file_name [MAX_FILE_NAME_LENGTH];
	char output_file_name [MAX_FILE_NAME_LENGTH];
	char temp_file_prefix [MAX_FILE_NAME_LENGTH];
	int is_BAM_input;

	unsigned long long all_reads;
	unsigned long long all_records;
	unsigned long long mapped_reads;

	int is_fragments_counted;
	int is_proppair_needed;
	int sort_buckets;
	int verbose;

	HashTable * split_fp_table;

} propMapped_context;





char * _PROPMAPPED_delete_tmp_prefix = NULL;
void PROPMAPPED_SIGINT_hook(int param)
{
	#ifdef MAKE_STANDALONE
	int xk1, last_slash = -1;
	if(_PROPMAPPED_delete_tmp_prefix != NULL)
	{
		char del2[MAX_FILE_NAME_LENGTH], del_suffix[MAX_FILE_NAME_LENGTH], del_name[MAX_FILE_NAME_LENGTH];
		SUBREADprintf("\n\nReceived a terminal signal. The temporary files were removed.\n");
		for(xk1=0; _PROPMAPPED_delete_tmp_prefix[xk1]; xk1++)
		{
			if(_PROPMAPPED_delete_tmp_prefix[xk1]=='/') last_slash = xk1;
			else if(_PROPMAPPED_delete_tmp_prefix[xk1]=='\\')
			{
				SUBREADprintf("The file name is unknown.\n");
				return;
			}
		}
		if(last_slash>=0)
		{
			memcpy(del2, _PROPMAPPED_delete_tmp_prefix, last_slash);
			del2[last_slash] = 0;
			strcpy(del_suffix , _PROPMAPPED_delete_tmp_prefix + last_slash + 1);
		}
		else
		{
			strcpy(del2,".");
			strcpy(del_suffix , _PROPMAPPED_delete_tmp_prefix);
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

	exit(param);
	#endif
}







int propMapped(propMapped_context * context) 
{

	SamBam_FILE * in_fp = SamBam_fopen(context -> input_file_name  ,context -> is_BAM_input?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM); 
	if(!in_fp)
	{
		SUBREADprintf("Unable to open file '%s'.\n", context -> input_file_name);
		return -1;
	}

	char * line_buffer = malloc(3000);

	while(1)
	{
		char * is_ret = SamBam_fgets(in_fp, line_buffer, 2999, 1);
		if(!is_ret) break;
		if(line_buffer[0]=='@') continue;

		char * tok_tmp;
		strtok_r (line_buffer,"\t", &tok_tmp);	// read name
		char * flags_str = strtok_r (NULL ,"\t", &tok_tmp);
		int flags = atoi(flags_str);

		context -> all_reads++;
		if( (flags & 4 ) == 0) 
			context -> mapped_reads ++;
	}


	SamBam_fclose(in_fp);
	return 0;
}

int write_result(propMapped_context * context)
{
	if(context -> output_file_name[0])
	{
		FILE * outfp = f_subr_open(context -> output_file_name, "a");
		#ifdef __MINGW32__
		fprintf(outfp, "%s,%I64u,%I64u,%f\n", context -> input_file_name, context -> all_reads, context -> mapped_reads,  context -> mapped_reads*1./context -> all_reads);
		#else
		fprintf(outfp, "%s,%llu,%llu,%f\n", context -> input_file_name, context -> all_reads, context -> mapped_reads,  context -> mapped_reads*1./context -> all_reads);
		#endif
		fclose(outfp);
	}
	char * objname = context -> is_fragments_counted? "fragment":"read";
	#ifdef __MINGW32__
	if(context -> verbose) SUBREADprintf("Finished. All records: %I64u; all %ss: %I64u; mapped %ss: %I64u; the mappability is %.2f%%\n", context->all_records, objname, context -> all_reads, objname, context -> mapped_reads, context -> mapped_reads*100./context -> all_reads);
	#else
	if(context -> verbose) SUBREADprintf("Finished. All records: %llu; all %ss: %llu; mapped %ss: %llu; the mappability is %.2f%%\n", context->all_records, objname, context -> all_reads, objname, context -> mapped_reads, context -> mapped_reads*100./context -> all_reads);
	#endif
	return 0;
}

static struct option propm_long_options[] =
{
	{"in",  required_argument, 0, 'i'},
	{"out",  required_argument, 0, 'o'},
	{0, 0, 0, 0}
};

void propMapped_usage()
{
	SUBREADprintf("\npropMapped v%s\n", SUBREAD_VERSION);
	SUBREADputs("");
	SUBREADputs("  Calculate the proportion of mapped reads/fragments.");
	SUBREADputs("");
	SUBREADputs("Usage:");
	SUBREADputs("");
	SUBREADputs("  ./prommapped [options] -i <file>");
	SUBREADputs("");
	SUBREADputs("Required arguments:");
	SUBREADputs("");
	SUBREADputs("  -i <string>  An input file containing read mapping results. Both SAM or BAM");
	SUBREADputs("               formats are supported.");
	SUBREADputs("");
	SUBREADputs("Optional arguments:");
	SUBREADputs("");
	SUBREADputs("  -o <string>  Name of output file including mapping statistics.");
	SUBREADputs("");
	SUBREADputs("  -f           If specified, fragments (read pairs) will be counted instead of");
	SUBREADputs("               individual reads. This option is only applicable for paired-end");
	SUBREADputs("               reads.");
	SUBREADputs("");
	SUBREADputs("  -p           If specified, only properly paired reads will be counted. This");
	SUBREADputs("               option is only applicable for paired-end reads.");
	SUBREADputs("");

}

int finalise_PE_split(propMapped_context * context)
{
	int bucket;
	KeyValuePair *cursor;

	for(bucket=0; bucket < context -> split_fp_table -> numOfBuckets; bucket++)
	{
		cursor = context -> split_fp_table -> bucketArray[bucket];
		while (1)
		{
			if (!cursor) break;
			FILE * fp = (FILE *) cursor -> value;
			fclose(fp);
			cursor = cursor->next;
		}
	}

	HashTableDestroy(context -> split_fp_table);
	
	return 0;
}

FILE * get_FP_by_read_name(propMapped_context * context, char * read_name)
{
	unsigned long long int hash_key = sort_SAM_hash(read_name);
	hash_key = hash_key % context -> sort_buckets;

	FILE * ret = HashTableGet(context -> split_fp_table, NULL+1+(int)(hash_key));
	if(!ret)
	{
		char fname [MAX_FILE_NAME_LENGTH+40];
		#ifdef __MINGW32__
		sprintf(fname, "%s-%I64u.bin", context->temp_file_prefix, hash_key);
		#else
		sprintf(fname, "%s-%llu.bin", context->temp_file_prefix, hash_key);
		#endif
		ret = f_subr_open(fname, "wb");
		HashTablePut(context -> split_fp_table, NULL+1+(int)(hash_key), ret);
	}

	return ret;
	
}

int add_read_flags(propMapped_context * context, FILE * fp, char * read_name, unsigned short flags)
{
	int x1;
	int rname_len = strlen(read_name);

	for(x1=rname_len-1; x1>0;x1--)
	{
		if(read_name[x1]=='/'){
			read_name[x1]=0;
			break;
		}
	}


	if(1&flags) if(!context->is_fragments_counted)
			{
				strcat(read_name, (flags&0x40)?"!!_1":"!!_2");
			}

	rname_len = strlen(read_name);
	if(rname_len>250)
		return -1;


	unsigned char rname_len_char = (unsigned char)rname_len;
	int wlen = fwrite(&rname_len_char,1,1,fp);
	if(wlen < 1) return -1;
	wlen = fwrite(read_name,rname_len,1,fp);
	if(wlen < 1) return -1;
	wlen = fwrite(&flags, sizeof(short), 1, fp);
	if(wlen < 1) return -1;
	return 0;
}

int init_PE_sambam(propMapped_context * context)
{
	char mac_rand[13];
	mac_or_rand_str(mac_rand);
	#ifdef MAKE_STANDALONE
	myrand_srand(time(NULL));
	#endif

	int x1;
	context->temp_file_prefix[0]=0;
	for(x1 = strlen(context->output_file_name) ; x1 >=0; x1--){
		if(context->output_file_name[x1]=='/'){
			memcpy(context->temp_file_prefix, context->output_file_name, x1);
			context->temp_file_prefix[x1] = 0;
		}
	}

	if(context->temp_file_prefix[0]==0) strcpy(context->temp_file_prefix, "./");

	sprintf(context->temp_file_prefix+strlen(context->temp_file_prefix), "/prpm-temp-sum-%06u-%s", getpid(), mac_rand);

	_PROPMAPPED_delete_tmp_prefix = context -> temp_file_prefix;
	signal (SIGTERM, PROPMAPPED_SIGINT_hook);
	signal (SIGINT, PROPMAPPED_SIGINT_hook);
	context -> split_fp_table = HashTableCreate(500);
	return 0;
}

int split_PE_sambam(propMapped_context * context)
{
	SamBam_FILE * in_fp = SamBam_fopen(context -> input_file_name  ,context -> is_BAM_input?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM); 
	if(!in_fp)
	{
		SUBREADprintf("Unable to open file '%s'.\n", context -> input_file_name);
		return -1;
	}

	char * line_buffer = malloc(3000);
	int is_error = 0;

	while(1)
	{
		char * is_ret = SamBam_fgets(in_fp, line_buffer, 2999, 1);
		if(!is_ret) break;
		if(line_buffer[0]=='@') continue;

		//SUBREADprintf("TELLO: %ld\n", ftello(in_fp -> os_file));
		char * tok_tmp;
		char * read_name = strtok_r (line_buffer,"\t", &tok_tmp);	// read name
		char * flags_str = strtok_r (NULL ,"\t", &tok_tmp);
		unsigned flags = atoi(flags_str);
		
		FILE * fp = get_FP_by_read_name(context , read_name);
		is_error = add_read_flags(context, fp, read_name, flags);
		if(is_error) break;
		context -> all_records++;
	}


	free(line_buffer);
	SamBam_fclose(in_fp);

	if(is_error) SUBREADprintf("ERROR: Unable to write into the temporary file. Please check the disk space in the output directory.");
	return is_error;
}

int prop_PE(propMapped_context * context)
{
	int bini;
	for(bini = 0; bini < context -> sort_buckets; bini++)
	{
		char fname[MAX_FILE_NAME_LENGTH+25];
		HashTable * rname_table = HashTableCreate(100000);
		HashTableSetKeyComparisonFunction(rname_table , fc_strcmp_chro);
		HashTableSetHashFunction(rname_table, fc_chro_hash);
		HashTableSetDeallocationFunctions(rname_table , free, NULL);

		sprintf(fname, "%s-%d.bin", context->temp_file_prefix, bini);
		FILE * fp = f_subr_open(fname, "rb");
		if(fp)
		{
				while(!feof(fp))
				{
					unsigned char read_len;
					unsigned short flags;
					char * read_name;
					int retvv = fread(&read_len,1,1,fp);
					if(retvv<1) break;
					if(feof(fp))break;
					read_name = malloc((int)read_len+1);
					retvv = fread(read_name,read_len,1,fp);
					if(retvv<1) break;
					read_name[read_len]=0;
					retvv = fread(&flags,1, sizeof(short), fp);
					if(retvv<1) break;

					int new_OK; 
					
					if(context -> is_fragments_counted){
						if(context -> is_proppair_needed )
							new_OK = (flags & 2)>0;
						else	new_OK=1;

						if(new_OK)
							new_OK = (!(flags & 4)) || ( (flags & 1) && !(flags&8) );
					}
					else
					{
						new_OK = ((flags & 4) == 0);
					}
					int old_OK = HashTableGet(rname_table, read_name) - NULL;
					if(old_OK)
					{
						old_OK--;

						// if the old read is unmapped or the old read is not properly mapped pair.
						if(new_OK && !old_OK )
							HashTablePut(rname_table, read_name, NULL+2);
					}
					else
						HashTablePut(rname_table, read_name, NULL + new_OK + 1);
				}

				fclose(fp);
				unlink(fname);


				int bucket;
				KeyValuePair *cursor;

				for(bucket=0; bucket < rname_table -> numOfBuckets; bucket++)
				{
					cursor = rname_table -> bucketArray[bucket];
					while (1)
					{
						if (!cursor) break;
						if(cursor -> value == NULL+2) context -> mapped_reads++;
						context -> all_reads++;
						cursor = cursor->next;
					}
				}
		}
		//SUBREADprintf("BIN %d  ALLR %llu  MAPPED %llu  ELES=%ld\n", bini, context -> all_reads, context -> mapped_reads, rname_table -> numOfElements);
		HashTableDestroy(rname_table);
	}
	return 0;
}

void ppm_warning_file_limit()
{
#ifndef __MINGW32__
	struct rlimit limit_st;
	getrlimit(RLIMIT_NOFILE, & limit_st);

	{
		if(min(limit_st.rlim_cur , limit_st.rlim_max) < 400) {
			SUBREADprintf("Your operation system does not allow a single process to open more then 400 files. You may need to change this setting by using a 'ulimit -n 500' command, or the program may crash.\n");
		}
	}
#endif
}



#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int propmapped(int argc, char ** argv)
#endif
{
	int ret = 0;
	int c;
	int option_index = 0;
	_PROPMAPPED_delete_tmp_prefix = NULL;

	propMapped_context * context = malloc(sizeof(propMapped_context));
	memset(context, 0 , sizeof(propMapped_context));

	optind = 0;
	opterr = 1;
	optopt = 63;

	context -> sort_buckets = 253;


	while((c = getopt_long (argc, argv, "Vi:o:bfph", propm_long_options, &option_index)) != -1)
	{
		switch(c){
			case 'i':
				strcpy(context -> input_file_name, optarg);
				break;
			case 'o':
				strcpy(context -> output_file_name, optarg);
				break;
			case 'f':
				context -> is_fragments_counted = 1;
				break;
			case 'V':
				context -> verbose = 1;
			case 'p':
				context -> is_proppair_needed = 1;
				break;
			case 0:	// long names
				break;
			default:
				propMapped_usage();
				return 0;

		}
	}		

	if(!context -> input_file_name[0])
	{
		propMapped_usage();
		return 0;
	}

	int is_bam = is_certainly_bam_file(context -> input_file_name, NULL, NULL);

	if(1==is_bam){
		context -> is_BAM_input = 1;
	}
	else if(is_bam < 0)
	{
		ret = -1;
		SUBREADprintf("Unable to open input file '%s' or the input file is empty.\n", context -> input_file_name);
	}
	if(context -> verbose) SUBREADprintf("The input file is opened as a %cAM file.\nThe %ss in the input file are being counted.\n", context -> is_BAM_input?'B':'S', context -> is_fragments_counted?"fragment":"read");

	ppm_warning_file_limit ();
	ret = ret || init_PE_sambam(context);
	ret = ret || split_PE_sambam(context);
	ret = ret || finalise_PE_split(context);
	ret = ret || prop_PE(context);
	ret = ret || write_result(context);

	free(context);

	return ret;
}


