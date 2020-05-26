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
	{"out",  required_argument, 0, 'o'},
	{"informat",  required_argument, 0, 'f'},
	{"outformat",  required_argument, 0, 'F'},
	{"sort",  required_argument, 0, 'S'},
	{0, 0, 0, 0}
};

void subtools_usage()
{
	SUBREADprintf("\nVersion %s\n\n", SUBREAD_VERSION);
	SUBREADputs("Usage:");
	SUBREADputs("");
	SUBREADputs(" ./subtools -i input_file -o output_file --informat <SAM|BAM> --outformat <SAM|BAM>  {--sort <byname>}");
	SUBREADputs("");
}

int main(int argc, char ** argv)
{
	int c;
	int option_index = 0;	

	optind = 1;
	opterr = 1;
	optopt = 63;

	int sort_needed = 0;
	char in_name[MAX_FILE_NAME_LENGTH];
	char out_name[MAX_FILE_NAME_LENGTH];
	int in_SAM = 1;
	int out_SAM = 1;

	in_name[0] = out_name[0]=0;

	while ((c = getopt_long (argc, argv, "i:o:f:F:S:", long_options, &option_index)) != -1)
	{
		switch(c){

			case 'S':
				sort_needed = 1;
				break;
			case 'f':
				in_SAM = (strcmp(optarg, "SAM")==0);
				break;
			case 'F':
				out_SAM = (strcmp(optarg, "SAM")==0);
				break;
			case 'i':
				strcpy(in_name, optarg);
				break;
			case 'o':
				strcpy(out_name, optarg);
				break;

			default:
				subtools_usage();
				return 0;
		}
	}


	if((!out_name[0])||(!in_name[0]))
	{
		subtools_usage();
		return 0;
	}

	if(in_SAM && out_SAM && !sort_needed)
	{
		SUBREADprintf("Subtools can convert between SAM and BAM files and/or pair up the alignment results by the read names.\n");
		SUBREADprintf("No operation was specified. The output file was not generated.\n");
		return 0;
	}

	if(strcmp(out_name, in_name)==0)
	{
		SUBREADprintf("ERROR: the input file and the output file are the same!\n");
		return 0;
	}

	SamBam_FILE * in_fp = SamBam_fopen(in_name,in_SAM?SAMBAM_FILE_SAM:SAMBAM_FILE_BAM); 
	if(!in_fp)
	{
		SUBREADprintf("Unable to open the input file: %s\n",in_name);
		return 1;
	}

	FILE * out_fp = NULL;
	SamBam_Writer out_writer;

	if(out_SAM)
	{
		out_fp = f_subr_open(out_name, "w");
	

		if(!out_fp)
		{
			SUBREADprintf("Unable to open the output file: %s\n",out_name);
			SamBam_fclose(in_fp);
			return 1;
		}
	}
	else
	{

		int swret =  SamBam_writer_create( & out_writer, out_name);
		if(swret)
		{
			SUBREADprintf("Unable to open the output file: %s\n",out_name);
			SamBam_fclose(in_fp);
			return 1;
		}
	}


	if(sort_needed)
	{	
		char fline[3000], temp_file_name[MAX_FILE_NAME_LENGTH], mac_rand[13];

		SAM_sort_writer writer;
		mac_or_rand_str(mac_rand);
		sprintf(temp_file_name, "./temp-subt-%06u-%s.sam", getpid(), mac_rand);
		int ret  =sort_SAM_create(&writer, temp_file_name, ".");
		if(ret)
		{
			SUBREADprintf("ERROR: temporary file '%s' is not able to be created.\n", temp_file_name);
		}

		unsigned long long int added_lines = 0;

		double t0 = miltime();
		while(1)
		{
			char * is_ret = SamBam_fgets(in_fp, fline, 2999, 0);
			if(!is_ret) break;
			int ret = sort_SAM_add_line(&writer, fline, strlen(fline));
			if(ret<0) 
			{
				SUBREADprintf("ERROR: read name is too long; check the input format.\n");
				break;
			}

			added_lines++;
			//if(added_lines>40000000) break;
		}
			//printf("N1=%llu\n",  writer.unpaired_reads);

		double t1 = miltime();
		sort_SAM_finalise(&writer);

		double t2 = miltime();

		SUBREADprintf("Loading time: %.2f, Sorting time: %.2f\n", t1-t0, t2-t1);

		if(writer.unpaired_reads && 0)
		SUBREADprintf("%llu reads were paired up.\n", writer.written_reads);

		SamBam_fclose(in_fp);

		strcpy(in_name, temp_file_name);
		in_SAM = 1;
		in_fp = SamBam_fopen(in_name,SAMBAM_FILE_SAM); 
	
	}

	int is_header = 1;

	while(1){
		char in_buff[3000];

		if(!SamBam_fgets(in_fp, in_buff, 2999, 1)) break;
		if(out_SAM)
			fputs(in_buff, out_fp);
		else {
			int tail_pos = strlen(in_buff)-1;

			if(in_buff[tail_pos] == '\n') in_buff[tail_pos] =0;
			if(in_buff[0]=='@')
				SamBam_writer_add_header(&out_writer, in_buff, 1);
			else
			{
				char * rname = NULL, *val_str = NULL;
				int flags = -1;
				char * chro = NULL;
				int pos = -1;
				int mapq = -1;
				char * cigar = NULL;
				char * mate_chro = NULL;
				int mate_pos = -1;
				int tlen = -1;
				char * read_text = NULL;
				char * qual_text = NULL;
				char * extra = NULL;
				char * toktmp;

				if(is_header && !out_SAM)
				{
					if(out_writer.chromosome_name_table->numOfElements<1)
						SUBREADprintf("WARNING: no chromosome was found from the header.\nThe BAM output will be useless.\n");
				}

				is_header= 0;

				rname = strtok_r(in_buff, "\t", &toktmp);
				val_str = strtok_r(NULL, "\t", &toktmp);
				if(val_str) flags = atoi(val_str);
				chro = strtok_r(NULL, "\t", &toktmp);
				val_str = strtok_r(NULL, "\t", &toktmp);
				if(val_str) pos = atoi(val_str);
				val_str = strtok_r(NULL, "\t", &toktmp);
				if(val_str) mapq = atoi(val_str);
				cigar = strtok_r(NULL, "\t", &toktmp);

				mate_chro = strtok_r(NULL, "\t", &toktmp);
				val_str = strtok_r(NULL, "\t", &toktmp);
				if(val_str) mate_pos = atoi(val_str);
				val_str = strtok_r(NULL, "\t", &toktmp);
				if(val_str) tlen = atoi(val_str);
				read_text = strtok_r(NULL, "\t", &toktmp);
				qual_text = strtok_r(NULL, "\t", &toktmp);
				extra = toktmp;


				if(4 & flags)
				{
					if(chro[0]!='*' && mate_chro[0]=='=')
						mate_chro = chro;

					chro = "*";
					pos = 0;
					tlen = 0;
				}
				if(8 & flags)
				{
					mate_chro = "*";
					mate_pos = 0;
					tlen = 0;
				}

				SamBam_writer_add_read(& out_writer, rname, flags, chro, pos, mapq, cigar, mate_chro, mate_pos, tlen, strlen(read_text), read_text, qual_text, extra);
			}
		}
	}


	if(out_SAM)
		fclose(out_fp);
	else
		SamBam_writer_close(&out_writer);

	SamBam_fclose(in_fp);
	if(sort_needed && (memcmp(in_name, "./temp-subt", 11)==0)) unlink(in_name);


	return 0;
}
