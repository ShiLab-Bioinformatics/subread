#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>

#include "subread.h"
#include "hashtable.h"
#include "core.h"
#include "HelperFunctions.h"


void parse_line(char * new_line, int * flags, char * chro_name, unsigned int * chro_pos, char * cigar)
{
	int tab_no = 0, read_cursor;
	int curr_line_len = strlen(new_line);
	int field_pnt = 0;
	chro_name[0]=0;
	cigar[0]=0;
	for(read_cursor=0; read_cursor<curr_line_len; read_cursor++)
	{
		char nch = new_line[read_cursor];
		if(nch == '\t')
		{
			tab_no++;
			if(tab_no > 6)break;
			field_pnt=0;
		}
		else
		{
			if(tab_no == 1)	// flags 
				(*flags) = (*flags) * 10 + (nch - '0');
			else if(tab_no == 2)
			{
				chro_name[field_pnt++]=nch;
				chro_name[field_pnt]=0;
			}
			else if(tab_no == 3)
				(*chro_pos) = (*chro_pos)*10 + (nch - '0');
			else if(tab_no == 5)
			{
				cigar[field_pnt++]=nch;
				cigar[field_pnt]=0;
			}

		}
	}
}

int load_junc_table(char * file_name, HashTable * tab, HashTable * edge_table_l, HashTable * edge_table_r, HashTable * bin_table)
{
	FILE * fp=f_subr_open(file_name, "r");
	char new_fl[200];
	if(!fp) return -1;
	while(1)
	{
		char * ll = fgets(new_fl, 199, fp);
		char * tmp_tok=NULL;
		if(!ll || strlen(ll)<4) break;

		strtok_r(ll, "\t", &tmp_tok);    //name
		char * chrostr = strtok_r(NULL, "\t", &tmp_tok);  //chro
		char * pos1str = strtok_r(NULL, "\t", &tmp_tok);  //pos1
		char * pos2str = strtok_r(NULL, "\t", &tmp_tok);  //pos2

		unsigned int pos1 = (unsigned int)atoll(pos1str); 
		unsigned int pos2 = (unsigned int)atoll(pos2str); 

		if(pos1==0||pos2==0)continue;

		char * chro_mem = malloc(30);
		chro_mem[0]=0;
		if(strlen(chrostr)<3)strcpy(chro_mem, "chr");
		strcat(chro_mem, chrostr);

		long long int pos_key;
		if(pos2<pos1){
			HashTablePut(edge_table_l, NULL+pos2, chro_mem);
			HashTablePut(edge_table_r, NULL+pos1, chro_mem);
			pos_key = (pos2*1LLU<<32)|pos1;
		}
		else{
			HashTablePut(edge_table_l, NULL+pos1, chro_mem);
			HashTablePut(edge_table_r, NULL+pos2, chro_mem);
			pos_key = (pos1*1LLU<<32)|pos2;
		}

		HashTablePut(bin_table, NULL+pos1/50,NULL+pos1);
		HashTablePut(bin_table, NULL+pos2/50,NULL+pos2);
		HashTablePut(bin_table, NULL+pos1/50-1,NULL+pos1);
		HashTablePut(bin_table, NULL+pos2/50-1,NULL+pos2);
		HashTablePut(bin_table, NULL+pos1/50+1,NULL+pos1);
		HashTablePut(bin_table, NULL+pos2/50+1,NULL+pos2);

		HashTablePut(tab, (void *)pos_key, chro_mem);

		//printf("NEW JUNC: %s: %u - %u\n", chro_mem, pos1, pos2);
	}
	fclose(fp);
	return 0;
}

int main(int argc, char ** argv)
{

	char * junc_table_file;
	HashTable * junc_table = HashTableCreate(337);
	HashTable * junc_table_l = HashTableCreate(337);
	HashTable * junc_table_r = HashTableCreate(337);
	HashTable * bin_table = HashTableCreate(9337);
	if(argc<2)
	{
		printf("Usage cat my.SAM | filterJunc junc.table | featureCounts -i STDIN -a myAnnot.txt \n");
		return -1;
	}

	junc_table_file = argv[1];
	if(load_junc_table(junc_table_file, junc_table, junc_table_l, junc_table_r, bin_table))
	{
		printf("Junction table not found!\n");
		return -1;
	}
	char * linebuf=malloc(5001);
	char * linebuf2=malloc(5001);
	while(1)
	{
		int support_found = 0, ii, jj, alternative_found = 0;
		char * new_line = fgets(linebuf, 5000, stdin);
		if(!new_line) break;

		if(new_line[0]=='@')
		{
			fputs(new_line, stdout);
			continue;
		}


		char * new_line2 = fgets(linebuf2, 5000, stdin);
		for(ii=0;ii<2;ii++)
		{
			int flag=0;
			char chro[20];
			unsigned int pos=0;
			char cigar[50];
			unsigned int start_points[FC_CIGAR_PARSER_ITEMS];
			unsigned int start_read_points[FC_CIGAR_PARSER_ITEMS];
			unsigned short section_lengths[FC_CIGAR_PARSER_ITEMS];
			
			parse_line(ii?new_line2:new_line, &flag, chro, &pos, cigar);
			int sections = RSubread_parse_CIGAR_string(cigar, start_points, start_read_points, section_lengths);
			for(jj=0; jj<sections-1; jj++)
			{
				unsigned int edge1 = start_points[jj] + pos + section_lengths[jj];
				unsigned int edge2 = start_points[jj+1] + pos - 1;

				unsigned long long int junc_table_key = (edge1*1LLU) <<32 | edge2;
				char * junc_chro=HashTableGet(junc_table, (void *)junc_table_key);

				if(junc_chro && strcmp(junc_chro, chro)==0)
				{
					support_found = 1;
				}
				else
				{
					char * junc_chro_l=HashTableGet(junc_table_l, NULL+edge1); 
					char * junc_chro_2=HashTableGet(junc_table_r, NULL+edge2); 
					if(		(junc_chro_l && strcmp(junc_chro_l, chro)==0)
						||	(junc_chro_2 && strcmp(junc_chro_2, chro)==0))
					alternative_found = 1;
					
				}
			}

			//if(support_found) break;

			if(!support_found && !alternative_found)
				for(jj=0; jj<sections; jj++)
				{
					unsigned int edge2 = start_points[jj] + pos + section_lengths[jj];
					unsigned int edge1 = start_points[jj] + pos;

					unsigned int real_edge1 = HashTableGet(bin_table, NULL+edge1/50)-NULL; 
					unsigned int real_edge2 = HashTableGet(bin_table, NULL+edge2/50)-NULL;

					if(real_edge1)
					{
						long long int diff = real_edge1;
						diff -= edge1;

						if(abs(diff)<section_lengths[jj]) 
							alternative_found=1;
						// this section should contain the read_edge, but it does not! (it was a section as a whole)
					}

					if(real_edge2)
					{
						long long int diff = real_edge2;
						diff -= edge1;

						if(abs(diff)<section_lengths[jj]) 
							alternative_found=1;
						// this section should contain the read_edge, but it does not! (it was a section as a whole)
					}


				}
		}

		if(alternative_found)
		//if(support_found)
		{
			fputs( new_line, stdout);
			fputs( new_line2, stdout);
		}
		//printf("READ: (%d)  %s,%u   %s\n\n", flag2, chro2, pos2, cigar2);

	}

	free(linebuf);
	free(linebuf2);
}
