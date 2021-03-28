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
  
  
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>


#ifdef MACOS

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <sys/sysctl.h>
#include <net/if.h>
#include <net/if_dl.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <stdlib.h>
#include <string.h>
#include <mach/mach.h>
#include <mach/vm_statistics.h>

#else

#ifndef __MINGW32__
#include <sys/ioctl.h>
#include <netinet/in.h>
#include <net/if.h>
#include <sys/sysinfo.h>
#endif
#include <sys/types.h>
#endif


#include "subread.h"
#include "input-files.h"
#include "seek-zlib.h"
#include "gene-algorithms.h"
#include "HelperFunctions.h"

size_t get_sys_mem_info(char * keyword){
	FILE * mfp = fopen("/proc/meminfo","r");
	if(mfp==NULL) return -1;
	char linebuf[1000];
	size_t ret = -1;
	while(1){
		char * rret = fgets(linebuf, 999, mfp);
		if(memcmp( keyword, linebuf, strlen(keyword) ) == 0 && strstr(linebuf," kB")) {
			ret=0;
			int ii ,state=0;
			for(ii=strlen(keyword);; ii++){
				//SUBREADprintf("CH[%d] = %d '%c' at state %d\n", ii, linebuf[ii], linebuf[ii], state);
				if(state == 0 && linebuf[ii]==' ') state = 1;
				if(state == 1 && linebuf[ii]!=' ') state = 2;
				if(state == 2 && linebuf[ii]==' ') state = 9999;

				if(state == 2 && !isdigit(linebuf[ii])){
					SUBREADprintf("WRONG MEMORY INFO '%s'\n", linebuf);
					ret = -1;
					break;
				}
				if(state == 2) ret = ret*10 + ( linebuf[ii] - '0' );
				if(state >= 9999) {
					ret *=1024;
					break;
				}
			}
		}
		if(!rret) break;
	}
	fclose(mfp);
	return ret;
}

char * get_short_fname(char * lname){
	char * ret = lname;

	int x1;
	for(x1 = strlen(lname)-2; x1>=0; x1--){
		if(lname [x1] == '/'||lname [x1] == '\\'){
			ret = lname + x1 + 1;
			break;
		}
	}
	return ret;
}



// This assumes the first part of Cigar has differet strandness to the main part of the cigar.
// Pos is the LAST WANTED BASE location before the first strand jump (split by 'b' or 'n').
// The first base in the read actually has a larger coordinate than Pos. 
// new_cigar has to be at least 100 bytes.
unsigned int reverse_cigar(unsigned int pos, char * cigar, char * new_cigar) {
	int cigar_cursor = 0;
	new_cigar[0]=0;
	unsigned int tmpi=0;
	int last_piece_end = 0;
	int last_sec_start = 0;
	unsigned int chro_pos = pos, this_section_start = pos, ret = pos;
	int is_positive_dir = 0;
	int read_cursor = 0;
	int section_no = 0;

	for(cigar_cursor = 0 ;  ; cigar_cursor++)
	{
		if( cigar [cigar_cursor] == 'n' ||  cigar [cigar_cursor] == 'b' ||  cigar [cigar_cursor] == 0)
		{
			int xk1, jmlen=0, nclen=strlen(new_cigar);
			char jump_mode [13];

			if(cigar [cigar_cursor] !=0)
			{
				sprintf(jump_mode, "%u%c", tmpi,  cigar [cigar_cursor] == 'b'?'n':'b');
				jmlen = strlen(jump_mode);
			}

			for(xk1=nclen-1;xk1>=0; xk1--)
				new_cigar[ xk1 +  last_piece_end + jmlen - last_sec_start ] = new_cigar[ xk1 ];
			new_cigar [nclen + jmlen + last_piece_end - last_sec_start ] = 0;

			memcpy(new_cigar , jump_mode, jmlen);
			memcpy(new_cigar + jmlen , cigar + last_sec_start, last_piece_end - last_sec_start);

			last_sec_start = cigar_cursor+1;

			if(is_positive_dir && cigar [cigar_cursor] !=0)
			{
				if(cigar [cigar_cursor] == 'b') chro_pos -= tmpi - read_cursor - 1;
				else	chro_pos += tmpi - read_cursor - 1;
			}
			if((!is_positive_dir) && cigar [cigar_cursor] !=0)
			{
				if(cigar [cigar_cursor] == 'b') chro_pos = this_section_start - tmpi - read_cursor - 1;
				else	chro_pos = this_section_start + tmpi - read_cursor - 1;
			}

			this_section_start = chro_pos;

			if(section_no == 0)
				ret = chro_pos;

			is_positive_dir = ! is_positive_dir;
			section_no++;
			tmpi=0;
		}
		else if(isalpha(cigar [cigar_cursor]))
		{
			if(cigar [cigar_cursor]=='M' || cigar [cigar_cursor] == 'S')
				read_cursor += tmpi;
			tmpi=0;
			last_piece_end = cigar_cursor+1;
		}
		else tmpi = tmpi*10 + (cigar [cigar_cursor] - '0');

		if(cigar [cigar_cursor] == 0)break;
	}

	SUBREADprintf("REV CIGAR: %s  =>  %s\n", cigar, new_cigar);
	return ret;
}

unsigned int find_left_end_cigar(unsigned int right_pos, char * cigar){
	int delta_from_right = 0;
	int cigar_cursor = 0;
	unsigned int tmpi = 0;
	while(1){
		int nch = cigar[cigar_cursor++];
		if(nch == 0) break;
		if(isdigit(nch)){
			tmpi = tmpi * 10 + nch - '0';
		}else{
			if(nch == 'M'||nch == 'D' || nch == 'N'){
				delta_from_right +=tmpi;
			}
			tmpi = 0;
		}
	}
	return right_pos - delta_from_right;
}


char contig_fasta_int2base(int v){
	if(v == 1) return 'A';
	if(v == 2) return 'T';
	if(v == 3) return 'G';
	if(v == 4) return 'C';
	return 'N';
}

int contig_fasta_base2int(char base){
	base = tolower(base);
	if((base) == 'a'){ return 1;}
	else if((base) == 't' || (base) == 'u'){ return 2;}
	else if((base) == 'g'){ return 3;}
	else if((base) == 'c'){ return 4;}
	else return 15 ;
}

int get_contig_fasta(fasta_contigs_t * tab, char * chro, unsigned int pos, int len, char * out_bases){
	unsigned int this_size = HashTableGet( tab -> size_table, chro ) - NULL;
	if(this_size > 0){
		if(this_size >= len && pos <= this_size - len){
			char * bin_block = HashTableGet(tab -> contig_table, chro );
			unsigned int bin_byte = pos / 2;
			int bin_bit = 4*(pos % 2), x1;

			for(x1 = 0 ;x1 < len; x1++)
			{
				int bin_int = (bin_block[bin_byte] >> bin_bit) & 0xf;
				if(bin_bit == 4) bin_byte++;
				bin_bit = (bin_bit == 4)?0:4;
				out_bases[x1] = contig_fasta_int2base(bin_int);
			}

			return 0;
		}
	} 
	return 1;
}

void destroy_contig_fasta(fasta_contigs_t * tab){
	HashTableDestroy( tab -> size_table );
	HashTableDestroy( tab -> contig_table );
}
int read_contig_fasta(fasta_contigs_t * tab, char * fname){
	autozip_fp fp;
	int rv = autozip_open(fname, &fp);
	if(rv>=0){
		tab -> contig_table = HashTableCreate(3943);
		tab -> size_table = HashTableCreate(3943);

		HashTableSetDeallocationFunctions(tab -> contig_table, free, free);
		HashTableSetDeallocationFunctions(tab -> size_table, NULL, NULL);

		HashTableSetKeyComparisonFunction(tab -> contig_table, fc_strcmp_chro);
		HashTableSetKeyComparisonFunction(tab -> size_table, fc_strcmp_chro);

		HashTableSetHashFunction(tab -> contig_table, fc_chro_hash);
		HashTableSetHashFunction(tab -> size_table, fc_chro_hash);

		char chro_name[MAX_CHROMOSOME_NAME_LEN];
		unsigned int inner_cursor = 0, current_bin_space = 0;
		int status = 0;
		char * bin_block = NULL;
		chro_name[0]=0;

		while(1){
			int nch = autozip_getch(&fp);
			if(nch<0)break;
			if(status == 0){
				assert(nch == '>');
				status = 1;
			}else if(status == 1){
				if(inner_cursor == 0){
					bin_block = calloc(sizeof(char),10000);
					current_bin_space = 10000;
				}
				if(nch == '|' || nch == ' ') status = 2;
				else if(nch == '\n'){
					status = 3;
					inner_cursor = 0;
				}else{
					chro_name[inner_cursor++] = nch;
					chro_name[inner_cursor] = 0;
				}
			}else if(status == 2){
				if(nch == '\n'){
					status = 3;
					inner_cursor = 0;
				}
			}else if(status == 3){
				if(nch == '>' || nch <= 0){
					char * mem_chro = malloc(strlen(chro_name)+1);
					strcpy(mem_chro, chro_name);
					HashTablePut(tab -> size_table , mem_chro, NULL + inner_cursor);
					HashTablePut(tab -> contig_table , mem_chro, bin_block);
		//			SUBREADprintf("Read '%s' : %u bases\n", chro_name, inner_cursor);
					inner_cursor = 0;
					status = 1;
					if(nch <= 0) break;
				}else if(nch != '\n'){
					int bin_bytes = inner_cursor / 2;
					int bin_bits = 4*(inner_cursor % 2);
					int base_int = contig_fasta_base2int(nch);
					if(bin_bytes >= current_bin_space){
						unsigned int new_bin_space = current_bin_space / 4 * 5;
						if(current_bin_space > 0xffff0000 /5 * 4){
							assert(0);
							return 1;
						}
						bin_block = realloc(bin_block, new_bin_space);
						memset(bin_block + current_bin_space, 0, new_bin_space - current_bin_space);
						current_bin_space = new_bin_space;
					}
					bin_block[bin_bytes] |= (base_int << bin_bits);
					inner_cursor++;
				}
			}
		}

		autozip_close(&fp);
		return 0;
	}
	return 1;
}

int RSubread_parse_CIGAR_Extra_string(int FLAG, char * MainChro, unsigned int MainPos, const char * CIGAR_Str, const char * Extra_Tags, int max_M, char ** Chros, unsigned int * Staring_Chro_Points, unsigned short * Section_Start_Read_Pos, unsigned short * Section_Length, int * is_junction_read){
	int ret = RSubread_parse_CIGAR_string(MainChro, MainPos, CIGAR_Str, max_M, Chros, Staring_Chro_Points, Section_Start_Read_Pos, Section_Length, is_junction_read);

	char read_main_strand = (((FLAG & 0x40)==0x40) == ((FLAG & 0x10) == 0x10 ))?'-':'+';
	int tag_cursor=0;
	//SUBREADprintf("EXTRA=%s\n", Extra_Tags);
	int status = PARSE_STATUS_TAGNAME;
	char tag_name[2], typechar=0;
	int tag_inner_cursor=0;

	char current_fusion_char[MAX_CHROMOSOME_NAME_LEN];
	unsigned int current_fusion_pos = 0;
	char current_fusion_strand = 0;
	char current_fusion_cigar[max_M * 15];
	current_fusion_cigar [0] =0;
	current_fusion_char [0]=0;

	while(1){
		int nch = Extra_Tags[tag_cursor];
		if(status == PARSE_STATUS_TAGNAME){
			tag_name[tag_inner_cursor++] = nch;
			if(tag_inner_cursor == 2){
				status = PARSE_STATUS_TAGTYPE;
				tag_cursor += 1;
				assert(Extra_Tags[tag_cursor] == ':');
			}
		}else if(status == PARSE_STATUS_TAGTYPE){
			typechar = nch;
			tag_cursor +=1;
			assert(Extra_Tags[tag_cursor] == ':');
			tag_inner_cursor = 0;
			status = PARSE_STATUS_TAGVALUE;
		}else if(status == PARSE_STATUS_TAGVALUE){
			if(nch == '\t' || nch == 0 || nch == '\n'){
				if(current_fusion_cigar[0] && current_fusion_char[0] && current_fusion_pos && current_fusion_strand){
					//SUBREADprintf("ENTER CALC:%s\n", current_fusion_char );
					unsigned int left_pos = current_fusion_pos;
					if(current_fusion_strand!=read_main_strand)
						left_pos = find_left_end_cigar(current_fusion_pos, current_fusion_cigar);
					ret += RSubread_parse_CIGAR_string(current_fusion_char, left_pos, current_fusion_cigar, max_M - ret, Chros + ret, Staring_Chro_Points+ ret, Section_Start_Read_Pos+ ret, Section_Length + ret, is_junction_read);

					current_fusion_pos = 0;
					current_fusion_strand = 0;
					current_fusion_cigar [0] =0;
					current_fusion_char [0]=0;
					//SUBREADprintf("EXIT CALC:%s\n", current_fusion_char );
				}

				tag_inner_cursor = 0;
				status = PARSE_STATUS_TAGNAME;
			}else{
				if(tag_name[0]=='C' && tag_name[1]=='C' && typechar == 'Z'){
					current_fusion_char[tag_inner_cursor++]=nch;
					current_fusion_char[tag_inner_cursor]=0;
				}else if(tag_name[0]=='C' && tag_name[1]=='G' && typechar == 'Z'){
					current_fusion_cigar[tag_inner_cursor++]=nch;
					current_fusion_cigar[tag_inner_cursor]=0;
				}else if(tag_name[0]=='C' && tag_name[1]=='P' && typechar == 'i'){
					current_fusion_pos = current_fusion_pos * 10 + (nch - '0');
				}else if(tag_name[0]=='C' && tag_name[1]=='T' && typechar == 'Z'){
					//SUBREADprintf("pos=%d %c -> %c\n", tag_cursor, current_fusion_strand, nch);
					current_fusion_strand = nch;
					//SUBREADprintf("spo=%d %c -> %c\n", tag_cursor, current_fusion_strand, nch);
				}
			}
		}

		if(nch == 0 || nch == '\n'){
			assert(status == PARSE_STATUS_TAGNAME);
			break;
		}

		tag_cursor++;
		//SUBREADprintf("CUR=%d [%s], c=%d\n", tag_cursor, Extra_Tags, Extra_Tags[tag_cursor]);
	}
	return ret;
}

int RSubread_parse_CIGAR_string(char * chro , unsigned int first_pos, const char * CIGAR_Str, int max_M, char ** Section_Chromosomes, unsigned int * Section_Start_Chro_Pos,unsigned short * Section_Start_Read_Pos, unsigned short * Section_Chro_Length, int * is_junction_read)
{
	unsigned int tmp_int=0;
	int cigar_cursor=0;
	unsigned short current_section_chro_len=0, current_section_start_read_pos = 0, read_cursor = 0;
	unsigned int chromosome_cursor=first_pos;
	int ret=0, is_first_S = 1;

	for(cigar_cursor=0; ; cigar_cursor++)
	{
		char ch = CIGAR_Str[cigar_cursor];

		if(ch >='0' && ch <= '9')
		{
			tmp_int=tmp_int*10+(ch - '0');
		}
		else
		{
			if(ch == 'S'){
				if(is_first_S) current_section_start_read_pos = tmp_int;
				read_cursor += tmp_int;
			}
			else if(ch == 'M' || ch == 'X' || ch == '=') {
				read_cursor += tmp_int;
				current_section_chro_len += tmp_int;
				chromosome_cursor += tmp_int;
			} else if(ch == 'N' || ch == 'D' || ch=='I' || ch == 0) {
				if('N' == ch)(*is_junction_read)=1;
				if(ret < max_M)
				{
					if(current_section_chro_len>0)
					{
						Section_Chromosomes[ret] = chro;
						Section_Start_Chro_Pos[ret] = chromosome_cursor - current_section_chro_len;
						Section_Start_Read_Pos[ret] = current_section_start_read_pos;
						Section_Chro_Length[ret] = current_section_chro_len;
						ret ++;
					}
				}else assert(0);
				current_section_chro_len = 0;
				if(ch == 'I') read_cursor += tmp_int;
				else if(ch == 'N' || ch == 'D') chromosome_cursor += tmp_int;
				current_section_start_read_pos = read_cursor;

				if(ch == 0) break;
			}
			//printf("C=%c, TV=%d, CC=%d, RC=%d\n", ch, tmp_int, chromosome_cursor, current_section_chro_len);
			tmp_int = 0;
			is_first_S = 0;
		}
		if(cigar_cursor>100) return -1;
	}

	return ret;
}

void display_sections(char * CIGAR_Str)
{
	//int is_junc=0;
	int Section_Start_Chro_Pos[FC_CIGAR_PARSER_ITEMS];
	unsigned short Section_Start_Read_Pos[FC_CIGAR_PARSER_ITEMS];
	unsigned short Section_Chro_Length[FC_CIGAR_PARSER_ITEMS];

	int retv = 0;//RSubread_parse_CIGAR_string(CIGAR_Str, Section_Start_Chro_Pos, Section_Start_Read_Pos, Section_Chro_Length, &is_junc);

	int x1;
	SUBREADprintf("Cigar=%s ; Sections=%d\n", CIGAR_Str, retv);
	for(x1=0; x1<retv; x1++)
	{
		SUBREADprintf("   Section #%d: chro_offset=%d, read_offset=%u  length=%u\n",x1, Section_Start_Chro_Pos[x1], Section_Start_Read_Pos[x1], Section_Chro_Length[x1]);
	}
	SUBREADprintf("\n");
	
}


#define GECV_STATE_BEFORE 10
#define GECV_STATE_NAME 20
#define GECV_STATE_GAP 30
#define GECV_STATE_VALUE 40
#define GECV_STATE_QVALUE 50
#define GECV_STATE_QV_END 60
#define GECV_STATE_ERROR 9999

int GTF_extra_column_istoken_chr(char c)
{
	return (isalpha(c)||isdigit(c)||c=='_');
}

int GTF_extra_column_value(const char * Extra_Col, const char * Target_Name, char * Target_Value, int TargVal_Size)
{
	int state = GECV_STATE_BEFORE;
	int col_cursor = 0, is_correct_name=0;
	char name_buffer[200];
	int name_cursor = 0, value_cursor=-1;

	while(1)
	{
		if(name_cursor>190) return -1;
		char nch = Extra_Col[col_cursor];
		if(nch == '\n' || nch == '\r') nch = 0;
		if(state == GECV_STATE_BEFORE)
		{
			if(GTF_extra_column_istoken_chr(nch))
			{
				name_buffer[0] = nch;
				name_cursor = 1;
				state = GECV_STATE_NAME;
			}
			else if(nch != ' ' && nch != 0)
			{
				state = GECV_STATE_ERROR;
			}
		}
		else if(state == GECV_STATE_NAME)
		{
			if(nch == ' ' || nch == '=')
			{
				state = GECV_STATE_GAP;
				name_buffer[name_cursor] = 0;
				is_correct_name = (strcmp(name_buffer , Target_Name) == 0);
				//printf("CORR=%d : '%s'\n", is_correct_name, name_buffer);
			}
			else if(nch == '"')
			{
				name_buffer[name_cursor] = 0;
				is_correct_name = (strcmp(name_buffer , Target_Name) == 0);
				state = GECV_STATE_QVALUE;
				if(is_correct_name)
					value_cursor = 0;
			}
			else if(GTF_extra_column_istoken_chr(nch))
				name_buffer[name_cursor++] = nch;
			else
			{
				state = GECV_STATE_ERROR;
				//printf("ERR2  : '%c'\n", nch);
			}
			
		}
		else if(state == GECV_STATE_GAP)
		{
			if(nch == '"')
			{
				state = GECV_STATE_QVALUE;
				if(is_correct_name)
					value_cursor = 0;
			}
			else if(nch != '=' && isgraph(nch))
			{
				state = GECV_STATE_VALUE;
				if(is_correct_name)
				{
					Target_Value[0]=nch;
					value_cursor = 1;
				}
			}
			else if(nch != ' ' && nch != '=')
				state = GECV_STATE_ERROR;
		}
		else if(state == GECV_STATE_VALUE)
		{
			if(nch == ';' || nch == 0)
			{
				state = GECV_STATE_BEFORE;
				if(is_correct_name)
				{
					Target_Value[value_cursor] = 0;
				}
				is_correct_name = 0;
			}
			else{
				if(value_cursor < TargVal_Size-1 && is_correct_name)
					Target_Value[value_cursor++] = nch;
			}
		}
		else if(state == GECV_STATE_QVALUE)
		{
			if(nch == '"')
			{
				state = GECV_STATE_QV_END;
				if(is_correct_name)
					Target_Value[value_cursor] = 0;
				is_correct_name = 0;
			}
			else
			{
				if(value_cursor < TargVal_Size-1 && is_correct_name)
				{
					if(nch !=' ' || value_cursor>0)
						Target_Value[value_cursor++] = nch;
				}
			}
		}
		else if(state == GECV_STATE_QV_END)
		{
			if(nch == ';' || nch == 0)
				state = GECV_STATE_BEFORE;
			else if(nch != ' ')
				state = GECV_STATE_ERROR;
				
		}

		if (GECV_STATE_ERROR == state){
			Target_Value[0]=0;
			return -1;
		}
		if (nch == 0)
		{
			if(state == GECV_STATE_BEFORE && value_cursor>0)
			{
				int x1;
				for(x1 = value_cursor-1; x1>=0; x1--)
				{
					if(Target_Value[x1] == ' '){
						value_cursor --;
						Target_Value[x1]=0;
					}
					else break;
				}

				if(value_cursor>0)
					return value_cursor;
			}
			Target_Value[0]=0;
			return -1;
		}
		col_cursor++;
	}
}


void hpl_test2_func()
{
	char * extra_column = " gene_id \"PC4-013  \"; 013=ABCD  ; PC4 =  CCXX  ";
	char * col_name = "gene_id";
	char col_val[100];

	int col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);

	col_name = "013";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);

	col_name = "PC4";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);


	col_name = "XXX";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);

	extra_column = "gene_id =   \"PC4-013  ;=\"  ;013 = AXXD ; PC4=x";
	col_name = "013";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);

	col_name = "gene_id";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);


	col_name = "PC4";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);




	extra_column = " gene_id\"  PC4-013  ;=  \"; XXX='123' ;013 :ABCD  ; PC4 =  CCXX=  ;";
	col_name = "013";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);


	col_name = "XXX";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);


	col_name = "PC4";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);

	col_name = "gene_id";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);



	extra_column = "gene_id \"653635\"; transcript_id \"TR:653635\";";
	col_name = "gene_id";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);





}

void hpl_test1_func()
{
	display_sections("");
	display_sections("*");
	display_sections("5S10M2D10M800N12M3I12M450N12M12D99M6S");
	display_sections("110M2I10M800N32M3I12M6S");
	display_sections("200S110M2I10M800N32M3I12M200N40M");
	display_sections("3M1663N61M1045N36M3D20M66N10M2D10M77N3M1663N61M1045N36M3D20M66N103M1663N61M1045N36M3D20M66N9M");
}

void test_bam_reader(){
	SamBam_FILE * bamfp = SamBam_fopen("/usr/local/work/liao/arena/Rsubread_Paper_OCT2016/RNAseq-SimHG38/Runs-100/STAR/STAR-Simulation-15M-DXC.bamAligned.out.bam", SAMBAM_FILE_BAM);
	assert(bamfp);

	while(1){
		char linebuf[2000];
		char * ret = SamBam_fgets(bamfp, linebuf, 1999, 1);
		if(!ret) break;
		SUBREADprintf("%s", linebuf); // linebuf has "\n" in the end.
	}
	SamBam_fclose(bamfp);
}



#ifdef HELPER_TEST_SAMBAM_READER
main() {
	//hpl_test1_func();
	test_bam_reader();
//	hpl_test1_func();
}
#endif

char *str_replace(char *orig, char *rep, char *with) {
    char *result; // the return string
    char *ins;    // the next insert point
    char *tmp;    // varies
    int len_rep;  // length of rep
    int len_with; // length of with
    int len_front; // distance between rep and end of last rep
    int count;    // number of replacements

    if (!orig)
        return NULL;
    if (!rep)
        rep = "";
    len_rep = strlen(rep);
    if (!with)
        with = "";
    len_with = strlen(with);

    ins = orig;
    for (count = 0; NULL != (tmp = strstr(ins, rep)); ++count) {
        ins = tmp + len_rep;
    }
    tmp = result = malloc(strlen(orig) + (len_with - len_rep) * count + 1);

    if (!result)
        return NULL;

    while (count--) {
        ins = strstr(orig, rep);
        len_front = ins - orig;
        tmp = strncpy(tmp, orig, len_front) + len_front;
        tmp = strcpy(tmp, with) + len_with;
        orig += len_front + len_rep; // move to next "end of rep"
    }
    strcpy(tmp, orig);
    return result;
}



// rule: the string is ABC123XXXXXX...
// This is the priroity:
// First, compare the letters part.
// Second, compare the pure numeric part.
// Third, compare the remainder.
int strcmp_number(char * s1, char * s2)
{
	int x1 = 0;
	int ret = 0;

	while(1)
	{
		char nch1 = s1[x1];
		char nch2 = s2[x1];

		if((!nch1) || !nch2){return nch2?1:(nch1?(-1):0);}
		if(isdigit(nch1) && isdigit(nch2))break;

		ret = nch1 - nch2;
		if(ret) return ret;
		x1++;
	}

	int v1 = 0, v2 = 0;
	while(1)
	{
		char nch1 = s1[x1];
		char nch2 = s2[x1];
		if((!nch1) || !nch2){
			if(nch1 || nch2)
				return nch2?(-1):1;
			break;
		}
		int is_chr1_digit = isdigit(nch1);
		int is_chr2_digit = isdigit(nch2);

		if(is_chr1_digit || is_chr2_digit)
		{
			if(is_chr1_digit && is_chr2_digit)
			{
				v1 = v1*10+(nch1-'0');
				v2 = v2*10+(nch2-'0');
			}
			else
			{
				ret = nch1 - nch2;
				return ret;
			}
		}
		else break;
		x1++;
	}

	if(v1==v2)
		return strcmp(s1+x1, s2+x1);	
	else
	{
		ret = v1 - v2;
		return ret;
	}
}



int mac_str(char * str_buff)
{
#if defined(FREEBSD) || defined(__MINGW32__)
	return 1;
#else
#ifdef MACOS
    int         mib[6], x1, ret = 1;
	size_t		len;
    char            *buf;
    unsigned char       *ptr;
    struct if_msghdr    *ifm;
    struct sockaddr_dl  *sdl;


	for(x1 = 0 ; x1 < 40; x1++)
    {
		mib[0] = CTL_NET;
		mib[1] = AF_ROUTE;
		mib[2] = 0;
		mib[3] = AF_LINK;
		mib[4] = NET_RT_IFLIST;
		mib[5] = x1;

		if (sysctl(mib, 6, NULL, &len, NULL, 0) >=0) {
			buf = malloc(len);
			if (sysctl(mib, 6, buf, &len, NULL, 0) >=0) {

				ifm = (struct if_msghdr *)buf;
				sdl = (struct sockaddr_dl *)(ifm + 1);
				ptr = (unsigned char *)LLADDR(sdl);

				if(sdl -> sdl_nlen < 1) continue;

				char * ifname = malloc(sdl -> sdl_nlen + 1);

				memcpy(ifname, sdl -> sdl_data, sdl -> sdl_nlen);
				ifname[sdl -> sdl_nlen] = 0;
				if(ifname[0]!='e'){
					free(ifname);
					continue;
				}
				free(ifname);

				sprintf(str_buff,"%02X%02X%02X%02X%02X%02X",  *ptr, *(ptr+1), *(ptr+2),
					*(ptr+3), *(ptr+4), *(ptr+5));
				ret = 0;
				break;
			}
			free(buf);
		}
	}
    return ret;
#else
#if defined(IFHWADDRLEN)
    struct ifreq ifr;
    struct ifconf ifc;
    char buf[1024];
    int success = 0;

    int sock = socket(AF_INET, SOCK_DGRAM, IPPROTO_IP);
    if (sock == -1) { /* handle error*/ };

    ifc.ifc_len = sizeof(buf);
    ifc.ifc_buf = buf;
    if (ioctl(sock, SIOCGIFCONF, &ifc) == -1) { /* handle error */ }

    struct ifreq* it = ifc.ifc_req;
    const struct ifreq* const end = it + (ifc.ifc_len / sizeof(struct ifreq));

    for (; it != end; ++it) {
        strcpy(ifr.ifr_name, it->ifr_name);
        if (ioctl(sock, SIOCGIFFLAGS, &ifr) == 0) {
            if (! (ifr.ifr_flags & IFF_LOOPBACK)) { // don't count loopback
                if (ioctl(sock, SIOCGIFHWADDR, &ifr) == 0) {
                      success = 1;
                      break;
                }
            }
        }
    }

    close(sock);

    unsigned char mac_address[6];

    if (success){
	memcpy(mac_address, ifr.ifr_hwaddr.sa_data, 6);
	    int x1;
	    for(x1 = 0; x1 < 6; x1++){
		 sprintf(str_buff+2*x1, "%02X",mac_address[x1]);
	    }
		return 0;
	}
#endif
	return 1;
#endif
#endif
}

int rand_str(char * str_buff){
	int ret = 1;
	FILE * fp = fopen("/dev/urandom","r");
	if(fp){
		int x1;
		for(x1=0; x1<6; x1++){
			sprintf(str_buff + 2*x1 , "%02X", fgetc(fp));
		}
		fclose(fp);
		ret = 0;
	}
	return ret;
}

int mathrand_str(char * str_buff){
	myrand_srand((int)(miltime()*100));
	int x1;
	for(x1 = 0; x1 < 6; x1++){
		sprintf(str_buff+2*x1, "%02X", myrand_rand() & 0xff );
	}
	return 0;
}

int mac_or_rand_str(char * str_buff){
	return mac_str(str_buff) && rand_str(str_buff) && mathrand_str(str_buff);
}

#define PI_LONG 3.1415926535897932384626434L

long double fast_fractorial(unsigned int x, long double * buff, int buflen){
	if(x<2) return 0;
	
	if(buff != NULL && x < buflen && buff[x]!=0){
		return buff[x];
	}
	long double ret;

	if(x<50){
		int x1;
		ret =0.L;
		for(x1 = 2 ; x1 <= x; x1++) ret += logl((long double)(x1));
	}else{
		ret = logl(x)*1.0L*x - 1.0L*x + 0.5L * logl(2.0L*PI_LONG* x) + 1.L/(12.L*x) - 1.L/(360.L* x*x*x) +  1.L/(1260.L* x*x*x*x*x) -  1.L/(1680.L*x*x*x*x*x*x*x);
	}
	if(buff != NULL && x < buflen) buff[x]=ret;
	return ret;
}


#define BUFF_4D 36

long double fast_freq( unsigned int tab[2][2] , long double * buff, int buflen){
	int x0 = -1;
	if(buff && buflen > BUFF_4D * BUFF_4D * BUFF_4D * BUFF_4D && tab[0][0] < BUFF_4D && tab[0][1] < BUFF_4D && tab[1][0] < BUFF_4D && tab[1][1] < BUFF_4D ){
		x0 = buflen + tab[0][0]* BUFF_4D * BUFF_4D * BUFF_4D + tab[0][1] * BUFF_4D * BUFF_4D + tab[1][0] * BUFF_4D + tab[1][1];
		if(buff[x0]!=0) return buff[x0];
		
	}
	long double ret = fast_fractorial(tab[0][0]+tab[0][1],buff,buflen)+fast_fractorial(tab[1][0]+tab[1][1],buff,buflen) + 
		fast_fractorial(tab[0][0]+tab[1][0],buff,buflen)+fast_fractorial(tab[0][1]+tab[1][1],buff,buflen) -
		fast_fractorial(tab[0][0],buff,buflen) - fast_fractorial(tab[0][1],buff,buflen) - 
		fast_fractorial(tab[1][0],buff,buflen) - fast_fractorial(tab[1][1],buff,buflen) - 
		fast_fractorial(tab[0][0] + tab[0][1] + tab[1][0] + tab[1][1],buff,buflen);
	if(x0>=0) buff[x0] = ret;
	return ret;
}



double fast_fisher_test_one_side(unsigned int a, unsigned int b, unsigned int c, unsigned int d, long double * buff, int buflen){
	unsigned int tab[2][2];
	long double P0, P_delta, ret;

	tab[0][0]=a;
	tab[0][1]=b;
	tab[1][0]=c;
	tab[1][1]=d;

	int x_a, y_a, x_min=-1, y_min=-1;
	unsigned int min_a = 0xffffffff;
	for(x_a = 0; x_a < 2; x_a++)
		for(y_a = 0; y_a < 2; y_a++){
			if(tab[x_a][y_a]<min_a){
				min_a = tab[x_a][y_a];
				x_min = x_a;
				y_min = y_a;
			}
		}
	P_delta = fast_freq(tab, buff, buflen);
	P0 = ret = exp(P_delta);
	//printf("P0=%LG\n", P0);
	if(min_a>0){
		unsigned int Qa = min_a;
		unsigned int Qb = tab[x_min][!y_min];
		unsigned int Qc = tab[!x_min][y_min];
		unsigned int Qd = tab[!x_min][!y_min];
		for(; ; ){
			min_a --;
			Qb++;Qc++;
			P_delta -= logl(Qb*Qc);
			P_delta += logl(Qa*Qd);
			Qa--;Qd--;
			ret += expl(P_delta);
		//	printf("%LG %LG %LG\n", ret, 1 - (ret - P0), expl(P_delta));
			if(min_a < 1) break;
		}
	}

	double ret1 = ret, ret2 = 1 - (ret - P0);

	if(min(ret1, ret2) < 0){
		return 0.0;
	}
	return  min(ret1, ret2) ;

}

int load_features_annotation(char * file_name, int file_type, char * gene_id_column, char * transcript_id_column, char * used_feature_type,
 void * context, int do_add_feature(char * gene_name, char * transcript_name, char * chro_name, unsigned int start, unsigned int end, int is_negative_strand, void * context)  ){
	char * file_line = malloc(MAX_LINE_LENGTH+1);
	int lineno = 0, is_GFF_txid_warned = 0, is_GFF_geneid_warned = 0, loaded_features = 0;
	autozip_fp afp;
	int aret = autozip_open(file_name, &afp);

	if(aret < 0){
		SUBREADprintf("Error: unable to open the annotation file : %s\n", file_name);
		return -1;
	}

	while(1){
		int is_tx_id_found = 0, is_gene_id_found = 0, is_negative_strand = -1;
		char * token_temp = NULL, * feature_name, *transcript_id = NULL, * chro_name = NULL;
		char feature_name_tmp[FEATURE_NAME_LENGTH], txid_tmp[FEATURE_NAME_LENGTH];
		feature_name = feature_name_tmp;
		
		unsigned int start = 0, end = 0;
		aret = autozip_gets(&afp, file_line, MAX_LINE_LENGTH);
		if(aret < 1) break;

		lineno++;
		if(is_comment_line(file_line, file_type, lineno-1))continue;

		if(file_type == FILE_TYPE_RSUBREAD)
		{
			feature_name = strtok_r(file_line,"\t",&token_temp);
			int feature_name_len = strlen(feature_name);
			if(feature_name_len > FEATURE_NAME_LENGTH) feature_name[FEATURE_NAME_LENGTH -1 ] = 0;
			
			chro_name = strtok_r(NULL,"\t", &token_temp);
			int chro_name_len = strlen(chro_name);
			if(chro_name_len > MAX_CHROMOSOME_NAME_LEN) chro_name[MAX_CHROMOSOME_NAME_LEN -1 ] = 0;

			char * start_ptr = strtok_r(NULL,"\t", &token_temp);
			char * end_ptr = strtok_r(NULL,"\t", &token_temp);

			if(start_ptr == NULL || end_ptr == NULL){
				SUBREADprintf("\nWarning: the format on the %d-th line is wrong.\n", lineno);
			}
			long long int tv1 = atoll(start_ptr);
			long long int tv2 = atoll(end_ptr);

			if( isdigit(start_ptr[0]) && isdigit(end_ptr[0]) ){
				if(strlen(start_ptr) > 10 || strlen(end_ptr) > 10 || tv1 > 0x7fffffff || tv2> 0x7fffffff){
					SUBREADprintf("\nError: Line %d contains a coordinate greater than 2^31.\n", lineno);
					return -2;
				}
			}else{
				SUBREADprintf("\nError: Line %d contains a format error. The expected annotation format is SAF.\n", lineno);
				return -2;
			}
			
			start = atoi(start_ptr);// start
			end = atoi(end_ptr);//end
			
			char * strand_str = strtok_r(NULL,"\t", &token_temp);
			if(strand_str == NULL)
				is_negative_strand = 0;
			else
				is_negative_strand = ('-' ==strand_str[0]);
				
			is_gene_id_found = 1;
			
		} else if(file_type == FILE_TYPE_GTF) {
			chro_name = strtok_r(file_line,"\t",&token_temp);
			strtok_r(NULL,"\t", &token_temp);// source
			char * feature_type = strtok_r(NULL,"\t", &token_temp);// feature_type
			
			if(strcmp(feature_type, used_feature_type)==0){
				char * start_ptr = strtok_r(NULL,"\t", &token_temp);
				char * end_ptr = strtok_r(NULL,"\t", &token_temp);
				

				if(start_ptr == NULL || end_ptr == NULL){
					SUBREADprintf("\nWarning: the format on the %d-th line is wrong.\n", lineno);
				}
				long long int tv1 = atoll(start_ptr);
				long long int tv2 = atoll(end_ptr);
				
				
				if( isdigit(start_ptr[0]) && isdigit(end_ptr[0]) ){
					if(strlen(start_ptr) > 10 || strlen(end_ptr) > 10 || tv1 > 0x7fffffff || tv2> 0x7fffffff){
						SUBREADprintf("\nError: Line %d contains a coordinate greater than 2^31.\n", lineno);
						return -2;
					}
				}else{
					SUBREADprintf("\nError: Line %d contains a format error. The expected annotation format is GTF/GFF.\n", lineno);
					return -2;
				}
				start = atoi(start_ptr);// start
				end = atoi(end_ptr);//end

				if(start < 1 || end<1 ||  start > 0x7fffffff || end > 0x7fffffff || start > end)
					SUBREADprintf("\nWarning: the feature on the %d-th line has zero coordinate or zero lengths\n\n", lineno);

					
				strtok_r(NULL,"\t", &token_temp);// score
				is_negative_strand = ('-' == (strtok_r(NULL,"\t", &token_temp)[0]));//strand
				strtok_r(NULL,"\t",&token_temp);	// "frame"
				char * extra_attrs = strtok_r(NULL,"\t",&token_temp);   // name_1 "val1"; name_2 "val2"; ...
				if(extra_attrs && (strlen(extra_attrs)>2)){
					int attr_val_len = GTF_extra_column_value(extra_attrs , gene_id_column , feature_name_tmp, FEATURE_NAME_LENGTH);
					if(attr_val_len>0) is_gene_id_found=1;

					if(transcript_id_column){
							transcript_id = txid_tmp;
							attr_val_len = GTF_extra_column_value(extra_attrs , transcript_id_column , txid_tmp, FEATURE_NAME_LENGTH);
							if(attr_val_len>0) is_tx_id_found=1;
							else transcript_id = NULL;
					}
				}

				if(!is_gene_id_found){
					if(!is_GFF_geneid_warned){
						int ext_att_len = strlen(extra_attrs);
						if(extra_attrs[ext_att_len-1] == '\n') extra_attrs[ext_att_len-1] =0;
						SUBREADprintf("\nERROR: failed to find the gene identifier attribute in the 9th column of the provided GTF file.\nThe specified gene identifier attribute is '%s'.\nAn example of attributes included in your GTF annotation is '%s'.\n\n",  gene_id_column, extra_attrs);
					}
					is_GFF_geneid_warned++;
				}
					
				if(transcript_id_column && !is_tx_id_found){
					if(!is_GFF_txid_warned){
						int ext_att_len = strlen(extra_attrs);
						if(extra_attrs[ext_att_len-1] == '\n') extra_attrs[ext_att_len-1] =0;
						SUBREADprintf("\nERROR: failed to find the transcript identifier attribute in the 9th column of the provided GTF file.\nThe specified transcript identifier attribute is '%s'.\nAn example of attributes included in your GTF annotation is '%s'.\n\n", transcript_id_column, extra_attrs);
					}
					is_GFF_txid_warned++;
				}
			}
		}
		
		if(is_gene_id_found){
			do_add_feature(feature_name, transcript_id, chro_name, start, end, is_negative_strand, context);
			loaded_features++;
		}
		
	}
	autozip_close(&afp);
	free(file_line);

	if(is_GFF_txid_warned || is_GFF_geneid_warned)return -2;
	if(loaded_features<1){
		SUBREADprintf("\nERROR: No feature was loaded from the annotation file. Please check if the annotation format was correctly specified, and also if the feature type was correctly specified if the annotation is in the GTF format.\n\n");
		return -2;
	}
	return loaded_features;
}

HashTable * load_alias_table(char * fname) {
	FILE * fp = f_subr_open(fname, "r");
	if(!fp)
	{
		print_in_box(80,0,0,"WARNING unable to open alias file '%s'", fname);
		return NULL;
	}

	char * fl = malloc(2000);

	HashTable * ret = HashTableCreate(1013);
	HashTableSetDeallocationFunctions(ret, free, free);
	HashTableSetKeyComparisonFunction(ret, fc_strcmp_chro);
	HashTableSetHashFunction(ret, fc_chro_hash);
	
	while (1)
	{
		char *ret_fl = fgets(fl, 1999, fp);
		if(!ret_fl) break;
		if(fl[0]=='#') continue;
		char * sam_chr = NULL;
		char * anno_chr = strtok_r(fl, ",", &sam_chr);
		if((!sam_chr)||(!anno_chr)) continue;

		sam_chr[strlen(sam_chr)-1]=0;
		if(sam_chr[strlen(sam_chr)-1]=='\r') sam_chr[strlen(sam_chr)-1]=0;
		char * anno_chr_buf = malloc(strlen(anno_chr)+1);
		strcpy(anno_chr_buf, anno_chr);
		char * sam_chr_buf = malloc(strlen(sam_chr)+1);
		strcpy(sam_chr_buf, sam_chr);
//		SUBREADprintf("ALIAS_LOAD : '%s' => '%s'\n", sam_chr_buf, anno_chr_buf);
		HashTablePut(ret, sam_chr_buf, anno_chr_buf);
	}

	fclose(fp);

	free(fl);
	return ret;
}

int rebuild_command_line(char ** lineptr, int argc, char ** argv){
	int linecap = 1000,c;
	*lineptr = malloc(linecap);
	**lineptr = 0;

	for(c = 0; c<argc;c++)
	{
		int cline = strlen(argv[c]);
		if(strlen(*lineptr) + 100 + cline > linecap)
		{
			linecap += cline+500;
			*lineptr = realloc(*lineptr, linecap);
		}
		sprintf((*lineptr) + strlen(*lineptr), "\"%s\" ", argv[c]);
	}

	return strlen(*lineptr);
}

// The following code was downloaded from
//   https://openwall.info/wiki/people/solar/software/public-domain-source-code/md5
//   on the 23rd of Nov, 2018. I (Yang LIAO) made minor changes to the code.
//
// It says "Public Domain" in the licence but the original author should be acknwoledged (only for expressing my gratitude):
//   Alexander Peslyak, better known as Solar Designer <solar at openwall.com>

/* Any 32-bit or wider unsigned integer data type will do */

/*
 * The basic HelpFuncMD5 functions.
 *
 * F and G are optimized compared to their RFC 1321 definitions for
 * architectures that lack an AND-NOT instruction, just like in Colin Plumb's
 * implementation.
 */
#define F(x, y, z)			((z) ^ ((x) & ((y) ^ (z))))
#define G(x, y, z)			((y) ^ ((z) & ((x) ^ (y))))
#define H(x, y, z)			(((x) ^ (y)) ^ (z))
#define H2(x, y, z)			((x) ^ ((y) ^ (z)))
#define I(x, y, z)			((y) ^ ((x) | ~(z)))
 
/*
 * The HelpFuncMD5 transformation for all four rounds.
 */
#define STEP(f, a, b, c, d, x, t, s) \
	(a) += f((b), (c), (d)) + (x) + (t); \
	(a) = (((a) << (s)) | (((a) & 0xffffffff) >> (32 - (s)))); \
	(a) += (b);
 
/*
 * SET reads 4 input bytes in little-endian byte order and stores them in a
 * properly aligned word in host byte order.
 *
 * The check for little-endian architectures that tolerate unaligned memory
 * accesses is just an optimization.  Nothing will break if it fails to detect
 * a suitable architecture.
 *
 * Unfortunately, this optimization may be a C strict aliasing rules violation
 * if the caller's data buffer has effective type that cannot be aliased by
 * HelpFuncMD5_u32plus.  In practice, this problem may occur if these HelpFuncMD5 routines are
 * inlined into a calling function, or with future and dangerously advanced
 * link-time optimizations.  For the time being, keeping these HelpFuncMD5 routines in
 * their own translation unit avoids the problem.
 */
#define SET(n) (*(HelpFuncMD5_u32plus *)&ptr[(n) * 4])
#define GET(n) SET(n)
/*
 * This processes one or more 64-byte data blocks, but does NOT update the bit
 * counters.  There are no alignment requirements.
 */
static const void *body(HelpFuncMD5_CTX *ctx, const void *data, unsigned long size)
{
	const unsigned char *ptr;
	HelpFuncMD5_u32plus a, b, c, d;
	HelpFuncMD5_u32plus saved_a, saved_b, saved_c, saved_d;
 
	ptr = (const unsigned char *)data;
 
	a = ctx->a;
	b = ctx->b;
	c = ctx->c;
	d = ctx->d;
 
	do {
		saved_a = a;
		saved_b = b;
		saved_c = c;
		saved_d = d;
 
/* Round 1 */
		STEP(F, a, b, c, d, SET(0), 0xd76aa478, 7)
		STEP(F, d, a, b, c, SET(1), 0xe8c7b756, 12)
		STEP(F, c, d, a, b, SET(2), 0x242070db, 17)
		STEP(F, b, c, d, a, SET(3), 0xc1bdceee, 22)
		STEP(F, a, b, c, d, SET(4), 0xf57c0faf, 7)
		STEP(F, d, a, b, c, SET(5), 0x4787c62a, 12)
		STEP(F, c, d, a, b, SET(6), 0xa8304613, 17)
		STEP(F, b, c, d, a, SET(7), 0xfd469501, 22)
		STEP(F, a, b, c, d, SET(8), 0x698098d8, 7)
		STEP(F, d, a, b, c, SET(9), 0x8b44f7af, 12)
		STEP(F, c, d, a, b, SET(10), 0xffff5bb1, 17)
		STEP(F, b, c, d, a, SET(11), 0x895cd7be, 22)
		STEP(F, a, b, c, d, SET(12), 0x6b901122, 7)
		STEP(F, d, a, b, c, SET(13), 0xfd987193, 12)
		STEP(F, c, d, a, b, SET(14), 0xa679438e, 17)
		STEP(F, b, c, d, a, SET(15), 0x49b40821, 22)
 
/* Round 2 */
		STEP(G, a, b, c, d, GET(1), 0xf61e2562, 5)
		STEP(G, d, a, b, c, GET(6), 0xc040b340, 9)
		STEP(G, c, d, a, b, GET(11), 0x265e5a51, 14)
		STEP(G, b, c, d, a, GET(0), 0xe9b6c7aa, 20)
		STEP(G, a, b, c, d, GET(5), 0xd62f105d, 5)
		STEP(G, d, a, b, c, GET(10), 0x02441453, 9)
		STEP(G, c, d, a, b, GET(15), 0xd8a1e681, 14)
		STEP(G, b, c, d, a, GET(4), 0xe7d3fbc8, 20)
		STEP(G, a, b, c, d, GET(9), 0x21e1cde6, 5)
		STEP(G, d, a, b, c, GET(14), 0xc33707d6, 9)
		STEP(G, c, d, a, b, GET(3), 0xf4d50d87, 14)
		STEP(G, b, c, d, a, GET(8), 0x455a14ed, 20)
		STEP(G, a, b, c, d, GET(13), 0xa9e3e905, 5)
		STEP(G, d, a, b, c, GET(2), 0xfcefa3f8, 9)
		STEP(G, c, d, a, b, GET(7), 0x676f02d9, 14)
		STEP(G, b, c, d, a, GET(12), 0x8d2a4c8a, 20)
 
/* Round 3 */
		STEP(H, a, b, c, d, GET(5), 0xfffa3942, 4)
		STEP(H2, d, a, b, c, GET(8), 0x8771f681, 11)
		STEP(H, c, d, a, b, GET(11), 0x6d9d6122, 16)
		STEP(H2, b, c, d, a, GET(14), 0xfde5380c, 23)
		STEP(H, a, b, c, d, GET(1), 0xa4beea44, 4)
		STEP(H2, d, a, b, c, GET(4), 0x4bdecfa9, 11)
		STEP(H, c, d, a, b, GET(7), 0xf6bb4b60, 16)
		STEP(H2, b, c, d, a, GET(10), 0xbebfbc70, 23)
		STEP(H, a, b, c, d, GET(13), 0x289b7ec6, 4)
		STEP(H2, d, a, b, c, GET(0), 0xeaa127fa, 11)
		STEP(H, c, d, a, b, GET(3), 0xd4ef3085, 16)
		STEP(H2, b, c, d, a, GET(6), 0x04881d05, 23)
		STEP(H, a, b, c, d, GET(9), 0xd9d4d039, 4)
		STEP(H2, d, a, b, c, GET(12), 0xe6db99e5, 11)
		STEP(H, c, d, a, b, GET(15), 0x1fa27cf8, 16)
		STEP(H2, b, c, d, a, GET(2), 0xc4ac5665, 23)
 
/* Round 4 */
		STEP(I, a, b, c, d, GET(0), 0xf4292244, 6)
		STEP(I, d, a, b, c, GET(7), 0x432aff97, 10)
		STEP(I, c, d, a, b, GET(14), 0xab9423a7, 15)
		STEP(I, b, c, d, a, GET(5), 0xfc93a039, 21)
		STEP(I, a, b, c, d, GET(12), 0x655b59c3, 6)
		STEP(I, d, a, b, c, GET(3), 0x8f0ccc92, 10)
		STEP(I, c, d, a, b, GET(10), 0xffeff47d, 15)
		STEP(I, b, c, d, a, GET(1), 0x85845dd1, 21)
		STEP(I, a, b, c, d, GET(8), 0x6fa87e4f, 6)
		STEP(I, d, a, b, c, GET(15), 0xfe2ce6e0, 10)
		STEP(I, c, d, a, b, GET(6), 0xa3014314, 15)
		STEP(I, b, c, d, a, GET(13), 0x4e0811a1, 21)
		STEP(I, a, b, c, d, GET(4), 0xf7537e82, 6)
		STEP(I, d, a, b, c, GET(11), 0xbd3af235, 10)
		STEP(I, c, d, a, b, GET(2), 0x2ad7d2bb, 15)
		STEP(I, b, c, d, a, GET(9), 0xeb86d391, 21)
 
		a += saved_a;
		b += saved_b;
		c += saved_c;
		d += saved_d;
 
		ptr += 64;
	} while (size -= 64);
 
	ctx->a = a;
	ctx->b = b;
	ctx->c = c;
	ctx->d = d;
 
	return ptr;
}
 
void HelpFuncMD5_Init(HelpFuncMD5_CTX *ctx)
{
	ctx->a = 0x67452301;
	ctx->b = 0xefcdab89;
	ctx->c = 0x98badcfe;
	ctx->d = 0x10325476;
 
	ctx->lo = 0;
	ctx->hi = 0;
}
 
void HelpFuncMD5_Update(HelpFuncMD5_CTX *ctx, const void *data, unsigned long size)
{
	HelpFuncMD5_u32plus saved_lo;
	unsigned long used, available;
 
	saved_lo = ctx->lo;
	if ((ctx->lo = (saved_lo + size) & 0x1fffffff) < saved_lo)
		ctx->hi++;
	ctx->hi += size >> 29;
 
	used = saved_lo & 0x3f;
 
	if (used) {
		available = 64 - used;
 
		if (size < available) {
			memcpy(&ctx->buffer[used], data, size);
			return;
		}
 
		memcpy(&ctx->buffer[used], data, available);
		data = (const unsigned char *)data + available;
		size -= available;
		body(ctx, ctx->buffer, 64);
	}
 
	if (size >= 64) {
		data = body(ctx, data, size & ~(unsigned long)0x3f);
		size &= 0x3f;
	}
 
	memcpy(ctx->buffer, data, size);
}
 
#define OUT(dst, src) \
	(dst)[0] = (unsigned char)(src); \
	(dst)[1] = (unsigned char)((src) >> 8); \
	(dst)[2] = (unsigned char)((src) >> 16); \
	(dst)[3] = (unsigned char)((src) >> 24);
 
void HelpFuncMD5_Final(unsigned char *result, HelpFuncMD5_CTX *ctx)
{
	unsigned long used, available;
 
	used = ctx->lo & 0x3f;
 
	ctx->buffer[used++] = 0x80;
 
	available = 64 - used;
 
	if (available < 8) {
		memset(&ctx->buffer[used], 0, available);
		body(ctx, ctx->buffer, 64);
		used = 0;
		available = 64;
	}
 
	memset(&ctx->buffer[used], 0, available - 8);
 
	ctx->lo <<= 3;
	OUT(&ctx->buffer[56], ctx->lo)
	OUT(&ctx->buffer[60], ctx->hi)
 
	body(ctx, ctx->buffer, 64);
 
	OUT(&result[0], ctx->a)
	OUT(&result[4], ctx->b)
	OUT(&result[8], ctx->c)
	OUT(&result[12], ctx->d)
 
	memset(ctx, 0, sizeof(*ctx));
}


void Helper_md5sum(char * plain_txt, int plain_len, unsigned char * bin_md5_buff){
	HelpFuncMD5_CTX ctx;
	HelpFuncMD5_Init(&ctx);
	HelpFuncMD5_Update(&ctx, plain_txt, plain_len);
	HelpFuncMD5_Final(bin_md5_buff, &ctx);
}

unsigned long long plain_txt_to_long_rand(char * plain_txt, int plain_len){
	unsigned char md5v[16];
	Helper_md5sum(plain_txt, plain_len, md5v);
	unsigned long long ret = 0;
	memcpy(&ret, md5v, sizeof(ret));
	return ret;
}

void md5txt(char *s){
	unsigned char md5v[16];
	unsigned long long randv;
	Helper_md5sum(s, strlen(s), md5v);
	int x;
	randv = plain_txt_to_long_rand(s, strlen(s));

	for(x=0;x<16;x++){
		SUBREADprintf("%02X", md5v[x]);
	}
	#ifdef __MINGW32__
	SUBREADprintf("\t'%s'\t%016I64X\t%I64u\t%.9f\n", s, randv, randv, randv*1./0xffffffffffffffffllu);
	#else
	SUBREADprintf("\t'%s'\t%016llX\t%llu\t%.9f\n", s, randv, randv, randv*1./0xffffffffffffffffllu);
	#endif
}

//#define TESTHelpermain main

int TESTHelpermain(){
	int xx;
	for(xx=0; xx<100; xx++){
		char xt[10];
		sprintf(xt, "%08d", xx);
		md5txt(xt);
	}
	
	return 0;
}


// The following code was retrieved from 
//   https://github.com/google/omaha/blob/master/third_party/lzma/files/C/Helper_Sha256.c
// Thank the authors -- Igor Pavlov and Wei Dai.

/* Crypto/Helper_Sha256.c -- SHA-256 Hash
2010-06-11 : Igor Pavlov : Public domain
This code is based on public domain code from Wei Dai's Crypto++ library. */
#define SHA256_DIGEST_SIZE 32

typedef struct
{
  unsigned int state[8];
  unsigned long long count;
  unsigned char buffer[64];
} CHelper_Sha256;

/* define it for speed optimization */
/* #define _SHA256_UNROLL */
/* #define _SHA256_UNROLL2 */

void Helper_Sha256_Init(CHelper_Sha256 *p)
{
  p->state[0] = 0x6a09e667;
  p->state[1] = 0xbb67ae85;
  p->state[2] = 0x3c6ef372;
  p->state[3] = 0xa54ff53a;
  p->state[4] = 0x510e527f;
  p->state[5] = 0x9b05688c;
  p->state[6] = 0x1f83d9ab;
  p->state[7] = 0x5be0cd19;
  p->count = 0;
}

#define rotlFixed(x, n) (((x) << (n)) | ((x) >> (32 - (n))))
#define rotrFixed(x, n) (((x) >> (n)) | ((x) << (32 - (n))))

#define S0(x) (rotrFixed(x, 2) ^ rotrFixed(x,13) ^ rotrFixed(x, 22))
#define S1(x) (rotrFixed(x, 6) ^ rotrFixed(x,11) ^ rotrFixed(x, 25))
#define s0(x) (rotrFixed(x, 7) ^ rotrFixed(x,18) ^ (x >> 3))
#define s1(x) (rotrFixed(x,17) ^ rotrFixed(x,19) ^ (x >> 10))

#define blk0(i) (W[i] = data[i])
#define blk2(i) (W[i&15] += s1(W[(i-2)&15]) + W[(i-7)&15] + s0(W[(i-15)&15]))

#define Ch(x,y,z) (z^(x&(y^z)))
#define Maj(x,y,z) ((x&y)|(z&(x|y)))

#define a(i) T[(0-(i))&7]
#define b(i) T[(1-(i))&7]
#define c(i) T[(2-(i))&7]
#define d(i) T[(3-(i))&7]
#define e(i) T[(4-(i))&7]
#define f(i) T[(5-(i))&7]
#define g(i) T[(6-(i))&7]
#define h(i) T[(7-(i))&7]


#ifdef _SHA256_UNROLL2

#define R(a,b,c,d,e,f,g,h, i) h += S1(e) + Ch(e,f,g) + K[i+j] + (j?blk2(i):blk0(i));\
  d += h; h += S0(a) + Maj(a, b, c)

#define RX_8(i) \
  R(a,b,c,d,e,f,g,h, i); \
  R(h,a,b,c,d,e,f,g, i+1); \
  R(g,h,a,b,c,d,e,f, i+2); \
  R(f,g,h,a,b,c,d,e, i+3); \
  R(e,f,g,h,a,b,c,d, i+4); \
  R(d,e,f,g,h,a,b,c, i+5); \
  R(c,d,e,f,g,h,a,b, i+6); \
  R(b,c,d,e,f,g,h,a, i+7)

#else

#define R(i) h(i) += S1(e(i)) + Ch(e(i),f(i),g(i)) + K[i+j] + (j?blk2(i):blk0(i));\
  d(i) += h(i); h(i) += S0(a(i)) + Maj(a(i), b(i), c(i))

#ifdef _SHA256_UNROLL

#define RX_8(i) R(i+0); R(i+1); R(i+2); R(i+3); R(i+4); R(i+5); R(i+6); R(i+7);

#endif

#endif

static const unsigned int K[64] = {
  0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
  0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
  0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
  0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
  0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
  0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
  0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
  0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
  0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
  0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
  0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
  0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
  0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
  0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
  0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
  0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

static void Helper_Sha256_Transform(unsigned int *state, const unsigned int *data)
{
  unsigned int W[16];
  unsigned j;
  #ifdef _SHA256_UNROLL2
  unsigned int a,b,c,d,e,f,g,h;
  a = state[0];
  b = state[1];
  c = state[2];
  d = state[3];
  e = state[4];
  f = state[5];
  g = state[6];
  h = state[7];
  #else
  unsigned int T[8];
  for (j = 0; j < 8; j++)
    T[j] = state[j];
  #endif

  for (j = 0; j < 64; j += 16)
  {
    #if defined(_SHA256_UNROLL) || defined(_SHA256_UNROLL2)
    RX_8(0); RX_8(8);
    #else
    unsigned i;
    for (i = 0; i < 16; i++) { R(i); }
    #endif
  }

  #ifdef _SHA256_UNROLL2
  state[0] += a;
  state[1] += b;
  state[2] += c;
  state[3] += d;
  state[4] += e;
  state[5] += f;
  state[6] += g;
  state[7] += h;
  #else
  for (j = 0; j < 8; j++)
    state[j] += T[j];
  #endif
  
  /* Wipe variables */
  /* memset(W, 0, sizeof(W)); */
  /* memset(T, 0, sizeof(T)); */
}

#undef S0
#undef S1
#undef s0
#undef s1

static void Helper_Sha256_WritecharBlock(CHelper_Sha256 *p)
{
  unsigned int data32[16];
  unsigned i;
  for (i = 0; i < 16; i++)
    data32[i] =
      ((unsigned int)(p->buffer[i * 4    ]) << 24) +
      ((unsigned int)(p->buffer[i * 4 + 1]) << 16) +
      ((unsigned int)(p->buffer[i * 4 + 2]) <<  8) +
      ((unsigned int)(p->buffer[i * 4 + 3]));
  Helper_Sha256_Transform(p->state, data32);
}

void Helper_Sha256_Update(CHelper_Sha256 *p, const char *data, size_t size)
{
  unsigned int curBufferPos = (unsigned int)p->count & 0x3F;
  while (size > 0)
  {
    p->buffer[curBufferPos++] = (char)*data;
    p->count++;

	data++;
    size--;
    if (curBufferPos == 64)
    {
      curBufferPos = 0;
      Helper_Sha256_WritecharBlock(p);
    }
  }
}

void Helper_Sha256_Final( unsigned char *digest , CHelper_Sha256 *p)
{
  unsigned long long lenInBits = (p->count << 3);
  unsigned int curBufferPos = (unsigned int)p->count & 0x3F;
  unsigned i;
  p->buffer[curBufferPos++] = 0x80;
  while (curBufferPos != (64 - 8))
  {
    curBufferPos &= 0x3F;
    if (curBufferPos == 0)
      Helper_Sha256_WritecharBlock(p);
    p->buffer[curBufferPos++] = 0;
  }
  for (i = 0; i < 8; i++)
  {
    p->buffer[curBufferPos++] = (char)(lenInBits >> 56);
    lenInBits <<= 8;
  }
  Helper_Sha256_WritecharBlock(p);

  for (i = 0; i < 8; i++)
  {
    *digest++ = (unsigned char)(p->state[i] >> 24);
    *digest++ = (unsigned char)(p->state[i] >> 16);
    *digest++ = (unsigned char)(p->state[i] >> 8);
    *digest++ = (unsigned char)(p->state[i]);
  }
}

void Helper_sha256sum(char * plain_txt, int plain_len, unsigned char * bin_md5_buff){
	CHelper_Sha256 ctx;
	Helper_Sha256_Init(&ctx);
	Helper_Sha256_Update(&ctx, plain_txt, plain_len);
	Helper_Sha256_Final(bin_md5_buff, &ctx);
}

void sha256txt(char *s){
	unsigned char sha256v[32];
	Helper_sha256sum(s, strlen(s), sha256v);
	int x;

	for(x=0;x<32;x++){
		SUBREADprintf("%02X", sha256v[x]);
	}
	SUBREADprintf("\t'%s'\n", s);
}

//#define  TEST256Helpermain main
int TEST256Helpermain(){
	int xx;
	for(xx=0; xx<100; xx++){
		char xt[10];
		sprintf(xt, "%08d", xx);
		sha256txt(xt);
	}
	return 0;
}



// The code of Helper_erfinv() was retreived from
//    https://github.com/antelopeusersgroup/antelope_contrib/blob/master/lib/location/libgenloc/erfinv.c
// I (Yang LIAO) adopted and modified it because it is under the New BSD license.
//
// Copyright (c) 2014 Indiana University
// All rights reserved. 
//
// Written by Prof. Gary L. Pavlis, Dept. of Geol. Sci.,
//   Indiana University, Bloomington, IN
//
// This software is licensed under the New BSD license: 
//
// Redistribution and use in source and binary forms,
// with or without modification, are permitted provided
// that the following conditions are met:
//
// Redistributions of source code must retain the above
// copyright notice, this list of conditions and the
// following disclaimer.
//
// Redistributions in binary form must reproduce the
// above copyright notice, this list of conditions and
// the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Neither the name of Indiana University nor
// the names of its contributors may be used to endorse
// or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
// CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
// THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
// IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.


#include <math.h>
#include <float.h>
#define MAXDOUBLE DBL_MAX

#define CENTRAL_RANGE 0.7
double Helper_erfinv( double y) {
        double x=0,z,num,dem; /*working variables */
        /* coefficients in rational expansion */
        double a[4]={ 0.886226899, -1.645349621,  0.914624893, -0.140543331};
        double b[4]={-2.118377725,  1.442710462, -0.329097515,  0.012229801};
        double c[4]={-1.970840454, -1.624906493,  3.429567803,  1.641345311};
        double d[2]={ 3.543889200,  1.637067800};
        if(fabs(y) > 1.0) return (atof("NaN"));  /* This needs IEEE constant*/
        if(fabs(y) == 1.0) return((copysign(1.0,y))*MAXDOUBLE); 
        if( fabs(y) <= CENTRAL_RANGE ) 
        {
                z = y*y;
                num = (((a[3]*z + a[2])*z + a[1])*z + a[0]);
                dem = ((((b[3]*z + b[2])*z + b[1])*z +b[0])*z + 1.0);
                x = y*num/dem;
        }
        else if( (fabs(y) > CENTRAL_RANGE) && (fabs(y) < 1.0) )
        {
                z = sqrt(-log((1.0-fabs(y))/2.0));
                num = ((c[3]*z + c[2])*z + c[1])*z + c[0];
                dem = (d[1]*z + d[0])*z + 1.0;
                x = (copysign(1.0,y))*num/dem;
        }
        /* Two steps of Newton-Raphson correction */
        x = x - (erf(x) - y)/( (2.0/sqrt(M_PI))*exp(-x*x));
        x = x - (erf(x) - y)/( (2.0/sqrt(M_PI))*exp(-x*x));

        return(x);
}


double inverse_sample_normal(double p){
	return Helper_erfinv(2*p-1)*1.4142135623730950488;
}

//#define TestNormalMain main
void TestNormalMain(){
	int i;
	for(i=0; i<=10; i++){
		double ii = i*1./10;
		double p_0 = inverse_sample_normal(ii);
		SUBREADprintf("p of %.1f = %.40f\n\n" , ii, p_0);
	}
}


#ifndef MAKE_STANDALONE

message_queue_t mt_message_queue;


int safeRprintf(char * fmt,...){
	va_list args;
	va_start(args , fmt);
	FILE * ofp = fopen("/tmp/del4-Rlog.txt", "a");
	//fprintf(ofp, "NEW LOG : %p \"%s\"\n", fmt, fmt);
	fprintf(ofp, fmt, args);
	fclose(ofp);
	va_end(args);
	return 0;
}

#endif

void msgqu_destroy(){
	#ifndef MAKE_STANDALONE
	ArrayListDestroy(mt_message_queue.message_queue);
	subread_destroy_lock(&mt_message_queue.queue_lock);
	#endif
}

void msgqu_notifyFinish(){
	#ifndef MAKE_STANDALONE
	subread_lock_occupy(&mt_message_queue.queue_lock);
	mt_message_queue.is_finished = 1;
	subread_lock_release(&mt_message_queue.queue_lock);
	#endif
}

void msgqu_init(int thread_mode){
	#ifndef MAKE_STANDALONE
	mt_message_queue.is_thread_mode = thread_mode;
	if(!thread_mode) return;
	mt_message_queue.is_finished = 0;
	mt_message_queue.message_queue = ArrayListCreate(100);
	ArrayListSetDeallocationFunction(mt_message_queue.message_queue,free);
	subread_init_lock(&mt_message_queue.queue_lock);
//	subread_init_lock(&mt_message_queue.queue_notifier);
	#endif
}

void msgqu_main_loop(){
	#ifndef MAKE_STANDALONE
	while(1){
		subread_lock_occupy(&mt_message_queue.queue_lock);
		while(0< mt_message_queue.message_queue->numOfElements){
			char * amsg = ArrayListShift(mt_message_queue.message_queue);
			Rprintf("%s", amsg);
			free(amsg);
		}
		if(mt_message_queue.is_finished) break;
		subread_lock_release(&mt_message_queue.queue_lock);
		usleep(40000);
	}
	#endif
}

#define MSGQU_LINE_SIZE 3000
void msgqu_printf(const char * fmt, ...){
	va_list args;
	va_start(args , fmt);
	char * obuf = malloc(MSGQU_LINE_SIZE+1);
	vsnprintf( obuf, MSGQU_LINE_SIZE, fmt, args );
	va_end(args);
	obuf[MSGQU_LINE_SIZE]=0;
	#ifdef MAKE_STANDALONE
		fputs(obuf, stderr);
		free(obuf);
	#else
		if(!mt_message_queue.is_thread_mode){
			Rprintf("%s", obuf);
			free(obuf);
			return;
		}

		subread_lock_occupy(&mt_message_queue.queue_lock);
		ArrayListPush(mt_message_queue.message_queue, obuf);
		subread_lock_release(&mt_message_queue.queue_lock);
		// it will be freed in the Main thread after printing.
	#endif
}

#ifdef MAKE_TEST_MYRAND_TEST

int main(){
	myrand_srand(time(NULL));
	int x;
	for(x = 0; x<100000; x++){
		int mv = myrand_rand();
		SUBREADprintf("%d\n", mv);
	}
	return 0;
}

#endif




// retrieved from https://github.com/kokke/tiny-TNbignum-c/
// The original code was in public domain

/* Functions for shifting number in-place. */
static void _lshift_one_bit(struct bn* a);
static void _rshift_one_bit(struct bn* a);
static void _lshift_word(struct bn* a, int nwords);
static void _rshift_word(struct bn* a, int nwords);



/* Public / Exported functions. */
void TNbignum_init(struct bn* n)
{
  require(n, "n is null");

  int i;
  for (i = 0; i < BN_ARRAY_SIZE; ++i)
  {
    n->array[i] = 0;
  }
}


void TNbignum_from_int(struct bn* n, DTYPE_TMP i)
{
  require(n, "n is null");

  TNbignum_init(n);

  /* Endianness issue if machine is not little-endian? */
#ifdef WORD_SIZE
 #if (WORD_SIZE == 1)
  n->array[0] = (i & 0x000000ff);
  n->array[1] = (i & 0x0000ff00) >> 8;
  n->array[2] = (i & 0x00ff0000) >> 16;
  n->array[3] = (i & 0xff000000) >> 24;
 #elif (WORD_SIZE == 2)
  n->array[0] = (i & 0x0000ffff);
  n->array[1] = (i & 0xffff0000) >> 16;
 #elif (WORD_SIZE == 4)
  n->array[0] = i;
  DTYPE_TMP num_32 = 32;
  DTYPE_TMP tmp = i >> num_32; /* bit-shift with U64 operands to force 64-bit results */
  n->array[1] = tmp;
 #endif
#endif
}


int TNbignum_to_int(struct bn* n)
{
  require(n, "n is null");

  int ret = 0;

  /* Endianness issue if machine is not little-endian? */
#if (WORD_SIZE == 1)
  ret += n->array[0];
  ret += n->array[1] << 8;
  ret += n->array[2] << 16;
  ret += n->array[3] << 24;  
#elif (WORD_SIZE == 2)
  ret += n->array[0];
  ret += n->array[1] << 16;
#elif (WORD_SIZE == 4)
  ret += n->array[0];
#endif

  return ret;
}


void TNbignum_from_string(struct bn* n, char* str, int nbytes)
{
  require(n, "n is null");
  require(str, "str is null");
  require(nbytes > 0, "nbytes must be positive");
  require((nbytes & 1) == 0, "string format must be in hex -> equal number of bytes");

  TNbignum_init(n);

  DTYPE tmp;                        /* DTYPE is defined in bn.h - uint{8,16,32,64}_t */
  int i = nbytes - (2 * WORD_SIZE); /* index into string */
  int j = 0;                        /* index into array */

  /* reading last hex-byte "MSB" from string first -> big endian */
  /* MSB ~= most significant byte / block ? :) */
  while (i >= 0)
  {
    tmp = 0;
    sscanf(&str[i], SSCANF_FORMAT_STR, &tmp);
    //printf("SCAN_IN %d : v=%u\n", i, tmp);
    n->array[j] = tmp;
    i -= (2 * WORD_SIZE); /* step WORD_SIZE hex-byte(s) back in the string. */
    j += 1;               /* step one element forward in the array. */
  }
}


void TNbignum_to_string(struct bn* n, char* str, int nbytes)
{
  require(n, "n is null");
  require(str, "str is null");
  require(nbytes > 0, "nbytes must be positive");
  require((nbytes & 1) == 0, "string format must be in hex -> equal number of bytes");

  int j = BN_ARRAY_SIZE - 1; /* index into array - reading "MSB" first -> big-endian */
  int i = 0;                 /* index into string representation. */

  /* reading last array-element "MSB" first -> big endian */
  while ((j >= 0) && (nbytes > (i + 1)))
  {
    sprintf(&str[i], SPRINTF_FORMAT_STR, n->array[j]);
    //printf("WRITE:%d %s\n" , i, str+i);
    i += (2 * WORD_SIZE); /* step WORD_SIZE hex-byte(s) forward in the string. */
    j -= 1;               /* step one element back in the array. */
  }

  /* Count leading zeros: */
  j = 0;
  while (str[j] == '0')
  {
    j += 1;
  }
 
  /* Move string j places ahead, effectively skipping leading zeros */ 
  for (i = 0; i < (nbytes - j); ++i)
  {
    str[i] = str[i + j];
  }

  /* Zero-terminate string */
  str[i] = 0;
}


void TNbignum_dec(struct bn* n)
{
  require(n, "n is null");

  DTYPE tmp; /* copy of n */
  DTYPE res;

  int i;
  for (i = 0; i < BN_ARRAY_SIZE; ++i)
  {
    tmp = n->array[i];
    res = tmp - 1;
    n->array[i] = res;

    if (!(res > tmp))
    {
      break;
    }
  }
}


void TNbignum_inc(struct bn* n)
{
  require(n, "n is null");

  DTYPE res;
  DTYPE_TMP tmp; /* copy of n */

  int i;
  for (i = 0; i < BN_ARRAY_SIZE; ++i)
  {
    tmp = n->array[i];
    res = tmp + 1;
    n->array[i] = res;

    if (res > tmp)
    {
      break;
    }
  }
}


void TNbignum_add(struct bn* a, struct bn* b, struct bn* c)
{
  require(a, "a is null");
  require(b, "b is null");
  require(c, "c is null");

  DTYPE_TMP tmp;
  int carry = 0;
  int i;
  for (i = 0; i < BN_ARRAY_SIZE; ++i)
  {
    tmp = (DTYPE_TMP)a->array[i] + b->array[i] + carry;
    carry = (tmp > MAX_VAL);
    c->array[i] = (tmp & MAX_VAL);
  }
}


void TNbignum_sub(struct bn* a, struct bn* b, struct bn* c)
{
  require(a, "a is null");
  require(b, "b is null");
  require(c, "c is null");

  DTYPE_TMP res;
  DTYPE_TMP tmp1;
  DTYPE_TMP tmp2;
  int borrow = 0;
  int i;
  for (i = 0; i < BN_ARRAY_SIZE; ++i)
  {
    tmp1 = (DTYPE_TMP)a->array[i] + (MAX_VAL + 1); /* + number_base */
    tmp2 = (DTYPE_TMP)b->array[i] + borrow;;
    res = (tmp1 - tmp2);
    c->array[i] = (DTYPE)(res & MAX_VAL); /* "modulo number_base" == "% (number_base - 1)" if number_base is 2^N */
    borrow = (res <= MAX_VAL);
  }
}


void TNbignum_mul(struct bn* a, struct bn* b, struct bn* c)
{
  require(a, "a is null");
  require(b, "b is null");
  require(c, "c is null");

  struct bn row;
  struct bn tmp;
  int i, j;

  TNbignum_init(c);

  for (i = 0; i < BN_ARRAY_SIZE; ++i)
  {
    TNbignum_init(&row);

    for (j = 0; j < BN_ARRAY_SIZE; ++j)
    {
      if (i + j < BN_ARRAY_SIZE)
      {
        TNbignum_init(&tmp);
        DTYPE_TMP intermediate = ((DTYPE_TMP)a->array[i] * (DTYPE_TMP)b->array[j]);
        TNbignum_from_int(&tmp, intermediate);
        _lshift_word(&tmp, i + j);
        TNbignum_add(&tmp, &row, &row);
      }
    }
    TNbignum_add(c, &row, c);
  }
}


void TNbignum_div(struct bn* a, struct bn* b, struct bn* c)
{
  require(a, "a is null");
  require(b, "b is null");
  require(c, "c is null");

  struct bn current;
  struct bn denom;
  struct bn tmp;

  TNbignum_from_int(&current, 1);               // int current = 1;
  TNbignum_assign(&denom, b);                   // denom = b
  TNbignum_assign(&tmp, a);                     // tmp   = a

  const DTYPE_TMP half_max = 1 + (DTYPE_TMP)(MAX_VAL / 2);
  int overflow = 0;
  while (TNbignum_cmp(&denom, a) != LARGER)     // while (denom <= a) {
  {
    if (denom.array[BN_ARRAY_SIZE - 1] >= half_max)
    {
      overflow = 1;
      break;
    }
    _lshift_one_bit(&current);                //   current <<= 1;
    _lshift_one_bit(&denom);                  //   denom <<= 1;
  }
  if (!overflow)
  {
    _rshift_one_bit(&denom);                  // denom >>= 1;
    _rshift_one_bit(&current);                // current >>= 1;
  }
  TNbignum_init(c);                             // int answer = 0;

  while (!TNbignum_is_zero(&current))           // while (current != 0)
  {
    if (TNbignum_cmp(&tmp, &denom) != SMALLER)  //   if (dividend >= denom)
    {
      TNbignum_sub(&tmp, &denom, &tmp);         //     dividend -= denom;
      TNbignum_or(c, &current, c);              //     answer |= current;
    }
    _rshift_one_bit(&current);                //   current >>= 1;
    _rshift_one_bit(&denom);                  //   denom >>= 1;
  }                                           // return answer;
}


void TNbignum_lshift(struct bn* a, struct bn* b, int nbits)
{
  require(a, "a is null");
  require(b, "b is null");
  require(nbits >= 0, "no negative shifts");

  TNbignum_assign(b, a);
  /* Handle shift in multiples of word-size */
  const int nbits_pr_word = (WORD_SIZE * 8);
  int nwords = nbits / nbits_pr_word;
  if (nwords != 0)
  {
    _lshift_word(b, nwords);
    nbits -= (nwords * nbits_pr_word);
  }

  if (nbits != 0)
  {
    int i;
    for (i = (BN_ARRAY_SIZE - 1); i > 0; --i)
    {
      b->array[i] = (b->array[i] << nbits) | (b->array[i - 1] >> ((8 * WORD_SIZE) - nbits));
    }
    b->array[i] <<= nbits;
  }
}


void TNbignum_rshift(struct bn* a, struct bn* b, int nbits)
{
  require(a, "a is null");
  require(b, "b is null");
  require(nbits >= 0, "no negative shifts");
  
  TNbignum_assign(b, a);
  /* Handle shift in multiples of word-size */
  const int nbits_pr_word = (WORD_SIZE * 8);
  int nwords = nbits / nbits_pr_word;
  if (nwords != 0)
  {
    _rshift_word(b, nwords);
    nbits -= (nwords * nbits_pr_word);
  }

  if (nbits != 0)
  {
    int i;
    for (i = 0; i < (BN_ARRAY_SIZE - 1); ++i)
    {
      b->array[i] = (b->array[i] >> nbits) | (b->array[i + 1] << ((8 * WORD_SIZE) - nbits));
    }
    b->array[i] >>= nbits;
  }
  
}


void TNbignum_mod(struct bn* a, struct bn* b, struct bn* c)
{
  /*
 *     Take divmod and throw away div part
 *       */
  require(a, "a is null");
  require(b, "b is null");
  require(c, "c is null");

  struct bn tmp;

  TNbignum_divmod(a,b,&tmp,c);
}

void TNbignum_divmod(struct bn* a, struct bn* b, struct bn* c, struct bn* d)
{
  /*
 *     Puts a%b in d
 *         and a/b in c
 *
 *             mod(a,b) = a - ((a / b) * b)
 *
 *                 example:
 *                       mod(8, 3) = 8 - ((8 / 3) * 3) = 2
 *                         */
  require(a, "a is null");
  require(b, "b is null");
  require(c, "c is null");

  struct bn tmp;

  /* c = (a / b) */
  TNbignum_div(a, b, c);

  /* tmp = (c * b) */
  TNbignum_mul(c, b, &tmp);

  /* c = a - tmp */
  TNbignum_sub(a, &tmp, d);
}


void TNbignum_and(struct bn* a, struct bn* b, struct bn* c)
{
  require(a, "a is null");
  require(b, "b is null");
  require(c, "c is null");

  int i;
  for (i = 0; i < BN_ARRAY_SIZE; ++i)
  {
    c->array[i] = (a->array[i] & b->array[i]);
  }
}


void TNbignum_or(struct bn* a, struct bn* b, struct bn* c)
{
  require(a, "a is null");
  require(b, "b is null");
  require(c, "c is null");

  int i;
  for (i = 0; i < BN_ARRAY_SIZE; ++i)
  {
    c->array[i] = (a->array[i] | b->array[i]);
  }
}


void TNbignum_xor(struct bn* a, struct bn* b, struct bn* c)
{
  require(a, "a is null");
  require(b, "b is null");
  require(c, "c is null");

  int i;
  for (i = 0; i < BN_ARRAY_SIZE; ++i)
  {
    c->array[i] = (a->array[i] ^ b->array[i]);
  }
}


int TNbignum_cmp(struct bn* a, struct bn* b)
{
  require(a, "a is null");
  require(b, "b is null");

  int i = BN_ARRAY_SIZE;
  do
  {
    i -= 1; /* Decrement first, to start with last array element */
    if (a->array[i] > b->array[i])
    {
      return LARGER;
    }
    else if (a->array[i] < b->array[i])
    {
      return SMALLER;
    }
  }
  while (i != 0);

  return EQUAL;
}


int TNbignum_is_zero(struct bn* n)
{
  require(n, "n is null");

  int i;
  for (i = 0; i < BN_ARRAY_SIZE; ++i)
  {
    if (n->array[i])
    {
      return 0;
    }
  }

  return 1;
}


void TNbignum_pow(struct bn* a, struct bn* b, struct bn* c)
{
  require(a, "a is null");
  require(b, "b is null");
  require(c, "c is null");

  struct bn tmp;

  TNbignum_init(c);

  if (TNbignum_cmp(b, c) == EQUAL)
  {
    /* Return 1 when exponent is 0 -- n^0 = 1 */
    TNbignum_inc(c);
  }
  else
  {
    /* Copy a -> tmp */
    TNbignum_assign(&tmp, a);

    TNbignum_dec(b);
 
    /* Begin summing products: */
    while (!TNbignum_is_zero(b))
    {

      /* c = tmp * tmp */
      TNbignum_mul(&tmp, a, c);
      /* Decrement b by one */
      TNbignum_dec(b);

      TNbignum_assign(&tmp, c);
    }

    /* c = tmp */
    TNbignum_assign(c, &tmp);
  }
}

void TNbignum_isqrt(struct bn *a, struct bn* b)
{
  require(a, "a is null");
  require(b, "b is null");

  struct bn low, high, mid, tmp;

  TNbignum_init(&low);
  TNbignum_assign(&high, a);
  TNbignum_rshift(&high, &mid, 1);
  TNbignum_inc(&mid);

  while (TNbignum_cmp(&high, &low) > 0) 
  {
    TNbignum_mul(&mid, &mid, &tmp);
    if (TNbignum_cmp(&tmp, a) > 0) 
    {
      TNbignum_assign(&high, &mid);
      TNbignum_dec(&high);
    }
    else 
    {
      TNbignum_assign(&low, &mid);
    }
    TNbignum_sub(&high,&low,&mid);
    _rshift_one_bit(&mid);
    TNbignum_add(&low,&mid,&mid);
    TNbignum_inc(&mid);
  }
  TNbignum_assign(b,&low);
}


void TNbignum_assign(struct bn* dst, struct bn* src)
{
  require(dst, "dst is null");
  require(src, "src is null");

  int i;
  for (i = 0; i < BN_ARRAY_SIZE; ++i)
  {
    dst->array[i] = src->array[i];
  }
}


/* Private / Static functions. */
static void _rshift_word(struct bn* a, int nwords)
{
  /* Naive method: */
  require(a, "a is null");
  require(nwords >= 0, "no negative shifts");

  int i;
  if (nwords >= BN_ARRAY_SIZE)
  {
    for (i = 0; i < BN_ARRAY_SIZE; ++i)
    {
      a->array[i] = 0;
    }
    return;
  }

  for (i = 0; i < BN_ARRAY_SIZE - nwords; ++i)
  {
    a->array[i] = a->array[i + nwords];
  }
  for (; i < BN_ARRAY_SIZE; ++i)
  {
    a->array[i] = 0;
  }
}


static void _lshift_word(struct bn* a, int nwords)
{
  require(a, "a is null");
  require(nwords >= 0, "no negative shifts");

  int i;
  /* Shift whole words */
  for (i = (BN_ARRAY_SIZE - 1); i >= nwords; --i)
  {
    a->array[i] = a->array[i - nwords];
  }
  /* Zero pad shifted words. */
  for (; i >= 0; --i)
  {
    a->array[i] = 0;
  }  
}


static void _lshift_one_bit(struct bn* a)
{
  require(a, "a is null");

  int i;
  for (i = (BN_ARRAY_SIZE - 1); i > 0; --i)
  {
    a->array[i] = (a->array[i] << 1) | (a->array[i - 1] >> ((8 * WORD_SIZE) - 1));
  }
  a->array[0] <<= 1;
}


static void _rshift_one_bit(struct bn* a)
{
  require(a, "a is null");

  int i;
  for (i = 0; i < (BN_ARRAY_SIZE - 1); ++i)
  {
    a->array[i] = (a->array[i] >> 1) | (a->array[i + 1] << ((8 * WORD_SIZE) - 1));
  }
  a->array[BN_ARRAY_SIZE - 1] >>= 1;
}

/*
void TNbignum_pow_then_mod( struct bn * m, int e, struct bn * modulus, struct bn * res ){
  int xx;
  struct bn res, remove_times, remove_value, tmpv;
  TNbignum_from_string(&res, "01", 2);
  for(xx = 0; xx < 31; xx++){
    TNbignum_mul( m, m, tmpv ); 
    TNbignum_div( tmpv, modulus, remove_times);
    TNbignum_mul( remove_times, modulus, remove_value);
    TNbignum_sub( tmpv , remove_value, m);
    e = e >> 1;
  }

}
*/

#ifdef HELPER_TEST_BIGNUM
void main(){
  char restxt[1030];
  struct bn n,e,res;
  int exponent = 35;
  restxt[0]=0;

  char * modulus = "AC97941B27989F9DFAC7FB03A951BF8D39CE60248D8FB5A2C4CDAFAE2431667DD5EBB46B6FF9BBEF5DE2B587CF6B06B5D63BB6B71B35C43FA5141F156A1AC77231FD5D916053CA3E5FEAD3AFB10877D1A5440119C6420A6C205758D01B5B75F5B420F9E99815820B389F0BFBA00B60C0A18612C268C0ED9B263B02DF526EED851581E0BDE4D46227DAA2EAA8D13EDF3476C20B30C3188FA3D0FBC1636C5477A744D0A042BA27C77E52CEFAE0B4005BF2E22001BA0C40A5898F3C8FE664969E72F208FDF8D7DC3D039690911E32DF6B4948280EAF67952D231907CAA1A3D1433D5A8E5325E3B849721871458E8EA5860E4EFFCCC61D82BE5EED6EAF9C6A292A2B";

  TNbignum_from_string(&n, modulus,strlen(modulus) );
  TNbignum_pow_then_mod(&m, exponent, &n, &res); // (m ^ e) mod N. E doesn't need to be a BIG NUMBER.
  TNbignum_to_string(&res, restxt, 1030);
  printf("RES=%s\n", restxt);
  

}

#endif

#ifdef HELPER_TEST_CRC32
void main(){
  long c32=0;
  c32 = crc32(c32, "Hello.",6);
  printf("CRC = %08lX\n", c32);

  c32=0;
  c32 = crc32(c32, "Hel",3);
  c32 = crc32(c32, "lo",2);
  c32 = crc32(c32, ".",1);
  c32 = crc32(c32, "",0);
  printf("CRC = %08lX\n", c32);
}

#endif

#ifdef __MINGW32__
#include <windows.h>
#endif


int get_free_total_mem(size_t * total, size_t * free_mem){

#ifdef FREEBSD
    return -1;
#endif

#ifdef __MINGW32__
	MEMORYSTATUSEX statex;
	statex.dwLength = sizeof (statex);
	GlobalMemoryStatusEx (&statex);
	(*total) = statex.ullTotalPhys;
	(*free_mem) = statex.ullAvailPhys;
	return 0;
#else
#ifdef MACOS
    mach_msg_type_number_t count = HOST_VM_INFO_COUNT;
    vm_statistics_data_t vmstat;
    int page_size = getpagesize();
    if(KERN_SUCCESS != host_statistics(mach_host_self(), HOST_VM_INFO, (host_info_t)&vmstat, &count))
        return -1;
    //printf("PSIZE=%d\nACT=%u; INACT=%u; FREE=%u\n", page_size, vmstat.active_count, vmstat.inactive_count, vmstat.free_count);
	size_t btlen = sizeof(*total);
    if(sysctl( (int[]) { CTL_HW, HW_MEMSIZE }, 2, total, &btlen, NULL, 0)) return -1;
    *free_mem = (vmstat.free_count + vmstat.inactive_count) * 1llu * page_size;
    return 0;
#else
    struct sysinfo sinf;
    sysinfo(&sinf);
    size_t cached_mem = get_sys_mem_info("Cached:");
    if(cached_mem<0)cached_mem=0;
    *free_mem = cached_mem + sinf.bufferram+sinf.freeram;
    *total = sinf.totalram;
    return 0;
#endif
#endif
}

void worker_master_mutex_init(worker_master_mutex_t * wmt, int all_workers){
	memset(wmt,0,sizeof(worker_master_mutex_t));
	wmt -> conds_worker_wait = malloc(sizeof(pthread_cond_t) * all_workers);
	wmt -> mutexs_worker_wait = malloc(sizeof(pthread_mutex_t) * all_workers);
	wmt -> mutex_with_master = calloc(sizeof(int), all_workers);
	wmt -> worker_is_working = calloc(sizeof(int), all_workers);

	wmt -> workers = all_workers;
	int x1;
	for(x1=0;x1<all_workers;x1++){
		pthread_cond_init(&wmt -> conds_worker_wait [x1], NULL);
		pthread_mutex_init(&wmt -> mutexs_worker_wait [x1], NULL);
	}
}

void worker_thread_start(worker_master_mutex_t * wmt, int worker_id){
	pthread_mutex_lock(&wmt->mutexs_worker_wait[worker_id]);
	wmt -> mutex_with_master[worker_id] = 0;
}

void worker_master_mutex_destroy(worker_master_mutex_t * wmt){
	int x1;
	for(x1=0;x1<wmt -> workers;x1++){
		pthread_mutex_destroy(&wmt -> mutexs_worker_wait [x1]);
		pthread_cond_destroy(&wmt -> conds_worker_wait [x1]);
	}
	free( wmt -> worker_is_working);
	free( wmt -> mutex_with_master);
	free( wmt -> conds_worker_wait);
	free( wmt -> mutexs_worker_wait);
}

int worker_wait_for_job(worker_master_mutex_t * wmt, int worker_id){
	pthread_mutex_trylock(wmt->mutexs_worker_wait + worker_id);
	wmt->worker_is_working[worker_id] = 0;
	while(1){
		pthread_cond_wait(&wmt->conds_worker_wait[worker_id], &wmt->mutexs_worker_wait[worker_id]);
		if(wmt -> all_terminate)
			pthread_mutex_unlock(&wmt->mutexs_worker_wait[worker_id]);
		if(wmt -> mutex_with_master[worker_id] == 0)break;
		// otherwise, it is a spurious wakeup
	}
	return wmt -> all_terminate;
}

void master_wait_for_job_done(worker_master_mutex_t * wmt, int worker_id){
	if(!wmt -> mutex_with_master[worker_id] ){
		while(1){
			pthread_mutex_lock(&wmt->mutexs_worker_wait[worker_id]);
			if(!wmt->worker_is_working[worker_id])break;
			// In a rare situation, the worker hasn't been scheduled after the last notification.
			// Namely, the master will notify twice the worker when the worker thread hasn't started at all. The second notification will cause no condition signaled.
			// This test is to prevent such a situation.
			pthread_mutex_unlock(&wmt->mutexs_worker_wait[worker_id]);
			usleep(50);
		}
		if(wmt->worker_is_working[worker_id])
			SUBREADprintf("ERROR 3: HOW CAN I HAVE THIS LOCK : %d?? TH_%d\n", wmt->worker_is_working[worker_id], worker_id);
	}
	wmt -> mutex_with_master[worker_id] = 1;
}

// master collects results between master_wait_for_job_done and master_notify_worker
void master_notify_worker(worker_master_mutex_t * wmt, int worker_id){
	if(!wmt -> mutex_with_master[worker_id])SUBREADprintf("ERROR 2: HOW CAN I NOT HAVE THE LOCK : %d ; TERM=%d\n", worker_id, wmt -> all_terminate);

	wmt -> worker_is_working[worker_id] = 1;
	wmt -> mutex_with_master[worker_id] = 0;
	pthread_cond_signal(&wmt->conds_worker_wait[worker_id]);
	pthread_mutex_unlock(&wmt->mutexs_worker_wait[worker_id]);
}

// this function must be called when the master thread has all the worker locks in control.
// ie all the workers should be in the "wait_for_job" status.
void terminate_workers(worker_master_mutex_t * wmt){
	int x1;
	wmt -> all_terminate = 1;
	for(x1=0;x1<wmt -> workers;x1++)
		master_notify_worker(wmt,x1);
}

void *windows_memmem(const void *haystack_start, size_t haystack_len, const void *needle_start, size_t needle_len)
{

    const unsigned char *haystack = (const unsigned char *) haystack_start;
    const unsigned char *needle = (const unsigned char *) needle_start;
    const unsigned char *h = NULL;
    const unsigned char *n = NULL;
    size_t x = needle_len;

    /* The first occurrence of the empty string is deemed to occur at
 *     the beginning of the string.  */
    if (needle_len == 0)
        return (void *) haystack_start;

    /* Sanity check, otherwise the loop might search through the whole
 *         memory.  */
     if (haystack_len < needle_len)
       return NULL;

    for (; *haystack && haystack_len--; haystack++) {

        x = needle_len;
        n = needle;
        h = haystack;

        if (haystack_len < needle_len)
            break;

        if ((*haystack != *needle) || ( *haystack + needle_len != *needle + needle_len))
            continue;

        for (; x ; h++ , n++) {
            x--;

            if (*h != *n) 
                break;

           if (x == 0)
            return (void *)haystack;
        }
    }

    return NULL;
}
