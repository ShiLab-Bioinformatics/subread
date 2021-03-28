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

/***************************************************************

	The SamBam_reg2bin function was derived from the BAM
    specification. (The SAM Format Specication Working
    Group, September 7, 2011)

  ***************************************************************/
  
  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "subread.h"
#include "core.h"
#include "HelperFunctions.h"
#include "gene-algorithms.h"
#include "sambam-file.h"
#include "input-files.h"
#define Z_DEFAULT_MEM_LEVEL  8

int SamBam_fetch_next_chunk(SamBam_FILE *fp)
{
	int x1, room =  SAMBAM_INPUT_STREAM_SIZE - ( fp -> input_binary_stream_write_ptr - fp -> input_binary_stream_read_ptr); 

	if(room < 65536)
		return -1;


	for(x1=0; x1 < fp->input_binary_stream_write_ptr - fp -> input_binary_stream_read_ptr; x1 ++)
	{
		fp -> input_binary_stream_buffer [x1] = fp -> input_binary_stream_buffer [x1 + fp->input_binary_stream_read_ptr - fp -> input_binary_stream_buffer_start_ptr];
	}
	fp -> input_binary_stream_buffer_start_ptr = fp->input_binary_stream_read_ptr;

	char * in_buff = malloc(65537);
	unsigned int real_len = 0;
	int ret, have = 0;
	
	while (1){
			int nchunk = 0;
			ret = PBam_get_next_zchunk(fp -> os_file, in_buff, 65536, & real_len);
			if(ret > 0)
				nchunk = SamBam_unzip(fp -> input_binary_stream_buffer + fp->input_binary_stream_write_ptr - fp -> input_binary_stream_read_ptr + have, 65536 , in_buff , ret, 0);
			else if(ret == -2){
				SUBREADputs("ERROR: BAM format is broken.");
				return -2;
			}

			//printf("RET=%d; CHK=%d\n", ret, nchunk);

			if(nchunk>0)
				have += nchunk; 
			if(have > 3000) break;
			if(feof(fp->os_file)){
				fp->is_eof=1;
				break;
			}
	}
	free(in_buff);

	fp -> input_binary_stream_write_ptr += have;

	return have;

}

int is_paired_end_BAM(char * fn){
	FILE * fp = fopen(fn, "rb");
	if(!fp) return 0;
	unsigned char h2[2];
	int retvv = fread(h2, 1,2,fp);
	if(retvv <2) return 0;
	int is_bam = (h2[0]==31) && (h2[1]==139);
	fclose(fp);

	SamBam_FILE * sambam_reader;
	char fline[3000];
	sambam_reader = SamBam_fopen(fn, is_bam?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM);
	while(1){
		char * is_ret = SamBam_fgets(sambam_reader, fline, 2999, 1);
		if(!is_ret) return 0;
		if(fline[0]!='@') break;
	}
	char * tok = NULL;
	strtok_r(fline, "\t", &tok);
	int flags = atoi(strtok_r(NULL, "\t", &tok));
	SamBam_fclose(sambam_reader);
	return flags & 1;
}
SamBam_FILE * SamBam_fopen(char * fname , int file_type)
{
	SamBam_FILE * ret = (SamBam_FILE *)malloc(sizeof(SamBam_FILE));
	memset(ret, 0, sizeof(SamBam_FILE));
	ret -> file_type = file_type;

	if(file_type ==SAMBAM_FILE_SAM) 
	{
		ret -> os_file = f_subr_open(fname, "rb");
		if(!ret -> os_file)
		{
			free(ret);
			return NULL;
		}

		ret -> header_length = 0;
		char nch0=-1, nchold=-1;

		while(!feof(ret -> os_file)){
			nchold = nch0;
			nch0=fgetc(ret -> os_file);

			if(nchold=='\n' && nch0!='@') break;

			if(nch0!='@' && ret -> header_length == 0){
				ret -> header_length = 0;
				break;
			}
			ret -> header_length++;
		}
		//SUBREADprintf("The SAM header has %d bytes\n", ret -> header_length);
		fseek(ret -> os_file,0,SEEK_SET);
	}
	else
	{
		ret -> os_file = f_subr_open(fname, "rb");
		if(ret -> os_file == NULL)
		{
			free(ret);
			return NULL;
		}
		unsigned char first_ch = fgetc(ret->os_file);
		unsigned char second_ch = fgetc(ret->os_file);

		if(first_ch!=31 || second_ch!=139)
		{
			free(ret);
			SUBREADprintf("Not a BAM file! %d %d\n", first_ch, second_ch);
			return NULL;
		}

		fseek(ret->os_file, 0, SEEK_SET);

		ret -> input_binary_stream_buffer = (char *)malloc(SAMBAM_INPUT_STREAM_SIZE);
		ret -> input_binary_stream_read_ptr = 0;
		ret -> input_binary_stream_write_ptr = 0;
		ret -> input_binary_stream_buffer_start_ptr = 0;

		ret -> bam_file_stage = BAM_FILE_STAGE_HEADER;
		ret -> is_eof = 0;
		
		SB_FETCH(ret);

		if(SB_EOF(ret))
		{
			free(ret->input_binary_stream_buffer);
			free(ret);
			SUBREADprintf("FEOF 0.\n");
			return NULL;

		}

		int magic_4 = 0;
		memcpy(&magic_4 , SB_READ(ret), 4);
		SB_RINC(ret, 4);

		if(magic_4 != 21840194) // this number is the four bytes of "BAM1"
		{
			free(ret->input_binary_stream_buffer);
			free(ret);
			SUBREADprintf("FEOF 4 == %d.\n", magic_4);
			return NULL;
		}


		int l_text = 0;
		memcpy(&l_text, SB_READ(ret), 4);
		SB_RINC(ret, 4);

		ret -> bam_file_next_section_start = ret -> input_binary_stream_read_ptr + l_text;
		ret -> header_length =  ret -> bam_file_next_section_start;
	}
	return ret;
}

char cigar_op_char(int ch)
{
	if(ch<9)
		return "MIDNSHP=X"[ch];
	else
	{
		SUBREADprintf("Unknwon opt was found in the CIGAR string: '%c'.\n", ch);
		return 'M';
	}
}

char read_int_char(int ch)
{
	assert(ch<16);
	return "=ACMGRSVTWYHKDBN"[ch];
}

void SamBam_fclose(SamBam_FILE * fp)
{
	if(fp->file_type==SAMBAM_FILE_SAM)
	{
		fclose(fp->os_file);
		free(fp);
	}
	else
	{
		fclose(fp->os_file);
		free(fp -> input_binary_stream_buffer);
		free(fp -> bam_chro_table);
		free(fp);
	}
}

int SamBam_feof(SamBam_FILE * fp)
{
	if(fp->file_type ==SAMBAM_FILE_SAM) 
		return feof(fp->os_file);
	else return SB_EOF(fp); 
}

void SamBam_read_ref_info(SamBam_FILE * ret)
{
	unsigned int ref_info_size = 0;
	ret ->bam_chro_table_size = 0;
	//printf("CKK0\n");

	SB_FETCH(ret);
	if(SB_EOF(ret))
		return;

	//printf("CKK1\n");

	memcpy(&ref_info_size, SB_READ(ret),4);
	SB_RINC(ret, 4);

	int xk1;
	ret -> bam_chro_table = malloc(sizeof(SamBam_Reference_Info) * ref_info_size);
	for(xk1=0;xk1<ref_info_size;xk1++)
	{
		SB_FETCH(ret);
	
		if(SB_EOF(ret))
			break;

		int ref_name_len = 0;
		memcpy(&ref_name_len, SB_READ(ret),4);
		SB_RINC(ret, 4);

		int ref_readin_len = min(ref_name_len, BAM_MAX_CHROMOSOME_NAME_LEN-1);
		int ref_skip_len = ref_name_len - ref_readin_len;

		memcpy(ret -> bam_chro_table[xk1].chro_name, SB_READ(ret), ref_readin_len);
		ret -> bam_chro_table[xk1].chro_name[ref_readin_len] = 0;
		SB_RINC(ret, ref_readin_len + ref_skip_len);

		memcpy(&(ret -> bam_chro_table[xk1].chro_length), SB_READ(ret),4);
		SB_RINC(ret, 4);

		//SUBREADprintf("CHRO[%d] : %s [%d]\n", xk1+1, ret -> bam_chro_table[xk1].chro_name , ret -> bam_chro_table[xk1].chro_length);
	}
	ret ->bam_chro_table_size = ref_info_size;
}

char * SamBam_fgets(SamBam_FILE * fp, char * buff , int buff_len, int seq_needed)
{
	if(fp->file_type==SAMBAM_FILE_SAM){
		char * ret = fgets(buff, buff_len, fp->os_file);
		int strlenbuff = 0;
		if(ret){
			strlenbuff = strlen(buff);
			///if(buff[0]!='@')
			//	if(strlenbuff < 100)SUBREADprintf("WRONG LOAD: '%s'\n", buff);
		}
		if(strlenbuff < 1 || ret == NULL) return NULL;
		else{
			if(ret[strlenbuff-1]!='\n')
			{
				while(1)
				{
					int ch = getc(fp->os_file);
					if(ch == EOF || ch == '\n')break;
				}
				ret[strlenbuff-1] = '\n';
			}
			if(fp -> is_paired_end < 10){
				if(buff[0]!='@'){
					int tabs = 0,x1=0, tmpi=0;
					for(x1 = 0; x1 < strlenbuff; x1++){
						if(buff[x1] == '\t'){
							if(tabs == 1){
		//						SUBREADprintf("TMPI_SAM = %d\n", tmpi);
								fp -> is_paired_end = 10 + (tmpi & 1);
								break;
							} else tabs ++;
						}else{
							if(tabs == 1)tmpi = tmpi * 10 + buff[x1]-'0';
						}
					}
				}
			}

			if(buff[0] == '@'){
				fp -> header_length = ftello(fp->os_file) + strlenbuff+1;
			}
			return ret;
		}
	}
	else
	{
		int xk1;
		// decrypt the BAM mess.
		if(fp-> bam_file_stage == BAM_FILE_STAGE_HEADER)
		{
			char nch;
			xk1=0;
			SB_FETCH(fp);
			if(SB_EOF(fp))
				return NULL;

			while(1)
			{
				if(fp -> input_binary_stream_read_ptr >= fp -> bam_file_next_section_start)
					break;

				SB_FETCH(fp);
				nch = *(SB_READ(fp));
				SB_RINC(fp,1);

				//printf("%c", nch);
				if(nch == '\r')continue;
				if(SB_EOF(fp)||nch == '\n' || nch <0) break;
				if(xk1 < buff_len-2)
				{
					buff[xk1]=nch;
					xk1++;
				}
			}

			buff[xk1]='\n';
			buff[xk1+1]=0;

			//SUBREADprintf("BUFF=%s\n========================================================================\n\n", buff);
			//printf("RL=%d , PTR %d >? RECORD_START %d\n\n\n", xk1, fp -> input_binary_stream_read_ptr , fp -> bam_file_next_section_start);

			if(fp -> input_binary_stream_read_ptr >= fp -> bam_file_next_section_start)
			{
				SamBam_read_ref_info(fp);
				fp -> bam_file_stage = BAM_FILE_STAGE_ALIGNMENT;
				fp -> header_length = fp-> input_binary_stream_read_ptr;
			}
			return buff;
		}
		else
		{
			SamBam_Alignment *aln = &fp->aln_buff;
			int chunk_ptr = 0;
			SB_FETCH(fp);
			if(SB_EOF(fp)) return NULL;

			fp -> is_paired_end = 10 + ((*(fp -> input_binary_stream_buffer + fp -> input_binary_stream_read_ptr - fp -> input_binary_stream_buffer_start_ptr + 18)) & 1);
			//SUBREADprintf("FLAG=%d\n",  *(fp -> input_binary_stream_buffer + fp -> input_binary_stream_read_ptr - fp -> input_binary_stream_buffer_start_ptr + 18) );
			int text_len = PBam_chunk_gets(SB_READ(fp) , &chunk_ptr, fp -> input_binary_stream_write_ptr - fp -> input_binary_stream_read_ptr , fp -> bam_chro_table, buff , buff_len, aln, seq_needed);
			SB_RINC(fp, chunk_ptr);

			if(text_len>0) return buff;
			return NULL;
		}
	}
}



int PBam_get_next_zchunk(FILE * bam_fp, char * buffer, int buffer_length, unsigned int * real_len)
{
	unsigned char ID1, ID2, CM, FLG;
	unsigned short XLEN;
	int BSIZE=-1, rlen, is_file_broken = 0, rrtv;

	if(feof(bam_fp)) return -1;

	rrtv = fread(&ID1, 1, 1, bam_fp);
	if(rrtv <1) return -1;
	rrtv = fread(&ID2, 1, 1, bam_fp);
	if(rrtv <1) return -1;
	rrtv = fread(&CM, 1, 1, bam_fp);
	if(rrtv <1) return -1;
	rrtv = fread(&FLG, 1, 1, bam_fp);
	if(rrtv <1) return -1;
	if(feof(bam_fp)) return -1;

	if(ID1!=31 || ID2!=139 || CM!=8 || FLG!=4)
	{
		//SUBREADprintf("4CHR = %d, %d, %d, %d\n", ID1, ID2, CM, FLG);
		return -1;
	}
	fseeko(bam_fp, 6, SEEK_CUR);
	rrtv = fread(&XLEN,1, 2, bam_fp );
	if(rrtv < 2) return -1;

	int XLEN_READ = 0;
	while(1)
	{
		unsigned char SI1, SI2;
		unsigned short SLEN, BSIZE_MID;
		
		rrtv = fread(&SI1, 1, 1, bam_fp);
		if(rrtv <1) return -1;
		rrtv = fread(&SI2, 1, 1, bam_fp);
		if(rrtv <1) return -1;
		rlen = fread(&SLEN, 2, 1, bam_fp);
		if(rlen < 1) is_file_broken = 1;

		if(SI1==66 && SI2== 67 && SLEN == 2)
		{
			rrtv = fread(&BSIZE_MID, 1,2 , bam_fp);
			if(rrtv <2) return -1;
			BSIZE = BSIZE_MID;
		}
		else	fseeko(bam_fp, SLEN, SEEK_CUR);
		XLEN_READ += SLEN + 4;
		if(XLEN_READ>=XLEN) break;
	}

	if(BSIZE>19)
	{
		int CDATA_LEN = BSIZE - XLEN - 19;
		int CDATA_READING = min(CDATA_LEN, buffer_length);
		rrtv = fread(buffer, 1, CDATA_READING, bam_fp);
		if(rrtv < CDATA_READING) return -1;
		if(CDATA_READING<CDATA_LEN){
			SUBREADputs("ERROR: buffer insufficient");
			return -1;
		}
		fseeko(bam_fp, 4, SEEK_CUR);
		rlen = fread(&real_len, 4, 1, bam_fp);//Input SIZE (length of uncompressed data);uint32_t
		if(rlen < 1) is_file_broken = 1;

	//	SUBREADprintf("read_data=%u\n", CDATA_LEN);
		if(is_file_broken){
			SUBREADputs("ERROR: the input BAM file is broken.");
		}
		return is_file_broken?-2:CDATA_READING;
	}
	else
		return -1;
}


// returns 0 if the header finished.
// returns 1 if the header is going on.
// returns -1 if error.
int PBam_chunk_headers(char * chunk, int *chunk_ptr, int chunk_len, SamBam_Reference_Info ** bam_chro_table, int * table_size, int * table_items, int * state, int * header_txt_remainder, int * reminder_byte_len)
{

	if((*state)  == 0)
	{
		unsigned int header_txt_len ;
		if(0!=memcmp("BAM\x1",chunk + (*chunk_ptr),4))
			return -1;
		(*chunk_ptr)+=4;	// MAGIC
		(*state) = 1;

		memcpy(&header_txt_len, chunk + (*chunk_ptr),4);
		(*chunk_ptr)+=4;	
		if(header_txt_len + 8 < chunk_len)
		{
			(* state) = 2;
			(*chunk_ptr) += header_txt_len;
		}
		else
		{
			(* state) = 1;
			(* header_txt_remainder) = header_txt_len - (chunk_len - 8); 
			return 1;
		} 
	}

	if((*state) == 1)
	{
		if((*header_txt_remainder)<chunk_len)
		{
			(*state) = 2;
			(*chunk_ptr) += (*header_txt_remainder);
		}
		else if((*header_txt_remainder)==chunk_len)
		{
			(*state) = 2;
			return 1;
		}
		else	
		{
			(* header_txt_remainder) -= (chunk_len);
			return 1;
		}
	}

	if((*state) == 2 || (*state == 3))
	{
		int chrs, remainder_chrs;
		if((*state)==2)
		{
			memcpy(&chrs, chunk + (*chunk_ptr),4); 
			(*chunk_ptr)+=4;

			remainder_chrs = chrs;
		}
		else	remainder_chrs = (*header_txt_remainder);

		while((*chunk_ptr) < chunk_len && remainder_chrs>0)
		{
			int chro_name_len;
			unsigned int chro_len;
			(*reminder_byte_len) = chunk_len - (*chunk_ptr);

			if( (*chunk_ptr) < chunk_len-4)
			{
				memcpy(&chro_name_len, chunk + (*chunk_ptr),4);
				(*chunk_ptr)+=4;
				if( (*chunk_ptr) <= chunk_len-chro_name_len-4)
				{
					char * chro_name = chunk + (*chunk_ptr);
					(*chunk_ptr)+=chro_name_len;
					memcpy(&chro_len, chunk + (*chunk_ptr),4);
					(*chunk_ptr)+=4;

					(*reminder_byte_len) =0;

					//todo: insert item
					if(0==(* table_items))
					{
						(*table_size) = 50;
						(*bam_chro_table) = malloc(sizeof(SamBam_Reference_Info)*50);
					}
					else if((*table_size) <= (* table_items))
					{
						(*table_size) *= 2;
						(*bam_chro_table) = realloc((*bam_chro_table),sizeof(SamBam_Reference_Info)*(*table_size));
					}

					SamBam_Reference_Info * new_event = (*bam_chro_table) + (* table_items);
					if(strlen(chro_name)>=BAM_MAX_CHROMOSOME_NAME_LEN) chro_name[BAM_MAX_CHROMOSOME_NAME_LEN-1]=0;
					strcpy(new_event->chro_name, chro_name);
					new_event -> chro_length = chro_len;

					(* table_items)++;
					//SUBREADprintf("CHRO %d/%d added\n", (* table_items),(remainder_chrs));
					remainder_chrs --;
				}
				else break;
			}
			else break;

		}

		if(remainder_chrs)
		{
			(*state) = 3;
			(*header_txt_remainder) = remainder_chrs;
			return 1;
		}
		else{
			(*state) = 4;
			return 0;
		}
	}
	return -1;
}

int convert_BAM_binary_to_SAM( SamBam_Reference_Info * chro_table, char * bam_bin, char * sam_txt ){
	int bin_len = 0;
	memcpy(&bin_len, bam_bin, 4);
	bin_len += 4;

	int sam_ptr = 0, tmpint = 0;
	sam_ptr += sprintf(sam_txt + sam_ptr, "%s\t", bam_bin+36);

	memcpy(&tmpint, bam_bin + 16 ,4);
	sam_ptr += sprintf(sam_txt + sam_ptr, "%d\t", (tmpint >> 16) & 0xffff);
	int cigar_opts = tmpint & 0xffff;

	memcpy(&tmpint, bam_bin + 4  ,4);
	int r1chro = tmpint;
	sam_ptr += sprintf(sam_txt + sam_ptr, "%s\t", tmpint<0?"*":chro_table[tmpint].chro_name);
	memcpy(&tmpint, bam_bin + 8  ,4);
	sam_ptr += sprintf(sam_txt + sam_ptr, "%d\t", tmpint+1);
	memcpy(&tmpint, bam_bin + 12 ,4);
	sam_ptr += sprintf(sam_txt + sam_ptr, "%d\t", (tmpint >> 8) & 0xff);
	int name_len = tmpint & 0xff;
	int cigar_i;
	for(cigar_i = 0; cigar_i < cigar_opts; cigar_i ++){
		unsigned int cigarint=0;
		memcpy(&cigarint, bam_bin + name_len + 36 + cigar_i * 4,4);
		sam_ptr += sprintf(sam_txt + sam_ptr, "%u%c", cigarint >> 4, "MIDNSHP=X"[cigarint&0xf]);
	}
	sam_ptr += sprintf(sam_txt + sam_ptr, "%s\t", cigar_i<1?"*":"");
	
	memcpy(&tmpint, bam_bin + 24, 4);
	//SUBREADprintf("CHRO_IDX=%d\n", tmpint);
	sam_ptr += sprintf(sam_txt + sam_ptr, "%s\t", tmpint<0?"*":((tmpint == r1chro)?"=":chro_table[tmpint].chro_name));
	
	memcpy(&tmpint, bam_bin + 28, 4);
	sam_ptr += sprintf(sam_txt + sam_ptr, "%d\t", tmpint+1);

	memcpy(&tmpint, bam_bin + 32, 4);
	sam_ptr += sprintf(sam_txt + sam_ptr, "%d\t", tmpint);

	int seq_len;
	memcpy(&seq_len, bam_bin + 20,4);
	int seqi, flex_ptr=name_len + 36 +  cigar_opts * 4;
	for(seqi=0; seqi<seq_len; seqi++){
		sam_txt[sam_ptr++]="=ACMGRSVTWYHKDBN"[ (bam_bin[flex_ptr] >> ( 4*!(seqi %2) )) & 15 ];
		if(seqi %2) flex_ptr++;
	}
	sam_txt[sam_ptr++]='\t';
	if(seqi %2) flex_ptr++;
	for(seqi=0; seqi<seq_len; seqi++){
		unsigned char nch = (unsigned char) bam_bin[flex_ptr++];
		if(nch!=0xff||seqi == 0)
				sam_txt[sam_ptr++]=nch==0xff?'*':(nch+33);
	}
	sam_txt[sam_ptr++]='\t';

	while(flex_ptr < bin_len){
		sam_txt[sam_ptr++]=bam_bin[flex_ptr++];
		sam_txt[sam_ptr++]=bam_bin[flex_ptr++];
		sam_txt[sam_ptr++]=':';
		char tagtype = bam_bin[flex_ptr++];

		if(tagtype == 'B'){
			char elemtype = bam_bin[flex_ptr++];
			int elem_no=0, is_int_type = 0, type_bytes = 0, is_signed = 0;
			memcpy(&elem_no, bam_bin + flex_ptr, 4);
			flex_ptr += 4;
			sam_txt[sam_ptr++]='B';
			sam_txt[sam_ptr++]=':';
			sam_txt[sam_ptr++]=elemtype;
			sam_txt[sam_ptr++]=',';

			if(elemtype == 'i' || elemtype == 'I'){
				is_int_type = 1;
				type_bytes = 4;
				is_signed = elemtype == 'i' ;
			}else if(elemtype == 's' || elemtype == 'S'){
				is_int_type = 1;
				type_bytes = 2;
				is_signed = elemtype == 's' ;
			}else if(elemtype == 'c' || elemtype == 'C'){
				is_int_type = 1;
				type_bytes = 1;
				is_signed = elemtype == 's' ;
			}else if(elemtype == 'f'){
				type_bytes = 4;
			}

			int elemi;
			for(elemi =0; elemi < elem_no; elemi++){
				if(is_int_type){
					int tagval = 0;
					memcpy(&tagval,  bam_bin + flex_ptr, type_bytes);
					long long printv = is_signed?tagval:( (unsigned int) tagval );
					if(elemtype == 'i') sam_ptr += sprintf(sam_txt + sam_ptr, "%d,", (int)printv);
					if(elemtype == 'I') sam_ptr += sprintf(sam_txt + sam_ptr, "%u,", (unsigned int)printv);

					if(elemtype == 's') sam_ptr += sprintf(sam_txt + sam_ptr, "%d,", (short)printv);
					if(elemtype == 'S') sam_ptr += sprintf(sam_txt + sam_ptr, "%u,", (unsigned short)printv);

					if(elemtype == 'c') sam_ptr += sprintf(sam_txt + sam_ptr, "%d,", (char)printv);
					if(elemtype == 'C') sam_ptr += sprintf(sam_txt + sam_ptr, "%u,", (unsigned char)printv);
				}else{
					float tagval = 0;
					memcpy(&tagval,  bam_bin + flex_ptr, type_bytes);
					sam_ptr += sprintf(sam_txt + sam_ptr, "%f,", tagval);
				}
				flex_ptr += type_bytes;
			}

			sam_txt[sam_ptr-1] = '\t';
			sam_txt[sam_ptr] = 0;
		}else{
				int is_int_type = 0, is_float_type = 0, type_bytes = 0, is_string_type = 0, is_char_type = 0, is_signed = 0;
				if(tagtype == 'i' || tagtype == 'I'){
					is_int_type = 1;
					type_bytes = 4;
					is_signed = tagtype == 'i' ;
				}else if(tagtype == 's' || tagtype == 'S'){
					is_int_type = 1;
					type_bytes = 2;
					is_signed = tagtype == 's' ;
				}else if(tagtype == 'c' || tagtype == 'C'){
					is_int_type = 1;
					type_bytes = 1;
					is_signed = tagtype == 's' ;
				}else if(tagtype == 'f'){
					is_float_type = 1;
					type_bytes = 4;
				}else if(tagtype == 'Z' || tagtype == 'H'){
					is_string_type = 1;
					while(bam_bin[flex_ptr+(type_bytes ++)]);
				}else if(tagtype == 'A'){
					is_char_type = 1;
					type_bytes = 1;
				}


				sam_txt[sam_ptr++]=is_int_type?'i':tagtype;
				sam_txt[sam_ptr++]=':';

				if(is_int_type){
					int tagval = 0;
					memcpy(&tagval,  bam_bin + flex_ptr, type_bytes);
					long long printv = is_signed?tagval:( (unsigned int) tagval );
					#ifdef __MINGW32__
					sam_ptr += sprintf(sam_txt + sam_ptr, "%I64d\t", printv);
					#else
					sam_ptr += sprintf(sam_txt + sam_ptr, "%lld\t", printv);
					#endif
				}else if(is_string_type){
					// type_bytes includes \0
					memcpy(sam_txt + sam_ptr, bam_bin + flex_ptr, type_bytes -1);
					sam_txt[ sam_ptr + type_bytes -1 ] = '\t';

					//sam_txt[ sam_ptr + type_bytes +1]=0;
					//SUBREADprintf("STR_LEN=%d\tSTR=%s\n", type_bytes-1, sam_txt + sam_ptr);
					sam_ptr += type_bytes;
				}else if(is_float_type){
					float tagval = 0;
					memcpy(&tagval,  bam_bin + flex_ptr, type_bytes);
					sam_ptr += sprintf(sam_txt + sam_ptr, "%f\t", tagval);
				}else if(is_char_type){
					sam_txt[ sam_ptr++ ] = bam_bin[flex_ptr];
					sam_txt[ sam_ptr++ ] = '\t';
				}
				flex_ptr += type_bytes;
		}
	}

	sam_txt[sam_ptr-1]=0; //last '\t' 
	return sam_ptr-1;
}

int PBam_chunk_gets(char * chunk, int *chunk_ptr, int chunk_limit, SamBam_Reference_Info * bam_chro_table, char * buff , int buff_len, SamBam_Alignment*aln, int seq_needed)
{
	int xk1;
	// decrypt the BAM mess.
	unsigned int block_size;
	if((*chunk_ptr) +4> chunk_limit) return -1;

	memcpy(&block_size, chunk+(*chunk_ptr), 4);
	//SUBREADprintf("PBSIZE=%u\n", block_size);
	(*chunk_ptr)+=4;
	unsigned int next_start = block_size+(*chunk_ptr);

	int ref_id;
	memcpy(&ref_id, chunk+(*chunk_ptr), 4);
	(*chunk_ptr)+=4;

	if(ref_id == -1) aln -> chro_name = NULL;
	else aln -> chro_name = bam_chro_table[ref_id].chro_name; 

	memcpy(&(aln -> chro_offset), chunk+(*chunk_ptr), 4);
	(*chunk_ptr)+=4;

	unsigned int comb1;
	memcpy(&comb1, chunk+(*chunk_ptr), 4);
	(*chunk_ptr)+=4;

	aln -> mapping_quality = 0xff & (comb1 >> 8);

	unsigned int comb2;
	memcpy(&comb2, chunk+(*chunk_ptr), 4);
	(*chunk_ptr)+=4;

	aln -> flags = 0xffff&(comb2 >> 16);

	unsigned int read_len;
	memcpy(&read_len, chunk+(*chunk_ptr), 4);

	(*chunk_ptr)+=4;

	unsigned int mate_ref_id;
	memcpy(&mate_ref_id, chunk+(*chunk_ptr), 4);
	(*chunk_ptr)+=4;

	if(mate_ref_id == -1) aln -> mate_chro_name = NULL;
	else aln -> mate_chro_name = bam_chro_table[mate_ref_id].chro_name; 

	memcpy(&(aln -> mate_chro_offset), chunk+(*chunk_ptr), 4);
	(*chunk_ptr)+=4;

	memcpy(&(aln -> templete_length), chunk+(*chunk_ptr), 4);
	(*chunk_ptr)+=4;

	int read_name_len = comb1 & 0xff;
	assert(read_name_len < BAM_MAX_READ_NAME_LEN);

	memcpy(aln -> read_name, chunk+(*chunk_ptr), read_name_len);
	aln -> read_name[read_name_len] = 0;
	(*chunk_ptr)+=read_name_len;

	int cigar_ops = comb2 & 0xffff;
	aln -> cigar[0]=0; 
	for(xk1=0; xk1<cigar_ops;xk1++)
	{
		char cigar_piece_buf[BAM_MAX_CIGAR_LEN];
		unsigned int cigar_piece;

		if((*chunk_ptr) +4 > chunk_limit) return -1;
		memcpy(&cigar_piece,  chunk+(*chunk_ptr),4);
		(*chunk_ptr)+=4;

		sprintf(cigar_piece_buf, "%u%c", cigar_piece>>4, cigar_op_char(cigar_piece&0xf));
		if(strlen(cigar_piece_buf)+strlen(aln->cigar)<BAM_MAX_CIGAR_LEN-1)
			strcat(aln->cigar, cigar_piece_buf);
		else
		{
			SUBREADprintf("WARNING: cigar string is too long to the buffer.\n");
			SUBREADprintf("Please only use the compressed BAM format.\n");
			assert(0);
			return -1;
		}
	}

	char read_2_seq = 0;
	int seq_qual_bytes = read_len + (read_len /2)+(read_len%2);

	if(seq_needed)
		memcpy( aln-> buff_for_seq, chunk+(*chunk_ptr), seq_qual_bytes);
	(*chunk_ptr) += seq_qual_bytes;

	char extra_tags [CORE_ADDITIONAL_INFO_LENGTH];
	extra_tags[0]=0;
	int extra_len = 0;
	while( (*chunk_ptr) < next_start)
	{
		char extag[2];
		char extype;
		int delta, need_tag = 1;
		memcpy(extag,  chunk+(*chunk_ptr), 2);
		extype = chunk[2+(*chunk_ptr)];
		(*chunk_ptr)+=3;
		//fprintf(stderr, "COL_EXTYPE: %c\n", extype);
		if(extype == 'Z' || extype == 'H')
		{
			delta = 0;
			// 'Z' columns are NULL-terminated.
			while(chunk[delta + (*chunk_ptr)]) delta++;
			delta += 1;
		}
		else if(extype == 'A' || extype == 'c' || extype=='C') delta=1;
		else if(extype == 'i' || extype=='I' || extype == 'f') delta=4;
		else if(extype == 's' || extype=='S') delta=2;
		else if(extype == 'B') 
		{
			extype = chunk[(*chunk_ptr)];
		//	fprintf(stderr, "B_EXTYPE: %c\n", extype);

			(*chunk_ptr)++;
			if(extype == 'A' || extype=='Z') delta=1;
			else if(extype == 'c' || extype=='C') delta=1;
			else if(extype == 'i' || extype=='I' || extype == 'f') delta=4;
			else if(extype == 's' || extype=='S') delta=2;
			else break;

			int array_len;
			need_tag = 0;
			memcpy(&array_len, chunk+(*chunk_ptr), 4);
			(*chunk_ptr)+=4;
			delta *= array_len;
		}
		else{
		//	fprintf(stderr, "NO_EXTYPE: %c\n", extype);
			break;
		}

		if(need_tag){
			if(extype == 'c' || extype=='C' || extype == 'i' || extype=='I' || extype == 's' || extype=='S'){
				int tmpi = 0;
				memcpy(&tmpi, chunk+(*chunk_ptr),delta);
				if(tmpi >= 0 && extra_len < CORE_ADDITIONAL_INFO_LENGTH - 18){
					int sret = sprintf(extra_tags + strlen(extra_tags), "\t%c%c:i:%d", extag[0], extag[1], tmpi);
					extra_len += sret;
				}
			}else if(extype == 'Z'){
				if(extra_len < CORE_ADDITIONAL_INFO_LENGTH - 7 - delta){
					sprintf(extra_tags + strlen(extra_tags), "\t%c%c:Z:", extag[0], extag[1]);
					extra_len += 6;
					*(extra_tags + strlen(extra_tags)+delta-1) = 0;
					memcpy(extra_tags + strlen(extra_tags), chunk + (*chunk_ptr), delta - 1);
					extra_len += delta - 1;
				}
			}else if(extype == 'A'){
				if(extra_len < CORE_ADDITIONAL_INFO_LENGTH - 8){
					int sret = sprintf(extra_tags + strlen(extra_tags), "\t%c%c:A:%c", extag[0], extag[1], *(chunk + *chunk_ptr) );
					extra_len += sret;
				}
			}
		}

		if((*chunk_ptr) + delta > chunk_limit) return -1;
		(*chunk_ptr)+=delta;
		
	}

	if(next_start > chunk_limit) return -1;
	(*chunk_ptr) = next_start;

	if(seq_needed)
	{
		for(xk1=0;xk1<read_len;xk1++)
		{
			if(xk1 %2 == 0){
				read_2_seq = aln-> buff_for_seq[xk1/2];
			}
			if(xk1 < BAM_MAX_READ_LEN)
				aln -> sequence[xk1] = read_int_char(0xf&(read_2_seq >> (xk1%2?0:4)));
		}
		aln -> sequence[min(BAM_MAX_READ_LEN-1,read_len)] = 0;
		if(read_len >= BAM_MAX_READ_LEN-1)
			SUBREADprintf("WARNING: read is too long to the buffer\n");

		
		for(xk1=0;xk1<read_len;xk1++)
		{
			read_2_seq = aln -> buff_for_seq[(read_len /2)+(read_len%2) + xk1] ;
			if(xk1 < BAM_MAX_READ_LEN)
				aln -> seq_quality[xk1] = 33+read_2_seq;
		}
		aln -> seq_quality[min(BAM_MAX_READ_LEN-1,read_len)] = 0;
		if(aln -> seq_quality[0]==' ')
			strcpy(aln -> seq_quality, "*");
	}
	else
	{
		aln -> sequence[0]='N';
		aln -> sequence[1]=0;
		aln -> seq_quality[0]='#';
		aln -> seq_quality[1]=0;
	}

	char * chro_name = "*";
	char * cigar = "*";
	unsigned int chro_offset = 0;

	if(aln -> chro_name){
		chro_name = aln -> chro_name;
		chro_offset = aln -> chro_offset+1;
		if(aln -> cigar[0])
			cigar = aln -> cigar;
	}

	char * mate_chro_name = "*";
	unsigned int mate_chro_offset = 0;
	if(aln -> mate_chro_name)
	{
		if(aln -> mate_chro_name == chro_name) mate_chro_name = "=";
		else
			mate_chro_name = aln -> mate_chro_name;
		mate_chro_offset = aln -> mate_chro_offset+1;
	}

	long long int templete_length = aln -> templete_length;


	//fprintf(stderr, "HN_TAG=%d\n", nh_val	);

	#ifdef __MINGW32__
	int plen = snprintf(buff, buff_len-1, "%s\t%u\t%s\t%u\t%d\t%s\t%s\t%u\t%I64d\t%s\t%s%s\n%c", aln -> read_name, aln -> flags , chro_name, chro_offset, aln -> mapping_quality, cigar, mate_chro_name, mate_chro_offset, templete_length, aln -> sequence , aln -> seq_quality, extra_tags, 0);
	#else
	int plen = snprintf(buff, buff_len-1, "%s\t%u\t%s\t%u\t%d\t%s\t%s\t%u\t%lld\t%s\t%s%s\n%c", aln -> read_name, aln -> flags , chro_name, chro_offset, aln -> mapping_quality, cigar, mate_chro_name, mate_chro_offset, templete_length, aln -> sequence , aln -> seq_quality, extra_tags, 0);
	#endif
	//fprintf(stderr,"%s", buff);

	return plen;
}


int PBum_load_header(FILE * bam_fp, SamBam_Reference_Info** chro_tab, char * remainder_reads_data , int * remainder_reads_data_len)
{
	char * CDATA = malloc(80010);
	char * PDATA = malloc(1000000);

	int chro_tab_size = 0, chro_tab_items = 0, chro_tab_state = 0, header_remainder = 0, remainder_byte_len = 0, bam_is_broken = 0; 
	z_stream strm;
	while(1)
	{
		unsigned int real_len = 0;
		int rlen = PBam_get_next_zchunk(bam_fp,CDATA,80000, & real_len);
		if(rlen<0){
			bam_is_broken = (rlen == -2);
			if(bam_is_broken){
				SUBREADprintf("BAM file format error.\n");
				free(CDATA);
				free(PDATA);
				return -1;
			}
			break;
		}

		strm.zalloc = Z_NULL;
		strm.zfree = Z_NULL;
		strm.opaque = Z_NULL;
		strm.avail_in = 0;
		strm.next_in = Z_NULL;
		int ret = inflateInit2(&strm, SAMBAM_GZIP_WINDOW_BITS);
		if (ret != Z_OK)
		{
			free(CDATA);
			free(PDATA);
			return -1;
		}
		strm.avail_in = (unsigned int)rlen;
		strm.next_in = (unsigned char *)CDATA;


		strm.avail_out = 1000000 - remainder_byte_len;
		strm.next_out = (unsigned char *)(PDATA + remainder_byte_len);
		ret = inflate(&strm, Z_FINISH);
		int have = 1000000 - strm.avail_out;
		int PDATA_ptr=0;

		inflateEnd(&strm);

		ret = PBam_chunk_headers(PDATA, &PDATA_ptr, have, chro_tab, &chro_tab_size, &chro_tab_items, &chro_tab_state, &header_remainder,&remainder_byte_len);
		memcpy(PDATA , PDATA + have - remainder_byte_len, remainder_byte_len);
		if(ret<0)
		{
			SUBREADprintf("Header error.\n");
			free(CDATA);
			free(PDATA);
			return -1;
		}
		else if(ret == 0)
		{
			//SUBREADprintf("Header loaded = %d\n", (chro_tab_items));
			remainder_byte_len=0;
		}
		if(chro_tab_state>3){
			if(remainder_reads_data && PDATA_ptr < have)
			{
				memcpy(remainder_reads_data , PDATA + PDATA_ptr, have - PDATA_ptr);
				(*remainder_reads_data_len) =  have - PDATA_ptr ;
			}
			break;
		}
	}
	free(CDATA);
	free(PDATA);
	return 0;
}

int test_bamview(int argc, char ** argv)
{
	if(argc>1)
	{
		SamBam_FILE * fp = SamBam_fopen(argv[1], SAMBAM_FILE_BAM);
		assert(fp);
		/*
		while(1)
		{
			char buf[3000];
			char * buf2 = SamBam_fgets(fp,buf, 3000);
			//SUBREADprintf(">>%s<<\n",buf);
			//if(buf2)
			//	fwrite(buf,strlen(buf), 1, stdout);
			//else break;
		}
		*/
		SamBam_fclose(fp);
	}
	return 0;
}

int SamBam_writer_create(SamBam_Writer * writer, char * BAM_fname, int threads, int sort_reads_by_coord, int is_temp_BAM_file, char * tmpfname)
{
	memset(writer, 0, sizeof(SamBam_Writer));

	if(BAM_fname)
	{
		writer -> bam_fp = f_subr_open(BAM_fname, "wb");
		if(sort_reads_by_coord){
			char tname[MAX_FILE_NAME_LENGTH];
			sprintf(tname, "%s.bai", BAM_fname);
			writer -> BAI_fp = f_subr_open(tname, "wb");
			worker_master_mutex_init(&writer->sorted_notifier, threads);
		}
		if(!writer -> bam_fp) return -1;
	}
	#ifdef MAKE_STANDALONE
	else
		writer -> bam_fp = stdout;
	#endif

	writer -> threads = threads;
	writer -> sort_reads_by_coord = sort_reads_by_coord;
	writer -> fastest_compression = is_temp_BAM_file;
	writer -> compressed_chunk_buffer = malloc(70000); 
	writer -> chunk_buffer = malloc(70000); 
	writer -> chunk_buffer_max_size = 70000;
	strcpy(writer -> tmpf_prefix, tmpfname);

	if(threads>= 2){
		int x1;
		writer -> threads_chunk_buffer = malloc(sizeof(char *) * threads) ;
		writer -> threads_chunk_buffer_compressed = malloc(sizeof(char *) * threads) ;
		writer -> threads_chunk_buffer_used = malloc(sizeof(long long) * threads);
		writer -> threads_output_stream = malloc(sizeof(z_stream) * threads);
		writer -> threads_chunk_buffer_max_size = malloc(sizeof(long long) * threads);
		memset(writer -> threads_chunk_buffer_used, 0, sizeof(long long) * threads);
		for(x1 = 0; x1 < threads ; x1++){
			writer -> threads_chunk_buffer [x1] = malloc(70000);
			writer -> threads_chunk_buffer_compressed [x1] = malloc(70000);
			writer -> threads_chunk_buffer_max_size[x1] = 70000;
		}
		subread_init_lock(&writer -> thread_bam_lock);
	}
	writer -> chromosome_name_table = HashTableCreate(1603);
	writer -> chromosome_id_table = HashTableCreate(1603);
	writer -> chromosome_len_table = HashTableCreate(1603);
	writer -> header_plain_text_buffer = malloc(100000000);
	writer -> header_plain_text_buffer_max = 100000000;
	writer -> header_plain_text_buffer_used = 0;

	//memset(writer -> header_plain_text_buffer , 0 , 100000000);
	HashTableSetHashFunction(writer -> chromosome_name_table , fc_chro_hash);
	HashTableSetKeyComparisonFunction(writer -> chromosome_name_table , fc_strcmp_chro);
	HashTableSetDeallocationFunctions(writer -> chromosome_name_table , free, NULL);

	return 0;
}

void SamBam_writer_chunk_header(SamBam_Writer * writer, int compressed_size)
{

	// the four magic characters
	fputc(31,  writer -> bam_fp);
	fputc(139,  writer -> bam_fp);
	fputc(8,  writer -> bam_fp);
	fputc(4,  writer -> bam_fp);

	time_t time_now = 0;
	fwrite(&time_now,4,1, writer -> bam_fp);

	int tmp_i;
	// Extra flags and OS
	fputc(0,  writer -> bam_fp);
	fputc(0xff,  writer -> bam_fp); 

	// Extra length
	tmp_i = 6;
	fwrite(&tmp_i,2,1, writer -> bam_fp);


	// SI1 and SI2 magic numbers, and SLEN
	fputc(66,  writer -> bam_fp);
	fputc(67,  writer -> bam_fp);
	tmp_i = 2;
	fwrite(&tmp_i,2,1, writer -> bam_fp);
	tmp_i = compressed_size + 19 + 6;
	fwrite(&tmp_i,2,1, writer -> bam_fp);
}

unsigned int SamBam_CRC32(char * dat, int len)
{
	unsigned int crc0 = crc32(0, NULL, 0);
	unsigned int ret = crc32(crc0, (unsigned char *)dat, len);
	return ret;
}

void SamBam_writer_add_chunk(SamBam_Writer * writer, int thread_no)
{
	int compressed_size ; 
	unsigned int CRC32;
	z_stream * this_stream = thread_no < 0 ? &writer ->output_stream:writer -> threads_output_stream + thread_no;
	long long * this_buffer_used = thread_no < 0 ? &writer ->chunk_buffer_used : writer -> threads_chunk_buffer_used + thread_no;
	char * this_buffer = thread_no < 0 ? writer ->chunk_buffer: writer -> threads_chunk_buffer[thread_no];
	char * this_compressed_chunk_buffer = thread_no < 0 ? writer ->compressed_chunk_buffer : writer -> threads_chunk_buffer_compressed[thread_no];

	if(*this_buffer_used < 1){
		// magic block with 0 byte data(EOF)
		subread_lock_occupy(&writer -> thread_bam_lock);
		fwrite( "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00",1,28,writer -> bam_fp);
		writer -> current_BAM_pos = ftello(writer -> bam_fp);
		subread_lock_release(&writer -> thread_bam_lock);
		return;
	}

	this_stream -> avail_out = 70000;
	this_stream -> avail_in = (* this_buffer_used);
	CRC32 = SamBam_CRC32(this_buffer , * this_buffer_used);

	this_stream -> zalloc = Z_NULL;
	this_stream -> zfree = Z_NULL;
	this_stream -> opaque = Z_NULL;

	deflateInit2(this_stream, writer -> fastest_compression? SAMBAM_COMPRESS_LEVEL_FASTEST : SAMBAM_COMPRESS_LEVEL_NORMAL, Z_DEFLATED,
		SAMBAM_GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
	
	this_stream -> next_in = (unsigned char *)this_buffer;
	this_stream -> next_out = (unsigned char *)this_compressed_chunk_buffer;

	deflate(this_stream, Z_FINISH);
	deflateEnd(this_stream);

	compressed_size = 70000 - this_stream -> avail_out;

	subread_lock_occupy(&writer -> thread_bam_lock);
	SamBam_writer_chunk_header(writer, compressed_size);
	int chunk_write_size = fwrite(this_compressed_chunk_buffer, 1, compressed_size, writer -> bam_fp);
	fwrite(&CRC32 , 4, 1, writer -> bam_fp);
	fwrite(this_buffer_used , 4, 1, writer -> bam_fp);
	writer -> current_BAM_pos = ftello(writer -> bam_fp);
	subread_lock_release(&writer -> thread_bam_lock);

	if(chunk_write_size < compressed_size){
		if(!writer -> is_internal_error)SUBREADputs("ERROR: no space left in the output directory.");
		writer -> is_internal_error = 1;
	}

	(* this_buffer_used) = 0;
}

double sambam_t1 = 0;

void SamBam_writer_write_header(SamBam_Writer * writer)
{
	int header_ptr=0, header_block_start = 0;
	while(header_ptr < writer->header_plain_text_buffer_used)
	{
		if(( header_ptr - header_block_start > 55000 || header_ptr >= writer->header_plain_text_buffer_used-1) && writer -> header_plain_text_buffer[header_ptr] == '\n')
		{
			writer -> chunk_buffer_used = 0;
			if(header_block_start == 0)	// the very first block
			{
				memcpy(writer -> chunk_buffer, "BAM\1",4);
				writer -> chunk_buffer_used  = 4;
				memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used, &writer -> header_plain_text_buffer_used, 4);
				writer -> chunk_buffer_used += 4;
		
			}

			memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used , writer -> header_plain_text_buffer + header_block_start, header_ptr - header_block_start+1);
			writer -> chunk_buffer_used +=  header_ptr - header_block_start + 1;
			SamBam_writer_add_chunk(writer, -1);
			header_block_start = header_ptr + 1;
		}
		header_ptr ++;
	}

	free(writer -> header_plain_text_buffer);
	writer -> header_plain_text_buffer = NULL;

	// reference sequences
	writer -> chunk_buffer_used = 0;
	memcpy(writer -> chunk_buffer, & writer -> chromosome_name_table -> numOfElements, 4);
	writer -> chunk_buffer_used = 4;

	for( header_ptr=0 ;  header_ptr < writer -> chromosome_name_table -> numOfElements ; header_ptr ++)
	{
		//printf("D=%d\n", writer -> chromosome_id_table -> numOfElements);
		char * chro_name = HashTableGet(writer -> chromosome_id_table, NULL + 1 + header_ptr);
		unsigned int chro_len = HashTableGet(writer -> chromosome_len_table, NULL + 1 + header_ptr) - NULL - 1;
		assert(chro_name);
		int chro_name_len = strlen(chro_name)+1;

		memcpy(writer -> chunk_buffer +  writer -> chunk_buffer_used , &chro_name_len, 4);
		writer -> chunk_buffer_used += 4;

		strcpy(writer -> chunk_buffer +  writer -> chunk_buffer_used , chro_name);
		writer -> chunk_buffer_used += chro_name_len;

		memcpy(writer -> chunk_buffer +  writer -> chunk_buffer_used , &chro_len, 4);
		writer -> chunk_buffer_used += 4;

		if(header_ptr ==  writer -> chromosome_name_table -> numOfElements - 1 || writer -> chunk_buffer_used > 55000)
		{
			SamBam_writer_add_chunk(writer, -1);
			writer -> chunk_buffer_used = 0;
		}
	}

}

void SamBam_writer_finalise_one_thread(SamBam_Writer * writer);
void SamBam_writer_sort_bins_to_BAM(SamBam_Writer * writer);

int SamBam_writer_close(SamBam_Writer * writer)
{
	if(writer -> writer_state == 0)	// no reads were added
	{
		if(writer -> header_plain_text_buffer)
			SamBam_writer_write_header(writer);
	} else SamBam_writer_finalise_one_thread(writer);

	if(writer -> sort_reads_by_coord){
		SamBam_writer_sort_bins_to_BAM(writer);
		worker_master_mutex_destroy(&writer->sorted_notifier);
	}
	
	writer -> chunk_buffer_used = 0;
	SamBam_writer_add_chunk(writer, -1);
//	fputc(0, writer -> bam_fp);

	writer -> output_stream.next_in= NULL;
	writer -> output_stream.avail_in= 0;
	writer -> output_stream.next_out= NULL;
	writer -> output_stream.avail_out= 0;

	free(writer -> chunk_buffer);
	free(writer -> compressed_chunk_buffer);
	if(writer -> threads >=2){
		int x1;
		for(x1 = 0 ; x1 < writer -> threads; x1 ++){
			free(writer -> threads_chunk_buffer[x1]);
			free(writer -> threads_chunk_buffer_compressed[x1]);
		}
		free(writer -> threads_output_stream);
		free(writer -> threads_chunk_buffer);
		free(writer -> threads_chunk_buffer_compressed);
		free(writer -> threads_chunk_buffer_used);
	}
	HashTableDestroy(writer -> chromosome_name_table);
	HashTableDestroy(writer -> chromosome_id_table);
	HashTableDestroy(writer -> chromosome_len_table);
	#ifdef MAKE_STANDALONE
	if(stdout != writer -> bam_fp)
	#endif
	fclose(writer -> bam_fp);
	if(writer -> BAI_fp!=NULL) fclose(writer -> BAI_fp);

	return 0;
}

int SamBam_writer_add_header(SamBam_Writer * writer, char * header_text, int add_chro)
{
	int new_text_len = strlen(header_text);

	if(writer -> header_plain_text_buffer_max <= writer -> header_plain_text_buffer_used + new_text_len + 1)
	{
		//return 0;
		writer -> header_plain_text_buffer_max *=2;
		writer -> header_plain_text_buffer = realloc(writer -> header_plain_text_buffer ,  writer -> header_plain_text_buffer_max);
		//printf("REAL: %d : %llX\n",writer -> header_plain_text_buffer_max, (long long ) writer -> header_plain_text_buffer);
	}

	strcpy(writer -> header_plain_text_buffer + writer -> header_plain_text_buffer_used, header_text);
	writer -> header_plain_text_buffer_used += new_text_len;
	strcpy(writer -> header_plain_text_buffer + writer -> header_plain_text_buffer_used, "\n");
	writer -> header_plain_text_buffer_used ++;

//	SUBREADprintf("ADJHEADERXCD BIN=%s\n", header_text);
	if(add_chro && memcmp(header_text, "@SQ",3)==0)
	{
		char * chro = NULL;
		int chro_len = -1;
		char * toktmp = NULL;
		char * ret_tmp = strtok_r(header_text, "\t", &toktmp);

		while(1){
			if(!ret_tmp) break;

			if(memcmp(ret_tmp,"SN:", 3)==0) chro = ret_tmp + 3;
			else if(memcmp(ret_tmp,"LN:", 3)==0) chro_len = atoi(ret_tmp + 3);

			ret_tmp = strtok_r(NULL, "\t", &toktmp);
		}

		if(chro && (chro_len>0))
			SamBam_writer_add_chromosome(writer, chro, chro_len, 0);
		
	}

	//if(writer -> header_plain_text_buffer_used %97==0) printf("MV=%d\n",writer -> header_plain_text_buffer_used);

	return 0;
}

int SamBam_writer_add_chromosome(SamBam_Writer * writer, char * chro_name, unsigned int chro_length, int add_header)
{
	unsigned int chro_id = writer -> chromosome_name_table -> numOfElements;

	//SUBREADprintf("ADJHEADER_CHRO %s of %u\n", chro_name, chro_length);
	//assert(strlen(chro_name) < 30);

	char * chro_name_space = malloc(strlen(chro_name)+1);
	strcpy(chro_name_space , chro_name);
	HashTablePut(writer -> chromosome_name_table, chro_name_space, NULL+1+chro_id);
	HashTablePut(writer -> chromosome_id_table, NULL+1+chro_id, chro_name_space);
	HashTablePut(writer -> chromosome_len_table, NULL+1+chro_id, NULL + 1 + chro_length);

	if(add_header)
	{
		char * line_buf = malloc(1000);
		snprintf(line_buf,999, "@SQ\tSN:%s\tLN:%u", chro_name , chro_length);
		SamBam_writer_add_header(writer, line_buf, 0);
		free(line_buf);
	}

	return 0;
}


int SamBam_compress_cigar(char * cigar, int * cigar_int, int * ret_coverage, int max_secs)
{
	int tmp_int=0;
	int cigar_cursor = 0, num_opt = 0;
	int coverage_len = 0;
	(* ret_coverage) = 0;

	if(cigar[0]=='*') return 0;
	
	while(1)
	{
		char nch = cigar[cigar_cursor++];
		if(!nch)break;
		if(isdigit(nch))
		{
			tmp_int = tmp_int*10+(nch-'0');
		}
		else
		{
			int int_opt=0;
			if(nch == 'M' || nch == 'N' || nch == 'D') coverage_len += tmp_int;
			//if(nch == 'M' ||nch == 'D' || nch == '=' || nch == 'X') coverage_len += tmp_int;
			for(; int_opt<8; int_opt++) if("MIDNSHP=X"[int_opt] == nch)break;
			cigar_int[num_opt ++] = (tmp_int << 4) | int_opt; 
			tmp_int = 0;
			//SUBREADprintf("CIGARCOM: %d-th is %c\n", num_opt, nch);
			if(num_opt>=max_secs)break;
		}
	}

	(*ret_coverage) = coverage_len;
	return num_opt;
}

void SamBam_read2bin(char * read_txt, char * read_bin)
{
	int bin_cursor = 0, txt_cursor = 0;

	while(1)
	{
		char nch = read_txt[txt_cursor++];
		if(!nch)break;
		int fourbit;
		for(fourbit=0;fourbit<15;fourbit++) if("=ACMGRSVTWYHKDBN"[fourbit] == nch)break;

		if(bin_cursor %2 == 0)  read_bin[bin_cursor/2] =  fourbit<<4;
		else read_bin[bin_cursor/2] |=  fourbit;

		bin_cursor++;
	}
}

int SamBam_compress_additional(char * additional_columns, char * bin)
{
	int col_cursor = 0 , col_len = strlen(additional_columns);
	int bin_cursor = 0;

	while(col_cursor<col_len)
	{
		if(col_cursor==0 || additional_columns[col_cursor]=='\t')
		{
			if(additional_columns[col_cursor]=='\t') col_cursor++;

			bin[bin_cursor] = additional_columns[col_cursor];
			bin[bin_cursor+1] = additional_columns[col_cursor+1];

			char datatype = additional_columns[col_cursor+3];
			if(datatype=='i' || datatype == 'f')
			{
				int dig_len =0;
				while(additional_columns[dig_len+col_cursor+5] != '\t' && additional_columns[dig_len+col_cursor+5]) dig_len++;
				int val = 0;
				float fval = 0;
				if(datatype=='i') val = atoi(additional_columns+col_cursor+5);
				else val = atof(additional_columns+col_cursor+5);

				bin[bin_cursor+2]=datatype;
				memcpy(bin+bin_cursor+3, (datatype=='i')? ((void *)&val):((void *)&fval),4);
				bin_cursor += 3 + 4;
				col_cursor += 5 + dig_len;
			}
			else if(datatype=='Z' || datatype == 'H')
			{
				bin[bin_cursor+2]=datatype;
				bin_cursor +=3;
				int str_len = 0;
				col_cursor +=5;
				while(additional_columns[str_len+col_cursor] != '\t' && additional_columns[str_len+col_cursor])
				{
					bin[bin_cursor + str_len] = additional_columns[str_len+col_cursor];
					str_len++;
					if(bin_cursor + str_len > 780) break;
				}

				bin[bin_cursor + str_len] =0;

				bin_cursor += str_len + 1;
				col_cursor += str_len;
			}
			else if(datatype=='A')
			{
				bin[bin_cursor+2]='A';
				bin[bin_cursor+3]=additional_columns[col_cursor+5];
				col_cursor += 6;
				bin_cursor += 4;
			}
			else if(datatype=='B')
				//array
			{
				char celltype = additional_columns[col_cursor+5];
				int * items = (int *)(&bin[bin_cursor+4]);

				bin[bin_cursor+2]='B';
				bin[bin_cursor+3]=celltype;
				bin_cursor += 4 + 4;
				col_cursor += 7;

				(*items) = 0;

				int last_cursor = col_cursor;
				while(1){
					if(additional_columns[col_cursor] == ',' || additional_columns[col_cursor] == '\t' || additional_columns[col_cursor] == 0)
					{ // add new item 

						char cell_buff [30];
						if((col_cursor - last_cursor) < 29)
						{
							memcpy(cell_buff, additional_columns + last_cursor, (col_cursor - last_cursor));
							cell_buff[(col_cursor - last_cursor)] = 0;
							int intv = 0; float fltv = 0;
							if(celltype == 'i')intv = atoi(cell_buff);							
							else fltv = atof(cell_buff);
							if(bin_cursor < 780){
								memcpy(bin + bin_cursor, (celltype == 'i')?(void *)&intv:(void *)&fltv, 4);
								bin_cursor += 4;
								(*items) ++;
							}
						}
						last_cursor = col_cursor+1;
					}
					if(additional_columns[col_cursor] == '\t' || additional_columns[col_cursor] == 0)
						break;

					col_cursor++;
					
				}
				
			}

			if(bin_cursor>750) break;
			continue;
		}
		
		col_cursor++;
	}
	return bin_cursor;
}

int SamBam_reg2bin(int beg, int end)
{
	--end;
	if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
	if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
	if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
	if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
	if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
	return 0;
}

void SamBam_writer_finish_header( SamBam_Writer * writer ){
	if(writer -> writer_state == 0){
		if(writer -> header_plain_text_buffer)
		SamBam_writer_write_header(writer);
		writer -> writer_state = 10;
	}
}

#define FC_MAX_CIGAR_SECTIONS 96

// The caller has to free the returned pointer.
// It returns NULL if no such field.
char * duplicate_TAB_record_field(char * rline, int fno, int to_end){
	int i, fldi = 0, start_pos = -1;
	if(fno >=1)
		for(i=0;rline[i] && rline[i] != '\n'; i++){
			int nch = rline[i];
			if(nch == '\t'){
				fldi++;
				if(fldi == fno) start_pos = i+1;
				if(fldi == fno+1) break;
			}
		}
	else{
		start_pos = rline[0]>0?0:-1;
		if(start_pos == 0) for(i=0;rline[i] && rline[i] != '\n' && rline[i] != '\t'; i++);
		else i=-1;
	}
	int end_pos =i;
	if(to_end){
		end_pos = strlen(rline);
		if(end_pos<1)  return NULL;

		if(rline[end_pos-1]=='\n') end_pos--;
	}

	if(start_pos<0 || end_pos <= start_pos) return NULL;
	char * ret = malloc(end_pos - start_pos +1);
	memcpy(ret, rline+start_pos, end_pos -start_pos);
	ret[end_pos -start_pos] = 0;
	return ret;
}

int SamBam_writer_add_read_bin(SamBam_Writer * writer, int thread_no, char * rbin, int committable){
	char * this_chunk_buffer=NULL;
	srInt_64 *this_chunk_buffer_used=NULL;
	if(thread_no >= 0){
		if(writer -> sort_reads_by_coord && writer -> threads_chunk_buffer_max_size[thread_no] < writer -> threads_chunk_buffer_used[thread_no] + 12000){
			writer -> threads_chunk_buffer_max_size[thread_no] = writer -> threads_chunk_buffer_max_size[thread_no] * 7/4;
			writer -> threads_chunk_buffer[ thread_no ] = realloc(writer -> threads_chunk_buffer[ thread_no ], writer -> threads_chunk_buffer_max_size[thread_no]);
		}
		this_chunk_buffer = writer -> threads_chunk_buffer[ thread_no ];
		this_chunk_buffer_used = writer -> threads_chunk_buffer_used + thread_no;
	}else{
		if(writer -> sort_reads_by_coord && writer -> chunk_buffer_max_size < writer -> chunk_buffer_used + 12000){
			//SUBREADprintf("0THR_REALLOCATE MEM : %d -> %d\n", writer -> chunk_buffer_max_size,  writer -> chunk_buffer_max_size * 7/4);
			writer -> chunk_buffer_max_size = writer -> chunk_buffer_max_size * 7/4;
			writer -> chunk_buffer = realloc(writer -> chunk_buffer, writer -> chunk_buffer_max_size);
		}
		this_chunk_buffer = writer -> chunk_buffer;
		this_chunk_buffer_used = &writer -> chunk_buffer_used;
	}
	//SUBREADprintf("WRTXBIN: PTR %p of USED %lld at read %p\n", this_chunk_buffer, *this_chunk_buffer_used, rbin);
	int reclen=0;
	memcpy(&reclen, rbin,4);
	memcpy(this_chunk_buffer+(*this_chunk_buffer_used), rbin, reclen+4);
	(*this_chunk_buffer_used)+= reclen+4;

	if((*this_chunk_buffer_used)>55000 && committable && !writer -> sort_reads_by_coord)
		SamBam_writer_add_chunk(writer, thread_no);
	//SUBREADprintf("WRTXBIN: FIN WTR %p of USED %lld\n", this_chunk_buffer, *this_chunk_buffer_used);
	return 0;
}
int SamBam_writer_add_read_line(SamBam_Writer * writer, int thread_no, char * rline, int committable){
	char * read_name, * flag_str, *chro_name, *chro_position_str, * mapping_quality_str, * cigar, * next_chro_name, *next_chro_position_str, *temp_len_str, *read_text, *qual_text, *additional_columns;
	read_name = duplicate_TAB_record_field(rline, 0,0);
	flag_str = duplicate_TAB_record_field(rline, 1,0);
	chro_name = duplicate_TAB_record_field(rline, 2,0);
	chro_position_str = duplicate_TAB_record_field(rline, 3,0);
	mapping_quality_str = duplicate_TAB_record_field(rline, 4,0);
	cigar = duplicate_TAB_record_field(rline, 5,0);
	next_chro_name = duplicate_TAB_record_field(rline, 6,0);
	next_chro_position_str = duplicate_TAB_record_field(rline, 7,0);
	temp_len_str = duplicate_TAB_record_field(rline, 8,0);
	read_text = duplicate_TAB_record_field(rline, 9,0);
	qual_text = duplicate_TAB_record_field(rline, 10,0);
	additional_columns = duplicate_TAB_record_field(rline, 11,1);
	if(qual_text==NULL){
		SUBREADprintf("FATAL ERROR : bad read format: %s, %s, %s, %s\n", read_name, flag_str, chro_name, rline);
		return -1;
	}

	SamBam_writer_add_read(writer, thread_no, read_name, atoi(flag_str), chro_name, atoi(chro_position_str), atoi(mapping_quality_str), cigar, next_chro_name, atoi(next_chro_position_str), atoi(temp_len_str), strlen(read_text), read_text, qual_text, additional_columns, committable);

	if(additional_columns) free(additional_columns);
	free(qual_text);
	free(read_text);
	free(temp_len_str);
	free(next_chro_position_str);
	free(next_chro_name);
	free(cigar);
	free(mapping_quality_str);
	free(chro_position_str);
	free(chro_name);
	free(flag_str);
	free(read_name);
	return 0;
}

int SamBam_writer_add_read(SamBam_Writer * writer, int thread_no, char * read_name, unsigned int flags, char * chro_name, unsigned int chro_position, int mapping_quality, char * cigar, char * next_chro_name, unsigned int next_chro_position, int temp_len, int read_len, char * read_text, char * qual_text, char * additional_columns, int committable)
{
	assert(writer -> writer_state!=0);
	if(!qual_text || !read_text)	
	{
		SUBREADprintf("ERROR: sam file is incomplete.\n");
		return 1;
	}

	char additional_bin[1000];
	int cigar_opts[FC_MAX_CIGAR_SECTIONS], xk1, cover_length = 0;
	int cigar_opt_len = SamBam_compress_cigar(cigar, cigar_opts, & cover_length, FC_MAX_CIGAR_SECTIONS);
	int read_name_len = 1+strlen(read_name) ;
	int additional_bin_len = 0;
	if(additional_columns) additional_bin_len = SamBam_compress_additional(additional_columns, additional_bin);
	int record_length = 4 + 4 + 4 + 4 +  /* l_seq: */ 4 + 4 + 4 + 4 + /* read_name:*/ read_name_len + cigar_opt_len * 4 + (read_len + 1) /2 + read_len + additional_bin_len;


	char * this_chunk_buffer;
	srInt_64 * this_chunk_buffer_used;

	if(thread_no >= 0){
		if(writer -> sort_reads_by_coord && writer -> threads_chunk_buffer_max_size[thread_no] < writer -> threads_chunk_buffer_used[thread_no] + 12000){
			writer -> threads_chunk_buffer_max_size[thread_no] = writer -> threads_chunk_buffer_max_size[thread_no] * 7/4;
			writer -> threads_chunk_buffer[ thread_no ] = realloc(writer -> threads_chunk_buffer[ thread_no ], writer -> threads_chunk_buffer_max_size[thread_no]);
		}
		this_chunk_buffer = writer -> threads_chunk_buffer[ thread_no ];
		this_chunk_buffer_used = writer -> threads_chunk_buffer_used + thread_no;
	}else{
		if(writer -> sort_reads_by_coord && writer -> chunk_buffer_max_size < writer -> chunk_buffer_used + 12000){
			//SUBREADprintf("REALLOCATE MEM : %d -> %d\n", writer -> chunk_buffer_max_size,  writer -> chunk_buffer_max_size * 7/4);
			writer -> chunk_buffer_max_size = writer -> chunk_buffer_max_size * 7/4;
			writer -> chunk_buffer = realloc(writer -> chunk_buffer, writer -> chunk_buffer_max_size);
		}
		this_chunk_buffer = writer -> chunk_buffer;
		this_chunk_buffer_used = &writer -> chunk_buffer_used;
	}

	memcpy(this_chunk_buffer + (*this_chunk_buffer_used) , & record_length , 4);
	(*this_chunk_buffer_used) += 4;

	int bin = SamBam_reg2bin(chro_position -1, chro_position-1+cover_length);
	//if(chro_position>0)SUBREADprintf("CREATE_BIN=%d\n", bin);

	int refID = HashTableGet(writer -> chromosome_name_table, chro_name) - NULL - 1; 
	int bin_mq_nl = (bin<<16) | (mapping_quality << 8) | read_name_len ;
	int fag_nc = (flags<<16) | cigar_opt_len;
	int nextRefID = -1;

	if(next_chro_name[0] != '*' && next_chro_name[0]!='=')
		nextRefID = HashTableGet(writer -> chromosome_name_table, next_chro_name) - NULL - 1;
	else if(next_chro_name[0] == '=')
		nextRefID = refID;

	
	chro_position--;
	next_chro_position--;

	memcpy(this_chunk_buffer + (*this_chunk_buffer_used) , & refID , 4);
	(*this_chunk_buffer_used) += 4;
	memcpy(this_chunk_buffer + (*this_chunk_buffer_used) , & chro_position , 4);
	(*this_chunk_buffer_used) += 4;
	memcpy(this_chunk_buffer + (*this_chunk_buffer_used) , & bin_mq_nl , 4);
	(*this_chunk_buffer_used) += 4;
	memcpy(this_chunk_buffer + (*this_chunk_buffer_used) , & fag_nc , 4);
	(*this_chunk_buffer_used) += 4;
	memcpy(this_chunk_buffer + (*this_chunk_buffer_used) , & read_len , 4);
	(*this_chunk_buffer_used) += 4;
	memcpy(this_chunk_buffer + (*this_chunk_buffer_used) , & nextRefID , 4);
	(*this_chunk_buffer_used) += 4;
	memcpy(this_chunk_buffer + (*this_chunk_buffer_used) , & next_chro_position , 4);
	(*this_chunk_buffer_used) += 4;
	memcpy(this_chunk_buffer + (*this_chunk_buffer_used) , & temp_len , 4);
	(*this_chunk_buffer_used) += 4;
	strcpy(this_chunk_buffer + (*this_chunk_buffer_used) , read_name);
	(*this_chunk_buffer_used) += read_name_len;
	memcpy(this_chunk_buffer + (*this_chunk_buffer_used) , cigar_opts, 4*cigar_opt_len);
	(*this_chunk_buffer_used) += 4*cigar_opt_len;
	SamBam_read2bin(read_text  , this_chunk_buffer + (*this_chunk_buffer_used));
	(*this_chunk_buffer_used) += (read_len + 1) /2; 
	memcpy(this_chunk_buffer + (*this_chunk_buffer_used), qual_text, read_len);
	for(xk1=0; xk1<read_len; xk1++)
		this_chunk_buffer[(*this_chunk_buffer_used)+xk1] -= 33;
	
	(*this_chunk_buffer_used) += read_len; 
	memcpy(this_chunk_buffer + (*this_chunk_buffer_used), additional_bin, additional_bin_len);
	(*this_chunk_buffer_used) += additional_bin_len;


	if((*this_chunk_buffer_used)>55000 && committable && !writer -> sort_reads_by_coord)
		SamBam_writer_add_chunk(writer, thread_no);

	return 0;
}

int SamBam_unzip(char * out, int outlen , char * in , int inlen, int sync_only)
{
	z_stream strm;
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	strm.avail_in = 0;
	strm.next_in = Z_NULL;
	int ret = inflateInit2(&strm, SAMBAM_GZIP_WINDOW_BITS);
	if (ret != Z_OK)
		return -1;

	strm.avail_in = (unsigned int)inlen;
	strm.next_in = (unsigned char *)in;

	strm.avail_out = outlen;
	strm.next_out = (unsigned char *)out;
	ret = inflate(&strm, sync_only?Z_SYNC_FLUSH:Z_FINISH);
	if(ret != Z_STREAM_END && ret!=0)
	{
		inflateEnd(&strm);
		SUBREADprintf("DATA ERROR! code=%d\n", ret);
		assert(0);
		return -1;
	}
	int have = outlen - strm.avail_out;

	inflateEnd(&strm);
	//SUBREADprintf("DECOMPRESS GENERATED=%d\n", have);

	return have;
}


int SamBam_writer_sort_buff_one_compare(void * Lbin, void * Rbin, ArrayList * me){
	unsigned long long * Lv = Lbin, *Rv = Rbin;
	if(*Lv > *Rv) return 1;
	else if(*Lv < *Rv) return -1;
	return 0;
}
// return total reads in bin.
int SamBam_writer_sort_buff_one_write(SamBam_Writer * writer, char * bin, int binlen, int thread_id){
	int bin_cursor = 0;

	ArrayList* sort_linear_pos = ArrayListCreate(1000000);
	ArrayListSetDeallocationFunction(sort_linear_pos,  free);

	int ii_reads = 0;
	while(bin_cursor < binlen){
		int this_block_len = 0;
		memcpy(&this_block_len, bin + bin_cursor, 4);

		char * key_binpos = malloc(12);
		memcpy(key_binpos,   bin + bin_cursor+8, 4);
		memcpy(key_binpos+4,   bin + bin_cursor+4, 4);
		memcpy(key_binpos+8, &bin_cursor, 4);

		ArrayListPush(sort_linear_pos, key_binpos);
		bin_cursor +=4 + this_block_len;

		ii_reads ++;
	}
	ArrayListSort(sort_linear_pos, SamBam_writer_sort_buff_one_compare);

	char * nbin = NULL;
	if(binlen >0 && binlen < 0x7fffffff)nbin=malloc(binlen);
	int nb_cursor = 0, xx, wlen=0;
	for(xx=0; xx<ii_reads;xx++){
		int * binpos = ArrayListGet(sort_linear_pos, xx);
		int block_len = 0;
		memcpy(&block_len, bin + binpos[2] , 4);
		assert(block_len < 10000);
		memcpy(nbin + nb_cursor, bin + binpos[2], 4+block_len);
		nb_cursor += 4+block_len;
	}
	assert(binlen == nb_cursor);
	memcpy(bin, nbin, binlen);
	ArrayListDestroy(sort_linear_pos);

	char tmpfname[MAX_FILE_NAME_LENGTH+40];
	if(writer -> threads>1) subread_lock_occupy(&writer -> thread_bam_lock);
	sprintf(tmpfname, "%s-%06d.sortedbin", writer -> tmpf_prefix, writer -> sorted_batch_id++);
	if(writer -> threads>1) subread_lock_release(&writer -> thread_bam_lock);
	FILE * tofp  = fopen(tmpfname, "wb");
	if(tofp){
		if(binlen>0)wlen = fwrite(nbin, binlen,1, tofp);
		fclose(tofp);
	}
	free(nbin);
	if(wlen < 1 && binlen>0) {
		SUBREADprintf("ERROR: no space (%d bytes) in the temp directory (%s).\nThe program cannot run properly.\n", binlen, tmpfname);
		writer ->is_internal_error = 1;
		return -1;
	}
	return ii_reads;
}

unsigned long long SamBam_writer_sort_bins_to_BAM_FP_pos(FILE * fp){
	int intvals[3];
	unsigned long long ret = 0;
	int rlen = fread(intvals, sizeof(int),3,fp);
	if(rlen>0 && intvals[0]<10000){
		int chro_no, poss;
		chro_no=intvals[1];
		poss=intvals[2];
		ret =((1LLU*chro_no) << 32)| poss;
		if(ret == SUBREAD_MAX_ULONGLONG) ret = ret- 10 ;
		fseek(fp, -12, SEEK_CUR);
		return ret;
	}
	return SUBREAD_MAX_ULONGLONG;
}


int SamBam_writer_calc_cigar_span(char * bin){
	int cops = 0, rname_len = 0, ret = 0;
	memcpy(&cops, bin+12, 4);
	memcpy(&rname_len, bin+8, 4);
	rname_len = rname_len & 0xff;
	cops = cops & 0xffff;
	int ii;
	for(ii = 0; ii < cops ; ii++){
		unsigned int copt = 0;
		memcpy(&copt, bin+32+rname_len+4*ii, 4);
		int copt_char = copt & 0xf;
		unsigned int copt_len = copt >> 4;
		if(copt_char == 0 || copt_char == 2 || copt_char == 3 || copt_char == 7 || copt_char == 8) ret += copt_len;
	}

	return ret;
}

#define MAX_ALLOWED_GAP_IN_BAI_CHUNK 5 // the number of blocks, not the file positions.

void SamBam_writer_sort_bins_to_BAM_test_bins(SamBam_Writer * writer, HashTable * bin_tab, ArrayList * bin_list, ArrayList * win16k_list, int block_len, void ***last_chunk_ptr, int chro_no){
	int inbin_pos = writer -> chunk_buffer_used - block_len ; // point to the byte AFTER "block_len" int.
	int pos=0, bin_mq_nl = 0, binno=0;

	memcpy(&pos, writer -> chunk_buffer + inbin_pos + 4, 4);
	memcpy(&bin_mq_nl, writer -> chunk_buffer + inbin_pos + 8,4);
	binno = bin_mq_nl>>16;

	int cigar_span = SamBam_writer_calc_cigar_span(writer -> chunk_buffer + inbin_pos);

	int this_w16_no = (pos + cigar_span) >>14;	// WIN is calculated on 0-based pos.
	unsigned long long this_Vpos = writer -> this_bam_block_no<<16 | (inbin_pos-4);

	// if this read is after the maximum coordinate in the win16k list: all elements before last one and this one starts at this read.
	if(this_w16_no >=win16k_list->numOfElements){
		int bbi;
		for(bbi = win16k_list->numOfElements; bbi <=this_w16_no; bbi++)
			ArrayListPush(win16k_list, NULL+ this_Vpos);
	}

	// a read only belongs to ONE bin which is the smallest bin that can FULLY cover the full length of this read.
	ArrayList * this_bin_chunks = HashTableGet(bin_tab, NULL+binno+1);
	if(NULL == this_bin_chunks){
		this_bin_chunks = ArrayListCreate(5);
		HashTablePut(bin_tab, NULL+binno+1, this_bin_chunks);
		ArrayListPush(bin_list, NULL+binno);
	}

	int found = 0;
	// a bin is not necessarily continuous. Say, a top-level bin only contains a few reads (most reads a in low-level bins), but their locations are everywhere.
	if(this_bin_chunks -> numOfElements > 0){
		//SUBREADprintf("RESLOCS : %ld == %ld\n", this_Vpos, this_bin_chunks -> elementList [ this_bin_chunks -> numOfElements - 1 ] - NULL);
		long long diff = this_Vpos >>16;
		diff -=(this_bin_chunks -> elementList [ this_bin_chunks -> numOfElements - 1] - NULL)>>16;
		if(diff < MAX_ALLOWED_GAP_IN_BAI_CHUNK){
			this_bin_chunks -> elementList [ this_bin_chunks -> numOfElements -1] = NULL+this_Vpos + block_len+4;
			*last_chunk_ptr = this_bin_chunks -> elementList + this_bin_chunks -> numOfElements -1;
			found = 1;
		}
	}

	// if the last chunk in this bin isn't good to be extended (too far from the file location of the new read), a new chunk is created.
	if(!found){
		ArrayListPush(this_bin_chunks, NULL + this_Vpos);
		ArrayListPush(this_bin_chunks, NULL + this_Vpos + block_len+4);
		*last_chunk_ptr = this_bin_chunks -> elementList + this_bin_chunks -> numOfElements -1;
	}
}

void * SamBam_writer_sorted_compress(void * vptr0){
	void ** vptr = vptr0;
	SamBam_Writer * writer = vptr[0];
	int this_thread_no = vptr[1]-NULL;
	subread_lock_t * initlock = vptr[2];

	worker_thread_start(&writer->sorted_notifier, this_thread_no);
	subread_lock_release(initlock);
	free(vptr0);
	struct SamBam_sorted_compressor_st * me = writer -> writer_threads + this_thread_no;

	while(1){
		int termed = worker_wait_for_job(&writer->sorted_notifier, this_thread_no);
		if(termed)break;

		me->CRC32_plain = SamBam_CRC32(me->plain_text,  me->text_size);
		me->strm.next_in = (unsigned char*)me->plain_text;
		me->strm.avail_in = me->text_size;
		me->strm.next_out = (unsigned char*)me->zipped_bin;
		me->strm.avail_out = 70000;
		int deret = deflate(&me->strm,Z_FINISH);
		if(deret >=0){
			deflateReset(&me->strm);
			me->bin_size = 70000 - me->strm.avail_out;
			me->last_job_done = 1;
		}else{
			SUBREADprintf("Error: cannot compress BAM block #%d , which is %llu, had %d => 70000 [ %d ] bytes , return = %d\n", this_thread_no, me -> bam_block_no, me->text_size, me->bin_size, deret);
		}
	}
	return NULL;
}

void SamBam_thread_wait_merge_write(SamBam_Writer * writer, int thread_no){
	master_wait_for_job_done(&writer -> sorted_notifier, thread_no);

	if(writer -> writer_threads[thread_no].last_job_done){
		srInt_64 fpos = ftello(writer -> bam_fp);
		HashTablePut(writer -> block_no_p1_to_vpos_tab, NULL+1+writer -> writer_threads[thread_no].bam_block_no, NULL+fpos);

		SamBam_writer_chunk_header(writer, writer -> writer_threads[thread_no].bin_size);
		int rlen = fwrite( writer -> writer_threads[thread_no].zipped_bin ,1, writer -> writer_threads[thread_no].bin_size, writer -> bam_fp);
		if(rlen != writer -> writer_threads[thread_no].bin_size){
			SUBREADprintf("ERROR: cannot write output files.\n");
			assert(0);
		}
		rlen = fwrite(&writer -> writer_threads[thread_no].CRC32_plain , 4, 1, writer -> bam_fp);
		rlen = fwrite(&writer -> writer_threads[thread_no].text_size , 4, 1, writer -> bam_fp);

		writer -> writer_threads[thread_no].text_size = 0;
		writer -> writer_threads[thread_no].bin_size = 0;
		writer -> writer_threads[thread_no].bam_block_no = -1;
		writer -> writer_threads[thread_no].last_job_done = 0;
	}
}

void SamBam_thread_set_new_job(SamBam_Writer * writer, int thread_no){
	memcpy(writer -> writer_threads[thread_no].plain_text, writer ->chunk_buffer, writer->chunk_buffer_used);
	writer -> writer_threads[thread_no].text_size = writer->chunk_buffer_used;
	writer -> writer_threads[thread_no].bam_block_no = writer -> this_bam_block_no;
	writer -> chunk_buffer_used = 0;
	master_notify_worker(&writer -> sorted_notifier, thread_no);
}

void SamBam_writer_submit_sorted_compressing_task(SamBam_Writer * writer){
	SamBam_thread_wait_merge_write(writer, writer -> sorted_compress_this_thread_no);
	SamBam_thread_set_new_job(writer, writer -> sorted_compress_this_thread_no);

	writer -> sorted_compress_this_thread_no++;
	if(writer -> sorted_compress_this_thread_no == writer -> threads) writer -> sorted_compress_this_thread_no=0;
	writer -> this_bam_block_no ++;
}

void SamBam_writer_sort_bins_to_BAM_write_1R(SamBam_Writer * writer, FILE * fp, HashTable * bin_tab, ArrayList * bin_list, ArrayList * win16k_list, int chro_no){
	int block_len=0;
	int rlen = fread(&block_len, 4,1,fp);
	//SUBREADprintf("WBIN=%d\n", block_len);
	if(rlen<1 || block_len >= 10000){
		SUBREADprintf("ERROR: sorted bin files are broken. RLEN=%d , BLKLEN=%d\n", rlen, block_len);
		assert(rlen >=1);
	}

	memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used, &block_len, 4);
	writer ->chunk_buffer_used += 4;
	rlen = fread(writer -> chunk_buffer + writer -> chunk_buffer_used, 1, block_len ,fp);
	if(rlen < block_len){
		SUBREADprintf("ERROR: sorted bin files are broken.\n");
		assert(rlen >=block_len);
	}
	writer -> chunk_buffer_used+= rlen;

	void ** last_chunk_ptr = NULL;
	SamBam_writer_sort_bins_to_BAM_test_bins(writer, bin_tab, bin_list, win16k_list, block_len, &last_chunk_ptr, chro_no);

	if(writer -> chunk_buffer_used>55000){
		assert(writer -> chunk_buffer_used< 65000);
		SamBam_writer_submit_sorted_compressing_task(writer);
	}
}

#define SAMBAM_MERGE_MAX_FPS (350) 

void SamBam_writer_one_thread_merge_sortedbins(SamBam_Writer * writer){
	int new_bins = 0;

	if(writer -> sorted_batch_id > SAMBAM_MERGE_MAX_FPS){
		int merge_i;
		for(merge_i = 0; merge_i < writer -> sorted_batch_id; merge_i+= SAMBAM_MERGE_MAX_FPS){
			char tfp[MAX_FILE_NAME_LENGTH+50];
			int this_size = min(SAMBAM_MERGE_MAX_FPS, writer -> sorted_batch_id - merge_i);
			if(this_size < 2) break;

			//SUBREADprintf("CREATE_MERGE %d\n", new_bins + writer -> sorted_batch_id);
			sprintf(tfp , "%s-%06d.sortedbin", writer -> tmpf_prefix, writer -> sorted_batch_id + (new_bins++) );
			FILE * outbinfp = fopen(tfp,"wb");

			FILE ** sb_fps = malloc(sizeof(FILE *) * this_size);
			unsigned long long * current_min_fps = malloc(sizeof(unsigned long long) * this_size);
			unsigned long long current_min = SUBREAD_MAX_ULONGLONG; 
			int current_min_batch = -1, bii;

			for(bii = 0; bii < this_size; bii++){
				char tfpx[MAX_FILE_NAME_LENGTH+50];
				current_min_fps[bii] = SUBREAD_MAX_ULONGLONG;

				sprintf(tfpx , "%s-%06d.sortedbin", writer -> tmpf_prefix, bii + merge_i);
				sb_fps[bii] = fopen(tfpx,"rb");
				current_min_fps[bii] = SamBam_writer_sort_bins_to_BAM_FP_pos(sb_fps[bii]);
				if(current_min_fps[bii] < SUBREAD_MAX_ULONGLONG && current_min_fps[bii] < current_min){
					current_min = current_min_fps[bii];
					current_min_batch = bii;
				}
			}

			char * mtmp = malloc(10000);
			while(1){
				if(current_min_batch>-1){
					int blk_len = 0;
					int rrlen = fread(&blk_len, 4, 1, sb_fps[current_min_batch]);
					if(rrlen <1) assert(rrlen >0);
					assert(blk_len < 10000);
					rrlen = fread(mtmp, blk_len, 1, sb_fps[current_min_batch]);
					if(rrlen <1) assert(rrlen >0);
					fwrite(&blk_len,4,1, outbinfp);
					fwrite(mtmp, blk_len,1, outbinfp);

					current_min_fps[current_min_batch] = SamBam_writer_sort_bins_to_BAM_FP_pos(sb_fps[current_min_batch]);
					current_min = SUBREAD_MAX_ULONGLONG;
					current_min_batch = -1;
				}else break;

				for(bii = 0; bii < this_size; bii++){
					if(current_min_fps[bii] < SUBREAD_MAX_ULONGLONG && current_min_fps[bii] < current_min){
						current_min = current_min_fps[bii];
						current_min_batch = bii;
					}
				}
			}
			free(mtmp);

			for(bii = 0; bii < this_size; bii++){
				fclose(sb_fps[bii]);
				sprintf(tfp , "%s-%06d.sortedbin", writer -> tmpf_prefix, merge_i + bii);
				//SUBREADprintf("    DEL_MERGE %d\n", bii);
				unlink(tfp);
			}
			fclose(outbinfp);
		}
		
		writer -> sorted_batch_id += new_bins;
	}
}


#define SAMBAM_renew_BAItabs	{ if(bin_chunks_table != NULL) { HashTableDestroy(bin_chunks_table); ArrayListDestroy(bins_list); ArrayListDestroy(win16k_list) ; } \
	bin_chunks_table = HashTableCreate(10000); bins_list = ArrayListCreate(10000); win16k_list = ArrayListCreate(10000);\
	HashTableSetDeallocationFunctions(bin_chunks_table, NULL, (void(*)(void *))ArrayListDestroy); }

#define SAMBAM_write_empty_chros(nSTART, nEND) {int iikk; for(iikk=(nSTART); iikk<(nEND);iikk++){ int wwwlen = fwrite("\0\0\0\0\0\0\0\0", 1,8,writer -> BAI_fp); if(wwwlen < 8){ assert(wwwlen==8);}} }

int level_min_binno[] = {0, 1, 9, 73, 585, 4681};

int SamBam_writer_merge_chunks_compare(void * vL, void * vR, ArrayList * me){
	long long * lL = vL;
	long long * lR = vR;
	if( (*lL) > (*lR)) return 1;
	if( (*lL) < (*lR)) return -1;
	return 0;
}

void SamBam_writer_merge_chunks(ArrayList * chs){
	ArrayList * intvs = ArrayListCreate(chs -> numOfElements/2);
	ArrayListSetDeallocationFunction(intvs, free);
	int i;
	for(i=0; i<chs -> numOfElements; i+=2){
		long long * ll = malloc(sizeof(long long)*2);
		ll[0] = ArrayListGet(chs, i)-NULL;
		ll[1] = ArrayListGet(chs, i+1)-NULL;
		ArrayListPush(intvs, ll);
	}
	chs -> numOfElements = 0;

	ArrayListSort(intvs, SamBam_writer_merge_chunks_compare);
	long long * ll0 =  ArrayListGet(intvs, 0);
	ArrayListPush(chs, NULL+ll0[0]);
	ArrayListPush(chs, NULL+ll0[1]);

	for(i=1; i<intvs-> numOfElements; i++){
		long long last_end = ArrayListGet(chs, chs -> numOfElements-1) -NULL;
		long long * ll_this = ArrayListGet(intvs, i);

		long long dist = ll_this[0]>>16;
		dist -= last_end>>16;

		if(dist < MAX_ALLOWED_GAP_IN_BAI_CHUNK)
			chs-> elementList[chs -> numOfElements-1] = NULL+max(ll_this[1],last_end);
		else{
			ArrayListPush(chs, NULL+ll_this[0]);
			ArrayListPush(chs, NULL+ll_this[1]);
		}
	}

	//SUBREADprintf("MERGE CHUNKS : %ld -> %ld\n", intvs-> numOfElements, chs -> numOfElements/2);
	ArrayListDestroy(intvs);
}

void SamBam_writer_optimize_bins_level(HashTable *bin_tab, ArrayList *bin_arr, HashTable * new_tab, ArrayList * new_arrs, int this_level){
	int i,j, my_min_binno = level_min_binno[this_level], my_parent_min_binno = -1, my_child_min_binno = 999999;
	if(this_level>0) my_parent_min_binno = level_min_binno[this_level-1];
	if(this_level<5) my_child_min_binno = level_min_binno[this_level+1];


	// copy all parents and ancestors and children and descendant bins into new table & list
	for(i=0; i< bin_arr -> numOfElements; i++){
		int binno = ArrayListGet(bin_arr, i)-NULL;
		if(binno < my_min_binno || binno >= my_child_min_binno){
			ArrayList * old_bin_arr = HashTableGet(bin_tab, NULL+1+binno);
			if(old_bin_arr -> numOfElements > 1){
				HashTablePut(new_tab, NULL+1+binno, ArrayListDuplicate(old_bin_arr));
				ArrayListPush(new_arrs, NULL+binno);
			}
		}
	}

	// scan this_level bins. If small enough -> put to parent bin and delete it.
	// Otherwise put it into new tab&list.
	for(i=0; i< bin_arr -> numOfElements; i++){
		int binno = ArrayListGet(bin_arr, i)-NULL;
		if(binno < my_min_binno || binno >= my_child_min_binno )continue;
		ArrayList * small_bin_arr = HashTableGet(bin_tab, NULL+1+binno);
		if(small_bin_arr -> numOfElements < 2) continue;

		long long min_start_Voff = SUBREAD_MAX_LONGLONG;
		long long max_end_Voff = -1;

		for(j =0 ; j< small_bin_arr -> numOfElements; j+=2){
			long long this_start = ArrayListGet(small_bin_arr, j)-NULL; 
			long long this_end = ArrayListGet(small_bin_arr, j+1)-NULL; 
			min_start_Voff = min(min_start_Voff, this_start);
			max_end_Voff = max(max_end_Voff, this_end);
		}

		long long dist = max_end_Voff >> 16;
		dist -= min_start_Voff >>16;
		if(dist < MAX_ALLOWED_GAP_IN_BAI_CHUNK) {
			int parent_binno = my_parent_min_binno +((binno - my_min_binno)>>3);
			ArrayList * parent_bin_arr = HashTableGet(new_tab, NULL+1+parent_binno);
			if(NULL == parent_bin_arr){
				parent_bin_arr = ArrayListCreate(10);
				HashTablePut(new_tab, NULL+1+parent_binno, parent_bin_arr);
				ArrayListPush(new_arrs, NULL+parent_binno);
			}
			for(j =0 ; j< small_bin_arr -> numOfElements; j++)
				ArrayListPush(parent_bin_arr, ArrayListGet(small_bin_arr, j));

		} else {
			HashTablePut(new_tab, NULL+1+binno, ArrayListDuplicate(small_bin_arr));
			ArrayListPush(new_arrs, NULL + binno);
		}
	}

	// look into all parent_level bins, reduce its size 
	for(i=0; i< new_arrs -> numOfElements; i++){
		int binno = ArrayListGet(new_arrs, i)-NULL;
		if(binno >=my_min_binno || binno < my_parent_min_binno) continue;
		ArrayList * large_bin_arr = HashTableGet(new_tab, NULL+1+binno);
		SamBam_writer_merge_chunks(large_bin_arr);
	}

	//SUBREADprintf("Tab size : %ld => %ld\n",bin_tab -> numOfElements, new_tab -> numOfElements);
	HashTableDestroy(bin_tab);
	ArrayListDestroy(bin_arr);
}

void SamBam_writer_optimize_bins(HashTable *bin_tab, ArrayList *bin_arr, HashTable ** new_tab, ArrayList ** new_arrs){
	int this_level;
	for(this_level = 5; this_level > 2; this_level --){
	//	SUBREADprintf("OPTIMIZING AT LEVEL %d\n", this_level);
		HashTable * new_bin_tab = HashTableCreate(2000);
		HashTableSetDeallocationFunctions(new_bin_tab, NULL, (void(*)(void *))ArrayListDestroy);
		ArrayList * new_bin_list = ArrayListCreate(2000);

		SamBam_writer_optimize_bins_level(bin_tab, bin_arr, new_bin_tab, new_bin_list, this_level);
		bin_tab = new_bin_tab;
		bin_arr = new_bin_list;
		(*new_tab) = new_bin_tab;
		(*new_arrs) = new_bin_list;
	}
}

void SamBam_write_sorted_thread_collect(SamBam_Writer * writer){
	int thread_i;
	if(writer -> chunk_buffer_used>0) SamBam_writer_submit_sorted_compressing_task(writer);
	for(thread_i = 0; thread_i < writer -> threads; thread_i++){
		//SUBREADprintf("FINISHED %d-th CHRO, MERGING THR %d\n", old_chro_no, writer->sorted_compress_this_thread_no );
		SamBam_thread_wait_merge_write(writer, writer->sorted_compress_this_thread_no);
		writer->sorted_compress_this_thread_no ++;
		if(writer->sorted_compress_this_thread_no == writer -> threads) writer->sorted_compress_this_thread_no = 0;
	}
}

#define SAMBAM_Block2Vpos(v) (v)=((HashTableGet(writer -> block_no_p1_to_vpos_tab, NULL+(v>>16)+1)-NULL )<<16 )|(v & 0xffff);
void SamBam_write_BAI_for_1chr(SamBam_Writer * writer, HashTable ** bin_chunks_table, ArrayList ** bins_list, ArrayList ** win16k_list){
	SamBam_write_sorted_thread_collect(writer);

	if(1){
		HashTable * new_bin_tab = NULL;
		ArrayList * new_bin_list = NULL;
		SamBam_writer_optimize_bins(*bin_chunks_table ,*bins_list,&new_bin_tab,&new_bin_list);
		*bin_chunks_table = new_bin_tab;
		*bins_list = new_bin_list;
	}

	int n_bins = (*bins_list)->numOfElements;
	int win16_nos = (*win16k_list)->numOfElements;
	//SUBREADprintf("WIN_ITEMS=%d\t\t\tOLD_CHRO=%d  NEW_CHRO=%d  CURR_BATCH=%d\n", win16_nos, old_chro_no, chro_no , current_min_batch);

	fwrite(&n_bins, 4, 1, writer -> BAI_fp);
	int bini;
	for(bini = 0; bini < n_bins; bini ++){
		int bin_no = ArrayListGet(*bins_list, bini)-NULL;
		ArrayList * chunks_list = HashTableGet(*bin_chunks_table, NULL+ bin_no+1);

		int n_chunks = chunks_list -> numOfElements/2;
		fwrite(&bin_no, 4, 1, writer -> BAI_fp);
		fwrite(&n_chunks, 4, 1, writer -> BAI_fp);
		int chunk_i;
		for(chunk_i = 0; chunk_i < 2* n_chunks; chunk_i +=2){
			long long int Voff_st = ArrayListGet(chunks_list, chunk_i) - NULL;
			long long int Voff_en = ArrayListGet(chunks_list, chunk_i+1) - NULL;
			SAMBAM_Block2Vpos(Voff_st);
			SAMBAM_Block2Vpos(Voff_en);
			fwrite(&Voff_st, 8, 1, writer -> BAI_fp);
			fwrite(&Voff_en, 8, 1, writer -> BAI_fp);
		}
	}

	fwrite(&win16_nos, 4, 1, writer -> BAI_fp);
	for(bini = 0; bini < win16_nos; bini++){
		long long int Voff = ArrayListGet(*win16k_list , bini) - NULL;
		SAMBAM_Block2Vpos(Voff);
		fwrite(&Voff, 8, 1, writer -> BAI_fp);
	}

}

#define SAMBAM_reset_sorting_writer { HashTableRemoveAll( writer -> block_no_p1_to_vpos_tab ); }

// MUST be called by THREAD 0!
void SamBam_writer_sort_bins_to_BAM(SamBam_Writer * writer){
	SamBam_writer_one_thread_merge_sortedbins(writer);

	FILE ** sb_fps = malloc(sizeof(FILE *) * writer -> sorted_batch_id);
	unsigned long long * current_min_fps = malloc(sizeof(unsigned long long) * writer -> sorted_batch_id);
	int bii;
	unsigned long long current_min = SUBREAD_MAX_ULONGLONG; 
	int current_min_batch = -1;
	writer -> chunk_buffer_used = 0;

	for(bii = 0; bii < writer -> sorted_batch_id; bii++){
		char tfp[MAX_FILE_NAME_LENGTH+40];
		current_min_fps[bii] = SUBREAD_MAX_ULONGLONG;

		sprintf(tfp , "%s-%06d.sortedbin", writer -> tmpf_prefix, bii);
		sb_fps[bii] = fopen(tfp,"rb");
		if(sb_fps[bii]!=NULL){
			current_min_fps[bii] = SamBam_writer_sort_bins_to_BAM_FP_pos(sb_fps[bii]);
			if(current_min_fps[bii] < SUBREAD_MAX_ULONGLONG && current_min_fps[bii] < current_min){
				current_min = current_min_fps[bii];
				current_min_batch = bii;
			}
		}
	}

	ArrayList * bins_list = NULL; //	[Bin_Number_0, Bin_Number_1, ...]
	ArrayList * win16k_list = NULL;	// [ Voffset_0, Voffset_1, ... ]
	HashTable * bin_chunks_table=NULL;
	SAMBAM_renew_BAItabs;

	int chro_no = (int)(current_min >> 32);
	int old_chro_no = -1, last_written_BAI_chro = -1;;

	int wlen = fwrite("BAI\1", 4, 1, writer -> BAI_fp);
	if(wlen < 1) assert(wlen >0);
	int thread_i, n_ref = writer -> chromosome_name_table -> numOfElements;
	subread_lock_t initlocks[writer -> threads];
	wlen = fwrite(&n_ref, 4, 1, writer -> BAI_fp);
	SAMBAM_write_empty_chros(0, chro_no);

	writer -> block_no_p1_to_vpos_tab = HashTableCreate(100000);
	writer -> writer_threads = calloc(sizeof(struct SamBam_sorted_compressor_st), writer -> threads);
	for(thread_i = 0; thread_i < writer -> threads; thread_i++){
		memset(&(writer -> writer_threads[thread_i].strm),0,sizeof(z_stream));
		deflateInit2(&(writer -> writer_threads[thread_i].strm), writer -> fastest_compression? SAMBAM_COMPRESS_LEVEL_FASTEST : SAMBAM_COMPRESS_LEVEL_NORMAL , Z_DEFLATED,
			SAMBAM_GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
		subread_init_lock(initlocks+thread_i);
		subread_lock_occupy(initlocks+thread_i);

		void ** thparam = malloc(sizeof(void*)*3);
		thparam[0] = writer;
		thparam[1] = NULL + thread_i;
		thparam[2] = initlocks+thread_i;
		pthread_create(&writer -> writer_threads[thread_i].thread_stub, NULL,SamBam_writer_sorted_compress,thparam);
	}
	for(thread_i = 0; thread_i < writer -> threads; thread_i++){
		subread_lock_occupy(initlocks+thread_i);
		subread_destroy_lock(initlocks+thread_i);
	}

	SAMBAM_reset_sorting_writer;
	while(1){ // GO THROUGH EACH READ and write into BAM
		//if(current_min_batch<0)SUBREADprintf("FINALTEST: OLD/NEW = %d , %d\n", old_chro_no ,chro_no);
		if(old_chro_no >=0 && chro_no != old_chro_no){ // if there is old_chro && if has to write. Note: old_chro_no<0 means we reached the "unmapped" part of the sorted temp files
			//SUBREADprintf("\n====== W1CHR %d =====\n", old_chro_no);
			SamBam_write_BAI_for_1chr(writer, &bin_chunks_table, &bins_list, &win16k_list);
			last_written_BAI_chro = old_chro_no;
			SAMBAM_write_empty_chros(old_chro_no +1, chro_no < 0?n_ref:chro_no);
			SAMBAM_renew_BAItabs;
			SAMBAM_reset_sorting_writer;
		}

		if(current_min_batch>-1){
			SamBam_writer_sort_bins_to_BAM_write_1R(writer, sb_fps[current_min_batch], bin_chunks_table, bins_list, win16k_list, chro_no);
			current_min_fps[current_min_batch] = SamBam_writer_sort_bins_to_BAM_FP_pos(sb_fps[current_min_batch]);
			current_min = SUBREAD_MAX_ULONGLONG;
			current_min_batch = -1;
			old_chro_no = chro_no;
		}else{
			SamBam_write_sorted_thread_collect(writer);
			break; // all finished. quit.
		}
		

		for(bii = 0; bii < writer -> sorted_batch_id; bii++){
			if(current_min_fps[bii] < SUBREAD_MAX_ULONGLONG && current_min_fps[bii] < current_min){
				current_min = current_min_fps[bii];
				current_min_batch = bii;
			}
		}
		chro_no = (int)(current_min >> 32);
	}

	for(bii = 0; bii < writer -> sorted_batch_id; bii++){
		if(!sb_fps[bii]) continue;

		char tfp[MAX_FILE_NAME_LENGTH+40];
		sprintf(tfp , "%s-%06d.sortedbin", writer -> tmpf_prefix, bii);
		fclose(sb_fps[bii]);
		unlink(tfp);
	}
	if(bin_chunks_table != NULL) {
		HashTableDestroy(bin_chunks_table);
		ArrayListDestroy(bins_list);
		ArrayListDestroy(win16k_list);
	}

	terminate_workers(&writer->sorted_notifier);
	for(thread_i = 0; thread_i < writer -> threads; thread_i++){
		pthread_join(writer -> writer_threads[thread_i].thread_stub, NULL);
		deflateEnd(&(writer -> writer_threads[thread_i].strm));
	}
	HashTableDestroy(writer -> block_no_p1_to_vpos_tab);
	
	free(writer->writer_threads);
	free(current_min_fps);
	free(sb_fps);
} 

void SamBam_writer_finalise_thread(SamBam_Writer * writer, int thread_id){
	if(writer -> threads < 2){
		if(writer -> sort_reads_by_coord){
			SamBam_writer_sort_buff_one_write(writer, writer ->chunk_buffer, writer -> chunk_buffer_used, -1);
			writer -> chunk_buffer_used= 0;
		}else{
			if(writer -> chunk_buffer_used) SamBam_writer_add_chunk(writer, thread_id);
		}
	}else{
		if(writer -> sort_reads_by_coord){
			SamBam_writer_sort_buff_one_write(writer, writer -> threads_chunk_buffer[thread_id], writer -> threads_chunk_buffer_used[thread_id], thread_id);
			writer -> threads_chunk_buffer_used [thread_id] = 0;
		}else{
			if(writer -> threads_chunk_buffer_used[thread_id]) SamBam_writer_add_chunk(writer, thread_id);
		}
	}
}

void SamBam_writer_finalise_one_thread(SamBam_Writer * writer){
	if(writer -> threads < 2){
		if(writer -> sort_reads_by_coord){
			if(writer -> chunk_buffer_used>0){
				SamBam_writer_sort_buff_one_write(writer, writer -> chunk_buffer, writer -> chunk_buffer_used, -1);
				writer -> chunk_buffer_used = 0;
			}
		}else{
			if(writer -> chunk_buffer_used)
				SamBam_writer_add_chunk(writer, -1);
		}
	}
}
