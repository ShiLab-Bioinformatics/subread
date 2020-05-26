#include <stdio.h>
#include <signal.h>
#include <dirent.h> 
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>

#ifndef __MINGW32__
#include <sys/resource.h>
#endif

#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>
#include <assert.h>

#include "LRMfile-io.h"
#include "LRMsorted-hashtable.h"

int LRMgenekey2int(char key []){
	int i;
	int ret;

	ret = 0;
	for (i=0; i<16; i++)
		ret |= (LRMbase2int(key[i]))<<(2*(15-i));

	return ret;
}

int LRMgeinput_open(const char * filename, LRMgene_input_t * input){
	int ret = 0;
	if(strlen(filename)>LRMMAX_FILENAME_LENGTH-2)
		return 1;

	strcpy(input->filename, filename);
	FILE * TMP_FP = fopen(filename, "rb");

	if(TMP_FP == NULL)	
		return 1;

	int id1, id2;
	id1 = fgetc(TMP_FP);
	id2 = fgetc(TMP_FP);

	if(id1 == 31 && id2 == 139) {
		fclose(TMP_FP);
		input->input_fp = malloc(sizeof(seekable_zfile_t));
		input->file_type = LRMGENE_INPUT_GZIP_FASTQ;
		ret = LRMseekgz_open(filename, input->input_fp);
	}else{
		input->file_type = LRMGENE_INPUT_FASTQ;
		input->input_fp = TMP_FP;
		fseek(input->input_fp, 0, SEEK_SET);
	}

	return ret;
}

void LRMgeinput_close(LRMgene_input_t * input){
		if(input -> file_type == LRMGENE_INPUT_GZIP_FASTQ)
		LRMseekgz_close((seekable_zfile_t * ) input->input_fp);
	else
		fclose((FILE*)input->input_fp);
}

char * LRM__converting_char_table = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNGNNNCNNNNNNNNNNNNAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNGNNNCNNNNNNNNNNNNAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN  ";

void LRMreverse_read(char * InBuff, int read_len){
	int i;

	for (i=0; i<read_len/2; i++)
	{
		int rll1 = read_len - 1 - i;
		unsigned char tmp = InBuff[rll1];

		InBuff[rll1] = *(LRM__converting_char_table+InBuff[i]);
		InBuff[i] = *(LRM__converting_char_table+tmp);

	}
	if(i*2 == read_len-1)
	{
		InBuff[i] = *(LRM__converting_char_table+InBuff[i]);
	}
}



int LRMgeinput_getc(LRMgene_input_t * input){
	if(input -> file_type == LRMGENE_INPUT_GZIP_FASTQ){
		return LRMseekgz_next_char((seekable_zfile_t*)input -> input_fp);
	}else{
		return fgetc((FILE*)input -> input_fp);
	}
}

int LRMgeinput_readline(LRMgene_input_t * input, int buf_len, char * linebuffer){
	int ret =0;
	while(1){
		char ch = LRMgeinput_getc(input);
		if (ch == '\n' || ch == EOF) break;
		if(ret < buf_len-1)
			linebuffer[ret++] = ch;
	}
	linebuffer[ret]=0;
	return ret;
}


#define SKIP_LINE { nch=' '; while(nch != EOF && nch != '\n') nch = LRMgeinput_getc(input); }

int LRMgeinput_next_read(LRMgene_input_t * input, char * read_name, char * read_string, char * quality_string){
	char nch;
	int ret;

	//READ NAME
	if (read_name == NULL){
		SKIP_LINE;
		if(nch == EOF) return -1;
	} else {
		int ATnch = LRMgeinput_getc(input);
		assert(ATnch == '@' || ATnch <0);
		if(ATnch < 0) return -1;

		int retv = LRMgeinput_readline(input, LRMMAX_READ_NAME_LEN, read_name);
		if(retv<=0) return -1;

		int cursor = 1;
		while(read_name[cursor])
		{
			if(read_name[cursor] == ' ' || read_name[cursor] == '\t')
			{
				read_name [cursor] = 0;
				break;	
			}
			cursor++;
		}
	}
	// READ LINE 
	ret = LRMgeinput_readline(input, LRMMAX_READ_LENGTH, read_string);

	// SKIP "+"
	do{
		nch = LRMgeinput_getc(input);
	} while( nch == '\n' );
	SKIP_LINE;

	// QUAL LINE 
	if (quality_string)
		LRMgeinput_readline(input, LRMMAX_READ_LENGTH, quality_string);
	else
		SKIP_LINE;

	return ret;
}

int LRMfetch_next_read(LRMcontext_t * context, LRMthread_context_t * thread_context, unsigned int *read_len, char * read_name, char * read_text, char * qual_text, unsigned int * read_no_in_chunk){
	int this_number = -1;
	int this_rlen = 0;

	LRMthread_lock(&context -> input_lock); 
	this_rlen = LRMgeinput_next_read(&context -> input_file, read_name, read_text, qual_text);
	if(this_rlen > 0 ){
		this_number = context -> processed_reads_in_chunk;
		context -> processed_reads_in_chunk ++;
	}else context -> input_exhausted = 1;

	LRMthread_lockrelease(&context -> input_lock); 

	if(this_rlen && this_number>=0)
	{
		*read_no_in_chunk = this_number;
		*read_len = this_rlen;
		return 0;
	}else{
		*read_no_in_chunk = -1;
		return 1;
	}
}

void LRMreverse_quality(char * InBuff, int read_len){
	int i;
	if(!InBuff) return;
	if(!InBuff[0]) return;
	for (i=0; i<read_len/2; i++)
	{
		char tmp;
		tmp = InBuff[i];
		InBuff[i] = InBuff[read_len -1-i];
		InBuff[read_len -1-i] = tmp;
	}
}

void LRMquality_64_to_33(char *qs){
	int i;
	for(i=0; qs[i]; i++){
		qs[i] -= (64-33);
	}
}

int LRMgenerate_bam_record_encode_cigar(LRMcontext_t * context, int * cigar_int, char * cigar, int * mapped_length, int read_text_len ){
	int tmp_int=0;
	int reported_read_len = 0;
	int cigar_cursor = 0, num_opt = 0;
	(*mapped_length) = 0;

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
			if(nch == 'M' ||nch == 'N'||nch == 'D') (*mapped_length) += tmp_int;
			if(nch == 'M' ||nch == 'I'||nch == 'S') reported_read_len += tmp_int;
			for(; int_opt<8; int_opt++) if("MIDNSHP=X"[int_opt] == nch)break;
			cigar_int[num_opt ++] = (tmp_int << 4) | int_opt; 
			tmp_int = 0;
			if(num_opt >= context -> max_cigars_in_read){
				int_opt = 4; //'S'
				tmp_int = read_text_len - reported_read_len;
				cigar_int[num_opt ++] = (tmp_int << 4) | int_opt;
				LRMprintf("CIGAR_TRIMMED : %d bases\n", tmp_int);
				break;
			}
		}
	}

	return num_opt;
}

int LRMreg2bin(int beg, int end){
	--end;
	if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
	if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
	if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
	if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
	if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
	return 0;
}

int LRMgenerate_bam_record_encode_read_qual(char * bin, char * read, char * qual, int rlen){
	int w_ptr = 0, xk1;
		
	for(xk1 = 0; xk1 < rlen; xk1++){
		int fourbit;
		for(fourbit=0;fourbit<15;fourbit++) if("=ACMGRSVTWYHKDBN"[fourbit] == read[xk1])break;
		if( xk1 % 2 == 0){ 
			fourbit = fourbit << 4;
			bin[w_ptr]=0;
		}

		bin[w_ptr] |= fourbit;
		
		if(xk1 % 2 == 1)w_ptr++;
	}
	if(rlen %2) w_ptr++;
	
	for(xk1=0; xk1<rlen; xk1++)
		bin[w_ptr+xk1] = qual[xk1]-33;
	return w_ptr+rlen;
}

int LRMgenerate_bam_record(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, char * target_ptr, int flags, unsigned int chro_pos,char *  chro_name, int map_quality, char * cigar, int mis_matched){
	int mapped_length = 0;
	
	int chro_number = LRMHashTableGet(context->sam_bam_chromosome_table, chro_name) - NULL;
	chro_number -=1;
	memcpy(target_ptr+4, &chro_number, 4);
	memcpy(target_ptr+8, &chro_pos, 4);
	
	int name_len = strlen(iteration_context -> read_name)+1;
	
	memcpy(target_ptr+20, &iteration_context->read_length, 4);
	int next_non = -1;
	memcpy(target_ptr+24, &next_non, 4);
	memcpy(target_ptr+28, &next_non, 4);
	memset(target_ptr+32, 0, 4);	//TLEN
	memcpy(target_ptr+36, iteration_context->read_name, name_len);
	int bin_ptr = 36+name_len;
	
	int cigar_opts = LRMgenerate_bam_record_encode_cigar(context, (int *)(target_ptr + bin_ptr), cigar, &mapped_length, iteration_context->read_length);
	int flag_nc = flags << 16 | cigar_opts;
	memcpy(target_ptr+16, &flag_nc, 4);
	bin_ptr += cigar_opts*4;
	
	
	int bin_mq_nl = (LRMreg2bin(chro_pos, chro_pos+mapped_length)<<16) | map_quality << 8 | name_len;
	memcpy(target_ptr+12, &bin_mq_nl, 4);
	
	//LRMprintf("READ TEXT %s : len=%d \n", iteration_context->read_name, iteration_context->read_length);
	bin_ptr += LRMgenerate_bam_record_encode_read_qual(target_ptr + bin_ptr, iteration_context->read_text, iteration_context -> qual_text, iteration_context->read_length);
	
	memcpy(target_ptr+bin_ptr, "NM",2);
	target_ptr[bin_ptr+2]='C';
	target_ptr[bin_ptr+3]=mis_matched;
	
	memcpy(target_ptr, &bin_ptr, 4); // block len not include itself.
	bin_ptr += 4;
	return bin_ptr;
}

unsigned int LRM_CRC32(char * dat, int len){
	unsigned int crc0 = crc32(0, NULL, 0);
	unsigned int ret = crc32(crc0, (unsigned char *)dat, len);
	return ret;
}

int LRMwrite_chunk_compress_bam_block(LRMcontext_t * context,  LRMthread_context_t * thread_context, char * bin_buf, char * bam_buf, int bin_len){
	char * compressed_buff = bam_buf + 18;
	
	int compressed_size ; 
	unsigned int CRC32;
	thread_context -> bam_file_output_stream.avail_out = 66600;
	thread_context -> bam_file_output_stream.avail_in = bin_len;
	CRC32 = LRM_CRC32(bin_buf , bin_len);
	
 	int Z_DEFAULT_MEM_LEVEL = 8;
	thread_context -> bam_file_output_stream.zalloc = Z_NULL;
	thread_context -> bam_file_output_stream.zfree = Z_NULL;
	thread_context -> bam_file_output_stream.opaque = Z_NULL;

	deflateInit2(&thread_context -> bam_file_output_stream, LRMBAM_COMPRESS_LEVEL, Z_DEFLATED,
		LRMGZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
	
	thread_context -> bam_file_output_stream.next_in = (unsigned char*) bin_buf;
	thread_context -> bam_file_output_stream.next_out = (unsigned char*) compressed_buff;

	deflate(&thread_context -> bam_file_output_stream, Z_FINISH);
	deflateEnd(&thread_context -> bam_file_output_stream);

	compressed_size = 66600 -thread_context -> bam_file_output_stream.avail_out;
	
	//LRMprintf("COMPRESS: %d -> %d\n", bin_len, compressed_size);
	
	bam_buf[0]=31;
	bam_buf[1]=(char)139;
	bam_buf[2]=8;
	bam_buf[3]=4;
	memset(bam_buf+4, 0, 5);
	bam_buf[9] = 0xff;	// OS
	
	int tmpi = 6;
	memcpy(bam_buf+10, &tmpi, 2); //XLSN
	bam_buf[12]=66; // SI1
	bam_buf[13]=67; // SI2
	tmpi = 2;
	memcpy(bam_buf+14, &tmpi, 2); //BSIZE
	tmpi = compressed_size + 19 + 6;
	memcpy(bam_buf+16, &tmpi, 2); //BSIZE
	
	memcpy(bam_buf+18+compressed_size, &CRC32, 4);
	memcpy(bam_buf+18+compressed_size+4, &bin_len, 4);
	return compressed_size + 26;
}

int LRMhash_strcmp(const void * s1, const void * s2){
	return strcmp(s1, s2);
}

srUInt_64 LRMhash_strhash(const void * sv){
	unsigned char *s = (unsigned char *)sv;
	unsigned long ret = 0;
	while(*s){
		ret ^= (ret << 6);
		ret ^= *s;
		s++;
	}
	return ret;
}

int LRMreadline(FILE * fp, char * buf, int buflen){
	int cr = 0;
	while(!feof(fp)){
		char c = fgetc(fp);
		if(c == '\r'){
			continue;
		}else if(c == '\n'){
			break;
		}else if(cr < buflen-1)
			buf[cr++]=c;
	}
	buf[cr]=0;
	return cr;
}

int LRMload_offsets(LRMcontext_t * context){
	char fn[LRMMAX_FILENAME_LENGTH+20];
	FILE * fp;
	int padding = 0;
	unsigned int last_end = 0;

	LRMgehash_load_option(context -> index_prefix, LRMSUBREAD_INDEX_OPTION_INDEX_PADDING , &padding);
	assert(padding>16);
	
	sprintf(fn, "%s.reads", context->index_prefix);
	fp = fopen(fn, "r");

	if(!fp){
		LRMprintf("file not found :%s\n", fn);
		return 1;
	}

	context -> current_index_padding = padding;
	int chromosome_no = 0;
	while (!feof(fp))
	{
		int i=0, step = 0, j=0;
		unsigned int this_offset =0;
		
		LRMreadline(fp, fn, LRMMAX_FILENAME_LENGTH-1);
		if (strlen(fn)<2)continue;
		char * chro_name = malloc(LRMMAX_CHROMOSOME_NAME_LEN);
		while (fn[i])
		{
			if (fn[i] == '\t') {
				fn[i]=0;
				this_offset = (unsigned int)atoll(fn);
				step = 1;
			}else if (step)	{
				if(j<LRMMAX_CHROMOSOME_NAME_LEN-1){
					chro_name[j++]=fn[i];
					chro_name[j]=0;
				}
			}
			i++;
		}

		LRMHashTablePut(context->sam_bam_chromosome_table, chro_name, NULL + chromosome_no + 1);
		LRMArrayListPush(context->sam_bam_chromosome_list, chro_name);

		LRMHashTablePut(context->chromosome_size_table, chro_name, NULL + (this_offset - last_end + 16 - context->current_index_padding * 2));
		LRMArrayListPush(context->chromosome_size_list, NULL + this_offset);

		chromosome_no++;
		last_end = this_offset;
	}

	fclose(fp);
	return 0;
}

#define LRMcheck_resize_buff 			if(new_need >= thread_context -> out_buff_capacity){\
				thread_context -> out_buff_capacity = max(thread_context -> out_buff_capacity*2, new_need );\
				thread_context -> out_SAMBAM_buffer=realloc(thread_context -> out_SAMBAM_buffer, thread_context -> out_buff_capacity);\
			}\

#define LRMBAM_COMPRESS_BLOCK 63000
#define LRMBAM_COMPRESS_TRIGGER 53000

void LRMsambam_write_header(LRMcontext_t * context, LRMthread_context_t * thread_context){
	if(context -> sam_bam_file_header_written) return;
	
	thread_context -> out_SAMBAM_buffer = malloc(3000000);
	thread_context -> out_buff_capacity = 3000000;
	thread_context -> out_buff_used = 0;

	int chro_no;
	if(!context -> is_SAM_output){
		memcpy(thread_context -> out_SAMBAM_buffer,"BAM\1",4);
		thread_context -> out_buff_used = 8;
	}
	
	for(chro_no = -1; chro_no < context -> sam_bam_chromosome_list->numOfElements + 2; chro_no++){
		char * header_line = malloc(10100);	
		int wrlen = 0;
		if(chro_no >=0 && chro_no < context -> sam_bam_chromosome_list->numOfElements){
			char * chro_name = LRMArrayListGet(context -> sam_bam_chromosome_list, chro_no);
			int chro_length = LRMHashTableGet(context -> chromosome_size_table, chro_name) - NULL;
						
			wrlen = sprintf(header_line, "@SQ\tSN:%s\tLN:%d\n",chro_name,chro_length);
		}else if(chro_no == -1){
			wrlen = sprintf(header_line, "@HD\tVN:1.0\tSO:unsorted\n");
		}else if(chro_no == context -> sam_bam_chromosome_list->numOfElements ){
			wrlen = sprintf(header_line, "@PG\tID:subread-long-read-mapping\tPN:subread-long-read-mapping\tCL:%s\n", context -> user_command_line);
		}
		if(context -> is_SAM_output){
			fwrite( header_line, 1, wrlen, context -> sam_bam_file);
		}else{
			int new_need = thread_context -> out_buff_used + wrlen +1;
			LRMcheck_resize_buff;
			memcpy(thread_context -> out_SAMBAM_buffer + thread_context -> out_buff_used,header_line,wrlen);
			thread_context -> out_buff_used += wrlen;
		}
		free(header_line);
	}
	
	int new_need = thread_context -> out_buff_used + 10;
	LRMcheck_resize_buff;

	int BAM_text_len = thread_context -> out_buff_used-8;
	memcpy(thread_context -> out_SAMBAM_buffer+4, &BAM_text_len, 4);
	memcpy(thread_context -> out_SAMBAM_buffer + thread_context -> out_buff_used, &context ->sam_bam_chromosome_list->numOfElements , 4);
	thread_context -> out_buff_used +=4;
	for(chro_no = 0; chro_no < context -> sam_bam_chromosome_list->numOfElements ; chro_no++){
		char * chro_name = LRMArrayListGet(context -> sam_bam_chromosome_list, chro_no);
		int chro_namelen = strlen(chro_name)+1;
		int new_need = thread_context -> out_buff_used + chro_namelen + 9;
		LRMcheck_resize_buff;

		memcpy(thread_context -> out_SAMBAM_buffer+thread_context -> out_buff_used, &chro_namelen, 4);
		thread_context -> out_buff_used +=4;
		memcpy(thread_context -> out_SAMBAM_buffer+thread_context -> out_buff_used, chro_name, chro_namelen);
		thread_context -> out_buff_used += chro_namelen;
		int chro_length = LRMHashTableGet(context -> chromosome_size_table, chro_name) - NULL;
		memcpy(thread_context -> out_SAMBAM_buffer+thread_context -> out_buff_used, &chro_length, 4);
		thread_context -> out_buff_used +=4;
	}
	
	if(!context -> is_SAM_output){
		LRMwrite_chunk_check_buffer_write(context, thread_context, 1);
		
	}
	context -> sam_bam_file_header_written = 1;
	free(thread_context -> out_SAMBAM_buffer);
}

void LRMbam_generate_tail_binary(LRMcontext_t * context, LRMthread_context_t * thread_context){
	if(context -> bam_file_tail_length>0)return;
	context -> bam_file_tail_length = LRMwrite_chunk_compress_bam_block(context, thread_context, "", context -> bam_file_tail_binary , 0);
}

int LRMwrite_chunk_add_buffered_output(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, int flags,  char * chro_name, unsigned int chro_pos, int map_quality, char * cigar, int mis_matched){
	char * target_ptr;
	int cigar_len = strlen(cigar);
	int rname_len = strlen(iteration_context->read_name);
	int read_binlen, actural_target_len;
	read_binlen = cigar_len+rname_len+2.5*iteration_context->read_length+500;
	
	if(read_binlen + thread_context->out_buff_used >= thread_context -> out_buff_capacity){
		thread_context -> out_buff_capacity = max( thread_context -> out_buff_capacity*1.3, read_binlen);
		thread_context -> out_SAMBAM_buffer=realloc(thread_context -> out_SAMBAM_buffer, thread_context -> out_buff_capacity);
	}
	
	target_ptr = thread_context -> out_SAMBAM_buffer +  thread_context->out_buff_used;
	if(context -> is_Phred_64) LRMquality_64_to_33(iteration_context->qual_text) ;
	
	if(context->is_SAM_output){
		actural_target_len = sprintf(target_ptr,"%s\t%d\t%s\t%u\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:%d\n", iteration_context -> read_name, flags, chro_name, chro_pos + 1, map_quality, cigar, iteration_context->read_text, iteration_context->qual_text, mis_matched);
	}else{
		actural_target_len = LRMgenerate_bam_record(context, thread_context, iteration_context, target_ptr, flags, chro_pos, chro_name, map_quality,  cigar, mis_matched);
	}
	thread_context -> out_buff_used+=actural_target_len;
	LRMwrite_chunk_check_buffer_write(context, thread_context, 0);
	return 0;
}
int LRMwrite_chunk_check_buffer_write(LRMcontext_t * context, LRMthread_context_t * thread_context, int force_write){
	if(force_write || thread_context -> out_buff_used > LRMBAM_COMPRESS_TRIGGER){
		if(!context->is_SAM_output){
			int write_cursor = 0;
			int compressed_cursor = 0;
			for(write_cursor = 0; write_cursor < thread_context -> out_buff_used; write_cursor += LRMBAM_COMPRESS_BLOCK){
				char compressed_data [66666];
				int bin_len = min(LRMBAM_COMPRESS_BLOCK,thread_context -> out_buff_used - write_cursor);
				int compressed_len = LRMwrite_chunk_compress_bam_block(context, thread_context, thread_context -> out_SAMBAM_buffer + write_cursor, compressed_data, bin_len);
				if(compressed_len > bin_len){	// very unlikely
					int new_needed_size = thread_context -> out_buff_used  + (compressed_len - bin_len);
					if(new_needed_size > thread_context -> out_buff_capacity){
						thread_context -> out_buff_capacity = new_needed_size;
						thread_context -> out_SAMBAM_buffer=realloc(thread_context -> out_SAMBAM_buffer, thread_context -> out_buff_capacity);
					}
					int x1;
					for(x1 = thread_context -> out_buff_used - 1; x1 >= write_cursor + bin_len; x1 --)
						thread_context -> out_SAMBAM_buffer[ x1 + (compressed_len - bin_len) ] = thread_context -> out_SAMBAM_buffer[ x1 ];
				}
				memcpy(thread_context -> out_SAMBAM_buffer + compressed_cursor, compressed_data, compressed_len);
				compressed_cursor += compressed_len;	
			}
			thread_context -> out_buff_used = compressed_cursor;
		}
		LRMthread_lock(&context->sam_bam_file_lock);
		fwrite(thread_context -> out_SAMBAM_buffer,1, thread_context -> out_buff_used, context -> sam_bam_file);
		LRMthread_lockrelease(&context->sam_bam_file_lock);
		thread_context -> out_buff_used = 0;
	}
	
	return 0;
}



void LRMpos2txt(LRMcontext_t * context, unsigned int linear, char * txt){
	char * chro_name = NULL;
	int pos;

	int retv = LRMlocate_gene_position(context, linear, &chro_name, &pos);
	if(chro_name && !retv)
		sprintf(txt, "%s:%d", chro_name, pos+1);
	else	strcpy(txt, "*");
}

int LRMlocate_gene_position(LRMcontext_t * context, unsigned int linear, char ** chro_name, int * pos) {
	int n = 0;
	int total_offsets = context -> chromosome_size_list -> numOfElements;
	int jump_ns = total_offsets/4;

	//LRMprintf("LINEAR = %u\n", linear);

	(*chro_name)=NULL;
	(*pos) = -1;
	while (jump_ns > 5)
	{
		while(n+jump_ns < total_offsets && LRMArrayListGet(context->chromosome_size_list, n + jump_ns) - NULL <= linear)
			n+=jump_ns;
		jump_ns /=4;
	}

	for (; n < total_offsets; n++)
	{
		//LRMprintf("TESTING %u\n",(LRMArrayListGet(context->chromosome_size_list, n) - NULL) );
		if ( (LRMArrayListGet(context->chromosome_size_list, n) - NULL)  > linear)
		{
			*pos = linear;
			if(n>0)(*pos) -= (LRMArrayListGet(context->chromosome_size_list, n-1)-NULL);
			if( (*pos) < context -> current_index_padding )	return 1;
			else (*pos) -= context -> current_index_padding;
			(*chro_name) = LRMArrayListGet(context -> sam_bam_chromosome_list, n);

			return 0;
		}
	}
	return -1;
}

int LRMlocate_chro_length(LRMcontext_t * context, unsigned int linear, char ** chro_name, long long * end_pos){
	int n = 0;
	int total_offsets = context -> chromosome_size_list -> numOfElements;
	int jump_ns = total_offsets/4;
	
	while (jump_ns > 5)
	{
		while(n+jump_ns < total_offsets && LRMArrayListGet(context->chromosome_size_list, n + jump_ns) - NULL <= linear)
			n+=jump_ns;
		jump_ns /=4;
	}

	for (; n < total_offsets; n++)
	{
		//LRMprintf("TESTING %u\n",(LRMArrayListGet(context->chromosome_size_list, n) - NULL) );
		if ( (LRMArrayListGet(context->chromosome_size_list, n) - NULL)  > linear) {
			(* chro_name) = LRMArrayListGet(context -> sam_bam_chromosome_list, n);
			(* end_pos) = LRMArrayListGet(context->chromosome_size_list, n)-NULL;
			if((* end_pos) > context -> current_index_padding) (* end_pos) -=  context -> current_index_padding ;
			return 0;
		}
	}
	return -1;
}
