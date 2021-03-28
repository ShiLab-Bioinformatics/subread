#include <assert.h> 
#include "LRMconfig.h"
#include "LRMseek-zlib.h"

#define SEEKGZ_INIT_TEXT_SIZE (1024*1024)
#define SEEKGZ_BINBUFF_SIZE (1*1024*1024)

unsigned long long LRMseekgz_ftello(seekable_zfile_t * fp){
	unsigned long long ret = ftello(fp -> gz_fp);
	ret -= fp -> stem.avail_in;
	return ret;
}

unsigned int LRMcrc_pos(char * bin, int len){
	unsigned int crc0 = crc32(0, NULL, 0);
	unsigned int CRC32 = crc32(crc0, (unsigned char *) bin, len);
	return CRC32;
}

void LRMseekgz_binreadmore(seekable_zfile_t * fp){
	if(feof(fp->gz_fp))return;

	assert(fp -> stem.avail_in >= 0);
	if(fp -> stem.avail_in < SEEKGZ_BINBUFF_SIZE / 2 ) {
		if(fp -> in_pointer > 0 && fp -> stem.avail_in > 0){
			int i;
			for(i = 0 ; i <  fp -> stem.avail_in ; i ++){
				fp -> current_chunk_bin[i] = fp -> current_chunk_bin[i + fp -> in_pointer];
			}
		}
		fp -> in_pointer = 0;

		int readlen = fread(fp -> current_chunk_bin + fp -> stem.avail_in, 1 , SEEKGZ_BINBUFF_SIZE - fp -> stem.avail_in , fp -> gz_fp);
		if(readlen>0)
			fp -> stem.avail_in += readlen;
		fp -> stem.next_in = (unsigned char *)fp -> current_chunk_bin;
		//SEEKZLIBprintf("READIN: %d,  POS: %llu,  CRC:%u\n", fp -> stem.avail_in , ftello(fp -> gz_fp) , crc_pos(fp -> current_chunk_bin , fp -> stem.avail_in));
	}
}

int LRMseekgz_bingetc(seekable_zfile_t * fp){
	LRMseekgz_binreadmore(fp);	
	int ret = -1;

	if(fp -> stem.avail_in > 0)
	{
		ret = fp -> current_chunk_bin [ fp -> in_pointer ++];
		fp -> stem.next_in = (unsigned char *)(fp -> current_chunk_bin + fp -> in_pointer);
		fp -> stem.avail_in --;
		if(ret<0) ret=256+ret;
	}
	return ret;
	
}

int LRMseekgz_skip_header(seekable_zfile_t * fp, int tail_size){
	int id1, id2;

	if(tail_size){
		for(id1=0; id1<tail_size; id1++)
			LRMseekgz_bingetc(fp);
	}
	id1 = LRMseekgz_bingetc(fp);
	id2 = LRMseekgz_bingetc(fp);

	if(id1 != 31 || id2 != 139){
		//SEEKZLIBprintf("header:%d,%d\n", id1, id2);
		return 1;
	}

	 LRMseekgz_bingetc(fp); // CM
	int FLG= LRMseekgz_bingetc(fp); // FLG 
	LRMseekgz_bingetc(fp);
	LRMseekgz_bingetc(fp);
	LRMseekgz_bingetc(fp);
	LRMseekgz_bingetc(fp);
	LRMseekgz_bingetc(fp); // XFL
	LRMseekgz_bingetc(fp); // OS

	//fprintf(stderr, "FLG=%d, XFL=%d\n" , FLG, XFL);

	if(FLG & (1<<2)){	// FEXT
		unsigned short XLEN=0;
		XLEN = LRMseekgz_bingetc(fp);
		XLEN += LRMseekgz_bingetc(fp)*256;
		for(; XLEN>0; XLEN--)
			LRMseekgz_bingetc(fp);
	}

	for(id1 = 3; id1 <=4; id1++){
		if(FLG & (1<<id1)){	// FNAME or FCOMMENT
			while(1){
				int namec = LRMseekgz_bingetc(fp);
				if(0==namec) break;
			}
		}
	}
	if(FLG & (1<<1)){	// FCRC
		LRMseekgz_bingetc(fp);
		LRMseekgz_bingetc(fp);
	}

	fp -> next_block_file_offset = LRMseekgz_ftello(fp);
	if(fp -> block_start_in_file_offset<1)
		fp -> block_start_in_file_offset = fp -> next_block_file_offset;
	fp -> next_block_file_bits = 0; 
	fp -> dict_window_used = 0;
	fp -> dict_window_pointer = 0;

	fp -> is_the_last_chunk = 2;
	return 0;
}

int LRMseekgz_decompress_next_chunk(seekable_zfile_t * fp);
int LRMseekgz_open(const char * fname, seekable_zfile_t * fp){
	memset(fp, 0, sizeof(seekable_zfile_t));
	fp -> gz_fp = fopen(fname, "rb");
	if(NULL==fp -> gz_fp)return -1;
	fp -> current_chunk_bin = malloc(SEEKGZ_BINBUFF_SIZE);
	fp -> current_chunk_txt = malloc(SEEKGZ_INIT_TEXT_SIZE);
	fp -> current_chunk_txt_size = SEEKGZ_INIT_TEXT_SIZE;
	//fp -> txt_buffer_size = SEEKGZ_INIT_TEXT_SIZE;

	fp -> stem.zalloc = Z_NULL;
	fp -> stem.zfree = Z_NULL;
	fp -> stem.opaque = Z_NULL;
	fp -> stem.avail_in = 0;
	fp -> stem.next_in = Z_NULL;
	
	int ret = LRMseekgz_skip_header(fp,0);
	if(ret) return 1;
	ret = inflateInit2(&(fp -> stem), -15);
	if(ret) return 1;
	return 0;
}

void LRMseekgz_tell(seekable_zfile_t * fp, seekable_position_t * pos){
	pos -> block_gzfile_offset = fp -> block_start_in_file_offset;
	pos -> block_gzfile_bits = fp -> block_start_in_file_bits;
	memcpy(pos -> dict_window, fp -> block_dict_window, fp -> block_dict_window_size); 
	pos -> block_dict_window_size = fp -> block_dict_window_size;
	pos -> in_block_text_offset = fp -> in_block_offset;
}

void LRMseekgz_seek(seekable_zfile_t * fp, seekable_position_t * pos){
	//#warning "COMMENT THIS LINE !!!!."
	//fprintf(stderr, "SEEK => %llu[%d] + %u ; WIN=%d CRC=%u\n", pos -> block_gzfile_offset, pos -> block_gzfile_bits, pos -> in_block_text_offset, pos -> block_dict_window_size, crc_pos( pos -> dict_window, pos -> block_dict_window_size));
	fseeko(fp->gz_fp, pos -> block_gzfile_offset - (pos -> block_gzfile_bits?1:0), SEEK_SET);

	if(Z_OK!=inflateReset(&fp->stem))
		SEEKZLIBprintf("FATAL: UNABLE TO INIT STREAM.\n\n\n");
	if(pos -> block_dict_window_size>0){
		if(pos -> block_gzfile_bits){
			char nch = fgetc(fp->gz_fp);
			//fprintf(stderr, "SEEK 2 FPPOS:%llu, NCH=%d\n", ftello(fp->gz_fp) , nch);
			inflatePrime(&fp->stem, pos -> block_gzfile_bits, nch>>(8-pos -> block_gzfile_bits));
		}
		if(Z_OK != inflateSetDictionary(&fp->stem, (unsigned char *)pos -> dict_window, pos -> block_dict_window_size))
			SEEKZLIBprintf("FATAL: UNABLE TO RESET STREAM.\n\n\n");
	}

	fp -> stem.avail_in = 0;
	fp -> in_pointer = 0;
	fp -> txt_buffer_used = 0;
	fp -> in_chunk_offset = 0;
	memcpy(fp -> block_dict_window, pos -> dict_window, pos -> block_dict_window_size);
	memcpy(fp -> dict_window, pos -> dict_window, pos -> block_dict_window_size);
	fp -> block_dict_window_size = fp -> dict_window_used = pos -> block_dict_window_size;
	fp -> dict_window_pointer = (pos -> block_dict_window_size<SEEKGZ_ZLIB_WINDOW_SIZE)?pos -> block_dict_window_size:0;
	fp -> in_block_offset = 0;
	fp -> block_start_in_file_offset = pos -> block_gzfile_offset;
	fp -> block_start_in_file_bits = pos -> block_gzfile_bits;

	unsigned int chunk_end_block_offset=0;
	while(1){
		LRMseekgz_decompress_next_chunk(fp);
		if(fp -> internal_error) break;
		chunk_end_block_offset += fp -> txt_buffer_used;

		if(chunk_end_block_offset >= pos -> in_block_text_offset){
			fp -> in_chunk_offset = fp -> txt_buffer_used - (chunk_end_block_offset - pos -> in_block_text_offset);
			fp -> in_block_offset = pos -> in_block_text_offset;
			break;
		}
		assert(chunk_end_block_offset < fp -> current_chunk_txt_size && !feof(fp->gz_fp));
		fp -> txt_buffer_used=0;
	}
}



int LRMseekgz_decompress_next_chunk(seekable_zfile_t * fp){
	unsigned int this_chunk_size = 0;
	int loaded_blocks = 0;
	while(1){
		LRMseekgz_binreadmore(fp);
		if(loaded_blocks > 0)
			//SEEKZLIBprintf("LOADED BLOCKS=%d\n", loaded_blocks);
		if(fp -> txt_buffer_used >=  fp -> current_chunk_txt_size * 7 / 8){
			//SEEKZLIBprintf("TRE ALLOC CHUNK_TXT: %d -> %d\n", fp -> current_chunk_txt_size, (int)(fp -> current_chunk_txt_size*1.5));
			fp -> current_chunk_txt_size *= 1.5;
			assert(fp -> current_chunk_txt_size < 3*512*1024*1024);
			fp -> current_chunk_txt = realloc(fp -> current_chunk_txt, fp -> current_chunk_txt_size );
		}

		fp -> stem.avail_out =  fp -> current_chunk_txt_size  - fp -> txt_buffer_used;
		int out_start = fp -> txt_buffer_used;
		fp -> stem.next_out = (unsigned char *)(fp -> current_chunk_txt + out_start);

		int inlen = fp -> stem.avail_in ;
		//fprintf(stderr,"INFLATING_0 : LEN=%u, CRC=%u\n", fp -> stem.avail_in , crc_pos( fp -> stem.next_in  , fp -> stem.avail_in ));

		int ret = inflate(&(fp -> stem), Z_BLOCK);
		int have = ( fp -> current_chunk_txt_size  - fp -> txt_buffer_used) - fp -> stem.avail_out;
		int is_chunk_end = 0;

		//#warning "COMMENT NEXT LINE!!!!!."
		//fprintf(stderr,"INFLATING: INLEN=%d , OLEN=%d, POS=%lld, RET=%d, TOOL=%s\n", inlen , have, LRMseekgz_ftello(fp), ret, zlibVersion());
		if(ret != Z_OK && ret != Z_STREAM_END){ //any error
			#ifdef __MINGW32__
			SEEKZLIBprintf("FATAL: INFLATE-ERROR=%d   POS=%I64d\n", ret, LRMseekgz_ftello(fp));
			#else
			SEEKZLIBprintf("FATAL: INFLATE-ERROR=%d   POS=%lld\n", ret, LRMseekgz_ftello(fp));
			#endif
			fp -> internal_error = 1;
			return -1;
		}

		fp -> in_pointer += inlen -  fp -> stem.avail_in ;

		if(have > 0){
			fp -> txt_buffer_used += have;
			int one_length = 0, one_src_start = 0, one_dst_start = 0;
			int two_length = 0, two_src_start = 0, two_dst_start = 0;
			int new_pntr = 0;
			if(have <= SEEKGZ_ZLIB_WINDOW_SIZE - fp -> dict_window_pointer){
				one_length = 0;
				two_src_start = out_start;
				two_dst_start = fp -> dict_window_pointer;
				two_length = have;
				new_pntr = two_dst_start + two_length;
			}else if(have > SEEKGZ_ZLIB_WINDOW_SIZE - fp -> dict_window_pointer && have <= SEEKGZ_ZLIB_WINDOW_SIZE){
				one_src_start = out_start + SEEKGZ_ZLIB_WINDOW_SIZE - fp -> dict_window_pointer;
				one_dst_start = 0;
				one_length = have - SEEKGZ_ZLIB_WINDOW_SIZE + fp -> dict_window_pointer;
				two_src_start = out_start;
				two_dst_start = fp -> dict_window_pointer;
				two_length = SEEKGZ_ZLIB_WINDOW_SIZE - fp -> dict_window_pointer;
				new_pntr = one_dst_start + one_length;
			}else{
				one_src_start = out_start + have - fp -> dict_window_pointer;
				one_dst_start = 0;
				one_length = fp -> dict_window_pointer;
				two_src_start = out_start + have - SEEKGZ_ZLIB_WINDOW_SIZE;
				two_dst_start = fp -> dict_window_pointer;
				two_length = SEEKGZ_ZLIB_WINDOW_SIZE - fp -> dict_window_pointer;
				new_pntr = fp -> dict_window_pointer;
			}

			if(one_length > 0)memcpy(fp -> dict_window + one_dst_start, fp -> current_chunk_txt + one_src_start, one_length);
			//fprintf(stderr,"CPY: %d -> %d [%d] ; PNTR=%d, NEWPNTR=%d, have=%d\n", two_src_start, two_dst_start, two_length, fp -> dict_window_pointer, new_pntr, have);
			memcpy(fp -> dict_window + two_dst_start, fp -> current_chunk_txt + two_src_start, two_length);
			fp -> dict_window_pointer = new_pntr;
			fp -> dict_window_used = SEEKZLIBmin(fp -> dict_window_used + have, SEEKGZ_ZLIB_WINDOW_SIZE);

			is_chunk_end = (fp -> stem.data_type & 128) && !(fp -> stem.data_type & 64);
			if(is_chunk_end){
				fp -> is_the_last_chunk = 1;
				unsigned long long file_pos_after_avail = LRMseekgz_ftello(fp);
				fp -> next_block_file_offset = file_pos_after_avail;
				fp -> next_block_file_bits = fp->stem.data_type & 7;
			}
			this_chunk_size += have;
		}
		loaded_blocks ++;

		if( 0 == fp -> stem.avail_in ) this_chunk_size = 0;

		if(Z_STREAM_END == ret || ((is_chunk_end || 0 == fp -> stem.avail_in) && fp -> txt_buffer_used >=10)){
			if(Z_STREAM_END == ret){
				LRMseekgz_skip_header(fp, 8); 
				inflateReset(&fp->stem);
			}
			break;
		}
	}
	return 0;
}

int LRMseekgz_next_char(seekable_zfile_t * fp){
	if(fp -> internal_error) return -1;
	while(fp -> in_chunk_offset >= fp -> txt_buffer_used){
		if(feof(fp -> gz_fp) && fp -> stem.avail_in < 10 )
			return EOF;
		else {
			fp -> txt_buffer_used = 0;
			fp -> in_chunk_offset = 0;
			int decompress_ret = LRMseekgz_decompress_next_chunk(fp);
			if(decompress_ret) return -1;
		}
	}
	fp -> in_block_offset ++;
	char retc = fp -> current_chunk_txt[fp -> in_chunk_offset++];

	if(fp -> is_the_last_chunk && fp -> in_chunk_offset == fp -> txt_buffer_used){
		fp -> in_block_offset = 0;
		fp -> block_start_in_file_offset = fp -> next_block_file_offset;
		fp -> block_start_in_file_bits   = fp -> next_block_file_bits;

		if(1 == fp -> is_the_last_chunk){
			fp -> block_dict_window_size  = fp -> dict_window_used;

			if(fp -> dict_window_used < SEEKGZ_ZLIB_WINDOW_SIZE)
				memcpy(fp -> block_dict_window , fp -> dict_window, fp -> dict_window_used);	
			else{
				memcpy(fp -> block_dict_window , fp -> dict_window + fp -> dict_window_pointer, SEEKGZ_ZLIB_WINDOW_SIZE - fp -> dict_window_pointer);
				memcpy(fp -> block_dict_window + SEEKGZ_ZLIB_WINDOW_SIZE - fp -> dict_window_pointer, fp -> dict_window, fp -> dict_window_pointer);
			}
		}else
			fp -> block_dict_window_size  = 0;

		fp -> is_the_last_chunk = 0;
	}

	return retc;
}

int LRMseekgz_gets(seekable_zfile_t * fp, char * buf, int buf_size){
	int i=0;
	buf[0]=0;
	while(1){
		if(i >= buf_size - 1){
			buf[i]=0;
			return i;
		}
		int nch = LRMseekgz_next_char(fp);
		if(nch<0 || nch == '\n'){
			if(i<1 && nch <0) return 0;
			buf[i] = '\n';
			buf[i+1]=0;
			return i+1;
		}else buf[i++]=nch;
	}
}

void LRMseekgz_close(seekable_zfile_t * fp){
	fclose(fp -> gz_fp);
	free(fp -> current_chunk_txt);
	free(fp -> current_chunk_bin);
}
