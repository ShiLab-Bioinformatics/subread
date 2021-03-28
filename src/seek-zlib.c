#include <assert.h> 
#include "core.h"
#include "seek-zlib.h"
#include "input-files.h"
#include "gene-algorithms.h"

#define SEEKGZ_INIT_TEXT_SIZE (1024*1024)
#define SEEKGZ_BINBUFF_SIZE (1*1024*1024)

unsigned long long seekgz_ftello(seekable_zfile_t * fp){
	unsigned long long ret = ftello(fp -> gz_fp);
	ret -= fp -> stem.avail_in;
	return ret;
}

unsigned int crc_pos(char * bin, int len){
	unsigned int crc0 = 0;//crc32(0, NULL, 0);
	unsigned int CRC32 = crc32(crc0, (unsigned char *) bin, len);
	return CRC32;
}

void seekgz_try_read_some_zipped_data(seekable_zfile_t * fp, int * is_eof){
	if(feof(fp->gz_fp)){
		if(is_eof) (*is_eof)=1;
		return;
	}

	assert(fp -> stem.avail_in >= 0);
	if(fp -> stem.avail_in < SEEKGZ_BINBUFF_SIZE / 2 ) {
		if(fp -> in_zipped_buff_read_ptr > 0 && fp -> stem.avail_in > 0){
			int i;
			for(i = 0 ; i < fp -> stem.avail_in ; i ++){
				fp -> in_zipped_buffer[i] = fp -> in_zipped_buffer[i + fp -> in_zipped_buff_read_ptr];
			}
		}
		fp -> in_zipped_buff_read_ptr = 0;

		int readlen = fread(fp -> in_zipped_buffer + fp -> stem.avail_in, 1 , SEEKGZ_BINBUFF_SIZE - fp -> stem.avail_in , fp -> gz_fp);
		if(readlen>0)
			fp -> stem.avail_in += readlen;
		fp -> stem.next_in = (unsigned char *)fp -> in_zipped_buffer;
	}
}

int seekgz_get_next_zipped_char(seekable_zfile_t * fp){
	seekgz_try_read_some_zipped_data(fp, NULL);	
	int ret = -1;

	if(fp -> stem.avail_in > 0)
	{
		ret = fp -> in_zipped_buffer [ fp -> in_zipped_buff_read_ptr ++];
		fp -> stem.next_in = (unsigned char *)(fp -> in_zipped_buffer + fp -> in_zipped_buff_read_ptr);
		fp -> stem.avail_in --;
		if(ret<0) ret=256+ret;
	}
	return ret;
	
}

int seekgz_skip_gzfile_header(seekable_zfile_t * fp, int tail_size){
	int id1, id2;

	if(tail_size){
		for(id1=0; id1<tail_size; id1++)
			seekgz_get_next_zipped_char(fp);
	}
	id1 = seekgz_get_next_zipped_char(fp);
	id2 = seekgz_get_next_zipped_char(fp);

	if(id1 != 31 || id2 != 139){
		//SUBREADprintf("WRONG header:%d,%d\n", id1, id2);
		return 1;
	}

	 seekgz_get_next_zipped_char(fp); // CM
	int FLG= seekgz_get_next_zipped_char(fp); // FLG 
	seekgz_get_next_zipped_char(fp);
	seekgz_get_next_zipped_char(fp);
	seekgz_get_next_zipped_char(fp);
	seekgz_get_next_zipped_char(fp);
	seekgz_get_next_zipped_char(fp); // XFL
	seekgz_get_next_zipped_char(fp); // OS

	//fprintf(stderr, "FLG=%d, XFL=%d\n" , FLG, XFL);

	if(FLG & (1<<2)){	// FEXT
		unsigned short XLEN=0;
		XLEN = seekgz_get_next_zipped_char(fp);
		XLEN += seekgz_get_next_zipped_char(fp)*256;
		for(; XLEN>0; XLEN--)
			seekgz_get_next_zipped_char(fp);
	}

	for(id1 = 3; id1 <=4; id1++){
		if(FLG & (1<<id1)){	// FNAME or FCOMMENT
			while(1){
				int namec = seekgz_get_next_zipped_char(fp);
				if(0==namec) break;
			}
		}
	}
	if(FLG & (1<<1)){	// FCRC
		seekgz_get_next_zipped_char(fp);
		seekgz_get_next_zipped_char(fp);
	}

	fp -> next_block_file_offset = seekgz_ftello(fp);
	fp -> next_block_file_bits = 0; 
	fp -> rolling_dict_window_used = 0;
	return 0;
}
int seekgz_load_more_blocks( seekable_zfile_t * fp , int bytes_ahead, subread_lock_t * read_lock );

int seekgz_open(const char * fname, seekable_zfile_t * fp, FILE * old_fp){
	memset(fp, 0, sizeof(seekable_zfile_t));
	if(old_fp) fp -> gz_fp = old_fp;
	else fp -> gz_fp = f_subr_open(fname, "rb");
	if(NULL==fp -> gz_fp)return -1;

	fp -> in_zipped_buffer = malloc(SEEKGZ_BINBUFF_SIZE);
	subread_init_lock(&fp -> write_lock);
	//fp -> txt_buffer_size = SEEKGZ_INIT_TEXT_SIZE;

	fp -> blocks_in_chain = 0;
	fp -> block_chain_current_no = 0; // no blocks in the chain
	fp -> rolling_dict_window_used = 0;
	fp -> current_block_txt_read_ptr = 0;

	fp -> stem.zalloc = Z_NULL;
	fp -> stem.zfree = Z_NULL;
	fp -> stem.opaque = Z_NULL;
	fp -> stem.avail_in = 0;
	fp -> stem.next_in = Z_NULL;
	if( old_fp ){
		fp -> stem.avail_in = 2;
		fp -> in_zipped_buffer[0] = 0x1f;
		fp -> in_zipped_buffer[1] = 0x8b;
	}
	
	int ret = seekgz_skip_gzfile_header(fp,0);
	if(ret) return 1;
	ret = inflateInit2(&(fp -> stem), -15);
	if(ret) return 1;

	fp -> next_block_file_offset = seekgz_ftello( fp );
	fp -> next_block_file_bits = 0;
	ret = seekgz_load_more_blocks(fp,  300000, NULL);
	//SUBREADprintf("Gzip : opened %s of %p\n", fname, fp);
	return ret;
}

void seekgz_tell(seekable_zfile_t * fp, seekable_position_t * pos){
	if(fp -> blocks_in_chain>0){
		seekable_decompressed_block_t * curr_block = fp -> block_rolling_chain+fp -> block_chain_current_no;

		pos -> block_gzfile_offset = curr_block -> block_start_in_file_offset;
		pos -> block_gzfile_bits = curr_block -> block_start_in_file_bits;
		memcpy(pos -> dict_window, curr_block -> block_dict_window, curr_block -> block_dict_window_size); 
		pos -> block_dict_window_size = curr_block -> block_dict_window_size;
		pos -> in_block_text_offset = fp -> current_block_txt_read_ptr;
	}else{
		pos -> block_gzfile_offset = 0xffffffffffffffffllu;
	}
}

int seekgz_load_1_block( seekable_zfile_t * fp , int empty_block_no);
void seekgz_seek(seekable_zfile_t * fp, seekable_position_t * pos){
	//SUBREADprintf("SEEK:gzip -> fpos=%llu, block_offset=%d\n", pos -> block_gzfile_offset, pos -> in_block_text_offset);
	//
	if(pos -> block_gzfile_offset== 0xffffffffffffffffllu){
		fp -> block_chain_current_no = 0;
		fp -> blocks_in_chain = 0;
		fp -> stem.avail_in = 0;
		fseeko(fp->gz_fp, 0, SEEK_END);
		fgetc(fp->gz_fp);
		return ;
	}
	fseeko(fp->gz_fp, pos -> block_gzfile_offset - (pos -> block_gzfile_bits?1:0), SEEK_SET);

	if(Z_OK!=inflateReset(&fp->stem))
		SUBREADprintf("FATAL: UNABLE TO INIT STREAM.\n\n\n");
	if(pos -> block_dict_window_size>0){
		if(pos -> block_gzfile_bits){
			char nch = fgetc(fp->gz_fp);
			//fprintf(stderr, "SEEK 2 FPPOS:%llu, NCH=%d\n", ftello(fp->gz_fp) , nch);
			inflatePrime(&fp->stem, pos -> block_gzfile_bits, nch>>(8-pos -> block_gzfile_bits));
		}
		if(Z_OK != inflateSetDictionary(&fp->stem, (unsigned char *)pos -> dict_window, pos -> block_dict_window_size))
			SUBREADprintf("FATAL: UNABLE TO RESET STREAM.\n\n\n");
	}

	fp -> stem.avail_in = 0;

	int ii;
	for(ii = 0; ii < fp -> blocks_in_chain ; ii++){
		int iv = ii + fp -> block_chain_current_no;
		if(iv >=SEEKGZ_CHAIN_BLOCKS_NO) iv-=SEEKGZ_CHAIN_BLOCKS_NO;
		free(fp -> block_rolling_chain[iv].block_txt);
		free(fp -> block_rolling_chain[iv].linebreak_positions);
	}

	fp -> block_chain_current_no = 0;
	fp -> blocks_in_chain = 0;
	fp -> in_zipped_buff_read_ptr = 0;
	fp -> current_block_txt_read_ptr = 0;

	memcpy(fp -> rolling_dict_window, pos -> dict_window, pos -> block_dict_window_size);
	fp -> rolling_dict_window_used= pos -> block_dict_window_size;
	fp -> next_block_file_offset = pos -> block_gzfile_offset;
	fp -> next_block_file_bits = pos -> block_gzfile_bits;
	seekgz_load_more_blocks(fp, 30000, NULL);
	fp -> current_block_txt_read_ptr = pos -> in_block_text_offset; 

	//SUBREADprintf("SEK: %d WIN=%d BLK (%d / %d), %d AVI\n", fp -> blocks_in_chain,  pos -> block_dict_window_size, fp -> current_block_txt_read_ptr , fp -> block_rolling_chain[0].block_txt_size, fp -> stem.avail_in);
}


void seekgz_find_linebreaks(seekable_zfile_t * fp, int empty_block_no){
	int lbks= 0, arrsize = 5000;
	unsigned int * arr = malloc(sizeof(int) * arrsize);

	int ii;
	for(ii = 0; ii < fp -> block_rolling_chain[empty_block_no].block_txt_size ; ii++){
		char nch = fp-> block_rolling_chain[empty_block_no].block_txt[ii];
		if(nch == '\n'){
			if(lbks >= arrsize){
				arrsize*=2;
				arr = realloc(arr, sizeof(int) * arrsize);
			}

			arr[lbks++] = ii;
		}
	}
	fp -> block_rolling_chain[empty_block_no].linebreak_positions = arr;
	fp -> block_rolling_chain[empty_block_no].linebreaks = lbks;
}

void seekgz_update_current_window(seekable_zfile_t * fp , char * txt, int tlen){
	int fp_start=0, new_start=0, cplen=0;

	int moveback = 0;
	if(tlen + fp -> rolling_dict_window_used > SEEKGZ_ZLIB_WINDOW_SIZE)
		moveback = tlen + fp -> rolling_dict_window_used - SEEKGZ_ZLIB_WINDOW_SIZE;
	if(moveback>0 && moveback < fp -> rolling_dict_window_used){
		int ii;
		for(ii = 0; ii <  fp -> rolling_dict_window_used - moveback; ii++){
			fp -> rolling_dict_window[ii] = fp -> rolling_dict_window[ii + moveback];
		}
	}

	fp -> rolling_dict_window_used = moveback < fp -> rolling_dict_window_used ? (fp -> rolling_dict_window_used - moveback):0; 

	if(tlen > SEEKGZ_ZLIB_WINDOW_SIZE){
		fp_start =0;
		new_start = tlen - SEEKGZ_ZLIB_WINDOW_SIZE;
		cplen = SEEKGZ_ZLIB_WINDOW_SIZE;
	}else{
		new_start = 0;
		fp_start = fp -> rolling_dict_window_used;
		cplen = tlen;
	}

	memcpy(fp -> rolling_dict_window + fp_start, txt + new_start, cplen);
	fp -> rolling_dict_window_used += cplen;
}

int seekgz_eof(seekable_zfile_t * fp){
	if(fp -> stem.avail_in < 1 && feof(fp->gz_fp))return 1;
	return 0;
}


#define SEEKGZ_ONE_BLOCK_TEXT_INIT_SIZE (256*1024)
// 1block ret < 0 : bad error
// 1block ret > 0 : data is loaded
// 1block ret == 0 : no data is loaded

// A block MUST have at least 10000 bytes otherwise multiple blocks are concatinated

#define MINIMUM_GZIP_BLOCK_SIZE 10000

int seekgz_load_1_block( seekable_zfile_t * fp , int empty_block_no){
	if(seekgz_eof(fp))return 0;

	char * out_txt = malloc(SEEKGZ_ONE_BLOCK_TEXT_INIT_SIZE);
	int out_txt_size = SEEKGZ_ONE_BLOCK_TEXT_INIT_SIZE;
	int out_txt_used = 0;

	memcpy(fp -> block_rolling_chain[empty_block_no].block_dict_window, fp -> rolling_dict_window, fp -> rolling_dict_window_used);
	fp -> block_rolling_chain[empty_block_no].block_dict_window_size = fp -> rolling_dict_window_used;
	//fp -> rolling_dict_window_used = 0;

	fp -> block_rolling_chain[empty_block_no].block_start_in_file_offset = fp -> next_block_file_offset;
	fp -> block_rolling_chain[empty_block_no].block_start_in_file_bits = fp -> next_block_file_bits;

	int is_data_error = 0;
	int is_eof = 0;
	while(1){
		int is_block_end = 0;
		int is_gzip_unit_end = 0;
		while(1){
			seekgz_try_read_some_zipped_data(fp, &is_eof);
			if(fp -> stem.avail_in < 1 && is_eof)break;

			if(out_txt_size < out_txt_used *4/3){
				out_txt_size *= 2;
				out_txt = realloc(out_txt,out_txt_size);
			}

			int old_avail_in = fp -> stem.avail_in;
			fp -> stem.avail_out = out_txt_size - out_txt_used;
			fp -> stem.next_out = (unsigned char *)out_txt + out_txt_used;

			void * old_next_in = fp -> stem.next_in;
				
				//SUBREADprintf("i  INFLATINHG: FILEPOS=%llu  IN_AVAIL=%d  OUT_AVAIL=%d\n", seekgz_ftello(fp), fp -> stem.avail_in, fp -> stem.avail_out);
			int ret_ifl = inflate(&(fp -> stem), Z_BLOCK);
			int have = (out_txt_size - out_txt_used) - fp -> stem.avail_out;
			if(ret_ifl != Z_OK && ret_ifl != Z_STREAM_END){
				SUBREADprintf("ERR  INFLATINHG: RET=%d BY IN %d zipped bytes  ==>  %d PLAIN TXT\n", ret_ifl, old_avail_in - fp -> stem.avail_in, have);
				SUBREADprintf("ERR  INFLATINHG: BEND=%d\n", ( fp -> stem.data_type & 128 )&& !(fp -> stem.data_type & 64));
			}

			//if(ret_ifl != Z_STREAM_END){
			//	SUBREADprintf("INF: STREAM_END\n");
			//}
			if(ret_ifl != Z_OK && ret_ifl != Z_STREAM_END){
				is_data_error = 1;
				break;
			}
			int zipped_data_used = (void *)fp -> stem.next_in - old_next_in;

			fp -> in_zipped_buff_read_ptr += zipped_data_used;
			if( ( fp -> stem.data_type & 128 )&& !(fp -> stem.data_type & 64)) is_block_end = 1;
			if(ret_ifl == Z_STREAM_END) is_gzip_unit_end = 1;

			//if(is_block_end)
			//	fp -> rolling_dict_window_used = 0;
			//else
			seekgz_update_current_window(fp, out_txt + out_txt_used, have);

			out_txt_used += have;

			if(is_block_end){
				fp -> next_block_file_offset = seekgz_ftello(fp);
				fp -> next_block_file_bits = fp->stem.data_type & 7;
			}

			if(is_gzip_unit_end){
				seekgz_skip_gzfile_header(fp, 8);
				inflateReset(&fp->stem);
			}

			if(is_block_end || is_gzip_unit_end) break;
		}
		if(out_txt_used >= MINIMUM_GZIP_BLOCK_SIZE || is_data_error ) break;
		if(fp -> stem.avail_in < 1 && is_eof)break;
	}

	if(is_data_error) return -1;

	//SUBREADprintf("o  END_1BLK: GOT %d PLAIN TXT i; first char: '%c'  fp==%p\n\n", out_txt_used, out_txt[0], fp);
	if(out_txt_used>0){
		fp -> block_rolling_chain[empty_block_no].block_txt = out_txt;
		fp -> block_rolling_chain[empty_block_no].block_txt_size = out_txt_used;
		seekgz_find_linebreaks(fp, empty_block_no);
		return out_txt_used;
	}else{
		free(out_txt);
		return 0;
	}
}

int seekgz_load_more_blocks( seekable_zfile_t * fp , int bytes_ahead, subread_lock_t * read_lock ){
	int loaded_bytes=0, ret=0;

	while(1){
		int empty_block_no = -1;
		subread_lock_occupy(&fp -> write_lock );
		if(read_lock) subread_lock_occupy(read_lock);
		if(fp -> blocks_in_chain < SEEKGZ_CHAIN_BLOCKS_NO){
			empty_block_no = fp-> block_chain_current_no + fp -> blocks_in_chain;
			if(empty_block_no >= SEEKGZ_CHAIN_BLOCKS_NO) empty_block_no -= SEEKGZ_CHAIN_BLOCKS_NO;
		}
		if(read_lock) subread_lock_release(read_lock);

		if(empty_block_no < 0 || (bytes_ahead>=0 && loaded_bytes >= bytes_ahead)){
			subread_lock_release(&fp -> write_lock );
			break;
		}
		int oneblk_ret = seekgz_load_1_block(fp, empty_block_no);

		//SUBREADprintf("LOAD_1BLK: RET=%d EMPTY_NO=%d\n", oneblk_ret, empty_block_no);
		// 1block ret < 0 : bad error
		// 1block ret > 0 : no data is loaded 
		// 1block ret == 0 : normal
		if(oneblk_ret < 0){
			subread_lock_release(&fp -> write_lock );
			ret = oneblk_ret;
			break;
		}else if(0 < oneblk_ret){
			if(read_lock) subread_lock_occupy(read_lock);	
			fp -> blocks_in_chain++;
			if(read_lock) subread_lock_release(read_lock);	
			loaded_bytes+= oneblk_ret;
		}else if(seekgz_eof(fp)){
			subread_lock_release(&fp -> write_lock );
			break;
		}
		subread_lock_release(&fp -> write_lock );
	}
	return ret;
}

int seqs = 0;

int seekgz_preload_buffer( seekable_zfile_t * fp , subread_lock_t * read_lock){
	int do_preload = 0;
	seqs+=1;// (unsigned int)(fp->rolling_dict_window[0]);
	if(( read_lock || !fp -> has_multi_thread_accessed) && fp -> blocks_in_chain < 3 && !seekgz_eof(fp) ) do_preload =1;
	else if(read_lock && fp -> blocks_in_chain < SEEKGZ_CHAIN_BLOCKS_NO ){
		if(seqs >= 2000) {
			do_preload = 1;
			seqs = 0;
		}
	}
	if(!fp -> has_multi_thread_accessed && read_lock) fp -> has_multi_thread_accessed =1;

	if(do_preload){
		int decompress_ret = seekgz_load_more_blocks(fp, -1, read_lock);
		return decompress_ret;
	}
	return 0;
}

// return : > 0  : N bytes loaded
//          == 0 : EOF
//          < 0  :  data error
int seekgz_gets(seekable_zfile_t * fp, char * buff, int buff_len){
	//if(fp -> blocks_in_chain<3)SUBREADprintf("GTS: %d BLK, %d AVI\n", fp -> blocks_in_chain, fp -> stem.avail_in);
	int line_write_ptr = 0, is_end_line = 0;
	if(fp->blocks_in_chain < 1 && seekgz_eof(fp)){
		return 0;
	}
	while(1){
		int consumed_bytes;

		assert(fp->blocks_in_chain>0); 
		seekable_decompressed_block_t *cblk = fp -> block_rolling_chain+fp -> block_chain_current_no;
		if( cblk -> linebreaks >0 && fp-> current_block_txt_read_ptr <= cblk->linebreak_positions[cblk -> linebreaks-1] ){
			if(fp-> current_block_txt_read_ptr <=cblk->linebreak_positions[0]){
				consumed_bytes = cblk->linebreak_positions[0] - fp-> current_block_txt_read_ptr +1;
			}else{
				int high_idx = cblk -> linebreaks - 1 , low_idx = 0, idx = -1;
				while(1){
					if(high_idx <= low_idx+1) { 
						idx = min(low_idx, high_idx); 
						break;
					}
					int mid_idx = (high_idx+low_idx)/2;
					if(cblk->linebreak_positions[mid_idx] > fp -> current_block_txt_read_ptr)
						high_idx = mid_idx;
					else if( cblk->linebreak_positions[mid_idx] < fp -> current_block_txt_read_ptr )
						low_idx = mid_idx;
					else{
						idx = high_idx;
						break;
					}
				}
				
				if(idx>0)idx--;
				while( cblk->linebreak_positions[idx+1] < fp -> current_block_txt_read_ptr) idx++;
				assert( idx <  cblk -> linebreaks-1 );
				consumed_bytes = cblk->linebreak_positions[idx+1] - fp-> current_block_txt_read_ptr +1;
			}
			is_end_line = 1;
		} else
			consumed_bytes = cblk -> block_txt_size - fp-> current_block_txt_read_ptr ;

		if(buff){
			int cp_bytes = min(consumed_bytes, buff_len - line_write_ptr);
			memcpy( buff + line_write_ptr , cblk-> block_txt + fp-> current_block_txt_read_ptr, cp_bytes );
			line_write_ptr += cp_bytes;
			if(buff)buff[line_write_ptr] = 0;
		}

		//SUBREADprintf("OUT_TXT  in BLC %d / %d %p  ENDL=%d CUNS_BYTES=%d / %d / %d : '''%s'''\n", fp -> block_chain_current_no, fp -> blocks_in_chain, fp, is_end_line, consumed_bytes,  fp-> current_block_txt_read_ptr, cblk -> block_txt_size ,  buff);

		fp-> current_block_txt_read_ptr += consumed_bytes;

		if( fp-> current_block_txt_read_ptr >= cblk -> block_txt_size ){
			//SUBREADprintf("CHAIN CONSUME_0 : %d ( - 1 ) %p\n", fp -> block_chain_current_no, cblk->block_txt);
			free(cblk->block_txt);
			//SUBREADprintf("CHAIN CONSUME_X : %d ( - 1 ) %p\n", fp -> block_chain_current_no, cblk->linebreak_positions);
			free(cblk->linebreak_positions);
			//SUBREADprintf("CHAIN CONSUME_Z : %d ( - 1 )\n", fp -> block_chain_current_no);
			fp-> current_block_txt_read_ptr = 0;
			fp -> block_chain_current_no++;
			if(fp -> block_chain_current_no>= SEEKGZ_CHAIN_BLOCKS_NO) fp -> block_chain_current_no=0;
			fp->blocks_in_chain --;
		}

		if(is_end_line) break;
	}
	return line_write_ptr;
}

int seekgz_next_int8(seekable_zfile_t * fp){
	if(fp -> blocks_in_chain<1){
		seekgz_load_more_blocks(fp, -1, NULL);
		if(fp -> blocks_in_chain<1) return -1;
	}
	seekable_decompressed_block_t *cblk = fp -> block_rolling_chain+fp -> block_chain_current_no;
	int ret = cblk->block_txt[fp-> current_block_txt_read_ptr];
	fp-> current_block_txt_read_ptr++;
	if(fp-> current_block_txt_read_ptr == cblk -> block_txt_size){
		free(cblk->block_txt);
		free(cblk->linebreak_positions);
		fp-> current_block_txt_read_ptr = 0;
		fp -> block_chain_current_no++;
		if(fp -> block_chain_current_no>= SEEKGZ_CHAIN_BLOCKS_NO) fp -> block_chain_current_no=0;
		fp->blocks_in_chain --;
	}
	return ret<0?256+ret:ret;
}


int seekgz_next_char(seekable_zfile_t * fp){ // MUST BE PROTECTED BY read_lock
	//SUBREADprintf("gen_next_char: block_in_chain=%d\n", fp -> blocks_in_chain);
	//if(fp -> blocks_in_chain<3)SUBREADprintf("1CH: %d BLK, %d AVI\n", fp -> blocks_in_chain, fp -> stem.avail_in);
	if(fp -> blocks_in_chain<1){
		return EOF;
	}
	seekable_decompressed_block_t *cblk = fp -> block_rolling_chain+fp -> block_chain_current_no;
	int ret = cblk->block_txt[fp-> current_block_txt_read_ptr];
	fp-> current_block_txt_read_ptr++;
	if(fp-> current_block_txt_read_ptr == cblk -> block_txt_size){
		free(cblk->block_txt);
		free(cblk->linebreak_positions);
		fp-> current_block_txt_read_ptr = 0;
		fp -> block_chain_current_no++;
		if(fp -> block_chain_current_no>= SEEKGZ_CHAIN_BLOCKS_NO) fp -> block_chain_current_no=0;
		fp->blocks_in_chain --;
	}
	return ret;
}

void seekgz_close(seekable_zfile_t * fp){
	int ii;
	fclose(fp -> gz_fp);
	free(fp -> in_zipped_buffer);
	for(ii = 0; ii < fp -> blocks_in_chain ; ii++){
		int iv = ii + fp -> block_chain_current_no;
		if(iv >=SEEKGZ_CHAIN_BLOCKS_NO) iv-=SEEKGZ_CHAIN_BLOCKS_NO;
		free(fp -> block_rolling_chain[iv].block_txt);
		free(fp -> block_rolling_chain[iv].linebreak_positions);
	}

	inflateEnd(&fp -> stem);
	subread_destroy_lock(&fp -> write_lock);
}

int autozip_open(const char * fname, autozip_fp * fp){
	int ret = -1;
	memset(fp, 0, sizeof(autozip_fp));
	strcpy(fp -> filename, fname);

	FILE * tstfp = fopen(fname,"rb");
	if(!tstfp) return ret;

	int cc1, cc2;
	cc1 = fgetc(tstfp);
	cc2 = fgetc(tstfp);

	if(cc2 == 0x8b && cc1 == 0x1f){
		//fclose(tstfp);
		fp -> is_plain = 0;
		int iret = seekgz_open(fname, &fp -> gz_fp, tstfp);
		if(iret >= 0) ret = 1;
		else ret = -1;
	}else{
		if(cc1 != EOF && cc2 != EOF){
			fp -> first_chars[0] = cc1;
			fp -> first_chars[1] = cc2;
			fp -> is_first_chars = 2;
		}
		fp -> plain_fp = tstfp;
		fp -> is_plain = 1;
		ret = 0;
	}

	return ret;
}

void autozip_close(autozip_fp * fp){
	if(fp -> is_plain) fclose(fp -> plain_fp);
	else
		seekgz_close(&fp -> gz_fp);
} 

int autozip_gets(autozip_fp * fp, char * buf, int buf_size){
	int ret = 0;
	if(fp -> is_plain) {
		int base0 = 0;
		if( fp -> is_first_chars ){
			buf[0]=fp -> first_chars[0];
			buf[1]=fp -> first_chars[1];
			base0 = 2;
			fp -> is_first_chars = 0;
		}
		buf[2]=0;

		char * retc = fgets(buf + base0, buf_size, fp -> plain_fp);
		if(retc == NULL && base0 == 0) ret = 0;
		else ret = strlen(buf);
	}else{
		seekgz_preload_buffer(&fp->gz_fp, NULL);
		ret = seekgz_gets(&fp -> gz_fp, buf, buf_size);
	}
	//SUBREADprintf("READBUF '%s'\n", buf);
	return ret;
}

int autozip_getch(autozip_fp * fp){
	int ret = 0;
	if(fp -> is_plain) {
		if( fp -> is_first_chars ) return fp -> first_chars[2-( fp -> is_first_chars -- ) ];
		ret = fgetc(fp -> plain_fp);
		if(ret==EOF)ret = -1;
	}else{
		ret = seekgz_next_int8(&fp -> gz_fp);
	}
	return ret;
}

void autozip_rewind(autozip_fp * fp){
	char fname [MAX_FILE_NAME_LENGTH+1];
	strcpy(fname, fp -> filename);

	autozip_close(fp);
	autozip_open(fname, fp);
}

void parallel_gzip_writer_init(parallel_gzip_writer_t * pzwtr, char * output_filename, int total_threads){
	memset(pzwtr, 0, sizeof(parallel_gzip_writer_t));
	pzwtr -> threads = total_threads;
	pzwtr -> thread_objs = calloc(sizeof(parallel_gzip_writer_thread_t), total_threads);
	
	pzwtr -> os_file = f_subr_open(output_filename, "wb");
	fputc(31, pzwtr -> os_file);
	fputc(139, pzwtr -> os_file);
	fputc(8, pzwtr -> os_file);
	fputc(0, pzwtr -> os_file);
	fputc(0, pzwtr -> os_file); // MTIME x 4
	fputc(0, pzwtr -> os_file);
	fputc(0, pzwtr -> os_file);
	fputc(0, pzwtr -> os_file);
	fputc(4, pzwtr -> os_file); // XFL : fastest
	fputc(255, pzwtr -> os_file); // OS=unknown
	int x1;
	for(x1=0; x1<total_threads; x1++){
		pzwtr -> thread_objs[x1].thread_no = x1;
		deflateInit2(&pzwtr -> thread_objs[x1].zipper, SAMBAM_COMPRESS_LEVEL_NORMAL, Z_DEFLATED, SAMBAM_GZIP_WINDOW_BITS, 8, Z_DEFAULT_STRATEGY);
	}
	pzwtr -> CRC32 = crc32(0, NULL, 0);
}

void parallel_gzip_writer_add_text(parallel_gzip_writer_t * pzwtr, char * text, int tlen, int thread_no){
	parallel_gzip_writer_thread_t *tho = pzwtr -> thread_objs + thread_no;
	if(tlen + tho -> in_buffer_used >= PARALLEL_GZIP_TXT_BUFFER_SIZE){
		SUBREADprintf("Insufficient gzip buffer.\n");
		return;
	}
	memcpy(tho -> in_buffer + tho -> in_buffer_used, text, tlen);
	tho -> in_buffer_used += tlen;
}

void parallel_gzip_zip_texts(parallel_gzip_writer_t * pzwtr, int thread_no, int eof_marker){
	parallel_gzip_writer_thread_t *tho = pzwtr -> thread_objs + thread_no;
	int write_txt_ptr = 0;

	tho -> out_buffer_used = 0;
	tho -> CRC32 = crc_pos(tho -> in_buffer, tho -> in_buffer_used);

	while(tho -> in_buffer_used - write_txt_ptr > 0 || eof_marker){
		tho -> zipper . next_in = (unsigned char*)tho -> in_buffer + write_txt_ptr;
		tho -> zipper . avail_in = tho -> in_buffer_used - write_txt_ptr;
		tho -> zipper . next_out = (unsigned char*)tho -> out_buffer + tho -> out_buffer_used;
		tho -> zipper . avail_out = PARALLEL_GZIP_ZIPPED_BUFFER_SIZE - tho -> out_buffer_used;
		int defret = deflate(&tho -> zipper, eof_marker?Z_FINISH:Z_FULL_FLUSH);
		int consumed_input = tho -> in_buffer_used - write_txt_ptr - tho -> zipper . avail_in;
		int generated_output = PARALLEL_GZIP_ZIPPED_BUFFER_SIZE - tho -> out_buffer_used - tho -> zipper.avail_out;

		if(defret == Z_OK || defret == Z_STREAM_END ){
			write_txt_ptr += consumed_input;
			tho -> out_buffer_used += generated_output;
		}else{
			SUBREADprintf("Cannot compress the zipped output: %d with in_len=%d, consumed=%d and out_aval=%d\n", defret, tho -> in_buffer_used, consumed_input, tho -> zipper.avail_out);
			break;
		}
		if(eof_marker)break;
	}

	tho -> plain_length = tho -> in_buffer_used;
	tho -> in_buffer_used =0;
}

// because we have to keep sync between three fastq files, the flush function has to be manually called three times at the same time point.
// otherwise R1, I2 and R2 files will have inconsistent read orders.
// the outer program has to check if any of the three in_buffers is full.
void parallel_gzip_writer_flush(parallel_gzip_writer_t * pzwtr, int thread_no){
	parallel_gzip_writer_thread_t *tho = pzwtr -> thread_objs + thread_no;

	if(tho -> out_buffer_used>0){
		int fret = fwrite(tho -> out_buffer, 1, tho -> out_buffer_used, pzwtr -> os_file);
		if(fret != tho ->out_buffer_used) SUBREADprintf("Cannot write the zipped output: %d\n", fret);

		if(tho -> plain_length){
			unsigned int CRCnew = crc32_combine( pzwtr -> CRC32, tho -> CRC32, tho -> plain_length );
			pzwtr -> plain_length += tho -> plain_length;
			pzwtr -> CRC32 = CRCnew;
		}
	}
	tho -> out_buffer_used =0;
	tho -> plain_length = 0;
}

void plgz_finish_in_buffers(parallel_gzip_writer_t * pzwtr){
	int x1;
	for(x1=0; x1<pzwtr -> threads; x1++){
		if(pzwtr -> thread_objs[x1].in_buffer_used<1) continue;
		parallel_gzip_zip_texts(pzwtr, x1, 0);
		parallel_gzip_writer_flush(pzwtr, x1);
	}
}

void parallel_gzip_writer_close(parallel_gzip_writer_t * pzwtr){
	int x1;
	plgz_finish_in_buffers(pzwtr);

	// write the last "Z_FINISH" 0-len block to tell gzip that the zipped chunk end
	pzwtr -> thread_objs[0].in_buffer_used = 0;
	parallel_gzip_zip_texts(pzwtr, 0, 1);
	parallel_gzip_writer_flush(pzwtr, 0);

	for(x1=0; x1<pzwtr -> threads; x1++)
		deflateEnd(&pzwtr -> thread_objs[x1].zipper);
	
	fwrite(&pzwtr -> CRC32,        4, 1, pzwtr -> os_file);
	fwrite(&pzwtr -> plain_length, 4, 1, pzwtr -> os_file); // write the first (lowest) 4 bytes in a 64-bit integer -- as gzip file format required original (uncompressed) input data modulo 2^32.
	fclose(pzwtr -> os_file);
	free(pzwtr -> thread_objs);
}


int parallel_gzip_writer_add_read_fqs_scRNA(parallel_gzip_writer_t**outfps, char * bambin, int thread_no){
	int reclen=0;
	parallel_gzip_writer_t * outR1fp = outfps[0];
	parallel_gzip_writer_t * outI1fp = outfps[1];
	parallel_gzip_writer_t * outR2fp = outfps[2];

	memcpy(&reclen, bambin,4);
	int flag = 0, l_seq = 0, l_read_name = 0, n_cigar_ops = 0;
	memcpy(&l_read_name, bambin+12,1);
	memcpy(&n_cigar_ops, bambin+16,2);
	memcpy(&flag, bambin+18,2);
	memcpy(&l_seq, bambin+20,4);

	int x1=0;

	parallel_gzip_writer_add_text(outR2fp,"@",1,thread_no);
	parallel_gzip_writer_add_text(outR1fp,"@",1,thread_no);
	parallel_gzip_writer_add_text(outI1fp,"@",1,thread_no);
	char * readname = bambin+36;
	parallel_gzip_writer_add_text(outR1fp,readname, 12,thread_no);
	parallel_gzip_writer_add_text(outR2fp,readname, 12,thread_no);
	parallel_gzip_writer_add_text(outI1fp,readname, 12,thread_no);
	parallel_gzip_writer_add_text(outR1fp,"\n",1,thread_no);
	parallel_gzip_writer_add_text(outR2fp,"\n",1,thread_no);
	parallel_gzip_writer_add_text(outI1fp,"\n",1,thread_no);

	char * R1seq = bambin+36+13;
	int R1len = 0;
	for(R1len=0; R1seq[R1len] && R1seq[R1len]!='|' ;R1len++);
	char * R1qual = R1seq + R1len + 1;
	parallel_gzip_writer_add_text(outR1fp,R1seq, R1len,thread_no);
	parallel_gzip_writer_add_text(outR1fp,"\n+\n",3,thread_no);
	parallel_gzip_writer_add_text(outR1fp,R1qual, R1len,thread_no);
	parallel_gzip_writer_add_text(outR1fp,"\n",1,thread_no);

	char * I1seq = R1qual + R1len + 1;
	int I1len = 0;
	for(I1len=0; I1seq[I1len] && I1seq[I1len]!='|' ;I1len++);
	char * I1qual = I1seq + I1len + 1;
	parallel_gzip_writer_add_text(outI1fp,I1seq, I1len,thread_no);
	parallel_gzip_writer_add_text(outI1fp,"\n+\n",3,thread_no);
	parallel_gzip_writer_add_text(outI1fp,I1qual, I1len,thread_no);
	parallel_gzip_writer_add_text(outI1fp,"\n",1,thread_no);

	char oseq[l_seq+1];
	int seqbase = 36+l_read_name+n_cigar_ops*4;
	for(x1=0; x1<l_seq;x1++)oseq[x1]="=ACMGRSVTWYHKDBN"[(bambin[ seqbase+(x1/2)] >> (((x1%2)?0:1) *4) )&0xf];
	oseq[l_seq]=0;
	if(flag & 16) reverse_read(oseq,l_seq, GENE_SPACE_BASE);
	parallel_gzip_writer_add_text(outR2fp, oseq, l_seq,thread_no);
	parallel_gzip_writer_add_text(outR2fp,"\n+\n",3,thread_no);

	seqbase = 36+l_read_name+n_cigar_ops*4 + ( l_seq+1 )/2;
	for(x1=0; x1<l_seq;x1++)oseq[x1]=33+bambin[ seqbase+x1];
	if(flag & 16)reverse_quality(oseq, l_seq);
	oseq[l_seq]=0;
	parallel_gzip_writer_add_text(outR2fp, oseq, l_seq,thread_no);
	parallel_gzip_writer_add_text(outR2fp,"\n",1,thread_no);
	return 0;
}

