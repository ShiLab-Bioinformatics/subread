#include <assert.h>
#include "seek-zlib.h"


int main(int argc, char ** argv){
	unsigned int tested_cell_total = atoi(argv[4]);
	unsigned int tested_inc ;
	tested_inc = 337;

	seekable_zfile_t * fp = malloc(sizeof(seekable_zfile_t));
	seekable_position_t * pos = malloc(tested_cell_total*sizeof(seekable_position_t));
	seekgz_open(argv[1], fp, NULL);

	char buf[1002];
	char should[tested_cell_total][1002];
	long long int text_pos[tested_cell_total];
	long long int alll = 0, marked = 0;
	long long int full_size = atoll(argv[2]);
	long long int tested_cell_no = 0;
	long long int step = full_size / tested_cell_total;

	unsigned int rand_seed = atoi(argv[3]);
	//step = 10;

	int write_cell = -1;

	long long int first_all = 0;
	while(1){
		int rl = seekgz_gets(fp, buf, 1000);
		unsigned int inchunk = fp -> in_chunk_offset;

		if(0){
			if(rl>92)
				fprintf(stderr, "LEN=%d; READ=%s", rl, buf);
			if(rl>52 && rl<91)
				fprintf(stderr, "LEN=%d; READ=%s", rl, buf);
		}
		if(rl<1) break;
		//fprintf(stdout, "%s", buf);
		alll += rl;

		if(write_cell >=0){
			strcpy(should[write_cell], buf);
			write_cell = -1;
		}

		if(1&& alll - rand_seed > tested_cell_no * step )
		{
			if(tested_cell_no == 0)
				first_all = alll;
			if(tested_cell_no < tested_cell_total){
				write_cell = ((1+tested_cell_no) * tested_inc) % tested_cell_total;
				seekgz_tell(fp, pos+write_cell);
				//assert(pos[write_cell].block_gzfile_offset > 10);
				text_pos[write_cell] = alll;
				//if(alll==925826012||alll==889250153){
				//if(inchunk == fp->in_block_offset + 1 || inchunk == fp->in_block_offset - 1 )
				//	fprintf(stderr, "MATCH: IN_BLOCK_OFFSET=%u IN_CHUNK_OFFSET=%u/%u IS_LAST=%d\n", fp->in_block_offset, inchunk, fp -> txt_buffer_used, fp->is_the_last_chunk);
				//}
				if(alll==344780683){
					char * quickview = malloc(100000);
					quickview[0]=0;
					//memcpy(quickview, fp -> block_dict_window, fp -> block_dict_window_size);
					quickview[fp -> block_dict_window_size]=0;
					fprintf(stderr, "----------------------------------------------------\n%s\n=====================================================\n\n", quickview);
					free(quickview);
				}
			}
			tested_cell_no++;
		}
	}

	fprintf(stderr, "TOTAL=%lld\n",alll);
	assert(tested_cell_no >= tested_cell_total);

	int i, valid=0;
	long long int last_all = 0;
	for(i=0;i< tested_cell_total ;i++){
		//assert(pos[i].block_gzfile_offset > 10);
		seekgz_seek(fp, pos+i);
		unsigned int inchunk = fp -> in_chunk_offset, chunk_size = fp -> txt_buffer_used;
		//fprintf(stderr, "JUMPTO=%u\n", pos[i].block_gzfile_offset);
		int rl = seekgz_gets(fp, buf, 1000);
		if(rl <= 0) break;
		if(strcmp(should[i], buf)!=0)
		{
			char * quickshow = malloc(1000000);
			quickshow[0]=0;
			//memcpy(quickshow , pos[i].dict_window, pos[i].block_dict_window_size);
			quickshow[ pos[i].block_dict_window_size ] = 0;

			fprintf(stderr, "=================================================\nMATCH:LEN=%d; TOTAL=%lld; GZFP=%lld; INBLOCK=%u; INCHUNK=%u/%u\t\tMATCH=%d; \nREAD=%s\nORGN=%s\n%s\n", rl, text_pos[i], pos[i].block_gzfile_offset, pos[i].in_block_text_offset, inchunk, chunk_size, strcmp(should[i], buf), buf, should[i], quickshow);
			free(quickshow);
		}
		else
			valid++;
		last_all = text_pos[i];
	}
	
	fprintf(stderr, "FINISHED size=%lld    first=%lld  [-1]=%lld   DOTS=%d/%u  rand=%u\n", full_size, first_all, last_all, valid,tested_cell_total , rand_seed);
}
