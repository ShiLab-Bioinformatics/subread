#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>

#include "subread.h"
#include "HelperFunctions.h"
#include "seek-zlib.h"
#include "gene-algorithms.h"
#include "input-files.h"
#include "input-blc.h"
#include "sambam-file.h"

#ifdef __MINGW32__
#define SR_PATH_SPLIT_STR "\\"
#else
#define SR_PATH_SPLIT_STR "/"
#endif

struct iBLC_scan_t{
	char out_format_string[MAX_FILE_NAME_LENGTH];
	char filter_format_string[MAX_FILE_NAME_LENGTH];
	int found_answer;
	int filter_is_gzipped;
	int bcl_is_gzipped;
	int reads_per_cluster;
	int read_lengths[INPUT_BLC_MAX_READS];
	int read_is_index[INPUT_BLC_MAX_READS];
};

int iBLC_guess_scan(struct iBLC_scan_t * scancon, char * data_dir ){
	DIR * this_level = opendir(data_dir);
	if(this_level == NULL) return -1;
	struct dirent *dp;
	int filter_found = 0, bcl_found = 0;
	char testfile_name[MAX_FILE_NAME_LENGTH];
	while ((dp = readdir (this_level)) != NULL) {
		#ifdef __MINGW32__
		if(1){
		#else
		if(dp -> d_type == DT_DIR && dp->d_name[0]!='.'){
		#endif
			strcpy(testfile_name,data_dir);
			strcat(testfile_name, SR_PATH_SPLIT_STR);
			strcat(testfile_name, dp->d_name);
			//SUBREADprintf("DIG: %s\n", testfile_name);
		#ifdef __MINGW32__
			if(strcmp( dp->d_name, ".") && strcmp( dp->d_name,".."))
				if(0==iBLC_guess_scan( scancon, testfile_name)) continue;
		#else
			if(iBLC_guess_scan( scancon, testfile_name))return -1;
		}else if(dp -> d_type == DT_REG){
		#endif
			//SUBREADprintf( "%s  %s  %p  %p\n" , data_dir, dp->d_name , strstr( dp->d_name , "0001.bcl." ) , strstr( dp->d_name , ".bci") );
			if(0==strcmp(dp->d_name, "RunInfo.xml")){
				if(scancon->reads_per_cluster > 0){
					SUBREADprintf("ERROR: the root directory contains multiple scRNA data sets.\n");
					return -1;
				}

				strcpy(testfile_name, data_dir);    
				strcat(testfile_name, SR_PATH_SPLIT_STR);
				strcat(testfile_name, dp->d_name);
				FILE *fp = fopen(testfile_name,"rb");
				if(NULL == fp){
					SUBREADprintf("ERROR: cannot open the run info file: %s\n", testfile_name);
				}
				while(1){
					char inbuf[MAX_READ_LENGTH];
					if(!fgets( inbuf, MAX_READ_LENGTH-1, fp))break;
					if(strstr( inbuf, "<Read Number=\"" )){
						char * rbuf=NULL;
						int my_index = -1, is_idx = -1, rlen = -1, ii=0;
						strtok_r(inbuf, "\"", &rbuf);

						while(rbuf){
							char * sec = strtok_r(NULL, "\"", &rbuf);
							if(!sec) break;
							//printf("SEC %d : %s\n", ii, sec);
							if(ii == 0) my_index = atoi(sec);
							if(ii == 2) rlen = atoi(sec);
							if(ii == 4) is_idx = sec[0]=='Y';
							ii++;
						}
						assert(INPUT_BLC_MAX_READS>my_index);
						if(my_index >0 && is_idx >=0 && rlen>0){
							scancon -> read_lengths[my_index-1]=rlen;
							scancon -> read_is_index[my_index-1]=is_idx;
							scancon -> reads_per_cluster = max(scancon -> reads_per_cluster, my_index);
						}else assert( my_index >0 && is_idx >=0 && rlen>0 );
					}
				}
				fclose(fp);
				if(scancon -> reads_per_cluster <1){
					SUBREADprintf("ERROR: the format of RunInfo.xml is unknown\n");
					return -1;
				}
			}
			if(0==memcmp(data_dir+ strlen(data_dir)-5, SR_PATH_SPLIT_STR "L001",5 ) && strstr( dp->d_name , "s_1.filter")){
				autozip_fp tfp;
				strcpy(testfile_name, data_dir);    
				strcat(testfile_name, SR_PATH_SPLIT_STR);
				strcat(testfile_name, dp->d_name);
				int resop = autozip_open(testfile_name, &tfp);
				if(0 <= resop){
					autozip_close(&tfp);
					char * gen_fmt = str_replace(dp->d_name , "s_1.filter", "s_%d.filter");
					char * gen_fmt2 = str_replace(data_dir , SR_PATH_SPLIT_STR "L001", SR_PATH_SPLIT_STR "L%03d");
					strcpy(scancon -> filter_format_string, gen_fmt2);
					strcat(scancon -> filter_format_string, SR_PATH_SPLIT_STR);
					strcat(scancon -> filter_format_string, gen_fmt);
					free(gen_fmt2);
					free(gen_fmt);
					filter_found = resop + 1;
				}
			}
			if(0==memcmp(data_dir+ strlen(data_dir)-5, SR_PATH_SPLIT_STR "L001",5 ) && strstr( dp->d_name , "0001.bcl." ) && !strstr( dp->d_name , ".bci") ){
				int tti;
				bcl_found = 1;
				char * gen_fmt = str_replace(dp->d_name , "0001.bcl.", "%04d.bcl.");
				
				for(tti = 0; tti<22; tti++){
					strcpy(testfile_name, data_dir);	
					strcat(testfile_name, SR_PATH_SPLIT_STR);
					sprintf(testfile_name+strlen(testfile_name), gen_fmt, 1, 2+tti);
					autozip_fp tfp;
					int resop = autozip_open(testfile_name, &tfp);
					//printf("%d === %s    %s\n", resop, gen_fmt, testfile_name);
					if(0<=resop){
						scancon -> bcl_is_gzipped = resop;
						autozip_close(&tfp);
					}else bcl_found=0;
				}
				if(bcl_found){
					char * gen_fmt2 = str_replace(data_dir , SR_PATH_SPLIT_STR "L001", SR_PATH_SPLIT_STR "L%03d");
					strcpy(scancon -> out_format_string, gen_fmt2);	
					free(gen_fmt2);
					strcat(scancon -> out_format_string, SR_PATH_SPLIT_STR);
					strcat(scancon -> out_format_string, gen_fmt);
				}

				free(gen_fmt);
			}
		}
	}

	if(bcl_found && filter_found){
		scancon -> found_answer = 1;
		scancon -> filter_is_gzipped = filter_found - 1;
	}

	closedir(this_level);
	return 0;
}

int iBLC_guess_format_string(char * data_dir, int * cluster_bases, char * format_string, char * filter_format, int * bcl_is_gzipped, int * filter_is_gzipped, int * read_lens, int * is_index){
	struct iBLC_scan_t sct;
	memset(&sct, 0, sizeof(sct));
	int tii = iBLC_guess_scan(&sct, data_dir);

	if(tii || ! sct.found_answer) return -1;
	strcpy(format_string, sct.out_format_string);
	strcpy(filter_format, sct.filter_format_string);
	*filter_is_gzipped = sct.filter_is_gzipped;
	*bcl_is_gzipped = sct.bcl_is_gzipped;
	*cluster_bases=0;

	for(tii=0; tii<sct.reads_per_cluster; tii++){
		if(sct.read_lengths[tii]<1) return -1;
		read_lens[tii] = sct.read_lengths[tii];
		is_index[tii] = sct.read_is_index[tii];
		(*cluster_bases) += sct.read_lengths[tii];
		read_lens[tii+1]=0;
	}
		
	return 0;
}

void iBLC_close_batch(input_BLC_t * blc_input){
	int ii;
	if(blc_input->is_EOF) return;
	if(NULL == blc_input -> bcl_gzip_fps && blc_input -> bcl_is_gzipped)return;
	if(NULL == blc_input -> bcl_fps && !blc_input -> bcl_is_gzipped)return;
	for(ii=0; ii < blc_input->total_bases_in_each_cluster; ii++){
		if(blc_input -> bcl_is_gzipped){
			seekgz_close(blc_input -> bcl_gzip_fps[ii]);
			free(blc_input -> bcl_gzip_fps[ii]);
			blc_input -> bcl_gzip_fps[ii] = NULL;
		}else{
			fclose(blc_input -> bcl_fps[ii]);
			blc_input -> bcl_fps[ii] = NULL;
		}
	}
	if(blc_input -> filter_is_gzipped){
		seekgz_close(blc_input -> filter_gzip_fp);
		free(blc_input -> filter_gzip_fp);
		blc_input -> filter_gzip_fp = NULL;
	}else{
		fclose(blc_input -> filter_fp);
		blc_input -> filter_fp = NULL;
	}

	if( blc_input -> bcl_is_gzipped ){
		free(blc_input -> bcl_gzip_fps);
		blc_input -> bcl_gzip_fps = NULL;
	}else{
		free(blc_input -> bcl_fps);
		blc_input -> bcl_fps = NULL;
	}
}

int iBLC_open_batch(input_BLC_t * blc_input ){
	char fname[MAX_FILE_NAME_LENGTH];
	iBLC_close_batch(blc_input);
	int fii, xx;
	blc_input -> is_EOF=1;

	if(blc_input -> bcl_gzip_fps == NULL) blc_input -> bcl_gzip_fps = calloc( sizeof(void *), blc_input -> total_bases_in_each_cluster ); // for both FILE** and seekgz **
	for(fii = 0; fii < blc_input -> total_bases_in_each_cluster; fii++){
		sprintf(fname, blc_input -> bcl_format_string, blc_input -> current_lane, fii+1);
		if(blc_input -> bcl_is_gzipped){
			blc_input -> bcl_gzip_fps[fii] = calloc( sizeof(seekable_zfile_t), 1);
			int rv = seekgz_open(fname, blc_input -> bcl_gzip_fps[fii], NULL);
			if(rv){
				SUBREADprintf("ERROR: Unable to open %s\n", fname);
				return -1;
			}
			for(xx = 0; xx < 4; xx++) seekgz_next_int8(blc_input -> bcl_gzip_fps[fii]); // skip the first 32-b integer
		}else{
			blc_input -> bcl_fps[fii] = fopen(fname, "rb");
			 if(NULL == blc_input -> bcl_fps[fii]){
				SUBREADprintf("ERROR: Unable to open %s\n", fname);
				return -1;
			}
			for(xx = 0; xx < 4; xx++) fgetc(blc_input -> bcl_fps[fii]); // skip the first 32-b integer
		}
	}

	sprintf(fname, blc_input -> filter_format_string, blc_input -> current_lane,blc_input -> current_lane);
	if(blc_input -> filter_is_gzipped){
		blc_input -> filter_gzip_fp = calloc( sizeof(seekable_zfile_t), 1);
		int rv = seekgz_open(fname, blc_input -> filter_gzip_fp, NULL);
		if(rv){
			SUBREADprintf("ERROR: Unable to open %s\n", fname);
			return -1;
		}
		for(xx = 0; xx < 12; xx++) seekgz_next_int8(blc_input -> filter_gzip_fp); // skip the 12-byte header
	}else{
		blc_input -> filter_fp = fopen(fname, "rb");
		if(NULL == blc_input -> filter_fp){
			SUBREADprintf("ERROR: Unable to open %s\n", fname);
			return -1;
		}
		for(xx = 0; xx < 12; xx++) fgetc(blc_input -> filter_fp); // skip the 12-byte header
	}
	blc_input -> is_EOF=0;
	return 0;
}

int iCache_open_batch( cache_BCL_t * cache_input){
	cache_input -> bcl_gzip_fps = calloc(sizeof(autozip_fp), cache_input -> total_bases_in_each_cluster);
	return 0;
}

int cacheBCL_go_chunk_end( cache_BCL_t * cache_input ){
	cache_input -> read_no_in_chunk= cache_input -> reads_available_in_chunk;
	return 0;
}
int cacheBCL_go_chunk_start( cache_BCL_t * cache_input ){
	cache_input -> read_no_in_chunk=0;
	return 0;
}
void cacheBCL_close(cache_BCL_t * cache_input){
	int x1;
	for(x1=0;x1<cache_input -> total_bases_in_each_cluster;x1++){
		if( cache_input -> bcl_gzip_fps[x1]. plain_fp || cache_input -> bcl_gzip_fps[x1].gz_fp.gz_fp )autozip_close(&cache_input -> bcl_gzip_fps[x1]);
		free(cache_input -> bcl_bin_cache[x1]);
	}
	free(cache_input -> bcl_gzip_fps);
	if(cache_input -> filter_fp. plain_fp || cache_input -> filter_fp.gz_fp.gz_fp )autozip_close(&cache_input -> filter_fp);
	free(cache_input -> lane_no_in_chunk);
	free(cache_input -> flt_bin_cache);
}

int cacheBCL_init( cache_BCL_t * cache_input, char * data_dir, int reads_in_chunk, int all_threads ){
	memset(cache_input, 0, sizeof( cache_BCL_t));
	subread_init_lock(&cache_input -> read_lock);
	int rv = iBLC_guess_format_string(data_dir, &cache_input -> total_bases_in_each_cluster, cache_input -> bcl_format_string, cache_input -> filter_format_string, &cache_input -> bcl_is_gzipped, &cache_input -> filter_is_gzipped, cache_input -> single_read_lengths, cache_input -> single_read_is_index);
	if(rv) return -1;
	cache_input -> current_lane = 1;
	cache_input -> reads_per_chunk = reads_in_chunk;
	cache_input -> bcl_bin_cache = malloc(sizeof(char*) * cache_input -> total_bases_in_each_cluster);
	int x1;
	for(x1 = 0; x1 < cache_input -> total_bases_in_each_cluster; x1++)
		cache_input -> bcl_bin_cache[x1] = malloc(reads_in_chunk);
	cache_input -> flt_bin_cache = malloc(reads_in_chunk*2);
	cache_input -> flt_bin_cache_size = reads_in_chunk*2;
	cache_input -> lane_no_in_chunk = malloc(reads_in_chunk);
	cache_input -> chunk_end_lane = 1;	// no "0th lane"
	cache_input -> all_threads = all_threads;
	return iCache_open_batch(cache_input)?1:0;
}


void iCache_close_one_fp( cache_BCL_t * cache_input, int bcl_no){
	autozip_fp * tfp = bcl_no<0?&cache_input-> filter_fp :(cache_input -> bcl_gzip_fps + bcl_no);
//	SUBREADprintf("CLOSE_AUTO: %d\n", bcl_no);
	autozip_close(tfp);
	memset(tfp,0,sizeof(autozip_fp));
}

int iCache_open_one_fp( cache_BCL_t * cache_input, int bcl_no, int lane_no){
	autozip_fp * tfp = bcl_no<0?&cache_input-> filter_fp :(cache_input -> bcl_gzip_fps + bcl_no);
	char fname[MAX_FILE_NAME_LENGTH+1];
	if(bcl_no <0)
		sprintf(fname,  cache_input -> filter_format_string, lane_no, lane_no);
	else
		sprintf(fname,  cache_input -> bcl_format_string, lane_no, bcl_no+1);

	int rv = autozip_open(fname, tfp);
//	SUBREADprintf("OPEN_AUTO %s = %d\n", fname, rv);
	if(rv>=0){
		int sk = bcl_no<0?12:4;
		for(; sk>0; sk--){
			autozip_getch(tfp);
			//SUBREADprintf("SKEP=%d\n", nch);
		}
	}else{
		memset(tfp,0,sizeof(autozip_fp));
	}
	return rv < 0;
}

// it returns the number of bytes loaded. 0 if no bytes are available.
int iCache_continuous_read_lanes( cache_BCL_t * cache_input, int bcl_no){
	autozip_fp * tfp = bcl_no<0?&cache_input-> filter_fp :(cache_input -> bcl_gzip_fps + bcl_no);
	int my_lane = cache_input -> chunk_start_lane ,wptr=0;
	char * wpt =  bcl_no<0?cache_input-> flt_bin_cache:cache_input -> bcl_bin_cache[bcl_no];
	int total_valid_reads = 0;
	int raw_ptr = 0;
	while(1){
		if(!(tfp -> plain_fp || tfp -> gz_fp.gz_fp)){
			int fpr = iCache_open_one_fp( cache_input, bcl_no, my_lane );
			if(fpr){
				if(bcl_no<0)cache_input -> last_chunk_in_cache = 1;
				break;
			}
		}
		while(1){
			int nch = autozip_getch( tfp );
			//if(bcl_no>=0)SUBREADprintf("NCH=%d\n", nch);
			if(nch>=0){
				if( bcl_no < 0 || cache_input -> flt_bin_cache [raw_ptr] ){
					if(bcl_no <0){
						if(nch>0) cache_input -> lane_no_in_chunk[total_valid_reads ++] = my_lane;
						if(wptr == cache_input ->  flt_bin_cache_size){
							cache_input ->  flt_bin_cache_size *= 1.6;
							wpt = cache_input-> flt_bin_cache = realloc( wpt, cache_input ->  flt_bin_cache_size );
						}
					}else total_valid_reads++;

					//if(wptr == 1000000)SUBREADprintf("BaseNo=%d; Nch=%d\n", bcl_no, nch);
					wpt[wptr++] = nch&0xff;
					if(total_valid_reads == cache_input -> reads_per_chunk) break;
				}
			}else{
//				SUBREADprintf("CONTINUOUSLY [%d-th base] BREAK : %d\n", bcl_no, nch);
				break;
			}
			raw_ptr ++;
		}
		if(total_valid_reads == cache_input -> reads_per_chunk) break;
		iCache_close_one_fp(cache_input, bcl_no);
		my_lane++;
	}
	if(bcl_no <0){
		cache_input -> reads_available_in_chunk = total_valid_reads;
		cache_input -> chunk_end_lane = my_lane;
	}
	//if(bcl_no<0) SUBREADprintf("CONTINUOUSLY [%d-th base] READ %d reads into the %d-th lanel fpntr = %ld\n", bcl_no,  total_valid_reads, my_lane, tfp -> plain_fp ? ftello(tfp -> plain_fp ): -1);
	return total_valid_reads;
}

void * iCache_decompress_chunk_1T(void * arg){
	cache_BCL_t * cache_input = arg;
	while (1){
		int my_bcl_no ;
		subread_lock_occupy(&cache_input -> read_lock);
		for(my_bcl_no =0; my_bcl_no<cache_input -> total_bases_in_each_cluster; my_bcl_no++){
			if(!cache_input -> bcl_no_is_used[my_bcl_no]) {
				cache_input -> bcl_no_is_used[my_bcl_no]=1;
				break;
			}
		}
		subread_lock_release(&cache_input -> read_lock);
		if(my_bcl_no>= cache_input -> total_bases_in_each_cluster) return NULL;

		iCache_continuous_read_lanes( cache_input, my_bcl_no );
	}
	return NULL;
}

int cacheBCL_next_chunk(cache_BCL_t * cache_input){
	int x1;
	//SUBREADprintf("BCL: READ_CHUNK_NEXT ( %d )\n", cache_input -> chunk_no+1);
	cache_input -> chunk_start_lane = cache_input -> chunk_end_lane;
	memset(cache_input -> bcl_no_is_used, 0 , sizeof(int)* MAX_READ_LENGTH);
	pthread_t * threads = malloc(sizeof(pthread_t)*cache_input -> all_threads);
	iCache_continuous_read_lanes( cache_input, -1 ); // read filtering binary

	for(x1=0; x1<cache_input -> all_threads; x1++)
		pthread_create(threads+x1, NULL, iCache_decompress_chunk_1T, cache_input);

	for(x1=0; x1<cache_input -> all_threads; x1++)
		pthread_join(threads[x1],NULL);
	free(threads);
	cache_input -> read_no_in_chunk = 0;
	cache_input -> chunk_no ++;
	return 0;
}


int iCache_copy_read(cache_BCL_t * cache_input, char * read_name, char * seq, char * qual, srInt_64 rno){
	int bii, idx_offset, base_offset;
	int * srii = cache_input -> single_read_lengths;

	#ifdef __MINGW32__
	sprintf(read_name, "R%011I64u:", rno);
	#else
	sprintf(read_name, "R%011llu:", rno);
	#endif
	idx_offset  = srii[0];
	base_offset = srii[1] + idx_offset;

	read_name[13+idx_offset]='|';
	read_name[14+2*idx_offset]='|';
	read_name[15+base_offset+idx_offset]='|';
	sprintf(read_name +16 +2*base_offset, "|@RgLater@L%03d" , cache_input -> lane_no_in_chunk[cache_input -> read_no_in_chunk]);

	for(bii = 0; bii < cache_input -> total_bases_in_each_cluster; bii++){
		int nch = cache_input -> bcl_bin_cache[bii][cache_input -> read_no_in_chunk];
		//if(rno == 1000000){
		//	SUBREADprintf("COPY: %d base NCH=%d\n", bii, nch);
		//}
		if(nch<0) nch+=256;
		int nbase = 'N';
		int nqual = '#';
		if(nch > 0){
			nbase="ACGT"[nch%4];
			nqual=33+(nch>>2);
		}
		if(nqual >= '/' && bii < srii[0] + srii[1]) nqual++;
		if(bii < srii[0]){
			read_name[13+bii] = nbase;
			read_name[14+idx_offset+bii]= nqual;
		}else if(bii < srii[0] + srii[1]){
			read_name[15+idx_offset+bii] = nbase;
			read_name[16+base_offset+bii]= nqual;
		}else{
			seq[bii - srii[0]-srii[1] ] = nbase;
			qual[bii - srii[0]-srii[1] ] = nqual;
		}
	}

	cache_input -> read_no_in_chunk++;
	//SUBREADprintf("GOT READ #%d ; ret=%d\n" , cache_input -> read_no_in_chunk, srii[2]);
	return srii[2];
}



int cacheBCL_next_read(cache_BCL_t * cache_input, char * read_name, char * seq, char * qual, srInt_64 * read_number_in_all){
	if(cache_input -> read_no_in_chunk >= cache_input -> reads_available_in_chunk){
		if(cache_input -> last_chunk_in_cache) return 0;
		cacheBCL_next_chunk(cache_input);
		if(cache_input -> read_no_in_chunk >= cache_input -> reads_available_in_chunk) // no reads are loaded in the previous step
			return 0;
	}

	srInt_64 rnumb =(cache_input -> chunk_no -1)*1ll * cache_input -> reads_per_chunk +(cache_input -> read_no_in_chunk);
	if(read_number_in_all) *read_number_in_all = rnumb;
	return iCache_copy_read(cache_input, read_name, seq, qual, rnumb);
}

int input_BLC_init( input_BLC_t * blc_input , char * data_dir ){
	memset(blc_input, 0, sizeof(input_BLC_t));
	subread_init_lock(&blc_input -> read_lock);

	int rv = iBLC_guess_format_string(data_dir, &blc_input -> total_bases_in_each_cluster, blc_input -> bcl_format_string, blc_input -> filter_format_string, &blc_input -> bcl_is_gzipped, &blc_input -> filter_is_gzipped, blc_input -> single_read_lengths, blc_input -> single_read_is_index);
	if(rv) return -1;

	blc_input -> current_lane = 1;

	return iBLC_open_batch(blc_input)?1:0;
}
// load the next read W/O switch lane. 
// Return 0 if EOF, -1 if error or bases if loaded correctly.
int iBLC_current_lane_next_read(input_BLC_t * blc_input, char * readname , char * read, char * qual){
	int bii, idx_offset, base_offset;

	#ifdef __MINGW32__
	sprintf(readname, "R%011I64u:", blc_input -> read_number +1);
	#else
	sprintf(readname, "R%011llu:", blc_input -> read_number +1);
	#endif

	{
		idx_offset = blc_input -> single_read_lengths[0];
		base_offset = idx_offset + blc_input -> single_read_lengths[1];
	}

	readname[13+idx_offset]='|';
	readname[14+2*idx_offset]='|';
	readname[15+base_offset+idx_offset]='|';
	sprintf(readname +16 +2*base_offset, "|L%03d" , blc_input -> current_lane);

	while(1){
		int fch = blc_input -> filter_is_gzipped? seekgz_next_int8(blc_input -> filter_gzip_fp) :fgetc(blc_input -> filter_fp);
		if(fch < 0) return 0;
		int baseii =0;
		for(bii =0; bii< blc_input -> total_bases_in_each_cluster; bii++){
			int nch = blc_input -> bcl_is_gzipped?seekgz_next_int8(blc_input -> bcl_gzip_fps[bii]):fgetc(blc_input -> bcl_fps[bii]), bv, qv;

			if(fch == 1){
				if(0==nch){
					bv='N';
					qv='#';
				}else{
					bv="ACGT"[nch%4];
					qv=33+(nch>>2);
					if(qv >= '/') qv++;
				}
				if(bii < idx_offset){
					assert(bv !=0 && qv !=0);
					readname[13+ bii]=bv;
					readname[14+ bii+idx_offset]=qv;
				}else if(bii < base_offset){
					assert(bv !=0 && qv !=0);
					readname[15+ bii+ idx_offset]=bv;
					readname[16+ bii+ base_offset]=qv;
				}else{
					read[baseii] = bv;
					qual[baseii] = qv;
					baseii++;
				}
			}
		}

		if(fch==1){
			blc_input -> read_number ++;
			return baseii;
		}
	}
}

int iBLC_inc_lane(input_BLC_t * blc_input){
	blc_input -> current_lane ++;
	return iBLC_open_batch(blc_input); // this function automatically closes BCL fps and FILTER fp.
}

// return : -1: error, 0: end of all files and all lanes, >0: actual read is loaded (return the read len). The read name is the combination of the short-end and the index-end.
// NOTE: this only works with scRNA protocol!!
int input_BLC_next_read(input_BLC_t * blc_input , char * readname, char * read, char * qual){
	int nextlane, rrv=0;
	if(blc_input->is_EOF) return 0;

	subread_lock_occupy(&blc_input -> read_lock);
	for(nextlane = 0; nextlane <2; nextlane++){
		int rv = iBLC_current_lane_next_read(blc_input, readname, read, qual);
		if(rv >0 || rv <0){
			rrv = rv;
			break;
		}
		if(rv ==0 && nextlane){
			rrv = 0;
			break;
		}
		if(nextlane>0){
			rrv = -1;
			break;
		}
		
		rv = iBLC_inc_lane(blc_input);
		if(rv){
			rrv = 0;
			break;
		}
	}
	subread_lock_release(&blc_input -> read_lock);
	return rrv;
}

int input_BLC_tell ( input_BLC_t * blc_input , input_BLC_pos_t * pos ){
	int xx1;
	memset(pos,0, sizeof(*pos));
	pos -> lane_id = blc_input -> current_lane;
	pos -> read_number = blc_input -> read_number;
	pos -> is_EOF = blc_input -> is_EOF;
	if(pos->is_EOF) return 0;
	if(blc_input->bcl_is_gzipped){
		pos -> pos_of_bclgzs = calloc(sizeof(void *) , blc_input -> total_bases_in_each_cluster);
		for(xx1=0; xx1<blc_input -> total_bases_in_each_cluster; xx1++){
			pos -> pos_of_bclgzs[xx1] = malloc(sizeof(seekable_position_t));
			seekgz_tell(blc_input->bcl_gzip_fps[xx1], pos -> pos_of_bclgzs[xx1]);
		}
	}else{
		pos -> pos_of_bcls = calloc(sizeof(srInt_64) , blc_input -> total_bases_in_each_cluster);
		for(xx1=0; xx1<blc_input -> total_bases_in_each_cluster; xx1++)
			pos -> pos_of_bcls[xx1] = ftello(blc_input->bcl_fps[xx1]);
	}


	if(blc_input->filter_is_gzipped){
		pos -> pos_of_filtergz = malloc(sizeof(seekable_position_t));
		seekgz_tell(blc_input->filter_gzip_fp, pos -> pos_of_filtergz);
	}else pos -> pos_of_filter = ftello(blc_input->filter_fp);

	return 0;
}

int input_BLC_seek( input_BLC_t * blc_input , input_BLC_pos_t * pos ){
	int xx1;
	blc_input -> read_number = pos -> read_number;
	if(pos -> is_EOF){
		iBLC_close_batch(blc_input);
		blc_input -> is_EOF = pos -> is_EOF;
		blc_input -> current_lane = pos -> lane_id;
		return 0;
	}

	if(pos -> lane_id != blc_input -> current_lane){
		blc_input -> current_lane = pos -> lane_id;
		iBLC_open_batch(blc_input);
	}

	for(xx1=0; xx1<blc_input -> total_bases_in_each_cluster; xx1++)
		if(blc_input->bcl_is_gzipped) seekgz_seek(blc_input->bcl_gzip_fps[xx1], pos -> pos_of_bclgzs[xx1]); 
		else fseeko(blc_input->bcl_fps[xx1], pos -> pos_of_bcls[xx1], SEEK_SET);

	if(blc_input->filter_is_gzipped) seekgz_seek(blc_input->filter_gzip_fp, pos -> pos_of_filtergz);
	else fseeko(blc_input->filter_fp, pos -> pos_of_filter, SEEK_SET);
	
	return 0;
}

void input_BLC_destroy_pos(input_BLC_t * blc_input , input_BLC_pos_t *pos){
	int xx1;
	for(xx1=0; xx1<blc_input -> total_bases_in_each_cluster; xx1++){
		if(blc_input->bcl_is_gzipped) free(pos -> pos_of_bclgzs[xx1]);
	}
	free((blc_input->bcl_is_gzipped? (void*)pos -> pos_of_bclgzs:(void*)pos -> pos_of_bcls));
}

void input_BLC_close(input_BLC_t * blc_input){
	iBLC_close_batch(blc_input);
	subread_destroy_lock(&blc_input -> read_lock);
}

void iBLC_free_sample_items(void * sample_arr){
	ArrayList * ar = (ArrayList * ) sample_arr;
	ArrayListDestroy(ar);
}

void iBLC_free_3tp(void * t){
	char **tt = (char**)t;
	free(tt[1]);
	free(tt);
}

#define SHEET_FORMAT_RAWDIR_INPUT 10
#define SHEET_FORMAT_FASTQ_INPUT 20
#define SHEET_FORMAT_BAM_INPUT 30
HashTable * input_BLC_parse_SampleSheet(char * fname){
	HashTable * ret = StringTableCreate(30);
	HashTableSetDeallocationFunctions(ret, free, iBLC_free_sample_items);
	FILE * fp = fopen(fname, "rb");
	if(fp==NULL) return NULL;
	char linebuf[MAX_FILE_NAME_LENGTH];
	int state = -1, file_format=-1, stat1_line = 0;
	while(!feof(fp)){
		char * gret = fgets(linebuf, MAX_FILE_NAME_LENGTH-1, fp);
		if(gret == NULL) break;
		int linelen = strlen(linebuf);
		if(linelen<5)continue;
		if(linebuf[linelen -1] == '\r' || linebuf[linelen -1] == '\n') linebuf[linelen -1] =0;
		if(linebuf[linelen -2] == '\r' || linebuf[linelen -2] == '\n') linebuf[linelen -2] =0;
		
		if(state < 0 && strstr(linebuf,"EMFileVersion,4")) state = 0;
		if(state == 1 && linebuf[0]=='[') state = 99999;
		if(state == 1){
			if(0==stat1_line){
				if(strstr( linebuf, "Lane")){
					file_format = SHEET_FORMAT_RAWDIR_INPUT;
					continue;
				}

				if(strstr( linebuf, "BAMFile")){
					file_format = SHEET_FORMAT_BAM_INPUT;
					continue;
				}

				if(strstr( linebuf, "BarcodeUMIFile")){
					file_format = SHEET_FORMAT_FASTQ_INPUT;
					continue;
				}
			}
			if(strlen(linebuf)>10)stat1_line ++;

			char * tokp=NULL;
			int lane_no = 1;
			char * sample_name = NULL;
			char * sample_index = NULL;

			if(file_format == SHEET_FORMAT_RAWDIR_INPUT){
				lane_no = atoi(strtok_r(linebuf, ",", &tokp));
				strtok_r(NULL, ",", &tokp);
				sample_name = strtok_r(NULL, ",", &tokp);
				sample_index = strdup(strtok_r(NULL, ",", &tokp));
			}else if(file_format == SHEET_FORMAT_FASTQ_INPUT){
				strtok_r(linebuf, ",", &tokp); // fastq file 1
				strtok_r(NULL, ",", &tokp); // fastq file 2
				sample_name = strtok_r(NULL, ",", &tokp); // sample name
			}else{
				strtok_r(linebuf, ",", &tokp); // bam file
				sample_name = strtok_r(NULL, ",", &tokp); // sample name
			}

			char ** entry = malloc(sizeof(void*)*3);
			entry[0] = NULL + lane_no;
			entry[1] = sample_index;
			entry[2] = NULL + stat1_line; // the N-th line in the file , 1 based

			ArrayList * arr = HashTableGet(ret, sample_name);
			if(NULL == arr){
				arr = ArrayListCreate(16);
				ArrayListSetDeallocationFunction(arr, iBLC_free_3tp);
				HashTablePut( ret, strdup(sample_name), arr );
			}
			//SUBREADprintf("PUT_SAMPLE=%s of lane_no %d\n", sample_name, lane_no);
			ArrayListPush(arr,entry);
		}
		if(state == 0 && strstr(linebuf,"ata]")){
			state = 1;
			stat1_line = 0;
		}
	}
	if(state <1){
		SUBREADprintf("ERROR: the sample sheet doesn't contain any sample.\n");
		return NULL;
	}
	return ret;
}

ArrayList * input_BLC_parse_CellBarcodes(char * fname){
	autozip_fp fp;
	int resop = autozip_open(fname, &fp);
	if(resop<0) return NULL;

	ArrayList * ret = ArrayListCreate(10000000);
	ArrayListSetDeallocationFunction(ret, free);

	while(1){
		char tmp_fl[MAX_BARCODE_LEN+1];
		int skr = autozip_gets(&fp, tmp_fl, MAX_BARCODE_LEN);
		if(skr<1) break;
		//if(ret -> numOfElements <40) SUBREADprintf("LOAD_CELL_BAR : %s , %d\n", tmp_fl, skr);
		int x1;
		if(tmp_fl[skr-1]=='\n') tmp_fl[skr-1]=0;
		for(x1=0; tmp_fl[x1]; x1++) if(!isalpha(tmp_fl[x1])){
			tmp_fl[x1]=0;
			break;
		}
		ArrayListPush(ret, strdup(tmp_fl));
	}

	autozip_close(&fp);
	return ret;
}

int is_ATGC(char c){
	return c=='A'||c=='C'||c=='G'||c=='T'||c=='N';
}

int hamming_dist_ATGC_max2(char* s1, char* s2 ){
	int xx,ret=0;
	for(xx=0;;xx++){
		char nch1 = s1[xx];
		char nch2 = s2[xx];
		if(is_ATGC(nch1) && is_ATGC(nch2)){
			ret += nch1==nch2;
			if(xx -ret >2) return 999;
		}else break;
	}
	return xx-ret;
}

int iCache_get_sample_id(ArrayList * sample_sheet_list, char * sbc, int read_laneno){
	int x1;

	//SUBREADprintf("TOTAL_SBC=%ld\n", global_context -> scRNA_sample_barcode_list -> numOfElements);
	for(x1=0; x1 < sample_sheet_list -> numOfElements ; x1++ ){
		char ** lane_and_barcode = ArrayListGet(sample_sheet_list, x1);
		int lane_no = lane_and_barcode[0]-(char*)NULL;
		if(read_laneno == lane_no){
			int sample_no = lane_and_barcode[1]-(char*)NULL;
			char * knownbar = lane_and_barcode[2];
			int hd = hamming_dist_ATGC_max2( sbc, knownbar );
			if(hd<=2) return sample_no;
		}
	}
	return -1;
}

void iCache_copy_sample_table_2_list(void* ky, void* va, HashTable* tab){
	ArrayList * cbclist = tab -> appendix1;
	ArrayList * hashed_arr = va ;

	srInt_64 xx1;
	for(xx1 =0; xx1< hashed_arr -> numOfElements; xx1++){
		char ** push_arr = malloc(sizeof(char*)*3);
		char ** sbc_lane_sample = ArrayListGet(hashed_arr, xx1);
		srInt_64 lane_sample_int = sbc_lane_sample[0]-(char*)NULL;

		ArrayListPush(cbclist, push_arr);
		push_arr[0] = NULL + lane_sample_int; 
		push_arr[1] = NULL + cbclist -> numOfElements+1;
		push_arr[2] = sbc_lane_sample[1]; // Sample Barcode
	}

}

#define IMPOSSIBLE_MEMORY_ADDRESS 0x5CAFEBABE0000000

int iCache_get_cell_no(HashTable * cell_barcode_table, ArrayList * cell_barcode_list, char * cell_barcode, int cell_barcode_length){
	char tmpc [MAX_READ_NAME_LEN];
	int xx1;
	ArrayList * ret=NULL;

	for(xx1=0;xx1<3;xx1++){
		int xx2;
		if(xx1==1) ret = ArrayListCreate(100);

		if(xx1>0){
			tmpc[0] = (xx1==2)?'S':'F';
			for(xx2=0; xx2<cell_barcode_length/2 ; xx2++)
				tmpc[1+xx2] = cell_barcode[2*xx2+xx1-1];
			tmpc[1+cell_barcode_length/2]=0;
		}else{
			memcpy(tmpc, cell_barcode, cell_barcode_length);
			tmpc[cell_barcode_length]=0;
		}

		void *xrawarr = HashTableGet(cell_barcode_table, tmpc);

		if(xx1 == 0){
			//if(xrawarr) SUBREADprintf("CAFE ? %p\n", xrawarr);
			srInt_64 xint = xrawarr - NULL;
			if(( xint & 0xFFFFFFFFF0000000llu)== IMPOSSIBLE_MEMORY_ADDRESS){
				int only_cell_id = xint - IMPOSSIBLE_MEMORY_ADDRESS;
				// no memory was allocated.
				return only_cell_id;
			}
		}else{
			ArrayList * rawarr = xrawarr;
			if(rawarr){
				int xx3,xx2, found;
				for(xx2=0; xx2<rawarr->numOfElements; xx2++){
					int bcno = ArrayListGet(rawarr, xx2)-NULL;
					found=0;
					for(xx3=0;xx3<ret -> numOfElements;xx3++){
						if(ArrayListGet(ret, xx3)==NULL+bcno){
							found=1;
							break;
						}
					}

					if(!found)ArrayListPush(ret, NULL+bcno);
				}
			}
		}
	}


	int tb1=-1;
	for(xx1=0; xx1<ret -> numOfElements; xx1++){
		int tbcn = ArrayListGet(ret,xx1)-NULL;
		char * known_cbc = ArrayListGet(cell_barcode_list, tbcn);
		int hc = hamming_dist_ATGC_max2( known_cbc, cell_barcode );

	//	cbc[16]=0; if(hc <=3)SUBREADprintf("TEST_CBC %s ~ %s = %d\n", known_cbc, cbc, hc);
		if(hc==1){
			tb1 = tbcn;
			break;
		}
	}
	//SUBREADprintf("CANDIDATE CELL BARCODES=%ld ; hit = %d\n", ret->numOfElements, tb1);
	ArrayListDestroy(ret);

	return tb1;

}

void iCache_delete_bcb_key(void *k){
	if(((k-NULL) & 0xFFFFFFFFF0000000llu)!= IMPOSSIBLE_MEMORY_ADDRESS) ArrayListDestroy(k);
}

int cacheBCL_quality_test(int input_mode, char * datadir, HashTable * sample_sheet_table, ArrayList * cell_barcode_list, int testing_reads, int * tested_reads, int * valid_sample_index, int * valid_cell_barcode);
#ifdef MAKE_TEST_QUALITY_TEST
int main(int argc,char ** argv){
#else
int do_R_try_cell_barcode_files(int argc, char ** argv){
#endif
	assert(argc>=6);
	int input_mode = GENE_INPUT_BCL;
	if(strcmp("fastq", argv[5])==0) input_mode = GENE_INPUT_SCRNA_FASTQ;
	if(strcmp("bam", argv[5])==0) input_mode = GENE_INPUT_SCRNA_BAM;
	int testing_reads = atoi(argv[4]);
	ArrayList * cell_barcode_list = input_BLC_parse_CellBarcodes(argv[3]);
	SUBREADprintf("Loaded %lld cell barcodes from %s\n", cell_barcode_list -> numOfElements, argv[3]);

	HashTable * sample_sheet = NULL;
	if(input_mode==GENE_INPUT_BCL) input_BLC_parse_SampleSheet(argv[2]);
	char * input_BCL_raw_dir = argv[1];

	int tested_reads=0, valid_sample_index=0, valid_cell_barcode=0;
	int rv = cacheBCL_quality_test(input_mode, input_BCL_raw_dir, sample_sheet, cell_barcode_list, testing_reads, &tested_reads, &valid_sample_index, &valid_cell_barcode );

#ifdef MAKE_TEST_QUALITY_TEST
	printf("Samples = %ld loaded\n", sample_sheet -> numOfElements);
	printf("Cell barcodes = %ld loaded\n", cell_barcode_list -> numOfElements);
	printf("QUALITY-TEST : rv=%d , tested = %d , good-sample = %d , good-cell = %d\n", rv, tested_reads, valid_sample_index, valid_cell_barcode);
#else
	argv[6]=NULL+rv;
	argv[7]=NULL+tested_reads;
	argv[8]=NULL+valid_sample_index;
	argv[9]=NULL+valid_cell_barcode;
#endif
	ArrayListDestroy(cell_barcode_list);
	return 0;
}

int cacheBCL_qualTest_BAMmode(char * datadir, int testing_reads, int known_cell_barcode_length, ArrayList * sample_sheet_list, ArrayList * cell_barcode_list, HashTable * cell_barcode_table, int * tested_reads, int * valid_sample_index, int * valid_cell_barcode){
	input_scBAM_t * scBAM_input = malloc(sizeof(input_scBAM_t));
	int ret = input_scBAM_init(scBAM_input, datadir);
	SUBREADprintf("cacheBCL_qualTest_RET_BAMmode = %d for %s\n", ret, datadir);
	if(ret)return ret;
	while(1){
		char base[MAX_READ_LENGTH], qual[MAX_READ_LENGTH], rname[MAX_READ_NAME_LEN];
		base[0]=qual[0]=rname[0]=0;
		ret = scBAM_next_read(scBAM_input, rname , base, qual);
		if(ret<=0)break;
		char *cell_barcode = NULL, *sample_barcode=NULL, *lane_str=NULL;

		int xx=0;
		char *testi;
		for(testi = rname+1; * testi; testi ++){
			if( * testi=='|'){
				xx++;
				if(xx == 1) {
					cell_barcode = testi +1;
				}else if(xx == 2) {
					sample_barcode = testi +1;
				}else if(xx == 4){
					lane_str = testi+1;
					break;
				}
			}
		}

		int cell_no = iCache_get_cell_no(cell_barcode_table, cell_barcode_list, cell_barcode, known_cell_barcode_length);
		if(cell_no>0) (*valid_cell_barcode)++;

		(*tested_reads)++;
		if((*tested_reads) >= testing_reads)break;
	}
	input_scBAM_close(scBAM_input);
	free(scBAM_input);
	return 0;
}

int cacheBCL_qualTest_FQmode(char * datadir, int testing_reads, int known_cell_barcode_length, ArrayList * sample_sheet_list, ArrayList * cell_barcode_list, HashTable * cell_barcode_table, int * tested_reads, int * valid_sample_index, int * valid_cell_barcode){
	input_mFQ_t fqs_input;
	int ret = input_mFQ_init_by_one_string(&fqs_input, datadir);
	if(ret)return ret;
	while(1){
		char base[MAX_READ_LENGTH], qual[MAX_READ_LENGTH], rname[MAX_READ_NAME_LEN];
		base[0]=qual[0]=rname[0]=0;
		ret = input_mFQ_next_read(&fqs_input, rname , base, qual);
		if(ret<=0)break;

		char *cell_barcode = NULL;
		//SUBREADprintf("RETV=%d, RN=%s\n" ,ret, rname);
		int xx=0;
		char *testi;
		for(testi = rname+1; * testi; testi ++){
			if( * testi=='|'){
				xx++;
				if(xx == 1) {
					cell_barcode = testi+1;
					break;
				}
			}
		}

		int cell_no = iCache_get_cell_no(cell_barcode_table, cell_barcode_list, cell_barcode, known_cell_barcode_length);
		if(cell_no>0) (*valid_cell_barcode)++;

		(*tested_reads)++;
		if((*tested_reads) >= testing_reads)break;
	}
	input_mFQ_close(&fqs_input);
	return 0;
}


int cacheBCL_qualTest_BCLmode(char * datadir, int testing_reads, int known_cell_barcode_length, ArrayList * sample_sheet_list, ArrayList * cell_barcode_list, HashTable * cell_barcode_table, int * tested_reads, int * valid_sample_index, int * valid_cell_barcode){
	cache_BCL_t blc_input;
	int orv = cacheBCL_init(&blc_input, datadir,testing_reads+1, 1);
	if(orv)return -1;
	while(1){
		char base[MAX_READ_LENGTH], qual[MAX_READ_LENGTH], rname[MAX_READ_NAME_LEN];
		srInt_64 readno = 0;
		base[0]=qual[0]=rname[0]=0;
		orv = cacheBCL_next_read(&blc_input, rname, base, qual, &readno);
		if(0==orv) break;

		char * testi, *lane_str = NULL, * sample_barcode = NULL, * cell_barcode, * umi_barcode; // cell_barcode MUST be 16-bp long, see https://community.10xgenomics.com/t5/Data-Sharing/Cell-barcode-and-UMI-with-linked-reads/td-p/68376
		cell_barcode = rname + 13;
		umi_barcode = cell_barcode + known_cell_barcode_length;
		int xx=0, laneno=0;
		for(testi = umi_barcode+1; * testi; testi ++){
			if( * testi=='|'){
				xx++;
				if(xx == 2) {
					sample_barcode = testi +1;
				}else if(xx == 4){
					lane_str = testi+1;
					break;
				}
			}
		}
		assert(xx ==4 && (*lane_str)=='L');
		for(testi = lane_str+1; *testi; testi++){
			assert(isdigit(*testi));
			laneno = laneno*10 + (*testi)-'0';
		}

		int sample_no = -1;
		if(sample_sheet_list->numOfElements)iCache_get_sample_id(sample_sheet_list, sample_barcode, laneno);
		int cell_no = iCache_get_cell_no(cell_barcode_table, cell_barcode_list, cell_barcode, known_cell_barcode_length);
		//printf("CELL_CALL %d = %s\n", cell_no, cell_barcode);
		if(sample_no>0) (*valid_sample_index)++;
		if(cell_no>0) (*valid_cell_barcode)++;

		(*tested_reads)++;
		if((*tested_reads) >= testing_reads)break;
	}

	cacheBCL_close(&blc_input);
	return 0;
}

int cacheBCL_quality_test(int input_mode, char * datadir, HashTable * sample_sheet_table, ArrayList * cell_barcode_list, int testing_reads, int * tested_reads, int * valid_sample_index, int * valid_cell_barcode){
	ArrayList * sample_sheet_list = ArrayListCreate(100);
	ArrayListSetDeallocationFunction(sample_sheet_list, free);
	if(sample_sheet_table){
		sample_sheet_table -> appendix1 = sample_sheet_list;
		HashTableIteration(sample_sheet_table, iCache_copy_sample_table_2_list);
	}

	HashTable * cell_barcode_table = StringTableCreate(1000000);
	HashTableSetDeallocationFunctions(cell_barcode_table, free, iCache_delete_bcb_key);

	char bctmp[MAX_READ_NAME_LEN];
	int xx1,xx2, known_cell_barcode_length=-1;
	for(xx1=0; xx1< cell_barcode_list -> numOfElements ; xx1++){
		char * bcbstr = ArrayListGet(cell_barcode_list, xx1);
		if(-1==known_cell_barcode_length) known_cell_barcode_length = strlen(bcbstr);
		else if(known_cell_barcode_length != strlen(bcbstr)){
			SUBREADprintf("ERROR: the cell barcodes have different lengths (%d!=%ld at %d). The program cannot process the cell barcodes.\n", known_cell_barcode_length, strlen(bcbstr),xx1);
			return -1;
		}
		HashTablePut(cell_barcode_table, strdup(bcbstr), NULL+IMPOSSIBLE_MEMORY_ADDRESS+xx1);
		//if(cell_barcode_table -> numOfElements % 1000000 < 10)printf("size=%ld\n", cell_barcode_table -> numOfElements);

		for(xx2=0; xx2<2; xx2++){
			bctmp[0] = xx2?'S':'F';
			int xx3;
			for(xx3 = 0; xx3< known_cell_barcode_length/2; xx3++)
				bctmp[xx3+1] = bcbstr[ xx3*2+xx2 ];
			bctmp[known_cell_barcode_length/2+1]=0;

			ArrayList * array_of_codes = HashTableGet(cell_barcode_table, bctmp);
			if(!array_of_codes){
				array_of_codes = ArrayListCreate(4);
				HashTablePut(cell_barcode_table, strdup(bctmp), array_of_codes);
			}
			ArrayListPush(array_of_codes, NULL+xx1);
		}
	}
	if(known_cell_barcode_length<0){
		SUBREADprintf("ERROR: cannot load any cell barcode from database\n");
		return -1;
	}

	int ret = -1;
	if(input_mode == GENE_INPUT_SCRNA_FASTQ)
		ret=cacheBCL_qualTest_FQmode(datadir, testing_reads, known_cell_barcode_length, sample_sheet_list, cell_barcode_list, cell_barcode_table,  tested_reads,  valid_sample_index,  valid_cell_barcode);
	else if(input_mode == GENE_INPUT_SCRNA_BAM)
		ret=cacheBCL_qualTest_BAMmode(datadir, testing_reads, known_cell_barcode_length, sample_sheet_list, cell_barcode_list, cell_barcode_table,  tested_reads,  valid_sample_index,  valid_cell_barcode);
	else
		ret=cacheBCL_qualTest_BCLmode(datadir, testing_reads, known_cell_barcode_length, sample_sheet_list, cell_barcode_list, cell_barcode_table,  tested_reads,  valid_sample_index,  valid_cell_barcode);
	ArrayListDestroy(sample_sheet_list);
	HashTableDestroy(cell_barcode_table);
	return ret;
}

#ifdef MAKE_TEST_SAMPLESHEET
int main(int argc, char ** argv){
	assert(argc>1);
	char * samplesheet = argv[1];
	HashTable * st = input_BLC_parse_SampleSheet(samplesheet);
	printf("P=%p\n", st);
}
#endif
#ifdef MAKE_TEST_ICACHE
int main(int argc, char ** argv){
	assert(argc>1);
	cache_BCL_t blc_input;
	int orv = cacheBCL_init(&blc_input, argv[1],5000000, 5), total_poses = 0;
	printf("orv=%d, bases=%d, filter_gzip=%d, data_gzip=%d\nBCL-pattern = %s\nFLT-pattern=%s\n", orv, blc_input.total_bases_in_each_cluster, blc_input.filter_is_gzipped, blc_input.bcl_is_gzipped, blc_input.bcl_format_string, blc_input.filter_format_string);

	while(1){
		char base[1000], qual[1000], rname[200];
		srInt_64 readno = 0;
		base[0]=qual[0]=rname[0]=0;
		orv = cacheBCL_next_read(&blc_input, rname, base, qual, &readno);
		assert(orv>=0);
		if(0==orv) break;
		if(1||readno%1000000==0)printf("%s %s %s\n", rname, base, qual);
	}

	cacheBCL_close(&blc_input);
}
#endif

#ifdef MAKE_TEST_IBLC
int main(int argc, char ** argv){
	assert(argc>1);
	input_BLC_t blc_input;
	input_BLC_pos_t poses[12];
	int orv = input_BLC_init(&blc_input, argv[1]), total_poses = 0;
	printf("orv=%d, bases=%d, filter_gzip=%d, data_gzip=%d\n", orv, blc_input.total_bases_in_each_cluster, blc_input.filter_is_gzipped, blc_input.bcl_is_gzipped);

	while(1){
		char base[1000], qual[1000], rname[200];
		base[0]=qual[0]=rname[0]=0;
		orv = input_BLC_next_read(&blc_input, rname, base, qual);
		assert(orv>=0);
		if(0==orv) break;
		if(blc_input.read_number%1000000==1)printf("%s %s %s\n", rname, base, qual);
		if(blc_input.read_number % 20000000 == 1){
			input_BLC_tell(&blc_input, poses + total_poses);
			total_poses++;
		}
		if(blc_input.read_number > 180000010) break;
	}

	int ii;
	for(ii = 0; ii < total_poses; ii++){
		srInt_64 jj;
		printf("\n\n========== SEEKING %d OF %d ================\n", ii,total_poses);
		input_BLC_seek( &blc_input, poses+ii );
		input_BLC_destroy_pos( &blc_input, poses+ii );
		for(jj = 0 ; jj < 3000011; jj++){
			char base[1000], qual[1000], rname[200];
			orv = input_BLC_next_read(&blc_input, rname, base, qual);
			assert(orv>=0);
			if(blc_input.read_number%1000000==2)printf("%s %s %s\n", rname, base, qual);
		}
	} 
	
	printf("END CORRECTLY  R=%llu\n", blc_input.read_number);
	input_BLC_close(&blc_input);
	return 0;
}
#endif

void scBAM_tell(input_scBAM_t * bam_input, input_scBAM_pos_t * pos){
	pos -> current_BAM_file_no = bam_input -> current_BAM_file_no;
	pos -> section_start_pos = bam_input -> section_start_pos;
	pos -> current_read_no = bam_input -> current_read_no;
	pos -> in_section_offset = bam_input -> in_section_offset;
}

int scBAM_next_char(input_scBAM_t * bam_input);
int scBAM_next_string(input_scBAM_t * bam_input, char * strbuff, int lenstr);

int scBAM_next_int(input_scBAM_t * bam_input, int *ret){
	*ret = 0;
	int xk1;
	for(xk1=0;xk1<4;xk1++){
		int nbyte = scBAM_next_char(bam_input);
		if(nbyte<0) return -1;
		(*ret) += nbyte <<(8*xk1);
	}
	return 0;
}

int scBAM_skip_bam_header(input_scBAM_t * bam_input){
	int ret = 0, tmpi = 0, nref=0, x1;
	ret = scBAM_next_int(bam_input, &tmpi);
	if(ret < 0 || tmpi != 0x014d4142)return -1; // 0x014d4142 = 'BAM\1' little-endian
	scBAM_next_int(bam_input, &tmpi); // header txt length
	while(tmpi--) scBAM_next_char(bam_input); // skip header txt
	scBAM_next_int(bam_input, &nref); // n_ref
	bam_input -> chro_table = calloc(sizeof(SamBam_Reference_Info), nref);
	SUBREADprintf("OPEN '%s' : %d refs\n",bam_input -> BAM_file_names[bam_input -> current_BAM_file_no], nref);
	for(x1=0;x1<nref;x1++){
		scBAM_next_int(bam_input, &tmpi); // l_name
		scBAM_next_string(bam_input, bam_input -> chro_table [x1].chro_name, tmpi);
		ret = scBAM_next_int(bam_input, (int*)&bam_input -> chro_table [x1].chro_length);
		if(ret < 0) return -1;
	}
	return 0;
}

int scBAM_inner_fopen(input_scBAM_t * bam_input){
	bam_input -> os_file = f_subr_open(bam_input -> BAM_file_names[bam_input -> current_BAM_file_no], "rb");
	if(!bam_input -> os_file)return -1;
	return scBAM_skip_bam_header(bam_input);
}

void scBAM_inner_fclose(input_scBAM_t * bam_input){
	free(bam_input -> chro_table);
	fclose(bam_input -> os_file);
}

int scBAM_rebuffer(input_scBAM_t * bam_input){
	int bin_len=0;
	while(1){
		if(bam_input -> current_BAM_file_no== bam_input -> total_BAM_files) return -1;
		if(feof(bam_input -> os_file)){
			scBAM_inner_fclose(bam_input);
			bam_input -> current_BAM_file_no++;
			if(bam_input -> current_BAM_file_no== bam_input -> total_BAM_files){
				return -1;
			}
			scBAM_inner_fopen(bam_input);
		}
		char zipped_bam_buf[66000];
		bam_input -> section_start_pos = ftello(bam_input -> os_file);
		int ziplen = PBam_get_next_zchunk(bam_input -> os_file, zipped_bam_buf, 66000, (unsigned int*)&bin_len);
		if(ziplen<1) return -1;
		bin_len = SamBam_unzip(bam_input -> section_buff, 65536, zipped_bam_buf, ziplen, 0);
		if(bin_len>0){
			bam_input -> section_bin_bytes = bin_len;
			bam_input -> in_section_offset = 0;
			break;
		}else if(bin_len<0) return -1;
	}
	return bin_len;
}

// negative : EOF
int scBAM_next_char(input_scBAM_t * bam_input){
	if(bam_input -> current_BAM_file_no== bam_input -> total_BAM_files) return -1;
	if(bam_input -> in_section_offset == bam_input -> section_bin_bytes)
		if(0>scBAM_rebuffer(bam_input))return -1;
	
	int ret = bam_input -> section_buff[bam_input -> in_section_offset++];
	if(ret<0)ret+=256;
	return ret;
}


void scBAM_seek(input_scBAM_t * bam_input, input_scBAM_pos_t * pos){
	if(bam_input -> current_BAM_file_no!= pos -> current_BAM_file_no){
		if(bam_input -> current_BAM_file_no < bam_input -> total_BAM_files) scBAM_inner_fclose(bam_input);
		bam_input -> current_BAM_file_no = pos -> current_BAM_file_no;
		if(bam_input -> current_BAM_file_no < bam_input -> total_BAM_files) scBAM_inner_fopen(bam_input);
	}
	if(bam_input -> current_BAM_file_no < bam_input -> total_BAM_files){
		bam_input -> section_start_pos = pos -> section_start_pos;
		fseeko(bam_input -> os_file, bam_input -> section_start_pos, SEEK_SET);
		scBAM_rebuffer(bam_input);
		bam_input -> current_read_no = pos -> current_read_no;
		bam_input -> in_section_offset = pos -> in_section_offset;
	}
}

int scBAM_next_string(input_scBAM_t * bam_input, char * strbuff, int lenstr){
	int x1=0;
	while(lenstr--){
		int nbyte = scBAM_next_char(bam_input);
		if(nbyte<0) return -1;
		strbuff[x1++]=nbyte;
	}
	return x1;
}


// binbuf >= FC_LONG_READ_RECORD_HARDLIMIT 
// the alignment block with binlen will be written into binbuf.
// negative : EOF
int scBAM_next_alignment_bin(input_scBAM_t * bam_input, char * binbuf){
	int reclen = 0;
	int ret = scBAM_next_int(bam_input, &reclen);
	if(ret<0)return -1;
	if(reclen < 36 || reclen > FC_LONG_READ_RECORD_HARDLIMIT -4)return -1;
	memcpy(binbuf, &reclen, 4);
	return scBAM_next_string(bam_input, binbuf +4, reclen);
}

int scBAM_next_read(input_scBAM_t * bam_input, char * read_name, char * seq, char * qual){
	int ret = scBAM_next_alignment_bin(bam_input, bam_input -> align_buff);
	if(ret<0)return -1;
	int binlen=0, l_read_name = 0, flag=0, n_cigar_op=0, l_seq=0, x1;
	memcpy(&binlen, bam_input -> align_buff, 4);
	memcpy(&l_read_name, bam_input -> align_buff+12,1);
	memcpy(&n_cigar_op, bam_input -> align_buff+16,2);
	memcpy(&flag, bam_input -> align_buff+18,2);
	memcpy(&l_seq, bam_input -> align_buff+20,4);
	int rname_build_ptr = l_read_name-1;
	strcpy(read_name, bam_input -> align_buff+36);
	char * seq_start = bam_input -> align_buff+36+l_read_name+4*n_cigar_op;
	for(x1=0;x1<l_seq;x1++)
		seq[x1]="=ACMGRSVTWYHKDBN"[(seq_start[x1/2] >> ( 4*!(x1%2) ) ) & 15];

	char * qual_start = seq_start + (l_seq+1)/2;
	memcpy(qual, qual_start, l_seq);
	for(x1=0;x1<l_seq;x1++) qual[x1]+=33;
	if(flag & 16){
		reverse_quality(qual, l_seq);
		reverse_read(seq, l_seq, GENE_SPACE_BASE);
	}
	qual[l_seq] = 0;
	char * extag_start = qual_start + l_seq;
	char * tag_str_val=NULL, tag_type = 0;

	for(x1=0;x1<5;x1++){
		char * tagname = NULL;
		if(x1==0)tagname ="CR";
		if(x1==1)tagname ="UR";

		if(x1==2)tagname ="CY";
		if(x1==3)tagname ="UY";

		if(x1==4)tagname ="RG";

		tag_type = 0;
		SAM_pairer_iterate_tags((unsigned char*)extag_start, bam_input -> align_buff + 4 + binlen - extag_start, tagname, &tag_type, &tag_str_val);
		if(tag_type!='Z')return -1;
		int tag_str_len = strlen(tag_str_val);
		if(x1==0||x1==2||x1==4)read_name[rname_build_ptr++] = '|';
		memcpy( read_name +  rname_build_ptr, tag_str_val, tag_str_len );
		rname_build_ptr += tag_str_len;
	}

	read_name[rname_build_ptr] = 0;
	return l_seq;
}

int input_scBAM_init(input_scBAM_t * bam_input, char * bam_fnames){
	char * ubamnames = strdup(bam_fnames);
	char * tpl1 = NULL;
	char * fnl1 = strtokmm(ubamnames, SCRNA_FASTA_SPLIT1, &tpl1);
	int no_file = 0;

	memset(bam_input, 0, sizeof(input_scBAM_t));
	while(fnl1){
		bam_input -> BAM_file_names[no_file++] = strdup(fnl1);
		fnl1 = strtokmm(NULL, SCRNA_FASTA_SPLIT1, &tpl1);
	}
	bam_input -> total_BAM_files = no_file;
	free(ubamnames);
	return scBAM_inner_fopen(bam_input);
}

void input_scBAM_close(input_scBAM_t * bam_input){
	int x1;
	for(x1=0; x1<bam_input -> total_BAM_files; x1++)free(bam_input ->BAM_file_names[x1]);
	if(bam_input -> current_BAM_file_no < bam_input -> total_BAM_files)scBAM_inner_fclose(bam_input);
}

void input_mFQ_fp_close(input_mFQ_t * fqs_input){
	if(fqs_input -> autofp1.filename[0]){
		autozip_close(&fqs_input -> autofp1);
		if(fqs_input->files2)autozip_close(&fqs_input -> autofp2);
		autozip_close(&fqs_input -> autofp3);
	}
	fqs_input -> autofp1.filename[0] = 0;
}

int input_mFQ_guess_lane_no(char * f1_name){
	int fnamelen = strlen(f1_name);
	char * lst = strstr(f1_name + fnamelen - 24, "_L0");
	if(lst){
		int lno = 0;
		lno += lst[3]-'0';
		lno = lno*10+( lst[4]-'0' );
	}
	return 999;
}

int input_mFQ_open_files(input_mFQ_t * fqs_input){
	fqs_input -> current_guessed_lane_no = input_mFQ_guess_lane_no(fqs_input->files1[fqs_input-> current_file_no]);
	int gzipped_ret = autozip_open(fqs_input->files1[fqs_input-> current_file_no],&fqs_input -> autofp1);
	if(fqs_input -> files2)gzipped_ret = gzipped_ret < 0 ?gzipped_ret:autozip_open(fqs_input->files2[fqs_input-> current_file_no],&fqs_input -> autofp2);
	gzipped_ret = gzipped_ret < 0?gzipped_ret:autozip_open(fqs_input->files3[fqs_input-> current_file_no],&fqs_input -> autofp3);
	return gzipped_ret < 0;
}

int input_mFQ_next_file(input_mFQ_t * fqs_input){
	input_mFQ_fp_close(fqs_input);
	fqs_input-> current_file_no ++;

	if(fqs_input-> current_file_no >= fqs_input-> total_files)return -1;
	return input_mFQ_open_files(fqs_input);
}

int input_mFQ_init_by_one_string(input_mFQ_t * fqs_input, char * three_paired_fqnames){
	int total_files = 0;
	char * fnames = strdup(three_paired_fqnames);
	char ** files1 = malloc(sizeof(char*) * MAX_SCRNA_FASTQ_FILES);
	char ** files2 = malloc(sizeof(char*) * MAX_SCRNA_FASTQ_FILES);
	char ** files3 = malloc(sizeof(char*) * MAX_SCRNA_FASTQ_FILES);

	char * tpl1 = NULL, * tpl2 = NULL;
	char * fnl1 = strtokmm(fnames, SCRNA_FASTA_SPLIT1, &tpl1);
	int no_file2 = 0;
	while(fnl1){
		char * fnl2 = strtokmm(fnl1, SCRNA_FASTA_SPLIT2, &tpl2);
		files1[total_files] = fnl2;
		fnl2 = strtokmm(NULL, SCRNA_FASTA_SPLIT2, &tpl2);
		files2[total_files] = fnl2;
		no_file2 = no_file2 || strlen(fnl2)<2;
		fnl2 = strtokmm(NULL, SCRNA_FASTA_SPLIT2, &tpl2);
		files3[total_files] = fnl2;
		fnl1 = strtokmm(NULL, SCRNA_FASTA_SPLIT1, &tpl1);
		total_files++;
	}

	int rv = input_mFQ_init(fqs_input, files1, no_file2?NULL:files2, files3, total_files);
	free(fnames);
	free(files1);
	free(files2);
	free(files3);
	return rv;

}

int input_mFQ_init( input_mFQ_t * fqs_input, char ** files1, char ** files2, char** files3, int total_files ){
	int x1;
	memset(fqs_input, 0, sizeof(input_mFQ_t));
	fqs_input->files1 = malloc(sizeof(char*)*total_files);
	fqs_input->files2 = files2?malloc(sizeof(char*)*total_files):NULL;
	fqs_input->files3 = malloc(sizeof(char*)*total_files);
	fqs_input->total_files = total_files;

	for(x1=0; x1<fqs_input -> total_files; x1++){
		fqs_input->files1[x1] = strdup(files1[x1]);
		if(files2)fqs_input->files2[x1] = strdup(files2[x1]);
		fqs_input->files3[x1] = strdup(files3[x1]);
	}
	fqs_input->current_file_no = 0;
	fqs_input->current_read_no = 0;
	return input_mFQ_open_files(fqs_input);
}
int input_mFQ_next_read(input_mFQ_t * fqs_input, char * readname , char * read, char * qual ){
	char tmpline [MAX_READ_NAME_LEN+1];
	int ret = -1,x1;
	if(fqs_input -> current_file_no == fqs_input -> total_files) return -1;
	while(1){
		ret = autozip_gets(&fqs_input -> autofp1, tmpline, MAX_READ_NAME_LEN);
		int write_ptr=0;
		if(ret==0){
			ret = input_mFQ_next_file(fqs_input);
			if(ret >=0) continue;
			return -1;
		} else if(ret<0) return -1;

		#ifdef __MINGW32__
		sprintf(readname, "R%011I64d", fqs_input -> current_read_no);
		#else
		sprintf(readname, "R%011lld", fqs_input -> current_read_no);
		#endif
		readname[12]='|';
		ret = autozip_gets(&fqs_input -> autofp1, readname+13, MAX_READ_NAME_LEN);
		readname[13+ret-1]='|';
		write_ptr = 13+ret;

		autozip_gets(&fqs_input -> autofp1, tmpline, MAX_READ_NAME_LEN);

		ret = autozip_gets(&fqs_input -> autofp1, readname+write_ptr, MAX_READ_NAME_LEN);
		ret --;
		for(x1=write_ptr; x1<write_ptr+ret; x1++) if(readname[x1]>='/') readname[x1]++;
		readname[write_ptr+ret]='|';
		write_ptr += ret;

		if(fqs_input->files2){
			ret = autozip_gets(&fqs_input -> autofp2, tmpline, MAX_READ_NAME_LEN);
			if(ret<=0) return -1;
			ret = autozip_gets(&fqs_input -> autofp2, readname+write_ptr, MAX_READ_NAME_LEN);
			ret --;
			readname[write_ptr+ret]='|';
			write_ptr += ret;

			autozip_gets(&fqs_input -> autofp2, tmpline, MAX_READ_NAME_LEN);
			ret = autozip_gets(&fqs_input -> autofp2, readname+write_ptr, MAX_READ_NAME_LEN);
			for(x1=write_ptr; x1<write_ptr+ret; x1++) if(readname[x1]>='/') readname[x1]++;
			ret --;
			readname[write_ptr+ret]=0;
		}else write_ptr+=sprintf(readname+write_ptr,"|input#%04d@L%03d", fqs_input -> current_file_no, fqs_input  -> current_guessed_lane_no);

		ret = autozip_gets(&fqs_input -> autofp3, tmpline, MAX_READ_NAME_LEN);
		if(ret<=0) return -1;
		ret = autozip_gets(&fqs_input -> autofp3, read, MAX_READ_LENGTH);
		ret --; // read length excludes "\n"
		read[ret]=0;
		autozip_gets(&fqs_input -> autofp3, tmpline, MAX_READ_NAME_LEN);
		autozip_gets(&fqs_input -> autofp3, qual, MAX_READ_LENGTH);
		qual[ret]=0;

		break;
	}
	fqs_input->current_read_no ++;
	return ret;
}


void input_mFQ_close(input_mFQ_t * fqs_input){
	int x1;

	input_mFQ_fp_close(fqs_input);
	for(x1=0; x1<fqs_input -> total_files; x1++){
		free(fqs_input->files1[x1]);
		if(fqs_input->files2)free(fqs_input->files2[x1]);
		free(fqs_input->files3[x1]);
	}
	free(fqs_input->files1);
	if(fqs_input->files2)free(fqs_input->files2);
	free(fqs_input->files3);
}

int input_mFQ_seek(input_mFQ_t * fqs_input, input_mFQ_pos_t * pos ){
	if(fqs_input -> current_file_no != pos -> current_file_no){
		if(fqs_input -> current_file_no < fqs_input -> total_files)input_mFQ_fp_close(fqs_input);
		fqs_input -> current_file_no = pos -> current_file_no;
		if(fqs_input -> current_file_no < fqs_input -> total_files)input_mFQ_open_files(fqs_input);
	}
	if(fqs_input -> current_file_no < fqs_input -> total_files){
		fqs_input -> current_read_no = pos -> current_read_no;
		if(fqs_input -> autofp1.is_plain){
			fseeko(fqs_input -> autofp1.plain_fp,  pos -> pos_file1, SEEK_SET);
			if(fqs_input -> files2)fseeko(fqs_input -> autofp2.plain_fp,  pos -> pos_file2, SEEK_SET);
			fseeko(fqs_input -> autofp3.plain_fp,  pos -> pos_file3, SEEK_SET);
		}else{
			seekgz_seek(&fqs_input -> autofp1.gz_fp,&pos -> zpos_file1);
			if(fqs_input -> files2)seekgz_seek(&fqs_input -> autofp2.gz_fp,&pos -> zpos_file2);
			seekgz_seek(&fqs_input -> autofp3.gz_fp,&pos -> zpos_file3);
		}
	}
	return 0;
}

int input_mFQ_tell(input_mFQ_t * fqs_input, input_mFQ_pos_t * pos ){
	memset(pos, 0, sizeof(input_mFQ_pos_t));
	pos -> current_file_no = fqs_input -> current_file_no;
	pos -> current_read_no = fqs_input -> current_read_no;

	if(fqs_input -> current_file_no < fqs_input -> total_files){
		if(fqs_input -> autofp1.is_plain){
			pos -> pos_file1 = ftello(fqs_input -> autofp1.plain_fp);
			if(fqs_input -> files2)pos -> pos_file2 = ftello(fqs_input -> autofp2.plain_fp);
			pos -> pos_file3 = ftello(fqs_input -> autofp3.plain_fp);
		}else{
			seekgz_tell(&fqs_input -> autofp1.gz_fp,&pos -> zpos_file1);
			if(fqs_input -> files2)seekgz_tell(&fqs_input -> autofp2.gz_fp,&pos -> zpos_file2);
			seekgz_tell(&fqs_input -> autofp3.gz_fp,&pos -> zpos_file3);
		}
	}
	return 0;
}
