#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#ifndef __MINGW32__
#include <sys/mman.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "core.h"
#include "core-bigtable.h"
#include "gene-algorithms.h"

#define TABLE_PRELOAD_SIZE 24

#define calc_offset(pair_NO)  ((pair_NO) * (1 + global_context -> input_reads.is_paired_end_reads) + is_second_read) * (sizeof(short) * 3 * global_context -> config.big_margin_record_size + (sizeof(mapping_result_t) + global_context -> config.do_breakpoint_detection*sizeof(subjunc_result_t)) * global_context -> config.multi_best_reads)


void bigtable_lock(global_context_t * global_context){
	subread_lock_occupy(&global_context -> bigtable_lock);
}

void bigtable_unlock(global_context_t * global_context){
	subread_lock_release(&global_context -> bigtable_lock);
}

unsigned int get_handling_thread_number(global_context_t * global_context , subread_read_number_t pair_number){

	/*
	if(global_context -> config.all_threads<2) return 0;

	for(ret = 0; ret < global_context -> config.all_threads; ret++){
		if(global_context -> input_reads.start_read_number_blocks[ret] > pair_number) return ret - 1;
	}
	return global_context -> config.all_threads - 1;
	*/return 0;
}

#define get_inner_pair(a,b) (b)

bigtable_cached_result_t * bigtable_retrieve_cache(global_context_t * global_context , thread_context_t * thread_context , subread_read_number_t pair_number, int is_second_read, int load_more);

void bigtable_readonly_result(global_context_t * global_context , thread_context_t * thread_context , subread_read_number_t pair_number, int result_number, int is_second_read, mapping_result_t * return_ptr, subjunc_result_t * return_junction_ptr){
	int rlen = -1;

	if(global_context -> bigtable_cache_file_fp){
		int loadjunc;
		long long inner_pair_number = get_inner_pair(global_context, pair_number);

		if(global_context -> bigtable_cache_file_loaded_fragments_begin >=0){
			bigtable_write_thread_cache(global_context);
			global_context -> bigtable_cache_file_loaded_fragments_begin = -1;
		}

		for(loadjunc = 0; loadjunc < 2; loadjunc++){
			if(loadjunc){ if(!return_junction_ptr) continue; }
			else{ if(!return_ptr) continue; }

			unsigned long long offset = calc_offset(inner_pair_number);
			offset += sizeof(short) * 3 * global_context -> config.big_margin_record_size ;
			if(loadjunc) offset += sizeof(mapping_result_t) * global_context -> config.multi_best_reads + sizeof(subjunc_result_t) * result_number;
			else	offset += sizeof(mapping_result_t) * result_number;
			fseeko(global_context -> bigtable_cache_file_fp , offset , SEEK_SET);

			void * write_ptr = return_ptr;
			if(loadjunc) write_ptr = return_junction_ptr;

			rlen = fread(write_ptr, loadjunc?sizeof(subjunc_result_t):sizeof(mapping_result_t), 1,  global_context -> bigtable_cache_file_fp);
			if(rlen <1)
				SUBREADprintf("UNABLE TO READ RESULT\n");
		}
	}else{
		bigtable_cached_result_t * rett = bigtable_retrieve_cache(global_context , thread_context , pair_number, is_second_read,0);

		int best_offset = result_number;
		if(return_ptr)memcpy(return_ptr, rett -> alignment_res + best_offset, sizeof(mapping_result_t));
		if(return_junction_ptr)memcpy(return_junction_ptr, rett -> subjunc_res + best_offset, sizeof(subjunc_result_t));
	}
}


int init_bigtable_results(global_context_t * global_context, int is_rewinding)
{

	if(global_context -> config.use_memory_buffer) {
		global_context -> bigtable_chunked_fragments = global_context -> config.reads_per_chunk+1;
		global_context -> bigtable_cache_size = global_context -> bigtable_chunked_fragments * (1+global_context -> input_reads.is_paired_end_reads);
	} else {
		global_context -> bigtable_chunked_fragments = 110llu*1024*1024;
		global_context -> bigtable_cache_size = (1+0*global_context -> config.all_threads) * global_context -> bigtable_chunked_fragments * (1+global_context -> input_reads.is_paired_end_reads);
	}


	//SUBREADprintf("reads_per_chunk = %u ; cached_single_reads = %u ; size of each read = %d + %d\n",  global_context -> config.reads_per_chunk,  global_context -> bigtable_cache_size, sizeof(mapping_result_t) , sizeof(subjunc_result_t));
	//
	if(is_rewinding){
		memset( global_context -> bigtable_cache_malloc_ptr, 0, sizeof(mapping_result_t) *  global_context -> config.multi_best_reads * global_context -> bigtable_cache_size );
		if(global_context -> config.do_breakpoint_detection)
			memset(global_context -> bigtable_cache_malloc_ptr_junc , 0, sizeof(subjunc_result_t) * global_context -> config.multi_best_reads *  global_context -> bigtable_cache_size);
	}else{
		global_context -> bigtable_cache = malloc(sizeof(bigtable_cached_result_t) * global_context -> bigtable_cache_size);
		global_context -> bigtable_cache_malloc_ptr = calloc( sizeof(mapping_result_t) ,  global_context -> config.multi_best_reads * global_context -> bigtable_cache_size );
		if(global_context -> config.do_breakpoint_detection)
			global_context -> bigtable_cache_malloc_ptr_junc = calloc(sizeof(subjunc_result_t)  , global_context -> config.multi_best_reads *  global_context -> bigtable_cache_size);
		assert(NULL != global_context -> bigtable_cache_malloc_ptr );
	}

	int xk1;
	if(!is_rewinding) for(xk1 = 0; xk1 < global_context -> bigtable_cache_size; xk1++){
		global_context -> bigtable_cache [xk1].alignment_res = (mapping_result_t*) global_context -> bigtable_cache_malloc_ptr + global_context -> config.multi_best_reads * xk1;

		if(global_context -> config.do_breakpoint_detection)
			global_context -> bigtable_cache [xk1].subjunc_res = (subjunc_result_t*)global_context -> bigtable_cache_malloc_ptr_junc + global_context -> config.multi_best_reads * xk1;
	}
	if(global_context->config.do_big_margin_filtering_for_junctions)for(xk1 = 0; xk1 < global_context -> bigtable_cache_size; xk1++)
		memset(global_context -> bigtable_cache [xk1].big_margin_data, 0, sizeof(global_context -> bigtable_cache [xk1].big_margin_data));

	subread_init_lock(&global_context -> bigtable_lock);

	assert(global_context -> config.use_memory_buffer);
	global_context -> bigtable_cache_file_fp = NULL;
	return 0;
}

mapping_result_t * _global_retrieve_alignment_ptr(global_context_t * global_context, subread_read_number_t pair_number, int is_second_read, int best_read_id){
	mapping_result_t * ret;
	bigtable_retrieve_result(global_context, NULL, pair_number, best_read_id, is_second_read, &ret, NULL);
	return ret;
}

subjunc_result_t * _global_retrieve_subjunc_ptr(global_context_t * global_context, subread_read_number_t pair_number, int is_second_read, int best_read_id){
	subjunc_result_t * ret;
	bigtable_retrieve_result(global_context, NULL, pair_number, best_read_id, is_second_read, NULL, &ret);
	return ret;
}


#define calc_file_location(pair_no) (pair_no)* (1 + global_context -> input_reads.is_paired_end_reads) *  (sizeof(short) * 3 * global_context -> config.big_margin_record_size + (sizeof(mapping_result_t) + global_context -> config.do_breakpoint_detection*sizeof(subjunc_result_t)) * global_context -> config.multi_best_reads)


void bigtable_write_thread_cache(global_context_t * global_context){
	assert(global_context -> bigtable_cache_file_fp == NULL);
	return;
}


void wait_occupied(global_context_t * global_context , unsigned long long old_begin){
	while(old_begin == global_context -> bigtable_cache_file_loaded_fragments_begin){
		int i, all_released = 1;
		for(i=0; i< global_context -> bigtable_chunked_fragments ; i++)
		{
			if(global_context -> bigtable_cache[i].status == CACHE_STATUS_OCCUPIED)
				all_released = 0;
		}
		if(all_released) break;
	}
}

bigtable_cached_result_t * bigtable_retrieve_cache(global_context_t * global_context , thread_context_t * thread_context , subread_read_number_t pair_number, int is_second_read, int load_more)
{
	long long inner_pair_number = get_inner_pair(global_context, pair_number);
	long long load_start_pair_no = inner_pair_number - inner_pair_number % global_context -> bigtable_chunked_fragments;
	bigtable_cached_result_t * ret_cache = global_context -> bigtable_cache + (inner_pair_number - load_start_pair_no)* (1+global_context -> input_reads.is_paired_end_reads) + is_second_read;
	return ret_cache;
}

int bigtable_retrieve_result(global_context_t * global_context , thread_context_t * thread_context , subread_read_number_t pair_number, int result_number, int is_second_read, mapping_result_t ** return_ptr, subjunc_result_t ** return_junction_ptr)
{
	bigtable_cached_result_t * cache = bigtable_retrieve_cache(global_context, thread_context, pair_number, is_second_read, 0);

	int best_offset = result_number;

	if(return_ptr)(*return_ptr)= cache -> alignment_res + best_offset;
	if(return_junction_ptr)(*return_junction_ptr)= cache -> subjunc_res + best_offset;

	return 0;
}


unsigned short * _global_retrieve_big_margin_ptr(global_context_t * global_context, subread_read_number_t pair_number, subread_read_number_t is_second_read){
	bigtable_cached_result_t * cache = bigtable_retrieve_cache(global_context, NULL,  pair_number, is_second_read, 0);
	return cache -> big_margin_data;
}


void bigtable_release_result(global_context_t * global_context , thread_context_t * thread_context , subread_read_number_t pair_number, int commit_change){
	assert(NULL == "File result cache is no longer used. Do not call this function.");
}


int finalise_bigtable_results(global_context_t * global_context){
	assert(global_context -> bigtable_cache_file_fp == NULL);
	free(global_context -> bigtable_cache_malloc_ptr);
	if(global_context -> config.do_breakpoint_detection) free(global_context -> bigtable_cache_malloc_ptr_junc);
	free(global_context -> bigtable_cache);
	return 0;
}

























void bktable_key_deallocator(void * key){
	free(key);
}

void bktable_bucket_deallocator(void * vbuk){
	bucketed_table_bucket_t * buk = vbuk;
	free(buk -> positions);
	free(buk -> details);
	free(buk);
}

void bktable_init(bucketed_table_t * tab, unsigned int maximum_interval_length, unsigned int expected_items){
	memset(tab, 0, sizeof(bucketed_table_t));
	tab -> expected_items = expected_items;
	tab -> maximum_interval_length = maximum_interval_length;
	tab -> entry_table = HashTableCreate(expected_items / 3);

	HashTableSetDeallocationFunctions(tab->entry_table, bktable_key_deallocator, bktable_bucket_deallocator);
	HashTableSetKeyComparisonFunction(tab->entry_table, fc_strcmp_chro);
	HashTableSetHashFunction(tab->entry_table, HashTableStringHashFunction);
}

void bktable_destroy(bucketed_table_t * tab){
	HashTableDestroy(tab -> entry_table);
}

void bktable_append(bucketed_table_t * tab, char * chro, unsigned int pos, void * detail){
	unsigned int twokeys[2], keyi;
	//assert(detail != NULL);
	if(detail == NULL) return;
	twokeys[0] = pos - pos % tab -> maximum_interval_length;
	if(twokeys[0] > tab -> maximum_interval_length)
		twokeys[1] = twokeys[0] - tab -> maximum_interval_length;
	else	twokeys[1] = 0xffffffff;

	for(keyi = 0; keyi < 2 ; keyi++)
	{
		char static_key [20 + MAX_CHROMOSOME_NAME_LEN];
		unsigned int curr_key_pos = twokeys[keyi];
		if(curr_key_pos == 0xffffffff) continue;

		sprintf(static_key, "%s:%u", chro, curr_key_pos);

		bucketed_table_bucket_t * had_items = HashTableGet(tab -> entry_table, static_key);

		//SUBREADprintf("INSERT STATIC_KEY [%d] = %s ; FOUND = %p ; ITEMS=%d\n", keyi, static_key, had_items, had_items?had_items->items:-1);

		if(NULL == had_items)
		{
			had_items = malloc(sizeof(bucketed_table_bucket_t));
			memset(had_items, 0, sizeof(bucketed_table_bucket_t));

			had_items -> capacity = BUCKETED_TABLE_INIT_ITEMS;
			had_items -> positions = malloc(sizeof(int) * BUCKETED_TABLE_INIT_ITEMS);
			had_items -> details  = malloc(sizeof(void *) *BUCKETED_TABLE_INIT_ITEMS);
			had_items -> keyed_bucket = curr_key_pos;
			had_items -> maximum_interval_length = tab -> maximum_interval_length;

			char * dynamic_key = malloc(strlen(static_key) + 1);
			strcpy(dynamic_key, static_key);
			HashTablePut(tab->entry_table, dynamic_key, had_items);
		}
		if(had_items -> capacity <= had_items -> items){
			had_items -> capacity = max(had_items -> capacity + 5, had_items -> capacity * 1.3);
			had_items -> positions = realloc(had_items -> positions, had_items -> capacity * sizeof(int));
			had_items -> details = realloc(had_items -> details, had_items -> capacity * sizeof(void *));
		}

		had_items -> positions[ had_items -> items ] = pos;
		had_items -> details[ had_items -> items ] = detail;
		had_items -> items++;
	}

	tab -> fragments ++;
}


void bktable_free_ptrs(void * bukey, void * buckv, HashTable * tab){
	int x1;
	bucketed_table_bucket_t * buck = buckv;
	for(x1 = 0; x1 < buck -> items; x1++)
	{
		if(buck->positions[x1] - buck->positions[x1] % buck -> maximum_interval_length == buck -> keyed_bucket)
		{
			//SUBREADprintf("FREE : %u ~ %u\n",  buck->positions[x1] , buck -> keyed_bucket);
			free(buck->details[x1]);
		}
	}
}

int bktable_lookup(bucketed_table_t * tab, char * chro, unsigned int start_pos, unsigned int interval_length, unsigned int * hit_pos_list, void ** hit_ptr_list, int max_hits){
	unsigned int my_key_pos;
	my_key_pos = start_pos - start_pos % tab -> maximum_interval_length;

	char static_key [20 + MAX_CHROMOSOME_NAME_LEN];
	sprintf(static_key, "%s:%u", chro, my_key_pos);

	bucketed_table_bucket_t * had_items = HashTableGet(tab -> entry_table, static_key);

	//if(strcmp(chro, "X")==0 && abs(start_pos - 21067381 - 500) < 10){
	//	SUBREADprintf("LOOK STATIC_KEY = %s, BUCK=%p, ITEMS=%d\n", static_key, had_items, had_items?had_items->items:-1);
	//}

	if(!had_items) // no bucket at all.
		return 0;

	int item_i, ret = 0;
	for(item_i = 0 ; item_i < had_items->items; item_i++){
		unsigned int potential_pos = had_items -> positions[item_i];
		//if(strcmp(chro, "X")==0 && abs(start_pos - 21067381 - 500) < 10){
		//	SUBREADprintf("POTENTIAL_POS=%u,  ACCEPT=%u ~ %u\n", potential_pos, start_pos, start_pos + interval_length);
		//} 
		if(potential_pos >= start_pos && potential_pos < start_pos + interval_length)
		{
			hit_pos_list[ret] = potential_pos;
			hit_ptr_list[ret] = had_items -> details[item_i];
			ret ++;
			if(ret >= max_hits) break;
		}
	}

	return ret;
}


void fraglist_init(fragment_list_t * list){
	memset(list, 0, sizeof(fragment_list_t));
	list -> capacity = FRAGMENT_LIST_INIT_ITEMS;
	list -> fragment_numbers = malloc(sizeof(subread_read_number_t) * list -> capacity);
}

void fraglist_destroy(fragment_list_t * list){
	free(list -> fragment_numbers);
}

void fraglist_append(fragment_list_t * list, subread_read_number_t fragment_number){
	if(list -> fragments >= list -> capacity){
		//#warning "============== COMMENT DEBUG INFO ====================="
		//SUBREADprintf("REALLOC_PRT_IN = %d , %p\n",list -> capacity , list -> fragment_numbers );
		list -> capacity = max(list -> capacity + 5, list -> capacity * 1.3);
		list -> fragment_numbers = realloc(list -> fragment_numbers, sizeof(subread_read_number_t) * list -> capacity);
	//	SUBREADprintf("REALLOC_PRT_OUT = %d , %p\n",list -> capacity , list -> fragment_numbers );
	}

	list -> fragment_numbers[ list -> fragments ++ ] = fragment_number;
}

