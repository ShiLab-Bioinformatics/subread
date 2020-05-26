#ifndef __CORE_BIGTABLE_H_
#define __CORE_BIGTABLE_H_

#include "subread.h"
#include "core.h"
#include "hashtable.h"
#include "gene-algorithms.h"

#define CACHE_STATUS_RELEASED 0
#define CACHE_STATUS_OCCUPIED 1

// This function creates an empty data structure for all results.
// The number of reads is unknown at this stage. 
int init_bigtable_results(global_context_t * global_context, int is_rewind);

// This function tries to retrieve the required result data structure into memory and set the return_ptr to the address of the data structure.
// Junction ptr can be NULL.
// This function returns ZERO if the record is available. It returns -1 if the record is unavailable. 
int bigtable_retrieve_result(global_context_t * global_context , thread_context_t * thread_context , subread_read_number_t pair_number, int result_number, int is_second_read, mapping_result_t ** return_ptr, subjunc_result_t ** return_junction_ptr);

// This function notifies the bigtable subsystem to save changes and deallocate the memory block if necessary.
// Junction ptr can be NULL.
// If the data has been changed, commit_change must be set to an non-ZERO value
void bigtable_release_result(global_context_t * global_context , thread_context_t * thread_context , subread_read_number_t pair_number, int commit_change);

// This function destroies the buffers and deletes the temporary MMAP files.
int finalise_bigtable_results(global_context_t * global_context);

void bigtable_readonly_result(global_context_t * global_context , thread_context_t * thread_context , subread_read_number_t pair_number, int result_number, int is_second_read, mapping_result_t * return_ptr, subjunc_result_t * return_junction_ptr);

void bktable_append(bucketed_table_t * tab, char * chro, unsigned int pos, void * detail);

int bktable_lookup(bucketed_table_t * tab, char * chro, unsigned int start_pos, unsigned int interval_length, unsigned int * hit_pos_list, void ** hit_ptr_list, int max_hits);

void bktable_init(bucketed_table_t * tab, unsigned int maximum_interval_length, unsigned int expected_items);

void bktable_destroy(bucketed_table_t * tab);

void bktable_free_ptrs(void * bkey, void * buckv, HashTable * tab);

void fraglist_init(fragment_list_t * list);

void fraglist_append(fragment_list_t * list, subread_read_number_t fragment_number);

void fraglist_destroy(fragment_list_t * list);

void bigtable_write_thread_cache(global_context_t * global_context);

#endif
