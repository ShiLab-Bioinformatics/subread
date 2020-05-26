#ifndef __LRMCHRO_EVENT_H_
#define __LRMCHRO_EVENT_H_

#include "LRMconfig.h"

int LRMevents_load_annot(LRMcontext_t * context);
int LRMevents_filtering(LRMcontext_t * context);
int LRMevents_reorder(LRMcontext_t * context);
int LRMevents_build_entries(LRMcontext_t * context);
int LRMevents_lookup(LRMcontext_t * context, unsigned int abs_pos, int event_type_masks, int search_large_side, LRMevent_t ** res);

int LRMchro_event_new(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context, LRMevent_t * new_event);

// high_penalty_for_creating_gap is for creating less CIGAR options (<65535)
int LRMindel_dynamic_search(LRMcontext_t* context, LRMthread_context_t* thread_context, int expected_indels, unsigned int last_correct_base_on_chro, char * corrected_read, int last_correct_base, int first_correct_base, unsigned int * total_mismatched, int high_penalty_for_creating_gap);
void LRMreverse_read_and_qual(LRMcontext_t * context, LRMthread_context_t * thread_context, LRMread_iteration_context_t * iteration_context);
int LRMdynamic_to_ends(LRMcontext_t* context, LRMthread_context_t* thread_context, LRMread_iteration_context_t * iteration_context, int last_mapped_in_read, unsigned int last_correct_base_on_chro, int search_to_3end);
int LRMdynamic_in_middle(LRMcontext_t* context, LRMthread_context_t* thread_context, LRMread_iteration_context_t * iteration_context, int last_correct_base, int first_correct_base ,  unsigned int last_correct_base_on_chro, int expected_indels);
#endif
