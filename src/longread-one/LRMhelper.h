#ifndef __LRM_HELPER_H_
#define __LRM_HELPER_H_

void LRMbasic_sort(void * arr, int items, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r));

void LRMmerge_sort(void * arr, int arr_size, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2));

#endif
