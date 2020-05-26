#include<stdio.h>
#include<string.h>
#include"LRMconfig.h"
#include"LRMhelper.h"

void LRMbasic_sort_run(void * arr, int start, int items, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r)){
	int i, j;
	for(i=start; i< start + items - 1; i++)
	{
		int min_j = i;
		for(j=i + 1; j< start + items; j++)
		{
			if(compare(arr, min_j, j) > 0)
				min_j = j;
		}
		if(i!=min_j)
			exchange(arr, i, min_j);
	}
}

void LRMbasic_sort(void * arr, int items, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r)){
	LRMbasic_sort_run(arr, 0, items, compare, exchange);
}


void LRMmerge_sort_run(void * arr, int start, int items, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2))
{
	if(items > 6)
	{
		int half_point = items/2;
		LRMmerge_sort_run(arr, start, half_point, compare, exchange, merge);
		LRMmerge_sort_run(arr, start + half_point, items - half_point, compare, exchange, merge);
		merge(arr, start, half_point, items - half_point);
	}
	else
	{
		LRMbasic_sort_run(arr, start, items, compare, exchange);
	}
}
void LRMmerge_sort(void * arr, int arr_size, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2))
{
	LRMmerge_sort_run(arr, 0, arr_size, compare, exchange, merge);
}
