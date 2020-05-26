/* Double-Click To Select Code */
 
#include<stdio.h>
 
void merge_sort_run(void * arr, int start, int items, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2));
void merge_sort(void * arr, int arr_size, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2));
void merge_ints(void * arr, int start, int items1, int items2);
int compare_ints(void* arr, int l, int r);
int exchange_ints(void* arr, int l, int r);

void main()
{

int arr[30];
int val [] = {9,1,2,9,6,7,8,9,1,2,3,9,4,1,2,3,4,6,1,9,3,1,4,5,3,2,5,4,2,4,8};
int x;

for(x=0; x<30;x++)arr[x]=val[x];

merge_sort(arr, 30, compare_ints, exchange_ints, merge_ints);

for(x=0; x<30;x++)printf("V[%d]=%d\n",x, arr[x]);

}
 
int exchange_ints(void* arr, int l, int r)
{
	int *arri = arr;
	int tm;
	tm=arri[l];
	arri[l]=arri[r];
	arri[r]=tm;
}

int compare_ints(void* arr, int l, int r)
{
	int * arri = arr;
	if(arri[l]==arri[r])return 0;
	if(arri[l]>arri[r])return 1;
	return -1;
}



void merge_ints(void * arr, int start, int items1, int items2)
{
	int r1, r2;
	int * arri = arr;
	r1=start;
	r2=start+items1;
	int * tmp = malloc(sizeof(int)*(items1+items2));
	int x;

	for(x=0; x<items1+items2; x++)
	{
		if((r1>= start+items1)||(r2<start+items1+items2 && arri[r1]>=arri[r2]))
		{
			tmp[x]=arri[r2];
			r2++;
		}else{
			tmp[x]=arri[r1];
			r1++;
		}
	}

	memcpy(arri+start , tmp, sizeof(int)*(items1+items2));
	free(tmp);
}


void merge_sort_run(void * arr, int start, int items, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2))
{
	if(items > 4)
	{
		int xx,half_point = items/2;

		merge_sort_run(arr, start, half_point, compare, exchange, merge);
		merge_sort_run(arr, start + half_point, items - half_point, compare, exchange, merge);
		merge(arr, start, half_point, items - half_point);
		printf("IN: %d-%d-%d\n", start, start + half_point, start+items);
		for(xx=start; xx < items+start; xx++)
			printf("INNER: %d\n", ((int *)arr)[xx]);
	}
	else
	{
		int i, j, xx;
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
		printf("RD: %d-%d\n", start,start+items);
		for(xx=start; xx < items+start; xx++)
			printf("INRED: %d\n", ((int *)arr)[xx]);
	}
}
void merge_sort(void * arr, int arr_size, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2))
{
	merge_sort_run(arr, 0, arr_size, compare, exchange, merge);
}
