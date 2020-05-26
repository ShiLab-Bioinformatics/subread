/***************************************************************

   The Subread software package is free software package: 
   you can redistribute it and/or modify it under the terms
   of the GNU General Public License as published by the 
   Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   Subread is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty
   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   
   See the GNU General Public License for more details.

   Authors: Drs Yang Liao and Wei Shi

  ***************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define NORMAL_INTS 35

int mergeIntervals(unsigned int * intervals, unsigned int * result_ints, int nints)
{

	unsigned int * stack_buffer_mem; 
	int stock_pointer = 0;
	stack_buffer_mem = result_ints;
	

	int ii,jj;

	for(ii=0; ii<nints; ii++)
	{
		unsigned int small_jj_value = 0xffffffff;
		int small_jj_id = -1;
		for(jj = ii + 1; jj < nints; jj++)
		{
			if(small_jj_value > intervals[jj*2])
			{
				small_jj_value = intervals[jj*2];
				small_jj_id = jj;
			}
		}

		if(small_jj_value < intervals[ii*2])
		{
			unsigned int tmpv;
			tmpv = intervals[ii*2];
			intervals[ii*2] = small_jj_value;
			intervals[small_jj_id * 2] = tmpv;

			tmpv = intervals[ii*2 + 1];
			intervals[ii*2 + 1] = intervals[small_jj_id * 2 + 1];
			intervals[small_jj_id * 2 + 1] = tmpv;
		}
	}

	stack_buffer_mem[stock_pointer*2] = intervals[0];
	stack_buffer_mem[stock_pointer*2+1] = intervals[1];
	stock_pointer++;
 
	for (ii = 1 ; ii < nints; ii++)
	{
		//int top_start = stack_buffer_mem[stock_pointer*2-2];
		int top_stop = stack_buffer_mem[stock_pointer*2-1];
 
		if (top_stop < intervals[ii*2])
		{
			stack_buffer_mem[stock_pointer * 2] = intervals[ii*2];
			stack_buffer_mem[stock_pointer * 2+1] = intervals[ii*2+1];
			stock_pointer++;
		}
		else if (top_stop < intervals[ii*2+1])
		{
			top_stop = intervals[ii*2+1];
			stack_buffer_mem[stock_pointer*2-1] = top_stop;
		}
		// Otherwise the ii-th interval is useless because it is enclosed in the top item
	}
 
	return stock_pointer;
}
 
#ifdef MAKE_TEST_INTERVL_MERGE
void print_gaps(unsigned int * gaps, int gapn)
{
	int ii;
	for(ii=0; ii<gapn; ii++)
		printf("{%d,%d}  ", gaps[ii*2], gaps[ii*2+1]);
	printf("\n");
}
#endif

#ifdef MAKE_TEST_INTERVL_MERGE
int main()
{
	unsigned int inbuff[30];
	unsigned int outbuff[30];

	unsigned int q1 [12] = 	{10,30,  20,40,  120,140,  70,80,  5,15,  90,130};
	unsigned int q2 [12] =	{ 100,120,  30,31 , 25,35,  110,115,  125,135,  30,40};
	unsigned int q3 [12] =	{ 30,40,  40,50,  20,29,  29,30,  3,100,  1,2};
	unsigned int q4 [12] =  {30,140,  40,50,  20,29,  29,30,  3,100,  1,3};
	unsigned int q5 [12] =  {90,100,  105,130,  135,150,  80,95,  135,151, 90,100};

	memcpy(inbuff, q1, 4*12);
	int gaps = mergeIntervals(inbuff, outbuff, 6);
	print_gaps(outbuff, gaps);

	memcpy(inbuff, q2, 4*12);
	gaps = mergeIntervals(inbuff, outbuff, 6);
	print_gaps(outbuff, gaps);

	memcpy(inbuff, q3, 4*12);
	gaps = mergeIntervals(inbuff, outbuff, 6);
	print_gaps(outbuff, gaps);

	memcpy(inbuff, q4, 4*12);
	gaps = mergeIntervals(inbuff, outbuff, 6);
	print_gaps(outbuff, gaps);

	memcpy(inbuff, q5, 4*12);
	gaps = mergeIntervals(inbuff, outbuff, 6);
	print_gaps(outbuff, gaps);

	return 0;
}
#endif
