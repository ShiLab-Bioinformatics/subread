#include <stdio.h>
#include <stdlib.h>

int main(int argc, char ** argv)
{
	int chromosomes = atoi(argv[1]);
	int x1, chro_offset;

	for(x1=0; x1<chromosomes; x1++)
	{
		if(x1) puts("");
		printf(">chr%d", x1+1);
		for(chro_offset = 0; chro_offset<500000; chro_offset++)
		{
			if(chro_offset % 70==0) puts("");
			char nch_i = myrand_rand() % 4 ;
			char nch = nch_i?(nch_i<2?'T':(nch_i < 3?'C':'G')):'A';
			
			putchar(nch);
		}
	}
	return 0;
}
