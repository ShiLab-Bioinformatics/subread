#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


char * __converting_char_table = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNGNNNCNNNNNNNNNNNNAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNGNNNCNNNNNNNNNNNNAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN  ";


void reverse_qual(int read_len, char * InBuff)
{
    int i;
    if(!InBuff) return;
    if(!InBuff[0]) return;
    for (i=0; i<read_len/2; i++)
    {
        char tmp;
        tmp = InBuff[i];
        InBuff[i] = InBuff[read_len -1-i];
        InBuff[read_len -1-i] = tmp;
    }
}


void reverse_read(int read_len, char * InBuff)
{
    int i;
    for (i=0; i<read_len/2; i++)
    {
        int rll1 = read_len - 1 - i;
        unsigned char tmp = InBuff[rll1];

        InBuff[rll1] = *(__converting_char_table+InBuff[i]);
        InBuff[i] = *(__converting_char_table+tmp);

    }
    if(i*2 == read_len-1)
    {
        InBuff[i] = *(__converting_char_table+InBuff[i]);
    }

}


int main(int argc, char ** argv)
{
	char fn[200];
	char linebuf[3000];

	sprintf(fn, "%s-A.fq", argv[1]);
	FILE * op1=f_subr_open(fn, "w");
	sprintf(fn, "%s-B.fq", argv[1]);
	FILE * op2=f_subr_open(fn, "w");

	while(1)
	{
		int flds, flags=0;
		char * sret = fgets(linebuf, 2999, stdin);
		char * readtxt=NULL, *qualtxt=NULL;
		if(!sret) break;
		if(linebuf[0]=='@') continue;

		char * readname = strtok(linebuf,"\t");
		for(flds=1; flds<11; flds++)
		{
			char * fldstr = strtok(NULL, "\t");
			if(flds==1) flags = atoi(fldstr);
			else if(flds == 9)
				readtxt = fldstr;
			else if(flds == 10)
				qualtxt = fldstr;
		}

		int rlen = strlen(readtxt);
		if(qualtxt[rlen]) qualtxt[rlen]=0;

		if(flags & 16)
		{
			reverse_read( strlen(readtxt), readtxt);
			reverse_qual( strlen(readtxt), qualtxt);
		}

		fprintf((flags & 0x80)? op2: op1, "@%s\n%s\n+\n%s\n", readname, readtxt, qualtxt);
	
	}

	fclose(op1);
	fclose(op2);
}
