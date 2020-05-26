#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "zlib.h"


main()
{
	char * tdata = "ABCABC\n";
	int tdata_len = 7, ret;

	z_stream strm;
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	strm.avail_in = 0;
	strm.next_in = Z_NULL;
	ret = deflateInit(&strm, 0);
	if (ret != Z_OK)printf("Ohh!\n");

	char * out_buff = malloc(999999);
	strm.avail_out = 99999;
	strm.next_out = out_buff;
	strm.next_in = tdata;
	strm.avail_in = tdata_len;
	ret = deflate(&strm, Z_FINISH);

	int have = 99999 - strm.avail_out;

	printf("RET=%d; LEN=%d\n",ret, have);

	FILE * ofp = fopen("tt.gz","wb");
	fwrite(out_buff,1,have,ofp);
	fclose(ofp);

	z_stream strmx;
	strmx.zalloc = Z_NULL;
	strmx.zfree = Z_NULL;
	strmx.opaque = Z_NULL;
	strmx.avail_in = 0;
	strmx.next_in = Z_NULL;
	ret = inflateInit(&strmx);
	if (ret != Z_OK)printf("Ohh!\n");
	char * in_buff = malloc(999999);

	strmx.avail_out = 99999;
	strmx.next_out = in_buff;
	strmx.next_in = out_buff;
	strmx.avail_in = have;
	ret = inflate(&strmx, Z_FINISH);
	have = 99999 - strmx.avail_out;

	
	printf("XRET=%d; LEN=%d; RES=%s\n",ret, have, in_buff);
	FILE * fbig = fopen("t.bin","rb");
	int fpos = 0;
	while(!feof(fbig))
	{
		int nch = fgetc(fbig);
		if(nch<0) break;
		in_buff[fpos++]=nch;
	}
	fclose(fbig);

	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	strm.avail_in = 0;
	strm.next_in = Z_NULL;
	
	ret = deflateInit(&strm, 1);
	if (ret != Z_OK)printf("Ohh!\n");

	strm.next_in = in_buff;
	strm.avail_in = fpos;
	strm.next_out = out_buff;
	strm.avail_out = 999999;

	ret = deflate(&strm, Z_FINISH);

	have = 999999 - strmx.avail_out;
	
	printf("XRET=%d; LEN=%d; RES=%s\n",ret, have, in_buff);
	FILE * fbigo = fopen("tt.bin.gz","wb");
	fwrite(out_buff, 1, have, fbigo);
	fclose(fbigo);
}
