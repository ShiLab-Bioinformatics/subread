#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>



#define MEGA 1024llu*1024

int main(){
	int fd = open("/usr/local/work/liao/arena/del4.mem", O_TRUNC | O_CREAT|O_WRONLY , 0600);
	long long int x;

	for(x=0; x<100*MEGA; x++){
		write(fd, &x, 4);
	}

	close(fd);

	fd = open("/usr/local/work/liao/arena/del4.mem", O_RDWR);

	void * fd_ptr = mmap(NULL, 400*MEGA, PROT_READ|PROT_WRITE, MAP_PRIVATE, fd, 0);
	assert(fd_ptr != MAP_FAILED);
	printf("MEMPTR = %08X\n", (unsigned int)(fd_ptr));

	int * int_ptr = (int *)fd_ptr;

	for(x=0; x<100*MEGA; x+=456){
		int myint = int_ptr[x];
		//printf("MYI=%d\n", myint);
	}
	for(x=0; x<100*MEGA; x+=2){
		int_ptr[x]=0x0a;
	}

	printf("MEMORY PREPARIED\n");
	sleep(100);
	munmap(fd_ptr, 400*MEGA);
	close(fd);
}
