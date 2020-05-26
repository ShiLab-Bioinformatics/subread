CC_EXEC = gcc
OPT_LEVEL = 3

include ../makefile.version
include make.version

CCFLAGS = -mtune=core2 ${MACOS} -O${OPT_LEVEL} -Wall  -DMAKE_FOR_EXON  -D MAKE_STANDALONE  -D_FILE_OFFSET_BITS=64 -DSUBREAD_VERSION=\"${SUBREAD_VERSION}\" ${WARNING_LEVEL} ${MINGW32}
LDFLAGS =  -lpthread -lz -lm -O${OPT_LEVEL} -DMAKE_FOR_EXON -D MAKE_STANDALONE # -DREPORT_ALL_THE_BEST
CC = ${CC_EXEC}  ${CCFLAGS}  -fmessage-length=0  -ggdb

ALL_LIBS=LRMsorted-hashtable LRMbase-index LRMchro-event LRMhelper LRMseek-zlib LRMfile-io LRMhashtable
ALL_OBJECTS=$(addsuffix .o, ${ALL_LIBS})
ALL_H=$(addsuffix .h, ${ALL_LIBS})
ALL_C=$(addsuffix .c, ${ALL_LIBS})


all: $(ALL_OBJECTS) LRM

clean:
	rm -f *.o LRM

LRM: longread-mapping.c  ${ALL_OBJECTS} 
	${CC} -o LRM longread-mapping.c   ${ALL_OBJECTS}  ${LDFLAGS}

$(ALL_OBJECTS): $(ALL_C) $(ALL_H)
	$(CC) -o $@ -c $(subst .o,.c,$@)
	
