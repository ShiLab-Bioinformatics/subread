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
#include <stdarg.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include "subread.h"
#include "sublog.h"

#define MINIMUM_LOG_LEVEL SUBLOG_LEVEL_INFO

void remove_ESC_effects(char * txt)
{
	int x1;
	int ocur = 0;
	int state = 0;
	int trimmed = 0;

//	return 0;

	for(x1=0;x1<1199;x1++)
	{
		if(!txt[x1])break;
		if((state == 0) && (txt[x1]==CHAR_ESC))
		{
			state = 1;
			trimmed = 1;
			continue;
		}
		if(state == 0){
				if(x1>ocur)
						txt[ocur] = txt[x1];
				ocur++;
		}
		
		if((state == 1) && (txt[x1]=='m'))
			state = 0;
	}
	if(trimmed)
		txt[ocur]=0;
}


int is_ESC_removed()
{
	#if defined(MAKE_STANDALONE) || defined(RUNNING_ENV)
	return !isatty(fileno(stderr));
	#else
	return 1;
	#endif

}


void sublog_printf(int stage, int level, const char * pattern, ...)
{
	va_list args;
	va_start(args , pattern);
	if(level<MINIMUM_LOG_LEVEL) return;

	int to_remove_ESC = 1;
	#ifndef __MINGW32__
	#if defined(MAKE_STANDALONE) || defined(RUNNING_ENV)
	to_remove_ESC = is_ESC_removed();
	#endif
	#endif	

	if(to_remove_ESC)
	{
		char * vsbuf=malloc(1200);

		vsnprintf(vsbuf, 1199, pattern , args);
		remove_ESC_effects(vsbuf);

		SUBREADprintf("%s\n",vsbuf);

		free(vsbuf);
	}
	else
	{
		#if defined(MAKE_STANDALONE) || defined(RUNNING_ENV)
		vfprintf(stderr, pattern , args);
		fputs("\n", stderr);

		fflush(stderr);
		#endif
	}
	va_end(args);
}

void sublog_fwrite(int stage, int level, const char * pattern, ...)
{
	va_list args;
	va_start(args , pattern);

	if(level<MINIMUM_LOG_LEVEL) return;

	int to_remove_ESC = 1;
	#if defined(MAKE_STANDALONE) || defined(RUNNING_ENV)
	to_remove_ESC = is_ESC_removed();
	#endif

	if(to_remove_ESC)
	{
		char * vsbuf=malloc(1200);

		vsnprintf(vsbuf, 1199, pattern , args);
		remove_ESC_effects(vsbuf);
		if(strlen(vsbuf)>0)
			SUBREADprintf("%s",vsbuf);
		free(vsbuf);
	}
	else
	{
		#if defined(MAKE_STANDALONE) || defined(RUNNING_ENV)
		vfprintf(stderr, pattern , args);

		fflush(stderr);
		#endif
	}
	va_end(args);

}

int sambamout_fprintf(FILE * fp, const char * pattern, ...)
{
	int ret;
	va_list args;
	va_start(args , pattern);

	//printf("FP=%llu\n", (long long)fp);
	#ifdef MAKE_STANDALONE
	if(fp == NULL) fp = stdout;
	#endif
	assert(fp);
	
	ret = vfprintf(fp, pattern , args);
	va_end(args);
	return ret;
}
