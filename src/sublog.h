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
  
  
#ifndef _SUBLOG_H_
#define _SUBLOG_H_

#define SUBLOG_STAGE_DEV1		16
#define SUBLOG_STAGE_DEV_ALPHA	256
#define SUBLOG_STAGE_DEV_BETA	4096
#define SUBLOG_STAGE_DEV_RC		65536
#define SUBLOG_STAGE_RELEASED	1048576	

#define SUBLOG_LEVEL_NIL		10
#define SUBLOG_LEVEL_DEBUG		20
#define SUBLOG_LEVEL_DETAILS	110
#define SUBLOG_LEVEL_INFO		120
#define SUBLOG_LEVEL_ABNORMAL	210
#define SUBLOG_LEVEL_WARNING	220
#define SUBLOG_LEVEL_ERROR		310
#define SUBLOG_LEVEL_FATAL		900

void sublog_printf(int stage, int level, const char * pattern, ...);
void sublog_fwrite(int stage, int level, const char * pattern, ...);
int sambamout_fprintf(FILE * fp, const char * pattern, ...);

#endif
