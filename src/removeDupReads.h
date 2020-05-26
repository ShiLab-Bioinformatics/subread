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
  
  
#ifndef __REP_READ_REMO_H__
#define __REP_READ_REMO_H__
#include "subread.h"


typedef struct{
	unsigned short * repeating_times;
} read_voting_table_t;


// Giving temp_location 'NULL' makes the function to use the current directory to store temporary files.
// Giving read_count '0' makas the function use its default value: read_count = 400,000,000, namely ~50MB memory is used to store the selection table (other parts of the program may use more memory)
// Note that this function generates ( 4 * all_reads ) bytes of temporary files in the current director or the directory specified in temp_location.

int repeated_read_removal(char * in_SAM_file, int threshold, char * out_SAM_file, char * temp_location, int threads);

#endif
