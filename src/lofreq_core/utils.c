/* -*- mode: c; tab-width: 4; c-basic-offset: 4;  indent-tabs-mode: nil -*- */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

/*********************************************************************
 *
 * Copyright (C) 2011, 2012 Genome Institute of Singapore
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 *********************************************************************/

int file_exists(char *fname) 
{
  /* from 
   * http://stackoverflow.com/questions/230062/whats-the-best-way-to-check-if-a-file-exists-in-c-cross-platform 
   */
  if (access(fname, F_OK) != -1) {
      return 1;
  } else {
      return 0;
  }
}

/* from http://www.anyexample.com/programming/c/how_to_load_file_into_memory_using_plain_ansi_c_language.xml
 *
 * returns file size (number of bytes) on success or negative number
 * on error
 */
int ae_load_file_to_memory(const char *filename, char **result) 
{ 
	int size = 0;
	FILE *f = fopen(filename, "rb");
	if (f == NULL) 
	{ 
		*result = NULL;
		return -1; // -1 means file opening fail 
	} 
	fseek(f, 0, SEEK_END);
	size = ftell(f);
	fseek(f, 0, SEEK_SET);
	*result = (char *)malloc(size+1);
    if (NULL==result)
        return -2;
	if (size != fread(*result, sizeof(char), size, f)) 
	{ 
		free(*result);
		return -3; /* -2 means file reading fail  */
	} 
	fclose(f);
	(*result)[size] = 0;
	return size;
}
