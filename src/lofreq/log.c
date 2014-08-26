/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */
/*********************************************************************
 *
 * Copyright (C) 2011-2014 Genome Institute of Singapore
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

#include "log.h"

int debug = 0;
int verbose = 0;

/* Taken from the Linux kernel source and slightly modified.
 */
int
vout(FILE *stream, const char *fmt, ...)
{                
     va_list args;
     int rc;
     
     va_start(args, fmt);
     rc = vfprintf(stream, fmt, args);
     va_end(args);
     return rc;
}

