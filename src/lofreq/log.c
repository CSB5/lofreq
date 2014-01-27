/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*- */

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

