/*********************************************************************
* The MIT License (MIT)
* 
* Copyright (c) 2013,2014 Genome Institute of Singapore
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation files
* (the "Software"), to deal in the Software without restriction,
* including without limitation the rights to use, copy, modify, merge,
* publish, distribute, sublicense, and/or sell copies of the Software,
* and to permit persons to whom the Software is furnished to do so,
* subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
* BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
* ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
* CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*
************************************************************************/

#ifndef LOG_H
#define LOG_H

#include <stdarg.h>
#include <stdio.h>

extern int debug;
extern int verbose;

int
vout(FILE *stream, const char *fmt, ...);

/* print only if debug is true*/
#define LOG_DEBUG(fmt, args...)    {if (debug) {(void)vout(stderr, "DEBUG(%s|%s): " fmt, __FILE__, __FUNCTION__, ## args);}}
/* print only if verbose is true*/
#define LOG_VERBOSE(fmt, args...)  {if (verbose) {(void)vout(stderr, fmt, ## args);}}
/* always warn to stderr */
#define LOG_WARN(fmt, args...)     (void)vout(stderr, "WARNING(%s|%s): " fmt, __FILE__, __FUNCTION__, ## args)
/* always print errors to stderr*/
#define LOG_ERROR(fmt, args...)    (void)vout(stderr, "ERROR(%s|%s:%d): " fmt, __FILE__, __FUNCTION__, __LINE__, ## args)
/* always print errors to stderr*/
#define LOG_FATAL(fmt, args...)    (void)vout(stderr, "FATAL(%s|%s:%d): " fmt, __FILE__, __FUNCTION__, __LINE__, ## args)
/* always print fixme's */
#define LOG_FIXME(fmt, args...)    (void)vout(stderr, "FIXME(%s|%s:%d): " fmt, __FILE__, __FUNCTION__, __LINE__, ## args)
/* always print notes */
#define LOG_NOTE(fmt, args...)     (void)vout(stderr, "NOTE(%s|%s:%d): " fmt, __FILE__, __FUNCTION__, __LINE__, ## args)

#endif
