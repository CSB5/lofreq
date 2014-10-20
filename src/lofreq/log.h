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
