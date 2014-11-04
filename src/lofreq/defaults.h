#ifndef LOFREQ_DEFAULTS_H
#define LOFREQ_DEFAULTS_H
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


#define SANGER_PHRED_MAX 93

/* mapping quality filters: applied to all reads. don't set too high as
 * this is a mapper dependent value
 * in case of BWA it's also dependent on the alignment command used.
 */
#define DEFAULT_MIN_MQ 0
#define DEFAULT_MAX_MQ 255

/* minimum base quality of any base below which they are skipped.
   note: GATK doesn't recalibrate BQ <=5 */
#define DEFAULT_MIN_BQ 6
/* minimum base quality for alt bases: below and they are skipped */
#define DEFAULT_MIN_ALT_BQ 6
#define DEFAULT_DEF_ALT_BQ 0
/* -1: ref median, 0: keep original, >0: replace with this value */

#define DEFAULT_MIN_JQ 0
/* minimum merged quality for alt bases  */
#define DEFAULT_MIN_ALT_JQ 0
#define DEFAULT_DEF_ALT_JQ 0
/* -1: ref median, 0: keep original, >0: replace with this value */

/* non match quality for source qual */
#define DEFAULT_DEF_NM_QUAL -1
 
/* coverage thresholds */
#define DEFAULT_MIN_COV 1
#define DEFAULT_MAX_PLP_DEPTH 1000000

#define DEFAULT_BAQ_ON 1

/* make lofreq blind to anything below this value */
#define DEFAULT_MIN_PLP_BQ 3
#define DEFAULT_MIN_PLP_IDQ 10

#define DEFAULT_SIG 0.01

/* ---------------------------------------------------------------------- */

/* Four nucleotides, with one consensus, makes three
   non-consensus bases */
#define NUM_NONCONS_BASES 3

#define SNVCALL_USE_BAQ     1
#define SNVCALL_USE_MQ      2
#define SNVCALL_USE_SQ      4
#define SNVCALL_CONS_AS_REF 8
/* indel alignment quality */
#define SNVCALL_USE_IDAQ      16


/* private tag for actual baq values: "l"ofreseq "b"ase-alignment */
#define BAQ_TAG "lb"


#ifndef MIN
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#endif
#ifndef MAX
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#endif


#define AI_TAG "ai"
#define AD_TAG "ad"

#endif
