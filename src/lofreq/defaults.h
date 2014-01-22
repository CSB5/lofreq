#ifndef DEFAULTS_H
#define DEFAULTS_H

/*
 * This file contains default settings of values that can be set by user on command line 
 *
 */


/* mapping quality filters: applied to all reads. don't set too high as
 * this is a mapper dependent value
 * in case of BWA it's also dependent on the alignment command used.
 */
#define DEFAULT_MIN_MQ 13
#define DEFAULT_MAX_MQ 255

/* minimum base quality of any base. note: GATK doesn't recalibrate BQ <=5 */
#define DEFAULT_MIN_BQ 6

/* minimum base quality for alt bases */
#define DEFAULT_MIN_ALT_BQ 20
#define DEFAULT_DEF_ALT_BQ 20
/* FIXME make arg FIXME might affect on src qual untested */
#define DEFAULT_MIN_ALT_MERGEDQ 20

/* non match quality for source qual */
#define DEFAULT_DEF_NM_QUAL 20
 
/* coverage thresholds */
#define DEFAULT_MIN_COV 1
#define DEFAULT_MAX_PLP_DEPTH 1000000

#define DEFAULT_BAQ_ON 1

#endif
