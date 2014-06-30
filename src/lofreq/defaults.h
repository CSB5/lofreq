#ifndef LOFREQ_DEFAULTS_H
#define LOFREQ_DEFAULTS_H

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
#define DEFAULT_DEF_NM_QUAL 20
 
/* coverage thresholds */
#define DEFAULT_MIN_COV 1
#define DEFAULT_MAX_PLP_DEPTH 1000000

#define DEFAULT_BAQ_ON 1

/* make lofreq blind to anything below this value */
#define DEFAULT_MIN_PLP_BQ 3

#define DEFAULT_SIG 0.01

/* ---------------------------------------------------------------------- */

/* Four nucleotides, with one consensus, makes three
   non-consensus bases */
#define NUM_NONCONS_BASES 3

#define SNVCALL_USE_BAQ     0x01
#define SNVCALL_USE_MQ      0x02
#define SNVCALL_USE_SQ      0x04
#define SNVCALL_CONS_AS_REF 0x08

/* private tag for actual baq values: "l"ofreseq "b"ase-alignment */
#define BAQ_TAG "lb"


#ifndef MIN
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#endif
#ifndef MAX
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#endif

#endif
