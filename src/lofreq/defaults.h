#ifndef LOFREQ_DEFAULTS_H
#define LOFREQ_DEFAULTS_H
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
#define DEFAULT_MIN_PLP_IDQ 0

#define DEFAULT_SIG 0.01

/* ---------------------------------------------------------------------- */

/* Four nucleotides, with one consensus, makes three
   non-consensus bases */
#define NUM_NONCONS_BASES 3

#define VARCALL_USE_BAQ     1
#define VARCALL_USE_MQ      2
#define VARCALL_USE_SQ      4
/* indel alignment quality */
#define VARCALL_USE_IDAQ      8


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

/* base insertion and deletion qualities. GATK uses BI and BD. 
 * GATKs BI & BD: "are per-base quantities which estimate
 * the probability that the next base in the read was
 * mis-incorporated or mis-deleted (due to slippage, for
 * example)". See
 * http://www.broadinstitute.org/gatk/guide/article?id=44
 * and
 * http2://gatkforums.broadinstitute.org/discussion/1619/baserecalibratorprintreads-bd-and-bi-flags
 *
 */
#define BI_TAG "BI"
#define BD_TAG "BD"

#endif
