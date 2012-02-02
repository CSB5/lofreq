#!/usr/bin/env python
"""Detection of very rare / low frequency variants.

In a first stage an expectation-maximization will be used to get error
probabilities for each possible base to base conversion. These can
then be used to predict a first set of SNPs. This set can be further
filtered using a quality based model.
"""


#--- standard library imports
#
from __future__ import division
import sys
import logging
import os
# optparse deprecated from Python 2.7 on
from optparse import OptionParser, OptionGroup, SUPPRESS_HELP
import itertools

#--- third-party imports
#
# /

#--- project specific imports
#
from lofreq import pileup
from lofreq import snp
from lofreq import em
from lofreq import qual
from lofreq.utils import count_bases


__author__ = "Andreas Wilm"
__version__ = "0.1-2012-01-30"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = ""
__license__ = ""
__credits__ = [""]
__status__ = ""


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

# significance threshold
DEFAULT_SIG_THRESH = 0.05
# quality filters/settings for quality based method
NONCONS_DEFAULT_QUAL = 20
NONCONS_FILTER_QUAL = 20
DEFAULT_IGN_BASES_BELOW_Q = 3
# EM settings
DEFAULT_EM_NUM_PARAM = 12
# FIXME: make user args
EM_TRAINING_SAMPLE_SIZE = 10000
EM_TRAINING_MIN_COVERAGE = 10


def read_exclude_pos_file(fexclude):
    """Parse file containing ranges of positions to exclude and return
    positions as list.
    """

    excl_pos = []
    fhandle = open(fexclude, 'r')
    for line in fhandle:
        if line.startswith('#'):
            continue
        if len(line.strip()) == 0:
            continue

        start = int(line.split()[0])
        end = int(line.split()[1])
        # ignore the rest

        assert start < end, ("Invalid position found in %s" % fexclude)

        excl_pos.extend(range(start, end))
    fhandle.close()

    return excl_pos



def write_snps(snp_list, fname, append=False):
    """Writes SNPs to file. Frontened to snp.write_snp_file
    """

    if append:
        mode = 'a'
        write_header = False
    else:
        mode = 'w'
        write_header = True

    LOG.info("Writing %d SNPs to %s" % (len(snp_list), fname))
    fhandle = open(fname, mode)
    snp.write_snp_file(fhandle, snp_list, write_header)
    fhandle.close()

    


def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog [Options]\n" \
            + "\n" + __doc__
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true", dest="debug",
                      help="enable debugging")
    parser.add_option("", "--log",
                      dest="logfile",
                      help="Optional: log to file")

    parser.add_option("-i", "--input",
                      dest="fpileup", # type="string|int|float"
                      default='-',
                      help="Pileup input. Will read from stdin (default) if '-' or not set at all."
                      " Tip: Use '-d 100000' to prevent sample depth filtering by samtools."
                      " Also consider using -B/-E to switch off/influence BAQ computation")
    parser.add_option("-e", "--exclude",
                      dest="fexclude", # type="string|int|float"
                      help="Optional: Exclude positions listed in this file"
                      " format is: start end [comment ...]"
                      " , with zero-based, half-open coordinates")
    parser.add_option("-o", "--out",
                      dest="fsnp", # type="string|int|float"
                      help="SNP output file")
    parser.add_option("", "--append",
                      dest="append", # type="string|int|float"
                      action="store_true",
                      help="Append to SNP file")
    
    parser.add_option("", "--em-only",
                      action="store_true", dest="skip_qual_stage",
                      help="Skip quality based stage,"
                      " i.e. terminate after first stage (EM)")
    parser.add_option("", "--qual-only",
                      action="store_true", dest="skip_em_stage",
                      help="Skip EM stage,"
                      " i.e. only use quality based predictions")
    
    parser.add_option("-b", "--bonf",
                      dest="bonf", type="int",
                      help="Bonferroni correction factor"
                      " (best to set to (seqlen-numexclpos)*3")
    parser.add_option("-s", "--sig-level",
                      dest="sig_thresh", type="float",
                      default=DEFAULT_SIG_THRESH,
                      help="Optional: p-value significance values"
                      " (default: %g)" % DEFAULT_SIG_THRESH)
    parser.add_option("-Q", "--ignore-bases-below-q",
                          dest="ign_bases_below_q", type="int",
                          default=DEFAULT_IGN_BASES_BELOW_Q,
                          help="Optional: remove any base below this base call quality threshold from pileup"
                          " (default: %d)" % DEFAULT_IGN_BASES_BELOW_Q)


    em_group = OptionGroup(parser, "Advanced Options for EM-based Stage", "")
    em_group.add_option('', '--num-param', dest='em_num_param',
                        choices = ['4', '12'],
                        default=DEFAULT_EM_NUM_PARAM,
                        help='Use 4- or 12-parameter model'
                        ' (default: %d)' % DEFAULT_EM_NUM_PARAM)
    em_group.add_option('', '--error-prob-file', dest='em_error_prob_file',
                        help='Read EM error probs from this file (skips training).'
                        ' General format is: "<ref-base> <snp-base-1> <eprob-1> ...".'
                        ' 4-parameter model needs only one line: "N A <eprob> C <eprob> <eprob> G T <eprob>".'
                        ' 12-parameter model needs one line for each nucleotide.')
    #em_group.add_option('', '--convergence',
    #                    dest='conv_epsilon', type=float, default=DEFAULT_CONVERGENCE_EPSILON,
    #                    help='Optional: difference value for convergence')
    parser.add_option_group(em_group)

    
    qual_group = OptionGroup(parser, "Advanced Options for Quality-Based Stage", "")
    qual_group.add_option("", "--noncons-default-qual",
                          dest="noncons_default_qual", type="int",
                          default=NONCONS_DEFAULT_QUAL,
                          help="Optional: Base call quality used for non-consensus bases"
                          " (default: %d)" % NONCONS_DEFAULT_QUAL)
    qual_group.add_option("", "--noncons-filter-qual",
                          dest="noncons_filter_qual", type="int",
                          default=NONCONS_FILTER_QUAL,
                          help="Optional: Non-consensus bases below this threshold will be filtered"
                          " (default: %d)" % NONCONS_FILTER_QUAL)
    parser.add_option_group(qual_group)

    parser.add_option("--test-sensitivity", help=SUPPRESS_HELP,
                      dest="test_sensitivity", action="store_true") 

    # no need for coverage filter option. snps on low coverage regions
    # are easily called and easily filtered downstrea.

    return parser


def test_sensitivity():
    """FIXME
    """
    from lofreq import utils

    bonf = 1

    lofreqnq = em.EmBasedSNPCaller(DEFAULT_EM_NUM_PARAM, bonf, DEFAULT_SIG_THRESH)
    lofreqnq.set_default_error_probs()

    lofreqq = qual.QualBasedSNPCaller(
        NONCONS_DEFAULT_QUAL, NONCONS_FILTER_QUAL, bonf, DEFAULT_SIG_THRESH)

    print "Testing default LoFreqQ/LofreqNQ detection limits on fake pileup" \
        " with varying coverage and uniform quality / error probability" \
        " (sign.threshold = %f)" % DEFAULT_SIG_THRESH

    refbase = 'A'
    snpbase = 'G'
    coverage_range = [10, 100, 1000, 10000]
    quality_range = [20, 25, 30, 35, 40]

    for q in quality_range:
        print "\tQ=%d" % q,
    print


    for cov in coverage_range:
        print "%d" % cov,
        for q in quality_range:

            num_noncons = 1
            while [ True ]:
                bases = (cov-num_noncons)*refbase + num_noncons*snpbase
                quals = len(bases)*[q]
                snps = lofreqq.call_snp_in_column(
                    666, bases, quals, refbase)
                if len(snps):
                    print "\t%d" % (num_noncons),
                    break
                num_noncons += 1
                if num_noncons == cov:
                    break

            num_noncons = 1
            # turn quality into uniform error probability
            prob = utils.phredqual_to_prob(q)
            lofreqnq.error_probs[refbase][snpbase] = prob/3.0

            while [ True ]:
                bases = (cov-num_noncons)*refbase + num_noncons*snpbase
                snps = lofreqnq.call_snp_in_column(
                    666, bases, refbase)
                if len(snps):
                    print "/%d" % (num_noncons),
                    break
                num_noncons += 1
                if num_noncons == cov:
                    break

        print


def main():
    """
    The main function
    """

    parser = cmdline_parser()
    (opts, args) = parser.parse_args()

    if len(args):
        parser.error("Unrecognized arguments found: %s." % (
            ' '.join(args)))
        sys.exit(1)

    if opts.test_sensitivity:
        test_sensitivity()
        sys.exit(0)

    if opts.logfile:
        hdlr = logging.FileHandler(opts.logfile)
        formatter = logging.Formatter(
            '%(levelname)s [%(asctime)s]: %(message)s')
        hdlr.setFormatter(formatter)
        LOG.addHandler(hdlr)

    if opts.verbose:
        LOG.setLevel(logging.INFO)
        pileup.LOG.setLevel(logging.INFO)
        em.LOG.setLevel(logging.INFO)
        qual.LOG.setLevel(logging.INFO)

    if opts.debug:
        LOG.setLevel(logging.DEBUG)
        pileup.LOG.setLevel(logging.DEBUG)
        em.LOG.setLevel(logging.DEBUG)
        qual.LOG.setLevel(logging.DEBUG)

    if not opts.fsnp:
        parser.error("SNP output file argument missing.")
        sys.exit(1)
    if os.path.exists(opts.fsnp) and not opts.append:
        LOG.fatal(
            "Cowardly refusing to overwrite already existing file '%s'.\n" % (
                opts.fsnp))
        sys.exit(1)

    if opts.fexclude and not os.path.exists(opts.fexclude):
        LOG.fatal("file '%s' does not exist.\n" % (opts.fexclude))
        sys.exit(1)

    # quality filter settings (for quality based method)
    #
    noncons_filter_qual = opts.noncons_filter_qual
    noncons_default_qual = opts.noncons_default_qual
    ign_bases_below_q = opts.ign_bases_below_q

    # em setting
    #
    if opts.em_num_param:
        em_num_param = int(opts.em_num_param)
    if opts.em_error_prob_file:
        if not os.path.exists(opts.em_error_prob_file):
            LOG.fatal("file '%s' does not exist.\n" % (opts.em_error_prob_file))
            sys.exit(1)

    if opts.fpileup:
        if opts.fpileup == '-':
            pileup_fhandle = sys.stdin
        else:
            if not os.path.exists(opts.fpileup):
                LOG.fatal("file '%s' does not exist.\n" % (opts.fpileup))
                sys.exit(1)
            else:
                pileup_fhandle = open(opts.fpileup, 'r')
    # a buffer for pushing back those pileup lines used up for training
    pileup_line_buffer = []


    # exclude positions
    #
    excl_pos = []
    if opts.fexclude:
        excl_pos = read_exclude_pos_file(opts.fexclude)
        LOG.info("Ignoring %d positions found in %s" % (
            len(excl_pos), opts.fexclude))
        LOG.debug("DEBUG: excl_pos = %s" % excl_pos)

    # pvalue threshold and correction settings (needs exclude
    # positions)
    #
    if not opts.bonf:
        LOG.fatal("Missing argument: Bonferroni factor")
        sys.exit(1)
    bonf_factor = opts.bonf

    sig_thresh = opts.sig_thresh

    LOG.info("Commandline (workdir %s): %s" % (os.getcwd(), ' '.join(sys.argv)))

    snpcaller_em = em.EmBasedSNPCaller(
        num_param = em_num_param,
        bonf_factor= bonf_factor,
        sig_thresh = sig_thresh)
    
    snpcaller_qual = qual.QualBasedSNPCaller(
        noncons_default_qual = noncons_default_qual,
        noncons_filter_qual = noncons_filter_qual,
        bonf_factor = bonf_factor,
        sig_thresh = sig_thresh)
    
    
    # ################################################################
    #
    # Stage 1: run EM training to get error probabilities for each
    # base turning into another bases (12 parameter model).
    #
    # ################################################################
    
    # Sample from total columns minus exclude positions. Take first
    # EM_TRAINING_SAMPLE_SIZE pileup columns (if possible given
    # start_pos and end_pos), with a non-ambigious reference and a
    # coverage of at least EM_TRAINING_MIN_COVERAGE
    
    
    # Get pileup data for EM training. Need base-counts and cons-bases
    #
    # FIXME how does this behave if we have a perfect pileup (no errors?)
    #
    if not opts.skip_em_stage and not opts.em_error_prob_file:
        cons_seq = []
        base_counts = []
        LOG.info("Processing pileup for EM training")
        num_lines = 0
        for line in pileup_fhandle:
            num_lines += 1
            pcol = pileup.PileupColumn(line)
            pcol.rem_ambiguities()
            pcol.rem_bases_below_qual(ign_bases_below_q)
            pileup_line_buffer.append(line)
            
            # note: not all columns will be present in pileup
            
            if pcol.coord in excl_pos:
                LOG.debug("Skipping col %d because of exclusion" % (pcol.coord+1))
                continue
            if pcol.ref_base not in 'ACGT':
                LOG.debug("Skipping col %d because of amibigous reference base %s" % (
                    pcol.coord+1, pcol.ref_base))
                continue
     
            (col_base_counts, dummy_cons_base_est) = count_bases(pcol.read_bases)
            if col_base_counts.has_key('N'):
                del col_base_counts['N']
     
            if sum(col_base_counts.values()) < EM_TRAINING_MIN_COVERAGE:
                continue
     
            base_counts.append(col_base_counts)
            cons_seq.append(pcol.ref_base)
     
            if len(base_counts) > EM_TRAINING_SAMPLE_SIZE:
                break
        
        if num_lines == 0:
            LOG.fatal("Pileup was empty. Will exit now.")
            sys.exit(1)
        if len(base_counts) < EM_TRAINING_SAMPLE_SIZE:
            LOG.warn("Insufficient data acquired from pileup for EM training.")
        LOG.info("Using %d columns with an avg. coverage of %d for EM training " % (
            len(base_counts),
            sum([sum(c.values()) for c in base_counts])/len(base_counts)))
            
        snpcaller_em.em_training(base_counts, cons_seq)
        LOG.info("EM training completed.")

        #LOG.critical("Saving provs to schmock.error_prob.")
        #snpcaller_em.save_error_probs("schmock.error_prob")


    elif opts.em_error_prob_file:
        LOG.info("Skipping EM training and using probs from %s instead." % (opts.em_error_prob_file))
        snpcaller_em.load_error_probs(opts.em_error_prob_file)



    # ################################################################
    #
    # Stage 2 & 3: Call SNPs based on EM training probabilities. Use
    # the quality based model on top of this.
    #
    # ################################################################

    snp_list = []
    LOG.info("Processing pileup for SNP calls")
    num_lines = 0
    num_ambigious_ref = 0
    for line in itertools.chain(pileup_line_buffer, pileup_fhandle):
        num_lines += 1
        # note: pileup_column_generator will ignore empty columns, i.e
        # it might skip some
        pcol = pileup.PileupColumn(line)
        #LOG.critical("Before filtering: bases = %s" % [(pcol.read_bases.count(b) ,b) for b in set(pcol.read_bases)])
        #LOG.critical("Before filtering: quals = %s" % [(pcol.base_quals.count(q), q) for q in sorted(set(pcol.base_quals))])
        pcol.rem_ambiguities()
        pcol.rem_bases_below_qual(ign_bases_below_q)
        #LOG.critical("After filtering: bases = %s" % [(pcol.read_bases.count(b), b) for b in set(pcol.read_bases)])
        #LOG.critical("After filtering: quals = %s" % [(pcol.base_quals.count(q), q) for q in sorted(set(pcol.base_quals))])

        if (pcol.coord+1) % 100000 == 0:
            LOG.info("Calling SNPs in column %d" % (pcol.coord+1))
            
        if pcol.coord in excl_pos:
            LOG.debug("Skipping col %d because of exclusion" % (pcol.coord+1))
            continue
        
        if pcol.ref_base not in 'ACGT':
            LOG.debug("Skipping col %d because of amibigous consensus %s" % (
                pcol.coord+1, pcol.ref_base))
            num_ambigious_ref += 1
            continue

        if not opts.skip_em_stage:
            new_snps = snpcaller_em.call_snp_in_column(
                pcol.coord, pcol.read_bases, pcol.ref_base)

        # If a SNP was predicted by EM or if EM was skipped, then
        # predict SNPs with the quality based method (possibly
        # overwriting old SNPs). Do this here and not in a separate
        # stage to avoid having to call mpileup yet again.
        #
        if not opts.skip_qual_stage:
            if opts.skip_em_stage or len(new_snps):
                new_snps = snpcaller_qual.call_snp_in_column(
                    pcol.coord, pcol.read_bases, pcol.base_quals, pcol.ref_base)

        if len(new_snps):
            snp_list.extend(new_snps)
    if num_lines == 0:
        LOG.fatal("Pileup was empty. Will exit now.")
        sys.exit(1)
    if num_ambigious_ref > 0:
        LOG.warn("%d positions skipped, because of amibigious reference in pileup." % num_ambigious_ref)
        
    write_snps(snp_list, opts.fsnp, opts.append)
    


if __name__ == "__main__":
    
    main()
    LOG.warn("IMPROVEMENT Write SNPs immediately.")
    LOG.warn("IMPROVEMENT Add support for vcf-output: quick and dirty or https://github.com/jamescasbon/PyVCF")
    LOG.warn("IMPROVEMENT Return pvalues as -log(pvalue)")
    LOG.info("Successful program exit")
