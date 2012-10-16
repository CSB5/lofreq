#!/usr/bin/env python
"""Ultra-sensitive detection of very rare / low frequency variants.

We recommended to recalibrate base-call quality-scores and realign
indels (both can be done with GATK).

NOTE:
* Reported p-values are not Bonferroni corrected.
* Use lofreq_filter.py for filtering the variant candidates produced
  here
* Variant positions with pvalue*sig-level<bonferroni are not
  reported.
"""


# Copyright (C) 2011, 2012 Genome Institute of Singapore
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.




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
USE_SCIPY = False
if USE_SCIPY:
    from scipy.stats import fisher_exact


#--- project specific imports
#
from lofreq import conf
from lofreq import utils
from lofreq import sam
from lofreq import em
from lofreq import qual
from lofreq import snp
from lofreq import simple_vcf
if not USE_SCIPY:
    from lofreq_ext import kt_fisher_exact


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011, 2012 Genome Institute of Singapore"
__license__ = "GPL2"


# global logger
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')




def process_pileup_line(pcol, lofreq_nq=None, lofreq_q=None, what='snp'):
    """Will use lofreq_nq and then lofreq_q in that order (if not
    None) to process a pileup column to call SNVs.

    return type is define by what, which can be snp or vcf    
    """

    assert lofreq_nq or lofreq_q
    assert what in ['snp', 'vcf']
    
    ret = []
    

    # return if we don't have a consensus (nothing to call against)
    if pcol.cons_base not in 'ACGT':
        LOG.info("Skipping col %d because of ambigious consensus %s" % (
                pcol.coord+1, pcol.cons_base))
        return ret


    # derive base counts and coverage
    #
    base_counts = pcol.get_all_base_counts(
        min_qual=3, keep_strand_info=False)
    # Even though NQ is quality agnostic, ignore the ones marked with
    # Illumina's Read Segment Indicator Q2, since those are not
    # supposed to be used according to Illumina
    del base_counts['N']
    # exit early if no bases found
    if sum(base_counts.values())==0:
        LOG.info("Zero coverage in col %d" % (pcol.coord+1))
        return ret


    # deal with "consensus variants"
    #
    # consensus snp, i.e. we determine consensus different to initial
    # pileup ref (ref_base). Always call against cons_base. If
    # ref_base is called then drop silently later.
    #
    cons_var_snp = None
    if pcol.cons_base != pcol.ref_base:
        # FIXME: should become vcf record immediately
        # SNP is deprecated and should be generated from vcf
        info_dict = dict()
        # reusing base_counts from above
        for (k, v) in base_counts.iteritems():
            info_dict["basecount-%s" % k] = v
        coverage = sum(base_counts.values())
        info_dict['coverage'] = coverage
        info_dict["type"] = 'consensus-var'
        cons_var_snp = snp.ExtSNP(pcol.coord, pcol.ref_base, pcol.cons_base,
                                  base_counts[pcol.cons_base]/coverage,
                                  info_dict)
        # add cons_var_snp to new_snps later, otherwise it gets
        # re-evaluated by lofreq-q

        
    new_snps = []
    # use lofreq-nq if activated
    if lofreq_nq:
        new_snps = lofreq_nq.call_snp_in_column(
            pcol.coord, base_counts, pcol.cons_base)
        for snpcall in new_snps:
            LOG.info("LoFreq-NQ SNV: %s." % (snpcall))
            
    # use lofreq-q if lofreq-nq was off or to re-evaluate its
    # predictions (overwriting the originals)
    #
    if lofreq_q:
        # if nq is off or if some snps were predicted, than
        # re-evaluate
        if not lofreq_nq or len(new_snps) != 0:
            base_qual_hist = pcol.get_base_and_qual_hist(
                keep_strand_info=False)
            del base_qual_hist['N']
            new_snps = lofreq_q.call_snp_in_column(
                pcol.coord, base_qual_hist, pcol.cons_base)
            for snpcall in new_snps:
                LOG.info("LoFreq-Q SNV: %s." % (snpcall))

                
    # ---
    #  predictions in all stages done. the rest is stupid formatting
    # ---

    
    for snpcall in new_snps:
        snpcall.info["type"] = 'low-freq-var'
        
    # merge in consensus variants (but drop if we called the ref_base
    # a snv)
    if cons_var_snp:
        new_snps = [s for s in new_snps
                    if s.variant != cons_var_snp.wildtype]            
        new_snps.append(cons_var_snp)


    ret = []    
    for snpcall in new_snps:
        # the calling routines don't know about chrom, so add it here
        snpcall.chrom  = pcol.chrom
        LOG.info("Final LoFreq SNV: %s." % (snpcall))

        # add phred quality
        #
        if snpcall.info['type'] == 'low-freq-var':
            snpcall.info['pvalue-phred'] = utils.prob_to_phredqual(
                snpcall.info['pvalue'])
        else:
            # consensus-vars
            snpcall.info['pvalue-phred'] = 'NA'

        # add strand-bias info
        #
        if lofreq_q:
            # Q
            ref_qf = lofreq_q.ign_bases_below_q
            var_qf = max(lofreq_q.ign_bases_below_q, 
                         lofreq_q.noncons_filter_qual)
        else:
            # NQ
            ref_qf = 0
            var_qf = 0

        ref_counts = pcol.get_counts_for_base(
            snpcall.wildtype, ref_qf, keep_strand_info=True)
        var_counts = pcol.get_counts_for_base(
            snpcall.variant, var_qf, keep_strand_info=True)
        # FIXME should become VCF record immediately
        # SNP is deprecated and should be generated from vcf
        add_strandbias_info(snpcall, ref_counts, var_counts)
        # report extra phred scaled pvalue


        if what == 'vcf':
            dp4 = (ref_counts[0], ref_counts[1], var_counts[0], var_counts[1])
            # FIXME any way to mark the cons-variation?
            vcf_record = gen_vcf_record(pcol, snpcall.wildtype, snpcall.variant, 
                                        snpcall.info['pvalue-phred'], snpcall.freq, dp4,
                                        snpcall.info['strandbias-pvalue-uncorr-phred'])
            ret.append(vcf_record)
        else:
            ret.append(snpcall)

    return ret



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
                      help="Log messages to file")

    plp_group = OptionGroup(parser, "Pileup Options", "")
    plp_group.add_option("-b", "--bam",
                      dest="fbam", # type="string|int|float"
                      help="BAM file containing your stringently mapped reads")
    plp_group.add_option("-f", "--reffa",
                      dest = "freffa", # type="string|int|float"
                      help = "Reference fasta file")
    choices = ['on', 'off', 'extended']
    plp_group.add_option("", "--baq",
                         dest = "baq",
                         choices = choices,
                         default=conf.DEFAULT_BAQ_SETTING,
                         help="Optional: BAQ setting for pileup."
                         " One of %s. Default %s"% (', '.join(choices), conf.DEFAULT_BAQ_SETTING))
    plp_group.add_option("-l", "--regions",
                         dest = "fregion_bed",
                         help = "Optional: bed file containing regions to limit analysis to."
                         " One of %s. Default %s"% (', '.join(choices), conf.DEFAULT_BAQ_SETTING))
    plp_group.add_option("-d", "--maxdepth",
                         dest = "max_depth",
                         type=int,
                         default = conf.DEFAULT_MAX_PLP_DEPTH,
                         help = "Optional: maximum depth. Default %d"% (conf.DEFAULT_MAX_PLP_DEPTH))
    parser.add_option_group(plp_group)

    
    output_group = OptionGroup(parser, "Output Options", "")
    output_group.add_option("-o", "--out",
                      dest="fsnp", default="-", # type="string|int|float"
                      help="Variant output file or '-' for stdout (default)")
    choices = ['snp', 'vcf']
    output_group.add_option("", "--format",
                      dest="outfmt", 
                      choices=choices,
                      default='snp',
                      help="Output format. One of: %s." % ', '.join(choices))
    parser.add_option_group(output_group)


    model_group = OptionGroup(parser, "Model Options", "")
    default = True
    model_group.add_option("", "--lofreq-q-on",
                      action="store_true", 
                      dest="lofreq_q_on",
                      default=default,
                      help="Activate quality-aware SNV calling (LoFreq-Q; default = %s)" % default)
    model_group.add_option("", "--lofreq-q-off",
                      action="store_false", 
                      dest="lofreq_q_on",
                      help="De-activate quality-aware SNV calling (LoFreq-Q)")
    
    default = False
    model_group.add_option("", "--lofreq-nq-on",
                      action="store_true",
                      dest="lofreq_nq_on",
                      default=default,
                      help="Activate quality-agnostic SNV calling (LoFreq-NQ; default = %s)" % default)
    model_group.add_option("", "--lofreq-nq-off",
                      action="store_false",
                      dest="lofreq_nq_on",
                      help="De-activate quality-agnostic SNV calling (LoFreq-NQ)")
    parser.add_option_group(model_group)
    
    filter_group = OptionGroup(parser, "Filtering Options", "")
    filter_group.add_option("-Q", "--ignore-bases-below-q",
                          dest="ign_bases_below_q", type="int",
                          default=conf.DEFAULT_IGN_BASES_BELOW_Q,
                          help="Remove any base below this base call quality threshold from pileup"
                          " (default: %d)" % conf.DEFAULT_IGN_BASES_BELOW_Q)
    filter_group.add_option("", "--bonf",
                            dest="bonf", 
                            default='auto',
                            help="Bonferroni correction factor."
                            " Set to an integer or 'auto'. 'auto' is default will use a stringent"
                            " and recommmended seqlen*3.")
    filter_group.add_option("-s", "--sig-level",
                      dest="sig_thresh", type="float",
                      default=conf.DEFAULT_SIG_THRESH,
                      help="p-value significance value"
                      " (default: %g)" % conf.DEFAULT_SIG_THRESH)
    parser.add_option_group(filter_group)


    em_group = OptionGroup(parser, "Advanced Options for quality-agnostic stage", "")
    em_group.add_option('', '--num-param', dest='em_num_param',
                        choices = ['4', '12'],
                        default=conf.DEFAULT_EM_NUM_PARAM,
                        help='Use 4- or 12-parameter model'
                        ' (default: %d)' % conf.DEFAULT_EM_NUM_PARAM)
    em_group.add_option('', '--error-prob-file', dest='em_error_prob_file',
                        help='Read EM error probs from this file (skips training).'
                        ' General format is: "<ref-base> <snp-base-1> <eprob-1> ...".'
                        ' 4-parameter model needs only one line: "N A <eprob> C <eprob> <eprob> G T <eprob>".'
                        ' 12-parameter model needs one line for each nucleotide.')
    #em_group.add_option('', '--convergence',
    #                    dest='conv_epsilon', type=float, default=conf.DEFAULT_CONVERGENCE_EPSILON,
    #                    help='Optional: difference value for convergence')
    parser.add_option_group(em_group)

    
    qual_group = OptionGroup(parser, "Advanced Options for quality-aware stage", "")
    qual_group.add_option("", "--noncons-default-qual",
                          dest="noncons_default_qual", type="int",
                          default=conf.NONCONS_DEFAULT_QUAL,
                          help="Base call quality used for non-consensus bases"
                          " (default: %d)" % conf.NONCONS_DEFAULT_QUAL)
    qual_group.add_option("", "--noncons-filter-qual",
                          dest="noncons_filter_qual", type="int",
                          default=conf.NONCONS_FILTER_QUAL,
                          help="Non-consensus bases below this threshold will be filtered"
                          " (default: %d)" % conf.NONCONS_FILTER_QUAL)
    parser.add_option_group(qual_group)


    # hidden options
    #
    parser.add_option("--test-sensitivity", help=SUPPRESS_HELP,
                      dest="test_sensitivity", action="store_true") 


    return parser



def test_sensitivity(mode):
    """Simple sensitivity test with a mockup column and uniform quality
    """

    bonf = 1
    
    assert mode in ['Q', 'NQ'], ("Mode must be one of Q or NQ")

    if mode == 'NQ':
        lofreqnq = em.EmBasedSNPCaller(
            conf.DEFAULT_EM_NUM_PARAM, bonf, conf.DEFAULT_SIG_THRESH)
        lofreqnq.set_default_error_probs()
    else:
        lofreqq = qual.QualBasedSNPCaller(
            conf.NONCONS_DEFAULT_QUAL, conf.NONCONS_FILTER_QUAL, 
            conf.DEFAULT_IGN_BASES_BELOW_Q,
            bonf, conf.DEFAULT_SIG_THRESH)

    print "Testing default LoFreq%s detection limits on fake pileup" \
        " with varying coverage and uniform quality / error probability" \
        " (sign.threshold = %f)" % (mode, conf.DEFAULT_SIG_THRESH)

    refbase = 'A'
    snpbase = 'G'
    coverage_range = [10, 50, 100, 500, 1000, 5000, 10000]
    #coverage_range = [10]
    quality_range = [20, 25, 30, 35, 40]
    #quality_range = [40]

    for q in quality_range:
        print "\tQ=%d" % q,
    print

    for cov in coverage_range:
        print "%d" % cov,
        for q in quality_range:

            num_noncons = 1

            if mode == 'Q':
                while [ True ]:
                    base_qual_hist = dict(zip(
                            ['A', 'C', 'G', 'T'],
                            [dict(), dict(), dict(), dict()]
                            ))
                    base_qual_hist[refbase][q] = cov-num_noncons
                    base_qual_hist[snpbase][q] = num_noncons
    
                    snps = lofreqq.call_snp_in_column(666, base_qual_hist, refbase)

                    if len(snps):
                        print "\t%d" % (num_noncons),
                        break
                    num_noncons += 1
                    if num_noncons == cov:
                        break

            else:
                # turn quality into uniform error probability
                prob = utils.phredqual_to_prob(q)
                lofreqnq.error_probs[refbase][snpbase] = prob

                while [ True ]:
                    base_counts = dict(zip(['A', 'C', 'G', 'T'], 4*[0]))
                    base_counts[refbase] = cov-num_noncons
                    base_counts[snpbase] = num_noncons
                    snps = lofreqnq.call_snp_in_column(666, base_counts, refbase)

                    if len(snps):
                        print "\t%d" % (num_noncons),
                        break
                    num_noncons += 1
                    if num_noncons == cov:
                        break    
        print

    

def add_strandbias_info(snpcall, ref_counts, var_counts):
    """strand bias test
    """

    try:
        if USE_SCIPY:
            # alternative : two-sided, less, greater
            # Which alternative hypothesis to the null hypothesis the test uses. Default is two-sided.
            (oddsratio, fisher_twotail_pvalue) = fisher_exact(
                [[ref_counts[0], ref_counts[1]],
                 [var_counts[0], var_counts[1]]])
        else:
            # expects two tuples as input
            # returns left, right, twotail pvalues
            (left_pvalue, right_pvalue, fisher_twotail_pvalue) = kt_fisher_exact(
                (ref_counts[0], ref_counts[1]),
                (var_counts[0], var_counts[1]))

    except ValueError:
        fisher_twotail_pvalue = -1

    # report extra phred scaled pvalue
    snpcall.info['strandbias-pvalue-uncorr'] = fisher_twotail_pvalue
    if fisher_twotail_pvalue == -1:
        snpcall.info['strandbias-pvalue-uncorr-phred'] = "NA"
    else:
        snpcall.info['strandbias-pvalue-uncorr-phred'] = utils.prob_to_phredqual(fisher_twotail_pvalue)

        

def gen_vcf_record(pcol, ref, alt, phredqual, freq, dp4=None, sb_phred=None):
    """
    Wrapper to simple_vcf.create_record
    """

    if dp4:
        assert len(dp4) == 4
        
    vcf_info_dict = {
    'AF':float("%.*f" % (len(str(pcol.coverage)), freq)),
    'DP':pcol.coverage, 
    }
    if dp4:
        vcf_info_dict['DP4'] = '%d,%d,%d,%d' % (dp4)
    if sb_phred != None:  # can be zero
        vcf_info_dict['SB'] = sb_phred
    if  phredqual == 'NA':
        phredqual = None
    vcf_record = simple_vcf.create_record(
        pcol.chrom, pcol.coord,
        None,
        ref, alt,
        phredqual,
        None,
        vcf_info_dict)
    
    return vcf_record
    



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
        test_sensitivity('NQ')
        test_sensitivity('Q')
        sys.exit(0)

        
    if opts.logfile:
        hdlr = logging.FileHandler(opts.logfile)
        formatter = logging.Formatter(
            '%(levelname)s [%(asctime)s]: %(message)s')
        hdlr.setFormatter(formatter)
        LOG.addHandler(hdlr)

    if opts.verbose:
        LOG.setLevel(logging.INFO)
        sam.LOG.setLevel(logging.INFO)
        em.LOG.setLevel(logging.INFO)
        qual.LOG.setLevel(logging.INFO)

    if opts.debug:
        LOG.setLevel(logging.DEBUG)
        sam.LOG.setLevel(logging.DEBUG)
        em.LOG.setLevel(logging.DEBUG)
        qual.LOG.setLevel(logging.DEBUG)

    if not opts.fsnp:
        parser.error("Variant output file argument missing.")
        sys.exit(1)
    if opts.fsnp != '-' and os.path.exists(opts.fsnp):
        LOG.fatal(
            "Cowardly refusing to overwrite already existing file '%s'.\n" % (
                opts.fsnp))
        sys.exit(1)

    if not opts.lofreq_nq_on and not opts.lofreq_q_on:
        parser.error("Nothing to do. LoFreq-Q and LoFreq-NQ both switched off.")
        sys.exit(1)


    # file check
    for (filename, descr, direction, mandatory) in [
            (opts.fsnp, "Output file", 'out', True),
            (opts.fbam, "BAM file", 'in', True),
            (opts.freffa, "Reference fasta file", 'in', True),
            (opts.fregion_bed, "Region BED-file", 'in', False)
            ]:

        if not mandatory and not filename:
            continue
            
        if not filename:
            parser.error("%s argument missing." % descr)
            sys.exit(1)
        if filename == '-':
            continue
        
        if direction == 'in' and not os.path.exists(filename):
            LOG.fatal(
                "file '%s' does not exist.\n" % filename)
            sys.exit(1)
        if direction == 'out' and os.path.exists(filename):
            LOG.fatal(
                "Refusing to overwrite existing file '%s'.\n" % filename)
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
            LOG.fatal("file '%s' does not exist.\n" % (
                opts.em_error_prob_file))
            sys.exit(1)

    # pileup
    #
    pileup_obj = sam.Pileup(opts.fbam, opts.freffa)
    pileup_gen = pileup_obj.generate_pileup(
        opts.baq, opts.max_depth, opts.fregion_bed)
    # a buffer for pushing back those pileup lines used up for training
    pileup_line_buffer = []
    

    # pvalue threshold and correction settings
    #
    sig_thresh = opts.sig_thresh
    if opts.bonf == "auto":
        bonf_factor = sam.auto_bonf_factor(
               opts.fbam, opts.fregion_bed)
        LOG.info("Will use an automatically determined Bonferroni"
                 " factor of %d." % (bonf_factor))
    else:
        try:
            bonf_factor = int(opts.bonf)
            assert opts.bonf > 0
        except:
            LOG.critical("Something's wrong with the Bonferroni factor you provided...")
            sys.exit(1)
    
    lofreq_nq = None
    if opts.lofreq_nq_on:
        lofreq_nq = em.EmBasedSNPCaller(
            num_param = em_num_param,
            bonf_factor= bonf_factor,
            sig_thresh = sig_thresh)
    
    lofreq_q = None
    if opts.lofreq_q_on:
        lofreq_q = qual.QualBasedSNPCaller(
            noncons_default_qual = noncons_default_qual,
            noncons_filter_qual = noncons_filter_qual,
            ign_bases_below_q =  ign_bases_below_q,
            bonf_factor = bonf_factor, sig_thresh = sig_thresh)
        
    LOG.info("Commandline (workdir %s): %s" % (
        os.getcwd(), ' '.join(sys.argv)))

    
    # ################################################################
    #
    # Stage 1 (LoFreq-NQ training): run EM training to get error
    # probabilities for each base turning into another bases (12
    # parameter model).
    #
    # ################################################################
    
    # Sample from total columns. Take first
    # conf.EM_TRAINING_SAMPLE_SIZE pileup columns (if possible given
    # start_pos and end_pos), with a non-ambigious reference and a
    # coverage of at least conf.EM_TRAINING_MIN_COVERAGE
    
    
    # Get pileup data for EM training. Need base-counts and cons-bases
    #
    if opts.lofreq_nq_on and not opts.em_error_prob_file:
        cons_seq = []
        base_counts = []
        LOG.info("Processing pileup for EM training")
        num_lines = 0
        for pcol in pileup_gen:
            # note: not all columns will be present in pileup
            num_lines += 1
            pileup_line_buffer.append(pcol)

            if pcol.cons_base not in 'ACGT':
                LOG.info("Skipping col %d because of ambigious consensus base %s" % (
                    pcol.coord+1, pcol.cons_base))
                continue
     
            # Even though NQ is quality agnostics, ignore the ones
            # marked with Illumina's Read Segment Indicator Q2, since
            # those are not supposed to be used according to Illumina
            col_base_counts = pcol.get_all_base_counts(
                min_qual=3, keep_strand_info=False)
            del col_base_counts['N']

            coverage = sum(col_base_counts.values())
            if coverage < conf.EM_TRAINING_MIN_COVERAGE:
                continue
     
            base_counts.append(col_base_counts)
            cons_seq.append(pcol.cons_base)
     
            if len(base_counts) >= conf.EM_TRAINING_SAMPLE_SIZE:
                break
        
        if num_lines == 0:
            LOG.fatal("Pileup was empty. Will exit now.")
            sys.exit(1)

        if len(base_counts) == 0:
            LOG.fatal(
                "No useable data data acquired from pileup for EM training. Exiting" % (
                    len(base_counts)))
            sys.exit(1)

        if len(base_counts) < conf.EM_TRAINING_SAMPLE_SIZE:
            LOG.warn(
                "Insufficient data (%d) acquired from pileup for EM training." % (
                    len(base_counts)))

        LOG.info("Using %d columns with an avg. coverage of %d for EM training " % (
            len(base_counts),
            sum([sum(c.values()) for c in base_counts])/len(base_counts)))
            
        lofreq_nq.em_training(base_counts, cons_seq)
        LOG.info("EM training completed.")

        # fprobs="schmock.error_prob."
        #LOG.critical("Saving provs to %s." % fprobs)
        #lofreq_nq.save_error_probs(fprobs)


    elif opts.em_error_prob_file:
        LOG.info(
            "Skipping EM training and using probs from %s instead." % (
                opts.em_error_prob_file))
        lofreq_nq.load_error_probs(opts.em_error_prob_file)


    # ################################################################
    #
    # Stage 2 & 3: Call SNPs based on EM training probabilities. Use
    # the quality based model (LoFreq-Q) on top of this.
    #
    # ################################################################

    if opts.fsnp == '-':
        fhout = sys.stdout
    else:
        fhout = open(opts.fsnp, 'w') 
        
    if opts.outfmt == 'snp':
        # No wan'a
        #snp.write_header(fhout)
        pass
    else:
        simple_vcf.write_header(fhout)


    num_lines = 0
    for pcol in itertools.chain(pileup_line_buffer, pileup_gen):
        # note: not all columns will be present in pileup
        num_lines += 1

        new_snps = process_pileup_line(
            pcol, lofreq_nq, lofreq_q, opts.outfmt)
        for snpcall in new_snps:
            if opts.outfmt == 'snp':
                snp.write_record(snpcall, fhout)
            else:
                simple_vcf.write_record(snpcall, fhout)                

    if num_lines == 0:
        LOG.fatal("Pileup was empty. Will exit now.")
        sys.exit(1)
        
    LOG.info("SNVs written to %s" % fhout.name)
    if fhout != sys.stdout:
        fhout.close()  

        
if __name__ == "__main__":
    main()
    LOG.info("Use lofreq_filter.py to post-process the just produced SNV calls.")    
    LOG.info("Successful program exit")

