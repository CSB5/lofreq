#!/usr/bin/env python
"""Plots detection limits for lofreq_snpcaller.py
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
from optparse import OptionParser, OptionGroup
from math import log

#--- third-party imports
#
#/

#--- project specific imports
#
from lofreq_ext import snpcaller_qual
from lofreq import conf
from lofreq import sam
from lofreq import qual


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011, 2012 Genome Institute of Singapore"
__license__ = "GPL2"


# global logger
LOG = logging.getLogger("")

logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')




def det_lim_for_pileup_line(pcol, lofreq_q):
    """Determines and reports the approximate SNV detection limit
    (range tuple) based on reference qualities in this pileup column
    
    Logic mainly stolen from process_pileup_line() and qual.py
    WARNING: code duplication/
    """

   
    # return if we don't have a consensus (nothing to call against)
    if pcol.cons_base not in 'ACGT':
        LOG.info("Skipping col %d because of ambigious consensus %s" % (
                pcol.coord+1, pcol.cons_base))
        return -1


    # derive base counts and coverage
    #
    base_counts = pcol.get_all_base_counts(
        min_qual=3, keep_strand_info=False)
    del base_counts['N']
    # exit early if no bases found
    if sum(base_counts.values())==0:
        LOG.info("Zero coverage in col %d" % (pcol.coord+1))
        return -1


    base_qual_hist = pcol.get_base_and_qual_hist(
        keep_strand_info=False)
    del base_qual_hist['N']

    #
    # The logic below is borrowed from qual.py
    #
    
    ref_base = pcol.cons_base
    noncons_bases = [b for b in base_qual_hist.keys() 
                     if b != ref_base]

    # only need ref_base quals and the number of noncons_bases (from
    # which we construct a fake qual list). number of noncons_bases
    # and therefore fake qual list should increase from 1 and go back
    # if exceeded to determine exact value


    # get list of consensus qualities but remove if quality is
    # below ign_bases_below_q
    #
    cons_quals = []
    for (q, c) in base_qual_hist[ref_base].iteritems():
        if q < lofreq_q.ign_bases_below_q:
            continue
        cons_quals.extend([q]*c)
    cons_count = len(cons_quals)

    step = int(log(cons_count))# FIXME
    if step == 0:
        step = 1
    last_step = 1
    for i in xrange(1, cons_count, step):
        noncons_counts = (i, 0, 0)
        base_quals = list(cons_quals) # copy
        base_quals.extend([lofreq_q.noncons_default_qual] * i)
        pvalues = snpcaller_qual(sorted(base_quals), noncons_counts,
                                 lofreq_q.bonf_factor, lofreq_q.sig_thresh)
        # pvalues only reported if pvalue * bonf < sig
        for (base, count, pvalue) in zip(noncons_bases, noncons_counts, pvalues):
            if pvalue * lofreq_q.bonf_factor < lofreq_q.sig_thresh:
                # i == count
                return (last_step, i)
        last_step = i


def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog [Options]\n" \
            + "\n" + __doc__
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      dest="verbose",
                      action="store_true", 
                      help="be verbose")
    parser.add_option("", "--debug",
                      dest="debug",
                      action="store_true",
                      help="enable debugging")

    parser.add_option("-o", "--output",
                      dest="out",
                      help="FIXME output")

    plp_group = OptionGroup(parser, "Pileup Options", "")
    plp_group.add_option("-b", "--bam",
                      dest="bam", # type="string|int|float"
                      help="BAM file with (stringently) mapped reads")
    plp_group.add_option("-f", "--reffa",
                      dest = "ref_fasta", # type="string|int|float"
                      help = "Reference fasta file")
    choices = ['on', 'off', 'extended']
    plp_group.add_option("", "--baq",
                         dest = "baq",
                         choices = choices,
                         default=conf.DEFAULT_BAQ_SETTING,
                         help="BAQ setting for pileup."
                         " One of: %s (default: %s)"% (
                             ', '.join(choices), conf.DEFAULT_BAQ_SETTING))
    plp_group.add_option("-j", "--join-mapq-and-baseq",
                         dest = "join_mapq_and_baseq",
                         action="store_true",
                         help = "Join mapping and base quality")
    plp_group.add_option("-l", "--regions",
                         dest = "region_bed",
                         help = "Optional: bed file containing regions"
                         " to limit analysis to.")
    plp_group.add_option("-d", "--maxdepth",
                         dest = "max_depth",
                         type=int,
                         default = conf.DEFAULT_MAX_PLP_DEPTH,
                         help = "Maximum pileup depth (default: %d)"% (
                             conf.DEFAULT_MAX_PLP_DEPTH))
    parser.add_option_group(plp_group)

    
    filter_group = OptionGroup(parser, "Filtering Options", "")
    filter_group.add_option("-Q", "--ignore-bases-below-q",
                          dest="ign_bases_below_q", type="int",
                          default=conf.DEFAULT_IGN_BASES_BELOW_Q,
                          help="Remove any base below this base call"
                          " quality threshold from pileup."
                          " This will also be reflected in reported freqs"
                          " (default: %d)" % conf.DEFAULT_IGN_BASES_BELOW_Q)
    filter_group.add_option("", "--bonf",
                            dest="bonf", 
                            default='auto',
                            help="Bonferroni correction factor."
                            " Set to an integer, 'auto' (default) or"
                            " 'auto-ign-zero-cov'. 'auto' will use a"
                             " stringent and recommmended seqlen*3."
                            " 'auto-ign-zero-cov' is the same as 'auto'"
                            " but ignores zero coverage columns")
    filter_group.add_option("-s", "--sig-level",
                      dest="sig_thresh", type="float",
                      default=conf.DEFAULT_SIG_THRESH,
                      help="p-value significance value"
                      " (default: %g)" % conf.DEFAULT_SIG_THRESH)
    parser.add_option_group(filter_group)

    
    qual_group = OptionGroup(parser, "Advanced Options for quality-aware stage (LoFreq-Q)", "")
    qual_group.add_option("", "--noncons-default-qual",
                          dest="noncons_default_qual", type="int",
                          default=conf.NONCONS_DEFAULT_QUAL,
                          help="Base call quality used for non-consensus bases"
                          " (default: %d)" % conf.NONCONS_DEFAULT_QUAL)
    qual_group.add_option("", "--noncons-filter-qual",
                          dest="noncons_filter_qual", type="int",
                          default=conf.NONCONS_FILTER_QUAL,
                          help="Non-consensus bases below this threshold"
                          " will be removed before SNV calling"
                          " (but go into reported freqs;"
                          " default: %d)" % conf.NONCONS_FILTER_QUAL)
    parser.add_option_group(qual_group)


    return parser





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


    if opts.verbose:
        LOG.setLevel(logging.INFO)

    if opts.debug:
        LOG.setLevel(logging.DEBUG)


    # file check
    for (filename, descr, direction, mandatory) in [
            (opts.out, "Output file", 'out', True),
            (opts.bam, "BAM file", 'in', True),
            (opts.ref_fasta, "Reference fasta file", 'in', True),
            (opts.region_bed, "Region BED-file", 'in', False)
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
                
            
    if opts.out == '-':
        fhout = sys.stdout
    else:
        fhout = open(opts.out, 'w') 
        

    # quality filter settings (for quality based method)
    #
    noncons_filter_qual = opts.noncons_filter_qual
    noncons_default_qual = opts.noncons_default_qual
    ign_bases_below_q = opts.ign_bases_below_q

    # pvalue threshold and correction settings
    #
    sig_thresh = opts.sig_thresh
    if opts.bonf == "auto":
        bonf_factor = sam.auto_bonf_factor_from_depth(
            opts.bam, opts.region_bed, opts.ign_bases_below_q)
        LOG.info("Will use an automatically determined (len*3) Bonferroni"
                 " factor of %d." % (bonf_factor))
    elif opts.bonf == 'auto-ign-zero-cov':
        bonf_factor = sam.auto_bonf_factor(
            opts.bam, opts.region_bed)
        LOG.info("Will use an automatically determined"
                 " (non-zero-cov-len*3) Bonferroni"
                 " factor of %d." % (bonf_factor))
    else:
        try:
            bonf_factor = int(opts.bonf)
            assert opts.bonf > 0
        except:
            LOG.critical("Something's wrong with the Bonferroni factor"
                         " you provided...")
            sys.exit(1)

    
    lofreq_q = qual.QualBasedSNPCaller(
        noncons_default_qual = noncons_default_qual,
        noncons_filter_qual = noncons_filter_qual,
        ign_bases_below_q =  ign_bases_below_q,
        bonf_factor = bonf_factor, sig_thresh = sig_thresh)

    # pileup
    #
    pileup_obj = sam.Pileup(opts.bam, opts.ref_fasta)
    pileup_gen = pileup_obj.generate_pileup(
        opts.baq, opts.max_depth, opts.region_bed, opts.join_mapq_and_baseq)
   
        
    fhout.write("#pos lower upper\n")
    for pcol in pileup_gen:
        # note: not all columns will be present in pileup
        detlim = det_lim_for_pileup_line(pcol, lofreq_q)
        fhout.write("%d %d %d\n" % (
            pcol.coord+1, detlim[0], detlim[1]))
    if fhout != sys.stdout:
        fhout.close()  

        
if __name__ == "__main__":
    main()
    LOG.critical("Detection limit-range based on ref-quals only!")
    LOG.info("Successful program exit")

