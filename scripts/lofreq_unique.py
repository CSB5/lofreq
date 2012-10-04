#!/usr/bin/env python
"""Test wether a SNV call only made in one sample (e.g. tumour
tissue), but not in another (e.g. normal/blood tissue) is 'unique' or
significant'. Done by performing a binomial test on the prospective
snp-base counts in the normal sample given the SNV-call frequency.

NOTE: assuming that the initial SNV calls were made with default
pileup parameters
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
import tempfile
# optparse deprecated from Python 2.7 on
from optparse import OptionParser, OptionGroup

#--- third-party imports
#
# default is to use our own binomial extension instead of introducing
# a scipy dependency. but it's great for testing/validation
#
USE_SCIPY = False
if USE_SCIPY:
    from scipy.stats import binom
    binom_cdf = binom.cdf


#--- project specific imports
#
from lofreq import pileup
#from lofreq.utils import count_bases
from lofreq import snp
if not USE_SCIPY:
    from lofreq_ext import binom_cdf
from lofreq import conf

# invocation of ipython on exceptions
#import sys, pdb
#from IPython.core import ultratb
#sys.excepthook = ultratb.FormattedTB(mode='Verbose',
#                                     color_scheme='Linux', call_pdb=1)

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011, 2012 Genome Institute of Singapore"
__license__ = "GPL2"


#DEFAULT_SAMTOOLS_ARGS = "-d 100000"


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

if USE_SCIPY:
    LOG.warn("Using scipy functions instead of internal ones. Should only be used for debugging.")

MYNAME = os.path.basename(sys.argv[0])

    

def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog [Options]\n" \
            + "\n" + __doc__
    parser = OptionParser(usage=usage)

    parser.add_option("-d", "--snpdiff",
                      dest="snpdiff_file", # type="string|int|float"
                      help="List of SNPs, predicted from other sample and not predicted for this mapping/BAM-file (e.g. somatic calls)")
    parser.add_option("-b", "--bam",
                      dest="bam_file", # type="string|int|float"
                      help="BAM file to check (i.e. given SNVs where unique to another sample")
    parser.add_option("-r", "--reffa",
                      dest="ref_fasta_file", # type="string|int|float"
                      help="Reference fasta file")
    parser.add_option("-o", "--out",
                      dest="out", # type="string|int|float"
                      default="-",
                      help="Output file ('-' for stdout, which is default)")


    opt_group = OptionGroup(parser, "Optional arguments", "")
    opt_group.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="be verbose")
    opt_group.add_option("", "--debug",
                      action="store_true", dest="debug",
                      help="enable debugging")

    #opt_group.add_option("-a", "--samtools_args",
    #                  dest="samtools_args", # type="string|int|float"
    #                  default=DEFAULT_SAMTOOLS_ARGS,
    #                  help="Optional: Extra Samtools args (default: %s). Should be the same as used for the original SNP calling" % (DEFAULT_SAMTOOLS_ARGS))
    opt_group.add_option("-q", "--ignore-bases-below-q",
                      dest="ign_bases_below_q", type="int",
                      default=conf.DEFAULT_IGN_BASES_BELOW_Q,
                      help="Optional: remove any base below this base call quality threshold from pileup"
                      " (default: %d). Should be the same as used for the original SNP calling" % conf.DEFAULT_IGN_BASES_BELOW_Q)
    opt_group.add_option("", "--noncons-filter-qual",
                      dest="noncons_filter_qual", type="int",
                      default=0, #conf.NONCONS_FILTER_QUAL,
                      help="Optional: Non-consensus bases below this threshold will be filtered"
                          " (default: %d). Should be the same as used for the original SNP calling" % 0)#conf.NONCONS_FILTER_QUAL)
    opt_group.add_option("", "--uniform-freq",
                      dest="uniform_freq", type="int",
                      help="Optional: Assume this uniform frequency [%] instead of original SNP frequency")
    opt_group.add_option("", "--freq-factor",
                      dest="freq_fac", type="float",
                      help="Optional: Apply this factor (0.0<x<=1.0) to original SNP frequency")
    opt_group.add_option("-s", "--sig-level",
                      dest="sig_thresh", type="float",
                      default=conf.DEFAULT_SIG_THRESH,
                      help="Optional: p-value significance value"
                      " (default: %g)" % conf.DEFAULT_SIG_THRESH)
    opt_group.add_option("-c", "--chrom",
                      dest="chrom", # type="string|int|float"
                      help="Deprecated: chromsome/sequence for pileup. Only if SNVs don't have chrom markup")

    parser.add_option_group(opt_group)
    
    return parser



def main():
    """
    The main function

    doctest for our own cdf function
    >>> from scipy.stats import binom
    >>> x = binom.cdf(10, 1000, 0.01)
    >>> from lofreq_ext import binom_cdf
    >>> y = binom_cdf(10, 1000, 0.01)
    >>> print '%.6f' % x == '%.6f' % y
    True
    
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

    for (in_file, descr) in [(opts.snpdiff_file, "SNP diff"), 
                             (opts.bam_file, "BAM"),
                             (opts.ref_fasta_file, "Reference fasta")]:
        if not in_file:
            parser.error("%s input file argument missing." % descr)
            sys.exit(1)
        if not os.path.exists(in_file):
            sys.stderr.write(
                "file '%s' does not exist.\n" % in_file)
            sys.exit(1)
            
    if opts.out == '-':
        fh_out = sys.stdout
    else:
        if os.path.exists(opts.out):
            sys.stderr.write(
                "Refusing to overwrite existing file '%s'.\n" % opts.out)
            sys.exit(1)
        fh_out = open(opts.out, 'w')
            
    if opts.uniform_freq or opts.uniform_freq == 0:
        if opts.uniform_freq < 1 or opts.uniform_freq > 100:
            LOG.fatal("Frequency out of valid range (1-100%)\n")
            sys.exit(1)
    if opts.freq_fac or opts.freq_fac == 0:
        if opts.freq_fac < 0.0 or opts.freq_fac > 1.0:
            LOG.fatal("Frequency factor out of valid range (0.0<x<=1.0)")
            sys.exit(1)
    if opts.freq_fac and opts.uniform_freq:
        LOG.fatal("Can't use both, uniform frequency and frequency factor")
        sys.exit(1)                    
    snps = snp.parse_snp_file(opts.snpdiff_file)
    LOG.info("Parsed %d SNPs from %s" % (len(snps), opts.snpdiff_file))

    # fix chromosome info if necessary
    if set([s.chrom for s in snps if s.chrom]) == 0:
        if opts.chrom:
            for s in snps:
               s.chrom = opts.chrom
        else:
            LOG.fatal("SNV file did not contain chromosome info and none given on command line")
            sys.exit(1)

    LOG.info("Removing any base with quality below %d and non-cons bases with quality below %d" % (
            opts.ign_bases_below_q, opts.noncons_filter_qual))

    # drop a bed-file only containing regions of interest
    #
    fh_bed_region = tempfile.NamedTemporaryFile(delete=False)
    bed_region_file = fh_bed_region.name
    LOG.info("Creating region file samtools: %s" % bed_region_file)
    for s in snps:
        fh_bed_region.write("%s\t%d\t%d\n" % (s.chrom, s.pos, s.pos+1))
    fh_bed_region.close()


    pileup_obj = pileup.Pileup(opts.bam_file, opts.ref_fasta_file)
    pileup_gen = pileup_obj.generate_pileup(
        conf.DEFAULT_BAQ_SETTING, conf.DEFAULT_MAX_PLP_DEPTH, bed_region_file)
    for pcol in pileup_gen:
        
        snp_candidates = [s for s in snps if
                          s.pos == pcol.coord and s.chrom == pcol.chrom]
        assert len(snp_candidates) != 0, (
            "Oups..pileup for %s:%d has no matching SNP" % (pcol.chrom, pcol.coord+1))

        for this_snp in snp_candidates:
            ref_count = sum(pcol.get_counts_for_base(
                this_snp.wildtype, opts.ign_bases_below_q))
            nonref_counts = dict()
            for base in pileup.VALID_BASES:
                if base == this_snp.wildtype or base == 'N':
                    continue
                nonref_counts[base] = sum(pcol.get_counts_for_base(
                        base, max(opts.ign_bases_below_q, opts.noncons_filter_qual)))
            cov = sum([ref_count, sum(nonref_counts.values())])
            alt_count = nonref_counts[this_snp.variant]

            if cov == 0:
                fh_out.write("%s: not necessarily unique (not rejected)\n" % (this_snp))
                continue

            if opts.uniform_freq:
                cmp_freq = opts.uniform_freq/100.0
            else:
                cmp_freq = this_snp.freq
                if opts.freq_fac:
                    cmp_freq = cmp_freq*opts.freq_fac
            LOG.info("Testing SNP at pos. %d %s>%s %f. Counts and freq in BAM ref=%d snp=%d snp-freq=%f" % (
                    pcol.coord+1, this_snp.wildtype, this_snp.variant, cmp_freq,
                    ref_count, alt_count, alt_count/float(cov)))

            pvalue = binom_cdf(alt_count+1, cov, cmp_freq)

            if pvalue < opts.sig_thresh:
                fh_out.write("%s: unique (number of 'SNP' bases significantly low)\n" % (this_snp.identifier()))
            else:
                fh_out.write("%s: not necessarily unique (not rejected)\n" % (this_snp.identifier()))

            snps.remove(this_snp)

            
    if len(snps):
        LOG.error("Oups...some SNVs left after processing. This shouldn't happen")
        

    #LOG.warn("Not deleting tmp file %s" % bed_region_file)
    os.unlink(bed_region_file)
    #if len(snps):
    #    LOG.warn("Unremoved SNPs left")
    if fh_out != sys.stdout:
        fh_out.close()
        
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
