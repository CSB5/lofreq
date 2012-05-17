#!/usr/bin/env python
"""FIXME: add doc string
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
import subprocess
# optparse deprecated from Python 2.7 on
from optparse import OptionParser

#--- third-party imports
#
# default is to use our own binomial extension instead of introducing
# a scipy dependency. but it's great for testing/validation
#
USE_SCIPY = True
if USE_SCIPY:
    from scipy import stats
    #from scipy.stats.distributions import binom


#--- project specific imports
#
from lofreq import pileup
#from lofreq.utils import count_bases
from lofreq import snp
if not USE_SCIPY:
    from lofreq_ext import binom_sf
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


DEFAULT_SAMTOOLS_ARGS = "-d 100000"


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

MYNAME = os.path.basename(sys.argv[0])



    

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
    parser.add_option("-d", "--snpdiff",
                      dest="snpdiff_file", # type="string|int|float"
                      help="List of SNPs, predicted from other sample and not predicted for this mapping (e.g. somatic calls)")
    parser.add_option("-b", "--bam",
                      dest="bam_file", # type="string|int|float"
                      help="BAM file to check")
    parser.add_option("-r", "--reffa",
                      dest="ref_fasta_file", # type="string|int|float"
                      help="Reference fasta file")
    parser.add_option("-c", "--chrom",
                      dest="chrom", # type="string|int|float"
                      help="Chromsome/sequence for pileup")
    parser.add_option("-a", "--samtools_args",
                      dest="samtools_args", # type="string|int|float"
                      default=DEFAULT_SAMTOOLS_ARGS,
                      help="Extra Samtools args (default: %s). Should be the same as used for the original SNP calling" % (DEFAULT_SAMTOOLS_ARGS))
    parser.add_option("-q", "--ignore-bases-below-q",
                      dest="ign_bases_below_q", type="int",
                      default=conf.DEFAULT_IGN_BASES_BELOW_Q,
                      help="Optional: remove any base below this base call quality threshold from pileup"
                      " (default: %d). Should be the same as used for the original SNP calling" % conf.DEFAULT_IGN_BASES_BELOW_Q)
    parser.add_option("", "--noncons-filter-qual",
                      dest="noncons_filter_qual", type="int",
                      default=0, #conf.NONCONS_FILTER_QUAL,
                      help="Optional: Non-consensus bases below this threshold will be filtered"
                          " (default: %d). Should be the same as used for the original SNP calling" % 0)#conf.NONCONS_FILTER_QUAL)
    parser.add_option("", "--uniform-freq",
                      dest="uniform_freq", type="int",
                      help="Optional: Assume this uniform frequency [%] instead of original SNP frequency")
    parser.add_option("", "--freq-factor",
                      dest="freq_fac", type="float",
                      help="Optional: Apply this factor (0.0<x<=1.0) to original SNP frequency")
    parser.add_option("-s", "--sig-level",
                      dest="sig_thresh", type="float",
                      default=conf.DEFAULT_SIG_THRESH,
                      help="Optional: p-value significance value"
                      " (default: %g)" % conf.DEFAULT_SIG_THRESH)

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

    if not opts.chrom:
        parser.error("Chromosome argument argument.")
        sys.exit(1)

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


    LOG.info("Removing any base with quality below %d and non-cons bases with quality below %d" % (
            opts.ign_bases_below_q, opts.noncons_filter_qual))
    
    fh_bed_region = tempfile.NamedTemporaryFile(delete=False)
    bed_region_file = fh_bed_region.name
    LOG.info("Creating region file samtools: %s" % bed_region_file)
    for s in snps:
        fh_bed_region.write("%s\t%d\n" % (opts.chrom, s.pos+1))
    fh_bed_region.close()

    cmd_args = ["samtools", "mpileup", "-l", bed_region_file, '-f', opts.ref_fasta_file]
    cmd_args.extend(opts.samtools_args.split())
    cmd_args.append(opts.bam_file)
    LOG.info("Calling: %s" % (' '.join(cmd_args)))
    process = subprocess.Popen(cmd_args, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    #(p_stdout, p_stderr) =  process.communicate()
    for line in process.stdout:
        pcol = pileup.PileupColumn(line)

        # FIXME make sure we remove used up candidates from the list
        # can use enum in list comphrehension to delete element later
        # Also need to report left over SNPs at the end as no coverage SNPs
        snp_candidates = [s for s in snps if s.pos == pcol.coord]
        assert len(snp_candidates) != 0, (
            "Oups..pileup for column %d has no matching SNP" % (pcol.coord+1))

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
                print "%s: not rejected (no coverage)" % (this_snp)
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

            if USE_SCIPY:
                pvalue = stats.binom.cdf(alt_count+1, cov, cmp_freq)
            else:
                # FIXME
                raise ValueError, ("Haven't coded up cdf support yet. Use scipy")
            if pvalue < opts.sig_thresh:
                print "%s: number of 'SNP' bases significantly low" % (this_snp.identifier())
            else:
                print "%s: not rejected" % (this_snp.identifier())
                
        # FIXME report left over SNPs at the end as no coverage SNPs
            
            
        

    #LOG.warn("Not deleting tmp file %s" % bed_region_file)
    os.unlink(bed_region_file)
    #if len(snps):
    #    LOG.warn("Unremoved SNPs left")
        
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
