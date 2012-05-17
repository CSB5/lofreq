#!/usr/bin/env python
"""Exact fisher test on strand bias on already called SNP positions.
P-Values are not Bonferroni corrected
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
import os, sys
import tempfile
import subprocess
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser

USE_SCIPY = False

#--- third-party imports
#
if USE_SCIPY:
    from scipy.stats import fisher_exact

    
#--- project specific imports
#
from lofreq import pileup
from lofreq import snp
from lofreq import conf
if not USE_SCIPY:
    from lofreq_ext import kt_fisher_exact


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011, 2012 Genome Institute of Singapore"
__license__ = "GPL2"


BASES = ['A', 'C', 'G', 'T', 'N']
DEFAULT_SAMTOOLS_ARGS = "-d 100000"


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')




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
    parser.add_option("-s", "--snp_file",
                      dest="snp_file", # type="string|int|float"
                      help="SNP file)")
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
                      default=conf.NONCONS_FILTER_QUAL,
                      help="Optional: Non-consensus bases below this threshold will be filtered"
                          " (default: %d). Should be the same as used for the original SNP calling" % conf.NONCONS_FILTER_QUAL)
    
    LOG.warning("What about sig-level/bonf-factor?")
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

    for (in_file, descr) in [(opts.snp_file, "SNP"), 
                             (opts.bam_file, "BAM"),
                             (opts.ref_fasta_file, "Reference fasta")]:
        if not in_file:
            parser.error("%s input file argument missing." % descr)
            sys.exit(1)
        if not os.path.exists(in_file):
            sys.stderr.write(
                "file '%s' does not exist.\n" % in_file)
            sys.exit(1)


    snps = snp.parse_snp_file(opts.snp_file)
    LOG.info("Parsed %d SNPs from %s" % (len(snps), opts.snp_file))

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
        if len(line.strip())==0:
            continue
        
        pcol = pileup.PileupColumn(line)
        #pcol.parse_line(line, keep_strand_info=True, delete_raw_values=False)
        #pcol.rem_ambiguities()
        #pcol.rem_bases_below_qual(opts.ign_bases_below_q)

        snp_candidates = [s for s in snps if s.pos == pcol.coord]
        assert len(snp_candidates) != 0, (
            "Oups..pileup for column %d has no matching SNP" % (pcol.coord+1))

        for this_snp in snp_candidates:

            # FIXME there is a function for this in main snpcaller

            ref_counts = pcol.get_counts_for_base(
                this_snp.wildtype, opts.ign_bases_below_q)
            snp_counts = pcol.get_counts_for_base(
                this_snp.variant, max(opts.noncons_filter_qual, opts.ign_bases_below_q))
    
            try:
                if USE_SCIPY:
                    # alternative : two-sided, less, greater
                    # Which alternative hypothesis to the null hypothesis the test uses. Default is two-sided.
                    (oddsratio, pvalue) = fisher_exact([[ref_counts[0], ref_counts[1]],
                                                        [snp_counts[0], snp_counts[1]]])
                else:
                    # expects two tuples as input
                    # returns left, right, twotail pvalues
                    (left_pv, right_pv, pvalue) = kt_fisher_exact((ref_counts[0], ref_counts[1]),
                                                                  (snp_counts[0], snp_counts[1]))
            except ValueError:
                pvalue = -1

            print "%s | counts: ref-fw/ref-rv %d/%d var-fw/var-rv %d/%d | strand bias pvalue %f" % (
                this_snp, ref_counts[0], ref_counts[1], snp_counts[0], snp_counts[0], pvalue)
        
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
