#!/usr/bin/env python
"""FIXME
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
from math import log
import subprocess
 
#--- third-party imports
#
# default is to use our own binomial extension instead of introducing
# a scipy dependency. but it's great for testing/validation
#
# FIXME remove dependency
from scipy import stats

#--- project specific imports
#
##from lofreq import sam
##from lofreq.utils import count_bases
##from lofreq import snp
##from lofreq import conf

# invocation of ipython on exceptions
if False:
    import sys, pdb
    from IPython.core import ultratb
    sys.excepthook = ultratb.FormattedTB(mode='Verbose',
                                         color_scheme='Linux', call_pdb=1)

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011, 2012 Genome Institute of Singapore"
__license__ = "GPL2"



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

    parser.add_option("-n", "--normal",
                      dest="normal_bam",
                      help="BAM file for normal samples")
    parser.add_option("-t", "--tumor",
                      dest="tumor_bam",
                      help="BAM file for tumor samples")
    parser.add_option("-r", "--ref-fa",
                      dest="ref_fa",
                      help="Reference fasta file")
    parser.add_option("-l", "--bed",
                      dest="bed",
                      help="Bed file listing regions to run on")
    # FIXME add script to derive automatically from WGS
    parser.add_option("-o", "--outdir",
                      dest="outdir",
                      help="Output directory (must exist)")


    opt_group = OptionGroup(parser, "Optional arguments", "")
    opt_group.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="be verbose")
    opt_group.add_option("", "--debug",
                      action="store_true", dest="debug",
                      help="enable debugging")
    opt_group.add_option("", "--combine-pvalues",
                      action="store_true",
                      dest="combine_pvalues",
                      help="Combine SNV call and binomial test p-values for final calls")
    DEFAULT=20
    opt_group.add_option("", "--min-snp-phred",
                      dest="min_snp_phred",
                      default=DEFAULT,
                      help="SNV Phred cutoff for filtering (default=%s)" % DEFAULT)
    DEFAULT=10
    opt_group.add_option("", "--min-cov",
                      dest="mincov",
                      default=DEFAULT,
                      help="SNV Phred cutoff for filtering (default=%s)" % DEFAULT)
    opt_group.add_option("", "--sbfilter",
                      dest="sbhbfilter",
                      help="Holm-Bonferroni strandbias filtering")

    parser.add_option_group(opt_group)
    
    return parser



def snv_call_wrapper(bam, ref_fa, bed, snv_raw, bonf='auto-ign-zero-cov'):
    """Wrapper for lofreq_snpvaller.py

    sUse stdout if snv_raw is None
    """

    
    if snv_raw and os.path.exists(snv_raw):
        LOG.warn("Reusing %s." % snv_raw)
        return snv_raw

    if snv_raw:
        LOG.fatal("lofreq2 call does not suuport output arguments at the moment")
        return None
    cmd_list = ['lofreq2', 'call',
                '--bonf', "%s" % bonf,
                '-f', ref_fa,
                '-l', bed,
                #'-o', snv_raw if snv_raw else "-",
                '--verbose', 
                bam]
    return cmd_list



def snv_filter_wrapper(snv_raw, snv_filt, 
                       minphred=None, mincov=None, sbfilter=None):
    """Wrapper for lofreq2 filter

    Use stdin if snv_raw is None.
    Use stdout if snv_filt is None.
    """

    #if snv_filt and os.path.exists(snv_filt):
    #    LOG.warn("Reusing %s" % snv_filt)
    #    return snv_filt

    if not any([minphred, mincov, sbfilter]):
        LOG.fatal("No filtering requested. Will just copy files")
        return None

    cmd_list = ['lofreq2', 'filter']
    if sbfilter:
        cmd_list.append('--strandbias-holmbonf')
        
    if mincov:
        cmd_list.extend(['--min-cov', "%d" % mincov])
        
    if minphred:
        cmd_list.extend(['--snp-phred', "%d" % minphred])
        
    cmd_list.extend(['-i', snv_raw if snv_raw else "-", 
                     '-o', snv_filt if snv_filt else "-", 
                     '-v'])
    return cmd_list



def main():
    """The main function
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

    for (in_file, descr) in [(opts.normal_bam, "Normal BAM"), 
                             (opts.tumor_bam, "Tumor BAM"),
                             (opts.ref_fa, "Reference fasta"),
                             (opts.bed, "Bed (regions)")]:
        if not in_file:
            parser.error("%s input file argument missing." % descr)
            sys.exit(1)
        if not os.path.exists(in_file):
            LOG.fatal("file '%s' does not exist.\n" % in_file)
            sys.exit(1)
            
    if not opts.outdir or not os.path.isdir(opts.outdir):
        LOG.fatal("output directory does not exist or argument missing")
        sys.exit(1)

    if opts.combine_pvalues:
        # should go to cmdline option
        LOG.error("FIXME Experimental code pvalue merging code."
                  " Don't use upstream filtering on this."
                  " Bonf correction must come downstream!")

    
    snv_call_cmd = snv_call_wrapper(opts.tumor_bam, opts.ref_fa, opts.bed, None)
    LOG.debug("snv_call_cmd = %s" % snv_call_cmd)
    basename = os.path.splitext(os.path.basename(opts.tumor_bam))[0] 
    basename =  os.path.join(opts.outdir, basename)
    snv_filt_out = basename + ".snp"
    snv_filt_cmd = snv_filter_wrapper(None, snv_filt_out,
                                      opts.min_snp_phred, opts.mincov, opts.sbhbfilter)
    LOG.debug("snv_filt_cmd = %s" % snv_filt_cmd)
    
    p1 = subprocess.Popen(snv_call_cmd, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(snv_filt_cmd, stdin=p1.stdout)
    p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.

    (stdoutdata, stderrdata) =  p2.communicate()
    import pdb; pdb.set_trace()

    LOG.critical("Now run: bgzip $vcf and tabix ${vcf}.gz")
    LOG.critical("Now run: vcf-isec -c tumour.vcf.gz normal.vcf.gz and save to diff.vcf")

    """
    pileup_obj = sam.Pileup(opts.bam_file, opts.ref_fasta_file)
    pileup_gen = pileup_obj.generate_pileup(
        conf.DEFAULT_BAQ_SETTING, conf.DEFAULT_MAX_PLP_DEPTH, bed_region_file)
    for pcol in pileup_gen:
        # in theory these pos should come in same order as input snps
        # but cols can be empty as well if they have no coverage
        if pos_map_aln_to_other:
            offset_pos = pos_map_other_to_aln[pcol.coord]
        else:
            offset_pos = pcol.coord

        if offset_pos == -1:
            LOG.fatal("INTERNAL ERROR: No alignment match for this SNP")
            sys.exit(1)

        snp_candidates = [s for s in snps if
                          s.pos == offset_pos and 
                          s.chrom == pcol.chrom]
        assert len(snp_candidates) != 0, (
            "Oops..pileup for %s:%d has no matching SNP" % (
                pcol.chrom, pcol.coord+1))

        for this_snp in snp_candidates:
            ref_count = sum(pcol.get_counts_for_base(
                this_snp.wildtype, opts.ign_bases_below_q))
            nonref_counts = dict()
            for base in sam.VALID_BASES:
                if base == this_snp.wildtype or base == 'N':
                    continue
                nonref_counts[base] = sum(pcol.get_counts_for_base(
                        base, 
                        max(opts.ign_bases_below_q, opts.noncons_filter_qual)))
            cov = sum([ref_count, sum(nonref_counts.values())])
            alt_count = nonref_counts[this_snp.variant]

            if cov == 0:
                LOG.info("%s: not necessarily unique"
                         " (zero coverage in other sample)\n" % (this_snp))
                continue

            if opts.uniform_freq:
                cmp_freq = opts.uniform_freq/100.0
            else:
                cmp_freq = this_snp.freq
                if opts.freq_fac:
                    cmp_freq = cmp_freq*opts.freq_fac
                    
            LOG.info("Testing potentially unique SNP %s %d %s>%s %f"
                     " (using freq %f for test)."
                     " Counts and freq in 'other' BAM:"
                     " #alt=%d cov=%d snp-freq=%f" % (
                    pcol.chrom,  pcol.coord+1, this_snp.wildtype, 
                    this_snp.variant, this_snp.freq, cmp_freq,
                    alt_count, cov, alt_count/float(cov)))

            pvalue = binom_cdf(alt_count, cov, cmp_freq)

            if opts.combine_pvalues and s.info['type'] == 'low-freq-var':
                # fisher's method
                \"""
                p://en.wikipedia.org/wiki/Fisher's_method
                --- 
                
                breseq-0.18b/src/c/breseq/polymorphism_statistics.r:
                
                combined_log = - 2* ( log(ks_test_p_value) + log(fisher_test_p_value) )
                Y$bias_p_value[i] = pchisq(combined_log, 2*2, lower.tail=F)
                Y$bias_e_value[i] = Y$bias_p_value[i] * total_length
                
                ---
                
                R:
                
                > ks_test_p_value = 0.05
                > fisher_test_p_value = 0.05
                > combined_log = - 2* ( log(ks_test_p_value) + log(fisher_test_p_value) )
                > combined_log
                [1] 11.98293
                > pchisq(combined_log, 2*2, lower.tail=F)
                [1] 0.01747866
                
                ---
                
                from scipy.stats import chi2
                chi2.sf(11.98293, 4
                0.017478654584055061
                
                # 1 - chi2.cdf(11.98293, 1)
                # Out[25]: 0.00053690098750680537
                \"""
                comb_log = - 2* (log(pvalue) + log(float(this_snp.info['pvalue'])))
                pvalue = chi2.sf(comb_log, 2*2)
                
            if pvalue < opts.sig_thresh:
                LOG.info("%s: unique (number of candidate SNP bases"
                         " significantly low)\n" % (
                             this_snp.identifier()))
                snp.write_record(this_snp, fh_out)
            else:
                LOG.info("%s: not necessarily unique"
                         " (not rejected)\n" % (this_snp.identifier()))

            snps.remove(this_snp)

            
    if len(snps):
        LOG.warn("%d SNVs left after processing. Pileup potentially"
                 " empty for these regions. Not deleting bed-file (%s)"
                 " to make debugging easier. Left SNVs are: %s" % (
                     len(snps), bed_region_file, ','.join([s.identifier() 
                                                           for s in snps])))
    else:
        #LOG.warn("Not deleting tmp file %s" % bed_region_file)
        os.unlink(bed_region_file)
        
    if fh_out != sys.stdout:
        fh_out.close()
        """
        
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
