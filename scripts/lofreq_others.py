#!/usr/bin/env python
"""Implementation of a few other low frequency SNV calling "methods"
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
import os
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser

#--- third-party imports
#
#/

#--- project specific imports
#
from lofreq import pileup
from lofreq import snp
from lofreq import simple_vcf
from lofreq import utils


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011, 2012 Genome Institute of Singapore"
__license__ = "GPL2"


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


BASES = ['A', 'C', 'G', 'T']


def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog [Options]\n" \
            + __doc__
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true",
                      dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true", 
                      dest="debug",
                      help="enable debugging")
    choices=['goto_2011', 'wright_2011']
    parser.add_option("-m", "--method",
                      dest="method", 
                      choices=choices,
                      help="Which method to use (one of %s)" % (', '.join(choices)))
    parser.add_option("-e", "--exclude",
                      dest="fexclude",
                      help="Exclude positions listed in this file"
                      " format is: start end [comment ...]"
                      " , with zero-based, half-open coordinates")
    parser.add_option("-i", "--in",
                      dest="pileup", 
                      default="-",
                      help="Input [m]pileup or '-' for stdin (default)")
    parser.add_option("-o", "--out",
                      dest="fsnp",
                      default="-",
                      help="Variant output file or '-' for stdout (default)")
    choices = ['snp', 'vcf']
    parser.add_option("", "--format",
                      dest="outfmt", 
                      choices=choices, 
                      default='snp',
                      help="Output format. One of: %s. SNP is chromsome agnostic!" % (
                          ', '.join(choices)))

    return parser



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



def cordey_2010(pcol):
    """"
    Cordey et al., 2010
    -------------------

    Quote:

    The MAQ consensus sequences, including SNP detection, were generated
    for all regions with a minimum coverage of three bases. We conducted a
    statistical analysis of the counts of the number of mapped A/C/G/T to
    extract potential SNP positions for each position on the reference
    sequence. Counts are used to determine the 95% confidence interval of
    the probability of observing A/C/G/T at the position while assuming
    the probabilities to follow a beta distribution. When the confidence
    interval of two of the bases probability is above a 5% threshold, the
    position is considered as a statistically significant SNP.
    """
    
    raise NotImplementedError, ("FIXME: not implemented")




def eckerle_2010(pcol):
    """
    Eckerle et al., 2010
    --------------------
    
    Quote
    
    Position profiles and SNPs from the read mapping were then generated
    for each sample. A position profile is defined for each sample as the
    observed nucleotide distribution at each position along the reference
    genome, based upon the reads that mapped for that sample. Position
    profiles for each sample were also generated separately for forward
    reads and reverse reads, as well as bi-directionally. Before
    generating the position profiles, several quality filters were applied
    to the data. First, the last five bases of each read were trimmed to
    remove low quality bases that tend to be more frequent toward the end
    of each read. Second, any bases from reads that overlapped with a PCR
    primer region were trimmed. Third, any base with a quality value ,29
    was removed from the position profiles. Finally, for the bidirectional
    position profiles, bases that had ,20x coverage in either read
    direction were excluded. For the determination of SNPs, two additional
    filters were imposed: 1) at least 20 combined reads from both
    directions must support the SNP, and 2) the SNP must occur in at least
    5% of the reads in each direction. After initial analysis, it was
    determined that reads with .3 mismatches should be excluded to further
    reduce the noise level attributable to mapping errors. The list of
    SNPs identified was unaffected by this additional filter, further
    supporting the idea that the extra mismatches could be attributed to
    random mapping errors.
    """
    
    raise NotImplementedError, ("FIXME: not implemented")



def wright_2011(pcol):
    """Wright et al., 2011
    -------------------

    Quote:
    
    we used the quality score of each nucleotide read to compute the
    average probability of a sequencing error, p_i, at each site i.
    Typical values of p_i are around 0.1%. Assuming sequencing errors
    to be independent, we computed the expected number of such errors
    as the mean of the binomial distribution B(x; p_i,n_i), where n_i
    is the coverage of site i. If the observed number of mismatches
    exceeded this expected number of errors in both runs, then we
    excluded the possibility of a sequencing error.

    On the other hand, we hypothesized that the probability that PCR
    errors in both runs independently generated identical base changes
    at the same site was very low. Based on values quoted for the
    enzymes used, we estimated that the error rate for the combined
    RT-PCR amplification process was 7.7 x 10^6 per base pair copied
    (2, 31, 38). We therefore defined polymorphic sites that could not
    be attributed to sequencing errors and at which both the
    most-common and second-most-common nucleotides were the same
    between the two runs as qualitatively validated sites.

    For each site in the set of qualitatively validated polymorphisms,
    we computed the 95% confidence interval (95% CI) for the
    polymorphism frequencies, using the binomial distribution
    described above. If the 95% confidence intervals from the runs
    overlapped, we defined the polymorphism frequency estimates from
    the two runs to be in quantitative agreement.


    NOTE: Wright et al. differentiate between three polymorphisms,
    which are basically represent three consecutive levels:

    1. simple polymorphisms: possibility of a sequencing errorr can be
    ruled out (see above for definition)

    2. qualitatively validated sites: 1 and most-common and
    second-most-common nucleotides were the same between the two runs

    3. quantitative agreement: 2 and if the 95% confidence intervals
    (for the polymorphism frequencies, using the binomial distribution
    described above) from the runs overlap
    
    We can only implement #1 here! They also used in-house alignment
    and read-filtering routines (and had replicates).
    """
    
    snps = []
    keep_strand_info = False
    base_qual_hist = pcol.get_base_and_qual_hist(keep_strand_info)
    # example: {'A': {20: 3, 22: 1}, 'C': {}, 'T': {}, 'G': {}, 'N': {}}
    # del base_qual_hist['N']

    q_sum = 0
    coverage = 0
    for base in base_qual_hist.keys():
        q_sum += sum([q*c for (q, c) in base_qual_hist[base].iteritems()])
        coverage += sum(base_qual_hist[base].values())   
    avg_qual = q_sum/float(coverage)
    avg_prob = utils.phredqual_to_prob(avg_qual)
    # the "mean of the binomial distribution B(x; p_i,n_i)" is simply
    # p_i*n_i !?? Why didn't they do a binomial test?
    bmean = avg_prob * coverage

    for snp_base in [b for b in BASES if b != pcol.ref_base]:
        snp_base_count =  sum(base_qual_hist[snp_base].values())
        if snp_base_count > bmean:
            #print "snp at pos %d %c>%c (counts=%d > bmean=%f)" % (
            #    pcol.coord+1, pcol.ref_base, snp_base, snp_base_count, bmean)
            info_dict = dict()
            info_dict['coverage'] = coverage
            new_snp = snp.ExtSNP(pcol.coord, pcol.ref_base, snp_base,
                                 snp_base_count/float(coverage),
                                 info_dict)
            snps.append(new_snp)
        #else:
        #    print "no snp at pos %d %c>%c (counts=%d <= bmean=%f)" % (
        #1        pcol.coord+1, pcol.ref_base, snp_base, snp_base_count, bmean)
                        
    return snps


    
def goto_2011(pcol):
    """Goto et al., 2011
    -----------------

    Quote
    
    To call a site heteroplasmic, we require the frequency of reads
    supporting a particular allele to be >=0.02 (to be conservative, we
    doubled the threshold from 0.01 to 0.02) on each strand and the
    quality of the base aligning to such a position to be >=30 on the
    phred scale (corresponding to an error probability of 0.001).

    
    NOTE: in their paper they had PCR replicates which might reduce
    the number of false positives (and true positives) 
    """

    # NOTE: The recipe in the paper is incomplete. For example, is the
    # filtering done before counting? Also, what's the final frequency
    # then?

    snps = []
    
    freq_tresh = 0.02

    # count all
    min_qual = 0
    all_fw_counts = all_rv_counts = 0
    keep_strand_info = True
    for (base, strand_counts) in pcol.get_all_base_counts(
            min_qual, keep_strand_info).iteritems():
        all_fw_counts += strand_counts[0]
        all_rv_counts += strand_counts[1]

    coverage = all_fw_counts + all_rv_counts
    if coverage == 0:
        return []
        
    # count possible SNP bases incl qual treshhold
    min_qual = 30
    for snp_base in [b for b in BASES if b != pcol.ref_base]:
        (snp_fw_counts, snp_rv_counts) = pcol.get_counts_for_base(
            snp_base, min_qual, keep_strand_info)

        snp_freq_fw = 0.0
        if all_fw_counts:
            snp_freq_fw = snp_fw_counts/float(all_fw_counts)
            
        snp_freq_rv = 0.0
        if all_rv_counts:
            snp_freq_rv = snp_rv_counts/float(all_rv_counts)
            
        if snp_freq_fw >= freq_tresh and snp_freq_rv >= freq_tresh:
            #print "SSSSSSSSNNNNPP at pos %d %c>%c" % (
            #pcol.coord+1, pcol.ref_base, snp_base)
            info_dict = dict()
            info_dict['coverage'] = coverage
            new_snp = snp.ExtSNP(pcol.coord, pcol.ref_base, snp_base,
                                 (snp_freq_fw+snp_freq_rv)/2.0,
                info_dict)
            snps.append(new_snp)
            
    return snps


    
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
        pileup.LOG.setLevel(logging.INFO)

    if opts.debug:
        LOG.setLevel(logging.DEBUG)
        pileup.LOG.setLevel(logging.DEBUG)

    if opts.fexclude and not os.path.exists(opts.fexclude):
        LOG.fatal("file '%s' does not exist.\n" % (opts.fexclude))
        sys.exit(1)
        
    if opts.pileup == "-":
        LOG.info("Reading from stdin")
        pileup_fh = sys.stdin
    else:
        LOG.info("Reading from %s" % opts.pileup)
        pileup_fh = open(opts.pileup, 'r')

    if not opts.fsnp:
        parser.error("Variant output file argument missing.")
        sys.exit(1)
    if opts.fsnp != '-' and os.path.exists(opts.fsnp):
        LOG.fatal(
            "Cowardly refusing to overwrite already existing file '%s'.\n" % (
                opts.fsnp))
        sys.exit(1)
        
    if opts.outfmt == 'vcf':
        LOG.fatal("vcf writing not yet supported")
        sys.exit(1)

    if not opts.method:
        parser.error("Missing method argument")
        
    # exclude positions
    #
    excl_pos = []
    if opts.fexclude:
        excl_pos = read_exclude_pos_file(opts.fexclude)
        LOG.info("Ignoring %d positions found in %s" % (
            len(excl_pos), opts.fexclude))
        LOG.debug("DEBUG: excl_pos = %s" % excl_pos)

    if opts.method == 'goto_2011':
        snp_call_func = goto_2011
    elif opts.method == 'wright_2011':
        snp_call_func = wright_2011
    else:
        raise NotImplementedError, (
            "Unknown method '%s'" % opts.method)
    LOG.info("Calling SNPs according to %s" % opts.method)

    
    snps = []
    num_lines = 0
    for line in pileup_fh:
        # note: not all columns will be present in pileup
        num_lines += 1
        #pileup_line_buffer.append(line)
        pcol = pileup.PileupColumn(line)

        # skip if in excl_pos or refbase ambigious/unknown
        if pcol.coord in excl_pos:
            LOG.debug("Skipping col %d because of exclusion" % (pcol.coord+1))
            continue
        if pcol.ref_base not in 'ACGT':
            LOG.info("Skipping col %d because of amibigous reference base %s" % (
                pcol.coord+1, pcol.ref_base))
            continue
        
        snps.extend(snp_call_func(pcol))
        
    if num_lines == 0:
        LOG.fatal("Pileup was empty. Will exit now.")
        sys.exit(1)

    if pileup_fh != sys.stdin:
        pileup_fh.close()



    # save
    #
    #    
    if opts.fsnp == '-':
        fhout = sys.stdout
    else:
        LOG.info("Saving %d SNPs to %s" % (len(snps), opts.fsnp))
        fhout = open(opts.fsnp, 'w')
    if opts.outfmt == 'snp':
        snp.write_header(fhout)
    else:
        simple_vcf.write_header(fhout)                
    for s in snps:
        if opts.outfmt == 'snp':
            snp.write_record(s, fhout)
        else:
            LOG.fatal("vcf writing not yet supported")
            #simple_vcf.write_record(vcf_record, fhout)
    if fhout != sys.stdout:
        fhout.close()  
                


if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
