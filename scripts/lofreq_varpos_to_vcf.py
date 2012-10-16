#!/usr/bin/env python
"""Create a vcf file from variant positions in a pileup. Mainly used
as input for GATK quality recalibration.

Make sure to use a reference when creating the pileup. Otherwise all
reference bases in the pileup will appear to be N, which makes them
unusable.
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
import datetime
from collections import namedtuple
# optparse deprecated from Python 2.7 on
from optparse import OptionParser

#--- third-party imports
#
# /

#--- project specific imports
#
from lofreq import sam

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011, 2012 Genome Institute of Singapore"
__license__ = "GPL2"


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

MYNAME = os.path.basename(sys.argv[0])


VarPos = namedtuple('VarPos',
                    ['chrom', 'pos', 'refbase', 'altbase'])



class VcfWriter(object):
    """
    Based on function in coverage_profile_to_vcf.py

    Quick & dirty...
    """
    
    default_id = "."
    default_qual = "."
    default_filter = "PASS"
    default_info = "."

    
    def __init__(self, fhandle = None):
        """FIXME
        """

        self.fhandle = fhandle


    def write_header(self, source=None, reference=None):
        """Print VCF header
        """

        self.fhandle.write("##fileformat=VCFv4.0\n")
        self.fhandle.write("##fileDate=%s\n" % (
            datetime.datetime.now().strftime("%Y%m%d")))
        if source:
            self.fhandle.write("##source=%s\n" % (source))
        if reference:
            self.fhandle.write("##reference=%s\n" % reference)
    
        self.fhandle.write("#%s\n" % '\t'.join(
            ["CHROM", "POS", "ID", "REF",
             "ALT", "QUAL", "FILTER", "INFO"]))



    def write_variants(self, var_pos_list):
        """Print given variants
        """
        
        #for var_pos in sorted(var_pos_list, key = lambda v: v.pos):
        
        for var_pos in var_pos_list:
            write_vals = [var_pos.chrom,
                          var_pos.pos+1,
                          self.default_id,
                          var_pos.refbase,
                          var_pos.altbase,
                          self.default_qual,
                          self.default_filter,
                          self.default_info]
            self.fhandle.write("%s\n" % ('\t'.join(
                [str(x) for x in write_vals])))
        

    

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
    parser.add_option("-o", "--output",
                      dest="fvcf", # type="string|int|float"
                      help="vcf output file name")
    parser.add_option("-t", "--threshold",
                      dest="threshold", type="float",
                      help="report variants if frequencey is bigger than this threshold")
    parser.add_option("-i", "--input",
                      dest="fpileup", # type="string|int|float",
                      default="-",
                      help="Pileup input (default: -, i.e. stdin)")

    
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
        sam.LOG.setLevel(logging.DEBUG)

    if not opts.threshold:
        parser.error("threshold argument missing\n")
        sys.exit(1)
    if opts.threshold > 1.0 or opts.threshold < 0.0:
        parser.error("threshold outside valid range\n")
        sys.exit(1)
    
    if not opts.fvcf:
        parser.error("vcf output file argument missing\n")
        sys.exit(1)
    if os.path.exists(opts.fvcf):
        sys.stderr.write(
            "Cowardly refusig to overwrite already existing file '%s'.\n" % (
                opts.fvcf))
        sys.exit(1)

    if opts.fpileup == '-':
        pileup_fhandle = sys.stdin
    else:
        if not os.path.exists(opts.fpileup):
            LOG.fatal("file '%s' does not exist.\n" % (opts.fpileup))
            sys.exit(1)
        else:
            pileup_fhandle = open(opts.fpileup, 'r')


    if opts.fvcf == '-':
        vcf_handle = sys.stdout
    else:
        vcf_handle = open(opts.fvcf, 'w')
    vcfwriter = VcfWriter(vcf_handle)
    vcfwriter.write_header(MYNAME)
    
    freq_threshold = opts.threshold
    num_skipped = 0
    num_total = 0
    for line in pileup_fhandle:
        pcol = sam.PileupColumn(line)
        num_total += 1

        # skip processing of columns with ambigious reference bases
        # hoping that gatk will skip them anyway
        if pcol.ref_base not in 'ACGT':
            LOG.info(
                "Skipping col %d because of ambigious reference base %s" % (
                    pcol.coord+1, pcol.ref_base))
            num_skipped += 1
            continue

        # count bases, ignore quality completely. GATK should be able to handle this
        base_counts = pcol.get_all_base_counts(min_qual=0)
        del base_counts['N']

        # sum over each base, and strand
        coverage = sum([sum(c) for c in base_counts.values()])
        if coverage == 0:
            continue

        for (alt_base, strand_counts) in base_counts.iteritems():
            count = sum(strand_counts)
            if alt_base == pcol.ref_base:
                continue
            freq = count/float(coverage)
            if freq > freq_threshold:
                LOG.info("Variant position above threshold: %s %d %s>%s %f" % (
                        pcol.chrom, pcol.coord, pcol.ref_base, alt_base, freq))
                # need: chrom pos ref alt[|alt..]
                var_pos = VarPos(pcol.chrom, pcol.coord,
                                 pcol.ref_base, alt_base)
                vcfwriter.write_variants([var_pos])

    if vcf_handle != sys.stdout:
        vcf_handle.close()
                
    if opts.fpileup != '-':
        pileup_fhandle.close()

    if num_total == 0:
        LOG.critical("No pileup information seen!?") 
        sys.exit(1)
    elif num_skipped/float(num_total) > 0.5:
        LOG.critical("%d%% of all pileup columns contained unusable information"
                     " (probably all refs were N?!)" % (num_skipped/float(num_total)*100)) 
        sys.exit(1)

        
if __name__ == "__main__":
    LOG.info("FIXME: use simple_vcf instead of VcfWriter")
    main()
    LOG.info("Successful program exit")


