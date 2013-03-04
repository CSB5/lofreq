#!/usr/bin/env python
"""Apply number of filters to given list of SNVs
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
from optparse import OptionParser, SUPPRESS_HELP
from collections import namedtuple



#--- third-party imports
#
#/

#--- project specific imports
#
from lofreq2 import vcf
from lofreq2 import multiple_testing
from lofreq2.utils import prob_to_phredqual, phredqual_to_prob

# invocation of ipython on exceptions
#import sys, pdb
#from IPython.core import ultratb
#sys.excepthook = ultratb.FormattedTB(mode='Verbose',
#                                     color_scheme='Linux', call_pdb=1)


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
      + "\n" + __doc__ + "NOTE: filters are applied in there order of appearance here\n" # FIXME
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true", dest="debug",
                      help="enable debugging")
    DEFAULT = "-"
    parser.add_option("-i", "--vcf_in",
                      dest="vcf_in",
                      default=DEFAULT,
                      help="Input vcf file (- for stdin). Default = %s)" % DEFAULT)
    DEFAULT = "-"
    parser.add_option("-o", "--outfile",
                      dest="vcf_out",
                      default=DEFAULT,
                      help="Output vcf file (- for stdout). Default = %s)" % DEFAULT)

    parser.add_option("", "--strandbias-bonf",
                      dest="strandbias_bonf", 
                      action="store_true",
                      help="Optional: Bonferroni corrected strand bias"
                      " value has to be 0.05 (applied first!)")
    parser.add_option("", "--strandbias-holmbonf",
                      dest="strandbias_holmbonf", 
                      action="store_true",
                      help="Optional: Holm-Bonferroni corrected strand"
                      " bias value has to be 0.05 (applied first!)")
    parser.add_option("", "--strandbias-phred",
                      dest="max_strandbias_phred", 
                      type='int',
                      help="Optional: Ignore SNVs with strandbias"
                      " phred-score above this value")
    parser.add_option("", "--min-freq",
                      dest="min_freq", 
                      type="float",
                      help="Optional: Ignore SNVs below this freq treshold")
    parser.add_option("", "--max-cov",
                      dest="max_cov", 
                      type='int',
                      #default=sys.maxint, type=int,
                      help="Optional: Ignore variant sites if coverage above this cap")
    parser.add_option("", "--min-cov",
                      dest="min_cov", 
                      type='int',
                      #default=1, type=int,
                      help="Optional: Ignore variants sites if coverage"
                      " is below this value")    
    parser.add_option("", "--snp-phred", 
                      dest="min_snp_phred",
                      type='int',
                      #default=1, type=int,
                      help="Optional: Ignore variants sites if SNP"
                      " phred score is below this (float) value")
    parser.add_option("", "--window-size",
                      dest="window_size",
                      type='int',
                      help="Optional: Ignore any variants if at least"
                      " one more was called within this window size")

    parser.add_option("--force", help=SUPPRESS_HELP,
                      dest="force_overwrite", action="store_true") 

    return parser



def main():

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

    for (in_file, descr) in [(opts.vcf_in, "VCF input file")]:
        if not in_file:
            parser.error("%s input file argument missing." % descr)
            sys.exit(1)
        if not os.path.exists(in_file) and in_file != "-":
            sys.stderr.write(
                "file '%s' does not exist.\n" % in_file)
            sys.exit(1)
            
    for (out_file, descr) in [(opts.vcf_out, "VCF output file")]:
        if not out_file:
            parser.error("%s output file argument missing." % descr)
            sys.exit(1)
        if os.path.exists(out_file) and out_file!="-":
            sys.stderr.write(
                "Cowardly refusing to overwrite existing output file '%s'.\n" % out_file)
            sys.exit(1)

    if opts.vcf_in == '-':
        vcf_reader = vcf.VCFReader(sys.stdin)
    else:
        vcf_reader = vcf.VCFReader(open(opts.vcf_in,'r'))
    snvs = [r for r in vcf_reader]
    LOG.info("Parsed %d SNVs from %s" % (len(snvs), opts.vcf_in))

    
    # list of tuples: first element is a lambda test with a snp as
    # input (keep if lambda returns 1), second one is description
    tests = []

    # strand-bias [holm]/bonf have to be first filters
    
    if opts.strandbias_bonf:
        assert opts.max_strandbias_phred == None and opts.strandbias_holmbonf == None, (
            "Can't filter strand bias twice")

        vcf_filter = vcf._Filter(
            id="sbb", 
            desc="Strand-bias filter on Bonferroni corrected p-values")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer
        vcf_info = vcf._Info(
            id="SBBC", num=1, type='Integer',
            desc="Strand-bias Bonferroni corrected")        
        vcf_reader.infos[vcf_info.id] = vcf_info
        
        pvals = [phredqual_to_prob(s.INFO['SB']) for s in snvs]
        corr_pvals = multiple_testing.Bonferroni(pvals).corrected_pvals
        for (cp, s) in zip(corr_pvals, snvs):
            s.INFO[vcf_info.id] = prob_to_phredqual(cp)
            # FIXME PHREDQUAL_TO_PROB in snspcaller.c has no limit
            if s.INFO[vcf_info.id] > sys.maxint:
                s.INFO[vcf_info.id] = sys.maxint
                
        tests.append((
            lambda s, f: f.id if s.INFO["SBBC"] > 13 else None,
            vcf_filter
            ))
    
    if opts.strandbias_holmbonf:
        assert opts.max_strandbias_phred == None and opts.strandbias_bonf == None, (
            "Can't filter strand bias twice")

        vcf_filter = vcf._Filter(
            id="sbh", 
            desc="Strand-bias filter on Holm-Bonferroni corrected p-values")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer
        vcf_info = vcf._Info(
            id="SBHC", num=1, type='Integer',
            desc="Strand-bias Holm-Bonferroni corrected")        
        vcf_reader.infos[vcf_info.id] = vcf_info
        
        pvals = [phredqual_to_prob(s.INFO['SB']) for s in snvs]
        corr_pvals = multiple_testing.HolmBonferroni(pvals).corrected_pvals
        for (cp, s) in zip(corr_pvals, snvs):
            s.INFO[vcf_info.id] = prob_to_phredqual(cp)
            # FIXME PHREDQUAL_TO_PROB in snspcaller.c has no limit
            if s.INFO[vcf_info.id] > sys.maxint:
                s.INFO[vcf_info.id] = sys.maxint
            
        tests.append((
            lambda s, f: f.id if s.INFO["SBHC"] > 13 else None,
            vcf_filter
            ))                
        
    if opts.max_strandbias_phred != None:
        vcf_filter = vcf._Filter(
            id="sbp", 
            desc="Phred-based strand-bias filter (max)")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

        tests.append((
            lambda s, f: f.id if float(s.INFO['SB']) > opts.max_strandbias_phred else None,
            vcf_filter
            ))
        
    if opts.min_freq != None:  
        vcf_filter = vcf._Filter(
            id="minaf", 
            desc="Minimum allele frequency")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

        tests.append((
            lambda s, f: f.id if s.INFO['AF'] < opts.min_freq else None,
            vcf_filter
            ))

    if opts.max_cov != None:  
        vcf_filter = vcf._Filter(
            id="maxcov", 
            desc="Maximum coverage")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

        tests.append((
            lambda s, f: f.id if s.INFO['DP'] > opts.max_cov else None,
            vcf_filter
            ))

    if opts.min_cov != None:  
        vcf_filter = vcf._Filter(
            id="mincov", 
            desc="Minimum coverage")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

        tests.append((
            lambda s, f: f.id if s.INFO['DP'] < opts.min_cov else None,
            vcf_filter
            ))

    if opts.min_snp_phred != None:  
        vcf_filter = vcf._Filter(
            id="minqual", 
            desc="Minimum SNV quality")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

        tests.append((
            lambda s, f: f.id if s.QUAL != '.' and s.QUAL < opts.min_snp_phred else None,
            vcf_filter
            ))

    LOG.error("TEST TEST TEST")
         
    # The actual filtering:
    #
    # LOG.info("Will perform the following tests: \n%s" % (
    # '\n'.join(["- %s" % t[1] for t in tests])))
    # can't this be done with map()?
    for (test_func, test_filt) in tests:
        for (i, s) in enumerate(snvs):
            f = test_func(s, test_filt)
            if f:
                # just s = s.__replace() can't work
                if s.FILTER == '.' or s.FILTER == 'PASS':
                    snvs[i] = s._replace(FILTER=f)
                else:
                    snvs[i] = s._replace(FILTER="%s,%s" % (s.FILTER, f))

    if opts.window_size != None:  
        raise NotImplementedError # FIXME

    LOG.error("FIXME test each filter")

    for (i, s) in enumerate(snvs):
        if s.FILTER == '.':
            snvs[i] = s._replace(FILTER="PASS")

    # FIXME LOG.info("%d SNVS survived all filters. Writing to %s" % (
    #    len(snvs), opts.vcf_out))
    if opts.vcf_out == '-':
        fh_out = sys.stdout
    else:
        fh_out = open(opts.vcf_out, 'w')

    vcf_writer = vcf.VCFWriter(sys.stdout)
    vcf_writer.meta_from_reader(vcf_reader)
    LOG.error("Add corrected values as info fields and also filters")
    vcf_writer.write(snvs)
    if fh_out != sys.stdout:
        fh_out.close()

    
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
