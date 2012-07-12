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
from optparse import OptionParser

#--- third-party imports
#
#/


#--- project specific imports
#
from lofreq import snp
from lofreq import multiple_testing

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
            + "\n" + __doc__
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true", dest="debug",
                      help="enable debugging")
    parser.add_option("-i", "--snp_infile",
                      dest="snp_infile",
                      help="SNP input file (- for stdin)")
    parser.add_option("-o", "--outfile",
                      dest="snp_outfile",
                      default="-",
                      help="SNP output file (- for stdout = default)")

    parser.add_option("", "--strandbias-bonf",
                      dest="strandbias_bonf", 
                      action="store_true",
                      help="Optional: Bonferroni corrected strand bias value has to be 0.05 (applied first!)")
    parser.add_option("", "--strandbias-holmbonf",
                      dest="strandbias_holmbonf", 
                      action="store_true",
                      help="Optional: Holm-Bonferroni corrected strand bias value has to be 0.05 (applied first!)")
    parser.add_option("", "--strandbias-phred",
                      dest="max_strandbias_phred", 
                      type='int',
                      help="Optional: Ignore SNPs with strandbias phred score above this value")
    parser.add_option("", "--min-freq",
                      dest="min_freq", 
                      type="float",
                      help="Optional: Ignore SNPs below this freq treshold")
    parser.add_option("", "--max-cov",
                      dest="max_cov", 
                      type='int',
                      #default=sys.maxint, type=int,
                      help="Optional: Ignore variant sites if coverage above this cap")
    parser.add_option("", "--min-cov",
                      dest="min_cov", 
                      type='int',
                      #default=1, type=int,
                      help="Optional: Ignore variants sites if coverage is below this value")    
    parser.add_option("", "--snp-phred", 
                      dest="min_snp_phred",
                      type='int',
                      #default=1, type=int,
                      help="Optional: Ignore variants sites if SNP phred score is below this (float) value")
    parser.add_option("", "--window-size",
                      dest="window_size",
                      type='int',
                      help="Optional: Ignore any variants if at least one more was called within this window size")

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

    for (in_file, descr) in [(opts.snp_infile, "SNP file")]:
        if not in_file:
            parser.error("%s input file argument missing." % descr)
            sys.exit(1)
        if not os.path.exists(in_file) and in_file != "-":
            sys.stderr.write(
                "file '%s' does not exist.\n" % in_file)
            sys.exit(1)
            
    for (out_file, descr) in [(opts.snp_outfile, "SNP output file")]:
        if not out_file:
            parser.error("%s output file argument missing." % descr)
            sys.exit(1)
        if os.path.exists(out_file) and out_file!="-":
            sys.stderr.write(
                "Cowardly refusing to overwrite existing output file '%s'.\n" % out_file)
            sys.exit(1)


    snps = snp.parse_snp_file(opts.snp_infile)
    LOG.info("Parsed %d SNPs from %s" % (len(snps), opts.snp_infile))

    
    # list of tuples: first element is a lambda test with a snp as
    # input (keep if lambda returns 1), second one is description
    tests = []

    # strand-bias [holm]/bonf have to be first filters
    
    if opts.strandbias_bonf:
        assert opts.max_strandbias_phred == None and opts.strandbias_holmbonf == None, (
            "Can't filter strand bias twice")
        pvals = [float(s.info['strandbias-pvalue-uncorr']) for s in snps]
        corr_pvals = multiple_testing.Bonferroni(pvals).corrected_pvals
        for (cp, s) in zip(corr_pvals, snps):
            s.info['strandbias-pvalue-corr'] = str(cp)
        tests.append((
            lambda s: float(s.info['strandbias-pvalue-corr']) > 0.05,
            "Bonferroni corrected strandbias phred-value"
            ))
        
    if opts.strandbias_holmbonf:
        assert opts.max_strandbias_phred == None and opts.strandbias_bonf == None, (
            "Can't filter strand bias twice")
        pvals = [float(s.info['strandbias-pvalue-uncorr']) for s in snps]
        corr_pvals = multiple_testing.HolmBonferroni(pvals).corrected_pvals
        for (cp, s) in zip(corr_pvals, snps):
            s.info['strandbias-pvalue-corr'] = str(cp)
        tests.append((
            lambda s: float(s.info['strandbias-pvalue-corr']) > 0.05,
            "Holm-Bonferroni corrected strandbias phred-value"
            ))                
        
    if opts.max_strandbias_phred != None:
        tests.append((
            lambda s: float(s.info['strandbias-pvalue-uncorr-phred']) <= opts.max_strandbias_phred,
            "maximum strandbias phred-value"
            ))
        
    if opts.min_freq != None:  
        tests.append((
            lambda s: s.freq >= opts.min_freq,
            "minimum frequency"
            ))

    if opts.max_cov != None:  
        tests.append((
            lambda s: int(s.info['coverage']) <= opts.max_cov,
            "maximum coverage"
            ))

    if opts.min_cov != None:  
        tests.append((
            lambda s: int(s.info['coverage']) >= opts.min_cov,
            "minimum coverage"
            ))

    if opts.min_snp_phred != None:  
        tests.append((
            lambda s: int(s.info['pvalue-phred']) >= opts.min_snp_phred,
            "minimum SNP phred"
            ))

                

    # The actual filtering:
    #
    # LOG.info("Will perform the following tests: \n%s" % (
    # '\n'.join(["- %s" % t[1] for t in tests])))
    for (test_lambda, test_descr) in tests:
        snps = [s for s in snps if test_lambda(s)]
        LOG.info("%d SNPs left after applying %s filter" % (
            len(snps), test_descr))


    if opts.window_size != None:  
        raise NotImplementedError # FIXME


        
    LOG.info("%d SNPs survived all filters. Writing to %s" % (
        len(snps), opts.snp_outfile))
    if opts.snp_outfile == '-':
        fh_out = sys.stdout
    else:
        fh_out = open(opts.snp_outfile, 'w')
    snp.write_snp_file(fh_out, snps)
    if fh_out != sys.stdout:
        fh_out.close()

    
if __name__ == "__main__":
    LOG.critical("WARN: phred value filtering might fail on type consensus-var (opposed to low-freq-var)")
    main()
    LOG.info("Successful program exit")
