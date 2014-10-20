#!/usr/bin/env python
"""FIXME:add-doc
"""

#
# FIXME:Copyright
#


#--- standard library imports
#
import sys
import os
import argparse
import logging
import itertools

#--- third-party imports
#
#/

#--- project specific imports
#
from lofreq2 import vcf


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2013 Genome Institute of Singapore"
__license__ = "GPL2"



#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')



def cmdline_parser():
    """FIXME:add-doc
    """

    # http://docs.python.org/dev/howto/argparse.html
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--verbose", 
                        action="store_true",
                        help="be verbose")
    parser.add_argument("--debug", 
                        action="store_true",
                        help="enable debugging")
    parser.add_argument("--ign-filtered", 
                        action="store_true",
                        help="only consider passed i.e. un-filtered variants")
    parser.add_argument("-v", "--vcf", 
                        required=True,
                        help="VCF file (gzip supported)")
    parser.add_argument("-b", "--bam", 
                        required=True,
                        help="BAM file")
    parser.add_argument("-l", "--readlen", 
                        required=True,
                        type=int,
                        help="read length")
    return parser



def main():
    """FIXME:add-doc
    """
    
    parser = cmdline_parser()
    args = parser.parse_args()
    
    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)

    for (in_file, descr) in [(args.vcf, "VCF file"),
                             (args.bam, "BAM file")]:
        if not os.path.exists(in_file): # and in_file != "-":
            LOG.error("%s '%s' does not exist.\n" % (descr, in_file))
            #parser.print_help()
            sys.exit(1)
            
    
    if args.vcf == '-':
        vcf_reader = vcf.VCFReader(sys.stdin)
    else:
        vcf_reader = vcf.VCFReader(open(args.vcf,'r'))
        
    snvs = dict()
    for vcf_rec in vcf_reader:
        if snvs.has_key(vcf_rec.CHROM):
            snvs[vcf_rec.CHROM].append(vcf_rec)
        else:
            snvs[vcf_rec.CHROM] = [vcf_rec]

    chroms = snvs.keys()
    for c in chroms:
        LOG.info("Parsed %d SNVs from %s for chrom %s" % (len(snvs[c]), args.vcf, c))

    for c in chroms:
        snvs[c] = sorted(snvs[c], key=lambda x: x.POS)
        # http://stackoverflow.com/questions/464864/python-code-to-pick-out-all-possible-combinations-from-a-list
        for (i, s1) in enumerate(snvs[c]):
            for (j, s2) in enumerate(snvs[c][i+1:]):
                dist = s2.POS - s1.POS
                assert dist>=0
                if dist > args.readlen: # that's why we need them sorted
                    break
                print "Possible pair within read length: chrom %s pos %d and %d" % (c, s1.POS, s2.POS)
                print "paired_nuc_freqs_from_bam.py -1 %d -2 %d -r %s -b %s" %  (s1.POS, s2.POS, c, args.bam)
    

    
if __name__ == "__main__":
    main()
    #LOG.info("Successful program exit")
