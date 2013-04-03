#!/usr/bin/env python
"""Perform set operations on two vcf-files (variant call format).
"""

#
# FIXME:Copyright
#

#--- standard library imports
#
from __future__ import division
import sys
import logging
import os
import argparse
import gzip

#--- third-party imports
#
#/

#--- project specific imports
#
# FIXME a hack to make not installed versions work as well
if True:
    import sys, os
    d = os.path.join(
        os.path.dirname(sys.argv[0]), '..')
    if os.path.exists(os.path.join(d, 'lofreq2')):
        sys.stderr.write("Adding local dir %s to PYTHONPATH\n" % d)
        sys.path.insert(0, d)

from lofreq2 import vcf

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



def get_vcfreader(vcffile):
    """gzip aware convenience wrapper for vcf.VCFReader
    """
    
    if vcffile[-3:] == '.gz':
        return vcf.VCFReader(gzip.open(vcffile, 'r'))
    else:
        return vcf.VCFReader(open(vcffile, 'r'))

        
def cmdline_parser():
    """FIXME:add-doc
    """

    # http://docs.python.org/dev/howto/argparse.html
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="be verbose")
    parser.add_argument("--debug", action="store_true",
                        help="enable debugging")
    parser.add_argument("--ign-filtered", action="store_true",
                        help="only consider passed i.e. un-filtered variants")
    parser.add_argument("-1", "--vcf1", 
                        help="1st vcf file (gzip supported)")
    parser.add_argument("-2", "--vcf2", 
                        help="2nd vcf file (gzip supported)")
    parser.add_argument("-a", "--action", choices=['intersect', 'complement'],
                        help="Set operation to perform. "
                        " intersect: vcf1 and vcf2."
                        " complement (rel.): vcf1 \ vcf2")
    parser.add_argument("-o", "--vcfout", default="-",
                        help="Output file or '-' for stdout (default)."
                        " Meta-data will be copied from vcf1")
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

    if not args.action:
        LOG.error("Missing action argument")
        #parser.print_help()
        sys.exit(1)

    for (in_file, descr) in [(args.vcf1, "1st vcf file"),
                             (args.vcf2, "2nd vcf file")]:
        if not in_file:
            LOG.error("%s input file argument missing." % descr)
            #parser.print_help()
            sys.exit(1)
        if not os.path.exists(in_file): # and in_file != "-":
            LOG.error("file '%s' does not exist.\n" % in_file)
            #parser.print_help()
            sys.exit(1)
            
    for (out_file, descr) in [(args.vcfout, "VCF output file")]:
        if not out_file:
            LOG.error("%s output file argument missing." % descr)
            #parser.print_help()
            sys.exit(1)
        if os.path.exists(out_file) and out_file!="-":
            LOG.error("Cowardly refusing to overwrite existing"
                      " output file '%s'.\n" % out_file)
            sys.exit(1)

    # ----------------------------------------------------------------
    # arg logic check done
    # ----------------------------------------------------------------
    
    vcf1_reader = get_vcfreader(args.vcf1)
    vcf2_reader = get_vcfreader(args.vcf2)

    if args.vcfout == '-':
        vcf_writer = vcf.VCFWriter(sys.stdout)
    else:
        fh_vcfout = open(args.vcfout, 'w')
        vcf_writer = vcf.VCFWriter(fh_vcfout)
    # meta-data copied from first vcf file
    vcf_writer.meta_from_reader(vcf1_reader)
    # vcf_writer.write(snvs) =
    # self.write_metainfo()
    # self.write_header()
    # for v in vars:
    #    self.write_rec(v)     
                                     
    LOG.warn("args.ign_filtered = %s" % args.ign_filtered)
    LOG.warn("args.action = %s" % args.action)

    
    # read B into memory and parse A one by one
    snvs2 = [r for r in vcf2_reader]
    # FIXME make dict
    # Record(CHROM='chr1', POS=334248, ID='.', REF='T', ALT=['C'], QUAL='.', FILTER='.', INFO={'SB': 0, 'DP4': [0, 0, 1, 0], 'CONSVAR': True, 'DP': 1, 'AF': 1.0}, FORMAT=None, samples=
    # assert len(ALT) == 1 
    # assert False if not rec.INFO.has_key(INDEL) else rec.INFO['INDEL']
    # warn if INDEL

    # relative complement : elements in A but not B

    import pdb; pdb.set_trace()
    for snv1 in vcf1_reader:
        LOG.critical("If snv2 present in snvs1 delete in snvs1")
    
    
#    if opts.vcf_out == '-':
#        fh_out = sys.stdout
#    else:
#        fh_out = open(opts.vcf_out, 'w')
#        
#    
#    if fh_out != sys.stdout:
#        fh_out.close()

    
if __name__ == "__main__":
    main()
    #LOG.info("Successful program exit")
