#!/usr/bin/env python
"""Perform 'set' operations on two vcf-files (variant call format).

Two SNV are regarded identical if their chromosome, their position and
their (ref and) alt base are identical. Note, this definition differs
from vcftools were the bases are ignored.

The VCF meta information will be copied from first file.

The exact SNV annotation will always be taken from SNV coming from
first file.
"""

#
# FIXME:Copyright
#

#--- standard library imports
#
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
try:
    import lofreq2_local
except:
    pass    

try:
    from lofreq2 import vcf
except:
    sys.stderr.write("FATAL: Couldn't LoFreq's vcf module."
                     " Are you sure your PYTHONPATH is set correctly?\n")
    sys.exit(1)
    
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


def key_for_var(var):
    """FIXME:add-doc"""

    return "%s %d %s %s" % (var.CHROM, var.POS, 
                            var.REF, ''.join(var.ALT))


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

    parser.add_argument("-v", "--verbose", 
                        action="store_true",
                        help="be verbose")
    parser.add_argument("--debug", 
                        action="store_true",
                        help="enable debugging")
    parser.add_argument("--ign-filtered", 
                        action="store_true",
                        help="only consider passed i.e. un-filtered variants")
    parser.add_argument("-1", "--vcf1", 
                        required=True,
                        help="1st vcf file (gzip supported)")
    parser.add_argument("-2", "--vcf2", 
                        required=True,
                        help="2nd vcf file (gzip supported)")
    parser.add_argument("-a", "--action", 
                        required=True, 
                        choices=['intersect', 'complement'],
                        help="Set operation to perform. "
                        " intersect: vcf1 and vcf2."
                        " complement (rel.): vcf1 \ vcf2")
    parser.add_argument("-o", "--vcfout", 
                        default="-",
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

    LOG.debug("args = %s" % args)

    # ----------------------------------------------------------------
    # arg logic check done
    # ----------------------------------------------------------------
    
    vcf1_reader = get_vcfreader(args.vcf1)
    vcf2_reader = get_vcfreader(args.vcf2)

    if args.vcfout == '-':
        fh_vcfout = sys.stdout
    else:
        fh_vcfout = open(args.vcfout, 'w')
    vcf_writer = vcf.VCFWriter(fh_vcfout)
    # meta-data copied from first vcf file
    vcf_writer.meta_from_reader(vcf1_reader)
    # FIXME should we add ourselve as source just like the vcftools folks do?
    vcf_writer.write_metainfo()
    vcf_writer.write_header()

    #
    # recipe:read B into memory and parse from A one by one
    #
    
    snvs2 = dict()
    for var in vcf2_reader:
        if args.ign_filtered and var.FILTER not in ['.', 'PASS']:
            continue
            
        assert len(var.ALT) == 1, (
            "Can't handle more then one alt base" 
            " (doesn't look like this file came from LoFreq)"
            " and therefore can't process: %s" % str(var))
        if var.INFO.has_key('INDEL'):
            assert not var.INFO['INDEL'], (
                "Can't handle indels and therefore can't process"
                " : %s" % str(var))
        
        k = key_for_var(var)
        assert not snvs2.has_key(k), (
            "I'm confused. Looks like I've already seen a SNV with"
            " key %s" % k)
        snvs2[k] = var        

    for var in vcf1_reader:
        if args.ign_filtered and var.FILTER not in ['.', 'PASS']:
            continue
        k = key_for_var(var)

        if args.action == 'complement':
            # relative complement : elements in A but not B
            if not snvs2.has_key(k):
                vcf_writer.write_rec(var)
            else:
                del snvs2[k]
        elif args.action == 'intersect':
            if snvs2.has_key(k):
                vcf_writer.write_rec(var)
        else:
            raise ValueError
            

    if fh_vcfout != sys.stdout:
        fh_vcfout.close()

    
if __name__ == "__main__":
    main()
    #LOG.info("Successful program exit")
