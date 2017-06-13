#!/usr/bin/env python
"""Experimental implementation of Mutect's "observed in control" AKA
alt_allele_in_normal filter:

From the Cibulskis (2013): Eliminate false positives in the tumor data
by look- ing at the control data (typically from the matched normal
sample) for evidence of the alternate allele beyond what is expected
from random sequencing error. A candidate is rejected if, in the
control data, there are (i) >= 2 observations of the alternate allele
or they represent >= 3% of the reads; and (ii) their sum of quality
scores is > 20.

Note, this only makes senseif you're working in similar coverage and
quality ranges.

"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "The MIT License"


#--- standard library imports
#
import sys
import logging
import os
import argparse
import gzip

#--- third-party imports
#
import pysam
import vcf


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')



FILTER_TAG = "alt_allele_in_normal"


def cmdline_parser():
    """Returns an argparse instance
    """

    # http://docs.python.org/dev/howto/argparse.html
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--verbose",
                        action="store_true",
                        help="Be verbose")
    parser.add_argument("--debug",
                        action="store_true",
                        help="Enable debugging")
    parser.add_argument("-b", "--bam",
                        required=True,
                        help="Normal BAM file")
    parser.add_argument("-v", "--vcfin",
                        required=True,
                        help="VCF file containing somatic variant"
                        " candidates to filter")
    parser.add_argument("-o", "--vcfout",
                        default = "-",
                        help="Output VCF")
    parser.add_argument("-p", "--pass-only",
                        action="store_true",
                        help="Don't print filtered variants")

    return parser


def skip_read(r):
    """Decide whether to skip a read

    FIXME identical copy in lofreq2_bias.py
    """
    
    skip_flags = [0x4, 0x100, 0x200, 0x400]
    skip = False
    # FIXME combine
    for f in skip_flags:
        if r.flag & f:
            return True
    return False


def main():
    """The main function
    """

    parser = cmdline_parser()
    args = parser.parse_args()

    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)

    assert os.path.exists(args.bam), (
        "BAM file %s does not exist" % args.bam)
    samfh = pysam.Samfile(args.bam)

    # setup vcf_reader
    #
    if args.vcfin[-3:] == '.gz':
        fh_in = gzip.open(args.vcfin)
        compressed = True
    else:
        compressed = False
        if args.vcfin == '-':
            fh_in = sys.stdin
        else:
            fh_in = open(args.vcfin)
    vcf_reader = vcf.VCFReader(fh_in, compressed)


    # setup vcf_writer
    #
    if args.vcfout == '-':
        fh_out = sys.stdout
    else:
        if os.path.exists(args.vcfout):
            LOG.fatal("Cowardly refusing to overwrite already existing"
                      " file %s" % (args.vcfout))
            sys.exit(1)

        if args.vcfout[-3:] == '.gz':
            fh_out = gzip.open(args.vcfout, 'w')
        else:
            fh_out = open(args.vcfout, 'w')

    # pyvcf needs template as arg to VCFWriter, whereas LoFreq's vcf
    # clone didn't
    vcf_writer = vcf.VCFWriter(fh_out, vcf_reader, lineterminator=os.linesep)
    #vcf_writer = vcf.VCFWriter(fh_out)
    #vcf_writer.meta_from_reader(vcf_reader)
    # FIXME should add filter description to header

    for (var_no, var) in enumerate(vcf_reader):
        if var_no % 500 == 1:
            LOG.info("Analyzing variant %d" % (var_no))

        if 'INDEL' in var.INFO:
            LOG.warn("Skipping indel %s:%d" % (var.CHROM, var.POS))
            continue
        if len(var.REF)>1 or len(var.ALT)>1:
            LOG.warn("Skipping indel (not tagged as such) %s:%d" % (
                var.CHROM, var.POS))
            continue


        reads = list(samfh.fetch(reference=var.CHROM,
                                 start=var.POS-1, end=var.POS))
        LOG.debug("%s %d: %d (unfiltered) reads covering position" % (
           var.CHROM, var.POS, len(reads)))

        ref_bquals = []
        alt_bquals = []

        # FIXME huge code overlap with lofreq2_bias.py
        for r in reads:

            if skip_read(r):
                continue
            
            # determine position on read for variant to then determine
            # the current base and its basequal
            #
            vpos_on_read = [vpos_on_read
                            for (vpos_on_read, vpos_on_ref) in r.aligned_pairs
                            if vpos_on_ref==var.POS-1]
            #if False:
            #    if len(vpos_on_read)!=1:
            #        #import pdb; pdb.set_trace()
            #        from IPython import embed; embed()
            assert len(vpos_on_read)==1
            vpos_on_read = vpos_on_read[0]
            if vpos_on_read == None:# skip deletions
                continue

            b = r.query[vpos_on_read]
            bq = ord(r.qqual[vpos_on_read])-33

            assert len(var.REF)==1 and len(var.ALT)==1
            if b.upper() == var.REF[0].upper():
                ref_bquals.append(bq)
            elif b.upper() == str(var.ALT[0]).upper():
                alt_bquals.append(bq)
            else:
                LOG.debug("Skipping non-ref-alt base %s at %s:%d" % (
                    b, var.CHROM, var.POS))
                continue

        # " A candidate is rejected if, in the control data, there are
        # (i) >= 2 observations of the alternate allele or they represent
        # >= 3% of the reads; and (ii) their sum of quality scores is >=
        # 20."
        # FIXME set filter var.INFO['AN'] = True
        print_this_var = True
        num_alt = len(alt_bquals)
        num_ref = len(ref_bquals)
        num_both = num_alt+num_ref
        if num_both==0:
            LOG.warn("No alt or ref bases for var %s" % var)
            print_this_var = True
        else:
            if (num_alt>=2 or num_alt/float(num_both)>=0.03) and sum(alt_bquals)>20:
                var.FILTER.append(FILTER_TAG)
                if args.pass_only:
                    print_this_var = False
        if print_this_var:
            # LoFreq's vcf clone called this write_rec()
            vcf_writer.write_record(var)

    if fh_in != sys.stdout:
        fh_in.close()
    if fh_out != sys.stdout:
        fh_out.close()

if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
