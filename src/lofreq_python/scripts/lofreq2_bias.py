#!/usr/bin/env python
"""Experimental implementation of various bias filters
"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "GPL2"


#--- standard library imports
#
import sys
import logging
import os
import argparse
import gzip
from math import log

#--- third-party imports
#
# FIXME get rid of deps
import pysam
from scipy.stats import mannwhitneyu
from scipy.stats import chi2

#--- project specific imports
#
# sets PATH so that local scripts/binary is used if present, i.e.
# stuff can be run without installing it
try:
    import lofreq2_local
except ImportError:
    pass
from lofreq_star.utils import prob_to_phredqual
from lofreq_star import vcf


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                                        format='%(levelname)s [%(asctime)s]: %(message)s')


ALPHA = 0.05


def fisher_comb(pv1, pv2):
    """
    Fisher's method for combining p-values
    
    See for example
    http://en.wikipedia.org/wiki/Fisher's_method
    and
    breseq-0.18b:polymorphism_statistics.r
    """
    
    if pv1==0 or pv2==0:
        # not sure if this is correct.
        # see also http://stats.stackexchange.com/questions/58537/fishers-method-when-p-value-0
        return 0.0
    
    comb_log = -2.0 * (log(pv1) + log(pv2))
    # http://stackoverflow.com/questions/11725115/p-value-from-chi-sq-test-statistic-in-python
    comb_pv = 1.0 - chi2.cdf(comb_log, 4)    
    return comb_pv



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
                        help="Input BAM file matching vcf")
    parser.add_argument("-i", "--vcfin",
                        required=True,
                        help="Input VCF file containing variants to filter")
    parser.add_argument("-o", "--vcfout",
                        default = "-",
                        help="Output VCF")
    default = 13
    parser.add_argument("--mq-filter",
                        dest="min_mq",
                        type=int,
                        default=default,
                        help="Ignore reads with mapping quality below this value (default=%d)" % default)
    default = 6
    parser.add_argument("--bq-filter",
                        dest="min_bq",
                        type=int,
                        default=default,
                        help="Ignore bases with quality below this value (default=%d)" % default)
    parser.add_argument("-a", "--use-orphan",
                        action="store_true",
                        help="Don't ignore orphan-reads / anomalous read-pairs")
    parser.add_argument("-p", "--pass-only",
                        action="store_true",
                        help="Don't print variants filtered here")

    return parser


def main():
    """The main function
    """
    
    parser = cmdline_parser()
    args = parser.parse_args()
    
    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)
        import pdb
        from IPython.core import ultratb
        sys.excepthook = ultratb.FormattedTB(mode='Verbose',
                                             color_scheme='Linux', call_pdb=1)
        
    
    assert os.path.exists(args.bam), (
        "BAM file %s does not exist" % args.bam)
    samfh = pysam.Samfile(args.bam)


    # setup vcf_reader
    # 
    if args.vcfin == '-':
        vcf_reader = vcf.VCFReader(sys.stdin)
    else:
        if args.vcfin[-3:] == '.gz':
            vcf_reader = vcf.VCFReader(gzip.open(args.vcfin))
        else:
            vcf_reader = vcf.VCFReader(open(args.vcfin))
            
    variants = [r for r in vcf_reader]
    LOG.info("Loaded %d variants. Will use this for Bonferroni correction" % len(variants))
    bonf = len(variants)
    
    # setup vcf_writer
    #
    if args.vcfout == '-':
        fh_out = sys.stdout
    else:
        if os.path.exists(args.vcfout):
            LOG.fatal("Cowardly refusing to overwrite already existing file %s" % (args.vcfout))
            sys.exit(1)
            
        if args.vcfout[-3:] == '.gz':
            fh_out = gzip.open(args.vcfout, 'w')
        else:
            fh_out = open(args.vcfout, 'w')
    vcf_writer = vcf.VCFWriter(fh_out)
    vcf_writer.meta_from_reader(vcf_reader)

    outvars = []
    for var in variants:
        if var.INFO.has_key('INDEL'):
            LOG.warn("Skipping unsupported indel variant %s:%d" % (var.CHROM, var.POS))
            
        reads = list(samfh.fetch(reference=var.CHROM,
                                 start=var.POS-1, end=var.POS))
        LOG.info("%s %d: %d (unfiltered) reads covering position" % (
           var.CHROM, var.POS, len(reads)))

        ref_mquals = []
        alt_mquals = []
        ref_bquals = []
        alt_bquals = []
        # only for PE
        #ref_isize = []
        #alt_isize = []
        # following two meant to test
        #alt_vpos = [] 
        #rlens = []
        
        for r in reads:

            skip_flags = [0x4, 0x100, 0x200, 0x400]
            # FIXME combine
            for f in skip_flags:
                if r.flag & f:
                    continue
                
            orphan = (r.flag & 0x1) and not (r.flag & 0x2)
            if orphan and not args.use_orphan:
                continue

            if r.mapq < args.min_mq:
                continue
        
            vpos_on_read = [vpos_on_read 
                            for (vpos_on_read, vpos_on_ref) in r.aligned_pairs 
                            if vpos_on_ref==var.POS-1]
            assert len(vpos_on_read)==1
            vpos_on_read = vpos_on_read[0]
            if vpos_on_read == None:# skip deletions
                continue

            #alt_vpos.append(vpos_on_read)
            #rlens.append(r.rlen)
            
            b = r.query[vpos_on_read]
            bq = ord(r.qqual[vpos_on_read])-33
            mq = r.mapq

            if bq < args.min_bq:
                continue
            
            assert len(var.REF)==1 and len(var.ALT)==1
            if b.upper() == var.REF[0].upper():
                ref_mquals.append(mq)
                ref_bquals.append(bq)
                #if not args.use_orphan:
                #    ref_isize.append(abs(r.tlen))
            elif b.upper() == var.ALT[0].upper():
                alt_mquals.append(mq)
                alt_bquals.append(bq)
                #if not args.use_orphan:
                #    alt_isize.append(abs(r.tlen))
            else:            
                LOG.debug("Skipping non-ref-alt base %s at %s:%d" % (b, var.CHROM, var.POS))
                continue
            
        LOG.debug("After filtering at %s:%d: %d ref mquals and %d alt mquals" % (
            var.CHROM, var.POS, len(ref_mquals), len(alt_mquals)))
        
        # mannwhitneyu fails if all values the same
        if len(set(ref_mquals).union(alt_mquals))==1 or len(ref_mquals)==0 or len(alt_mquals)==0:
            m_pv = 1.0
        else:
            ustat = mannwhitneyu(ref_mquals, alt_mquals)
            m_pv = ustat[1]
            
        # same for bqs
        if len(set(ref_bquals).union(alt_bquals))==1 or len(ref_bquals)==0 or len(alt_bquals)==0:
            b_pv = 1.0
        else:
            ustat = mannwhitneyu(ref_bquals, alt_bquals)
            b_pv = ustat[1]

        # same for isize-qs
        #if len(ref_isize) and len(alt_isize):
        #    if len(set(ref_isize).union(alt_isize))==1:
        #        i_pv = 1
        #    else:
        #        ustat = mannwhitneyu(ref_isize, alt_isize)
        #        i_pv = ustat[1]
        #else:
        #    i_pv = 1
        
        c_pv = fisher_comb(m_pv, b_pv)
            
        #import pdb; pdb.set_trace()
        LOG.info("%s %d: mb %f bb %f cb %f" % (var.CHROM, var.POS, m_pv, b_pv, c_pv))

        var.INFO['MB'] = prob_to_phredqual(m_pv)
        var.INFO['BB'] = prob_to_phredqual(b_pv)
        #var.INFO['IB'] = prob_to_phredqual(i_pv)
        var.INFO['CB'] = prob_to_phredqual(c_pv)
        
        keep = True

        #import pdb; pdb.set_trace()
        for (pv, ftag) in [(m_pv, 'MB'), 
                           (b_pv, 'BB'),
                           (c_pv, 'CB'),
                           #(i_pv, 'IB')
                           ]:
            if pv*bonf < ALPHA:
                if var.FILTER == '.' or var.FILTER == 'PASS':
                    var = var._replace(FILTER=ftag)
                else:
                    var = var._replace(FILTER="%s;%s" % (var.FILTER, ftag))
                if args.pass_only:
                    keep = False
        if keep:
            outvars.append(var)

        #    LOG.info("Filtered: %s %d:\t%.8f\t%.8f" % (var.CHROM, var.POS, m_pv, b_pv))
        #else:
        #    LOG.info("Keeping: %s %d:\t%.8f\t%.8f" % (var.CHROM, var.POS, m_pv, b_pv))
        

    vcf_writer.write(outvars)
    if fh_out != sys.stdout:
        fh_out.close()
                        
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
