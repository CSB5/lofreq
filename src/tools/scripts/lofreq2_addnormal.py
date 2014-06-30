#!/usr/bin/env python
"""Complement VCF with simple pileup info from BAM
"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "GPL2"


# --- standard library imports
#
import sys
import os
import argparse
import logging
import copy
from collections import OrderedDict, namedtuple
import csv
import gzip

#--- third-party imports
#
import vcf
import pysam

#--- project specific imports
#
# /


# global logger
#
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


def median(x):
    """compute median of provided list"""

    if not len(x):
        return None
    # http://stackoverflow.com/questions/10482339/how-to-find-median/10482422#10482422 answer by user3100512
    return sorted(x)[len(x)//2]


def cmdline_parser():
    """
    creates an OptionParser instance
    """

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--verbose",
                      action="store_true",
                      dest="verbose",
                      help="be verbose")
    parser.add_argument("--debug",
                      action="store_true",
                      dest="debug",
                      help="enable debugging")
    parser.add_argument("-i", "--vcf-in",
                      dest="vcf_in",
                      required=True,
                      help="Input vcf file listing somatic variants"
                      " (gzip supported; - for stdin).")
    default = "-"
    parser.add_argument("-o", "--vcf-out",
                      dest="vcf_out",
                      default=default,
                      help="Output vcf file (gzip supported; - for stdout;"
                      " default: %s)." % default)
    parser.add_argument("-b", "--bam",
                        dest="bam",
                        required=True,
                        help="BAM file, e.g. from normal sample")
    return parser


def add_plp_to_vcf(vcf_in, vcf_out, bam):
    """process each var in vcf_reader and add plp info from sam_fh,
    writing to vcf_writer (both pyvcf classes)
    """

    assert os.path.exists(bam)
    sam_fh = pysam.Samfile(bam)

    # set up vcf_reader
    #
    if vcf_in == '-':
        vcf_reader = vcf.VCFReader(sys.stdin)
    else:
        assert os.path.exists(vcf_in)
        vcf_reader = vcf.VCFReader(filename=vcf_in)

        LOG.critical("Need to add 'FORMAT samples...' to #CHROM line")
    
    vcf_reader.formats['TD'] = vcf.parser._Format(
        id='TD', num=1, type='Integer', desc='Total depth')
    vcf_reader.formats['NR'] = vcf.parser._Format(
        id='NR', num=1, type='Integer', desc='Number of reference bases')
    vcf_reader.formats['NA'] = vcf.parser._Format(
        id='NA', num=1, type='Integer', desc='Number of alternate bases')
    vcf_reader.formats['OR'] = vcf.parser._Format(
        id='OR', num=1, type='Integer', desc='Number of orphan reads supporting reference bases')
    vcf_reader.formats['OA'] = vcf.parser._Format(
        id='OA', num=1, type='Integer', desc='Number of orphan reads supporting alternate bases')
    # FIXME num=3 for this and following entries?
    vcf_reader.formats['BR'] = vcf.parser._Format(
        id='BR', num=None, type='Integer', desc='Minimum, median and maximum base-qualities for reference bases')
    vcf_reader.formats['BA'] = vcf.parser._Format(
        id='BA', num=None, type='Integer', desc='Minimum, median and maximum base-qualities for alternate bases')
    vcf_reader.formats['MR'] = vcf.parser._Format(
        id='MR', num=None, type='Integer', desc='Minimum, median and maximum mapping-qualities for reference bases')
    vcf_reader.formats['MA'] = vcf.parser._Format(
        id='MA', num=None, type='Integer', desc='Minimum, median and maximum mapping-qualities for alternate bases')
    
        
    # set up vcf_writer/fh_out
    #
    if vcf_out == '-':
        fh_out = sys.stdout
    else:
        assert not os.path.exists(vcf_out)
        if vcf_out[-3:] == ".gz":            
            fh_out = gzip.open(vcf_out, 'wb')
        else:
            fh_out = open(vcf_out, 'wb')
    #vcf_writer = vcf.VCFWriter(sys.stdout, vcf_reader,
    #                           lineterminator=os.linesep)
    vcf_out = csv.writer(fh_out, delimiter='\t')

    import pdb; pdb.set_trace()
    # using csv for writing similar to vcf_melt script in pyvcf
    formats = vcf_reader.formats.keys()
    infos = vcf_reader.infos.keys()
    header = ["SAMPLE"] + formats + ['FILTER', 'CHROM', 'POS', 'REF', 'ALT', 'ID'] \
      + ['info.' + x for x in infos]
    vcf_out.writerow(header)

    
    # process each variant
    #
    for var in vcf_reader:
        # no support for indels
        if var.INFO.has_key('INDEL') or len(var.REF) > 1 or len(var.ALT) > 1:
            LOG.warn("Skipping unsupported variant) %s:%d:%s" % (
                var.CHROM, var.POS, var.REF))
            continue

        # pilvcf module works with unit offset. pysam with zero offset
        for plp_col in sam_fh.pileup(var.CHROM, var.POS-1, var.POS):
            # pileup() extracts all reads overlapping that region.
            # only look at the one of interest
            if plp_col.pos != var.POS-1:
                continue

            bqs = {'ref': [], 'alt': []}
            mqs = {'ref': [], 'alt': []}
            num_orphans = {'ref': 0, 'alt': 0}

            for plp_read in plp_col.pileups:
                aln_read = plp_read.alignment
                # most minimal filtering
                if aln_read.is_unmapped or aln_read.is_secondary or \
                   aln_read.is_qcfail or aln_read.is_duplicate:
                    continue

                if aln_read.is_paired and aln_read.mate_is_unmapped:
                    assert not aln_read.is_unmapped
                    is_orphan = True
                else:
                    is_orphan = False

                base = aln_read.seq[plp_read.qpos]
                mq = aln_read.mapq
                bq = ord(aln_read.qual[plp_read.qpos])-33

                if base == var.REF:
                    k = 'ref'
                elif base == var.ALT[0]:
                    k = 'alt'
                else:
                    continue

                bqs[k].append(bq)
                mqs[k].append(mq)
                if is_orphan:
                    num_orphans[k] += 1

            (min_bqs, median_bqs, max_bqs) = (
                {'ref': -1, 'alt': -1},
                {'ref': -1, 'alt': -1},
                {'ref': -1, 'alt': -1})
            (min_mqs, median_mqs, max_mqs) = (
                {'ref': -1, 'alt': -1},
                {'ref': -1, 'alt': -1},
                {'ref': -1, 'alt': -1})

            for k in ['ref', 'alt']:
                if len(bqs[k]):
                    (min_bqs[k], median_bqs[k], max_bqs[k]) = (
                        min(bqs[k]), median(bqs[k]), max(bqs[k]))
                if len(mqs[k]):
                    (min_mqs[k], median_mqs[k], max_mqs[k]) = (
                        min(mqs[k]), median(mqs[k]), max(mqs[k]))
        LOG.critical("TD:%d NR=%d NA=%d OR=%d OA=%d BR=%d,%d,%d BA=%d,%d,%d MR=%d,%d,%d MA=%d,%d,%d" % (
            plp_col.n, len(bqs['ref']), len(bqs['alt']),
            num_orphans['ref'], num_orphans['alt'],
            min_bqs['ref'], median_bqs['ref'], max_bqs['ref'],
            min_bqs['alt'], median_bqs['alt'], max_bqs['alt'],
            min_mqs['ref'], median_mqs['ref'], max_mqs['ref'],
            min_mqs['alt'], median_mqs['alt'], max_mqs['alt']))

        # there is no way to edit/add format fields in current
        # versions of pyvcf. see discussion here
        # https://github.com/jamescasbon/PyVCF/issues/82 for patches
        # and workarounds

        
# # FIXME
#        if not var.FORMAT:
#            var.FORMAT = ""
#        call_data = dict()
        add_vals = OrderedDict()
        for (fmt_key, val) in [
                ('TD', plp_col.n),
                ('NR', len(bqs['ref'])),
                ('NA', len(bqs['alt'])),
                ('OR', num_orphans['ref']),
                ('OA', num_orphans['alt']),
                ('BR', "%d,%d,%d" % (min_bqs['ref'], median_bqs['ref'], max_bqs['ref'])),
                ('BA', "%d,%d,%d" % (min_bqs['alt'], median_bqs['alt'], max_bqs['alt'])),
                ('MR', "%d,%d,%d" % (min_mqs['ref'], median_mqs['ref'], max_mqs['ref'])),
                ('MA', "%d,%d,%d" % (min_mqs['alt'], median_mqs['alt'], max_mqs['alt']))]:
            add_vals[fmt_key] = val
            
        #new_var = copy.deepcopy(var)
        #new_var.samples.appendc(ollections.namedtuple('CallData', add_vals.keys()))
        #new_var.samples[sx].data = new_record.samples[sx].data._make(new_vals)
        
        import pdb; pdb.set_trace()
        # finally set CallData
        var.samples.append(namedtuple('Call', namedtuple('CallData', add_vals.keys())))

#            # need to add FORMAT sample sample
#            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#            # >>> vars[0].samples
#            # [Call(sample=SPIKEIN, CallData(GT=0/1))]
#            call_data[fmt_key] = val
#            
#        # undocumented add_format always preprends ":"
#        var.FORMAT += ':'.join(call_data.keys())
#
#        samp_fmt = make_calldata_tuple(samp_fmt.split(':'))
#
#        # FIXME need to see vars[0].FORMAT
#        ## parser.py
#        ## create a call object
#        ## call = _Call(site, name, samp_fmt(*sampdat))
#        call = vcf.model._Call(var, os.path.basename(bam), call_data)
#        import pdb; pdb.set_trace()
#        #var.samples.append(call)

        vcf_writer.write_record(var)
        
    vcf_writer.close()


def main():
    """main function
    """

    parser = cmdline_parser()
    args = parser.parse_args()

    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)

    for (in_file, descr) in [(args.bam, "BAM"),
                             (args.vcf_in, "VCF input")]:
        if not in_file:
            parser.error("%s file argument missing." % descr)
            sys.exit(1)
        if not os.path.exists(in_file) and in_file != "-":
            LOG.fatal("file '%s' does not exist.\n" % in_file)
            sys.exit(1)

    for (out_file, descr) in [(args.vcf_out, "VCF output")]:
        if not out_file:
            parser.error("%s output file argument missing." % descr)
            sys.exit(1)
        if os.path.exists(out_file) and out_file != "-":
            LOG.fatal("Cowardly refusing to overwrite existing"
                      " output file '%s'.\n" % out_file)
            sys.exit(1)

    add_plp_to_vcf(args.vcf_in, args.vcf_out, args.bam)


    
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")

    # FIXME add filter info?

    # FIXME: add format
    ##FORMAT=<ID=TD,Number=1,Type=Integer,Description="Total depth">
    ##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Number of reference bases">
    ##FORMAT=<ID=NA,Number=1,Type=Integer,Description="Number of alternate bases">
    ##FORMAT=<ID=OR,Number=1,Type=Integer,Description="Number of orphan reads supporting reference bases">
    ##FORMAT=<ID=OA,Number=1,Type=Integer,Description="Number of orphan reads supporting alternate bases">
    ##FORMAT=<ID=BR,Number=.,Type=Integer,Description="Minimum, median and maximum base-qualities for reference bases">
    ##FORMAT=<ID=BA,Number=.,Type=Integer,Description="Minimum, median and maximum base-qualities for alternate bases">
    ##FORMAT=<ID=MR,Number=.,Type=Integer,Description="Minimum, median and maximum mapping-qualities for reference bases">
    ##FORMAT=<ID=MA,Number=.,Type=Integer,Description="Minimum, median and maximum mapping-qualities for alternate bases">
    #
    # mutect format:
    ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=BQ,Number=A,Type=Float,Description="Average base quality for reads supporting alleles">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=FA,Number=A,Type=Float,Description="Allele fraction of the alternate allele with regard to reference">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
    ##FORMAT=<ID=SS,Number=1,Type=Integer,Description="Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown">
    #
    # format values
    # GT:AD:BQ:DP:FA  0/1:6,2:35:8:0.250      0:7,0:.:6:0.00
    # FORMAT  synthetic.challenge.set4.tumour synthetic.challenge.set4.normal
