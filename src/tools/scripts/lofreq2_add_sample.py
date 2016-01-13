#!/usr/bin/env python
"""Complement VCF with simple pileup info from BAM files
"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "The MIT License"


# --- standard library imports
#
import sys
import os
import argparse
import logging
from collections import OrderedDict, namedtuple
import csv
import gzip

#--- third-party imports
#
import pysam
assert [int(x) for x in pysam.__version__.split('.')] >= [0, 8, 2]

#--- project specific imports
#
# /


# global logger
#
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


Variant = namedtuple('Variant',
                     ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
# all fields except POS (int) are strings and values are preserved as-is

Format = namedtuple('Format',
                    ['id', 'num', 'type', 'descr'])

def median(data):
    """compute median of provided list"""

    if not len(data):
        return None
    # http://stackoverflow.com/questions/10482339/how-to-find-median/10482422#10482422 answer by user3100512
    return sorted(data)[len(data)//2]


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
                        dest="bams", nargs="*",
                        required=True,
                        help="BAM files, e.g. normal and tumor bam")
    return parser


def fmt_to_line(fmt):
    """convert format class to vcf line"""
    return "##FORMAT=<ID=%s,Number=%s,Type=%s,Description=\"%s\">" % (
        fmt.id, fmt.num, fmt.type, fmt.descr)


def gen_formats():
    """Must be in sync with gen_plp_data
    """

    formats = OrderedDict()
    for (fid, num_str, type_str, descr) in [
            ('DP', '1', 'Integer', 'Read depth at this position for this sample'),# standard
            ('NR', '1', 'Integer', 'Number of reference bases'),
            ('NA', '1', 'Integer', 'Number of alternate bases'),
            ('OR', '1', 'Integer', 'Number of orphan reads supporting reference bases'),
            ('OA', '1', 'Integer', 'Number of orphan reads supporting alternate bases'),
            ('BR', '3', 'Integer', 'Minimum, median and maximum base-qualities for reference bases'),
            ('BA', '3', 'Integer', 'Minimum, median and maximum base-qualities for alternate bases'),
            ('MR', '3', 'Integer', 'Minimum, median and maximum mapping-qualities for reference bases'),
            ('MA', '3', 'Integer', 'Minimum, median and maximum mapping-qualities for alternate bases')]:
        formats[fid] = Format(id=fid, num=num_str, type=type_str, descr=descr)

    return formats


def gen_plp_data(sam_fh, var):
    """generate data must be in sync with gen_formats()
    """

    for plp_col in sam_fh.pileup(var.CHROM, var.POS-1, var.POS):
        # pileup() extracts all reads overlapping that region.
        # only look at the one of interest
        if plp_col.pos != var.POS-1:
            continue

        cov = plp_col.n
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

            base = aln_read.seq[plp_read.query_position]
            mq = aln_read.mapq
            bq = ord(aln_read.qual[plp_read.query_position])-33

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

    sample_data = OrderedDict()
    for (fmt_key, val) in [
            ('DP', "%d" % cov),
            ('NR', "%d" % len(bqs['ref'])),
            ('NA', "%d" % len(bqs['alt'])),
            ('OR', "%d" % num_orphans['ref']),
            ('OA', "%d" % num_orphans['alt']),
            ('BR', "%d,%d,%d" % (min_bqs['ref'], median_bqs['ref'], max_bqs['ref'])),
            ('BA', "%d,%d,%d" % (min_bqs['alt'], median_bqs['alt'], max_bqs['alt'])),
            ('MR', "%d,%d,%d" % (min_mqs['ref'], median_mqs['ref'], max_mqs['ref'])),
            ('MA', "%d,%d,%d" % (min_mqs['alt'], median_mqs['alt'], max_mqs['alt']))]:
        sample_data[fmt_key] = val

    return sample_data


def add_plp_to_vcf(vcf_in, vcf_out, bam_files):
    """process each var in vcf_in and add plp info from sam_fh,
    writing to vcf_out. is no way to edit/add format fields in current
    versions of pyvcf (as of 2014-06-30). see discussion here
    https://github.com/jamescasbon/PyVCF/issues/82 for patches and
    workarounds. chose to use csv module instead for simplicity
    """

    assert all([os.path.exists(b) for b in bam_files])

    # set up vcf_reader
    #
    if vcf_in == '-':
        fh_in = sys.stdin
    else:
        assert os.path.exists(vcf_in)
        if vcf_in[-3:] == ".gz":
            fh_in = gzip.open(vcf_in, 'rb')
        else:
            fh_in = open(vcf_in, 'rb')
    vcf_reader = csv.reader(fh_in, delimiter='\t')


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
    vcf_writer = csv.writer(fh_out, delimiter='\t',
                            quotechar='', quoting=csv.QUOTE_NONE,
                            lineterminator=os.linesep)

    formats = gen_formats()

    for row in vcf_reader:

        if row[0].startswith('#'):
            if row[0] == "#CHROM":
                assert len(row) == 8, (
                    "variant incomplete or FORMAT column already exists")

                # before writing header, add our format description.
                for fmt in formats.values():
                    vcf_writer.writerow([fmt_to_line(fmt)])

                row.append("FORMAT")

                for bam in bam_files:
                    row.append(os.path.basename(bam))

            vcf_writer.writerow(row)

        else:
            assert len(row) == 8, (
                "variant incomplete or FORMAT column already exists")
            var = Variant._make([row[0], int(row[1]), row[2], row[3],
                                 row[4], row[5], row[6], row[7]])

            # no support for indels
            if 'INDEL' in var.INFO.split(';') or len(var.REF) > 1 or len(var.ALT) > 1:
                LOG.warn("Skipping unsupported variant) %s:%d:%s" % (
                    var.CHROM, var.POS, var.REF))
                continue

            row.append(':'.join(formats.keys()))
            for bam in bam_files:
                assert os.path.exists(bam)
                sam_fh = pysam.AlignmentFile(bam)

                sample_data = gen_plp_data(sam_fh, var)
                assert sample_data.keys() == formats.keys(), (
                    "sample keys (%s) != format keys (%s)" % (sample_data.keys(), formats.keys()))
                row.append(':'.join(sample_data.values()))
            vcf_writer.writerow(row)


    if fh_in != sys.stdin:
        fh_in.close()
    if fh_out != sys.stdout:
        fh_out.close()


def main():
    """main function
    """

    parser = cmdline_parser()
    args = parser.parse_args()

    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)

    for (in_file, descr) in [#(args.bam, "BAM"),
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

    add_plp_to_vcf(args.vcf_in, args.vcf_out, args.bams)



if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
