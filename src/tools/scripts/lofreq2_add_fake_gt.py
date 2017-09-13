#!/usr/bin/env python
"""Complement VCF with unknown genotype for sample.
"""


# --- standard library imports
#
import sys
import os
import argparse
import logging
import csv
import gzip

#--- third-party imports
#
#/

#--- project specific imports
#
# /


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2016 Genome Institute of Singapore"
__license__ = "The MIT License"


# global logger
#
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


FORMAT_HEADER = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'


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
    parser.add_argument("-s", "--samples",
                        required=True, nargs='+',
                        help="Sample name/s")
    return parser


def add_fake_gt(vcf_in, vcf_out, sample_names):
    """Add fake genotype to header and variants"""

    assert len(set(sample_names)) == len(sample_names), ("Duplicate sample names found")
    
    # set up vcf_reader
    #
    if vcf_in == '-':
        fh_in = sys.stdin
    else:
        assert os.path.exists(vcf_in)
        if vcf_in[-3:] == ".gz":
            fh_in = gzip.open(vcf_in, 'rt')
        else:
            fh_in = open(vcf_in, 'rt')
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

    has_our_format_in_header = False
    for row in vcf_reader:
        # modify row if needed and finally print

        if row[0].startswith('##'):
            if row[0].startswith('##FORMAT'):
                if row[0] == FORMAT_HEADER:
                    has_our_format_in_header = True
                else:
                    LOG.fatal("Incompatible, pre-existing format definition found. Exiting")
                    raise ValueError(row)
            # don't touch header

        elif row[0].startswith('#CHROM'):
            # insert genotype format line
            if not has_our_format_in_header:
                extrarow = ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
                vcf_writer.writerow(extrarow)

            if not "FORMAT" in row:
                row.append("FORMAT")
            for name in sample_names:
                row.append(name)
        else:
            # add format and sample columns
            assert len(row) >= 8, (
                "variant incomplete or FORMAT column already exists")
            # Add GT column if not present
            if len(row) > 8:
                assert row[8] == 'GT'
            else:
                row.append("GT")
            # Add fake GT
            for _ in sample_names:
                row.append(".")

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
            LOG.fatal("file '%s' does not exist.\n", in_file)
            sys.exit(1)

    for (out_file, descr) in [(args.vcf_out, "VCF output")]:
        if not out_file:
            parser.error("%s output file argument missing." % descr)
            sys.exit(1)
        if os.path.exists(out_file) and out_file != "-":
            LOG.fatal("Cowardly refusing to overwrite existing"
                      " output file '%s'.\n", out_file)
            sys.exit(1)

    add_fake_gt(args.vcf_in, args.vcf_out, args.samples)



if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
