#!/usr/bin/env python
"""TMP: Get ignore positions (those which would be called cons-vars) from pileup, and list them (don't do anything else)
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
# /


#--- project specific imports
#
from lofreq import pileup

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011, 2012 Genome Institute of Singapore"
__license__ = "GPL2"


# global logger
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
    parser.add_option("-i", "--in",
                      dest="fpileup", # type="string|int|float"
                      help="Pileup input. Will read from stdin if '-' or not set at all."
                      " Tip: Use '-d 100000' to prevent sample depth filtering by samtools."
                      " Also consider using -B/-E to influence BAQ computation")
    parser.add_option("-o", "--out",
                      dest="fbed",
                      help="Output bed file with ignore positions or '-' for stdout")

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


    pileup_fh = sys.stdin
    if not opts.fpileup:
        LOG.warn("No pileup file/stream specified. Will try stdin\n")
    else:
        if opts.fpileup != '-':
            if not os.path.exists(opts.fpileup):
                LOG.fatal("file '%s' does not exist.\n" % (opts.fpileup))
                sys.exit(1)
            else:
                pileup_fh = open(opts.fpileup, 'r')

    bed_fh = sys.stdout
    if not opts.fbed:
        LOG.warn("No output file/stream specified. Will use stdout\n")
    else:
        if opts.fbed != '-' and os.path.exists(opts.fbed):
            LOG.fatal("Cowardly refusing to overwrite existing output file %s\n" % opts.fbed)
            sys.exit(1)
        else:
            bed_fh = open(opts.fbed, 'w')

    num_lines = 0
    for line in pileup_fh:

        # note: not all columns will be present in pileup
        num_lines += 1

        pcol = pileup.PileupColumn(line)            
        #if (pcol.coord+1) % 100000 == 0:
        #    LOG.info("Still alive...calling SNPs in column %d" % (pcol.coord+1))


        if pcol.cons_base not in 'ACGT':
            LOG.info("Skipping col %d because of ambigious consensus %s" % (
                    pcol.coord+1, pcol.cons_base))
            continue

        # Even though NQ is quality agnostics, ignore the ones
        # marked with Illumina's Read Segment Indicator Q2, since
        # those are not supposed to be used according to Illumina
        base_counts = pcol.get_all_base_counts(
            min_qual=3, keep_strand_info=False)
        del base_counts['N']
        if sum(base_counts.values())==0:
            LOG.info("Zero coverage in col %d" % (pcol.coord+1))
            continue
        
        # consensus snp, i.e. we determine consensus different to
        # initial pileup ref. always call against cons_base. if
        # ref_base is called then drop silently later.
        if pcol.cons_base != pcol.ref_base:
            LOG.debug("Listing position %s:%d because ref %s != cons %s (base counts = %s)" % (
                pcol.chrom, pcol.coord+1, 
                pcol.ref_base, pcol.cons_base,
                base_counts))
            # bed output (zero-based, exclusive)
            bed_fh.write("%s\t%d\t%d\n" % (pcol.chrom, pcol.coord, pcol.coord+1))
                

    if num_lines == 0:
        LOG.fatal("Pileup was empty. Will exit now.")
        sys.exit(1)

    if pileup_fh != sys.stdin:
        pileup_fh.close()
    if bed_fh != sys.stdout:
        bed_fh.close()
        


if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
