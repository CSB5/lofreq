#!/usr/bin/env python
"""Create summary of base calls per pileup colum given by stdin filtered
at different quality levels
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
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser

#--- third-party imports
#
#/

#--- project specific imports
#
from lofreq import sam
from lofreq import conf

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011, 2012 Genome Institute of Singapore"
__license__ = "GPL2"


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

BASES = ['A', 'C', 'G', 'T', 'N']


def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog [Options]\n" \
            + "\n" + __doc__
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      dest="verbose",
                      action="store_true",
                      help="be verbose")
    parser.add_option("", "--debug",
                      dest="debug",
                      action="store_true",
                      help="enable debugging")
    parser.add_option("", "--orig-samtools",
                      dest="orig_samtools",
                      action="store_true",
                      help="Pileup comes from samtools instead of lofreq_samtools")
    parser.add_option("-i", "--input",
                      dest="pileup",
                      default="-",
                      help="Pileup (- for stdin). Use of lofreq_samtools is recommended."
                      " If you use samtools you need --orig-samtools as well")
    parser.add_option("-q", "--qual",
                      dest="qual",
                      default=conf.DEFAULT_IGN_BASES_BELOW_Q,
                      type="int",
                      help="Quality threshold (default: %d)" % (
                          conf.DEFAULT_IGN_BASES_BELOW_Q))
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

    if not opts.qual and opts.qual != 0:
        parser.error("Missing quality threshold argument.")
        sys.exit(1)

    if opts.verbose:
        LOG.setLevel(logging.INFO)
        sam.LOG.setLevel(logging.INFO)

    if opts.debug:
        LOG.setLevel(logging.DEBUG)
        sam.LOG.setLevel(logging.DEBUG)


    if opts.pileup == "-":
        LOG.info("Reading from stdin")
        fh = sys.stdin
    else:
        LOG.info("Reading from %s" % opts.pileup)
        fh = open(opts.pileup, 'r')

    qual = int(opts.qual)


    for line in fh:
        if len(line.strip())==0:
            continue

        if opts.orig_samtools:
            pcol = sam.PileupColumn(line)
        else:
            pcol = sam.LoFreqPileupColumn(line)

        base_counts = pcol.get_all_base_counts(qual)
        print "coord\tbases\ttotal\tfw\trv"
        for base in sorted(base_counts.keys()):
            print "%d\t%s\t%d\t%d\t%d" % (pcol.coord+1, base,
                                          sum(base_counts[base]),
                                          base_counts[base][0],
                                          base_counts[base][1])

        if opts.orig_samtools:
            print "coord\tins\tdel\trstart\trend"
            print "%d\t%d\t%d\t%d\t%d" % (
                pcol.coord+1,
                pcol.num_ins_events,
                pcol.num_del_events,
                pcol.num_read_starts,
                pcol.num_read_ends)
        else:
            print "coord\tins\tins-len\tdel\tdel-len\trstart\trend"
            print "%d\t%d\t%.1f\t%d\t%.1f\t%d\t%d" % (
                pcol.coord+1,
                pcol.num_ins_events,
                pcol.avg_ins_len,
                pcol.num_del_events,
                pcol.avg_del_len,
                pcol.num_read_starts,
                pcol.num_read_ends)

        
    if fh != sys.stdin:
        fh.close()
    



if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
