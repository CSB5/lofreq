#!/usr/bin/env python
"""Determine Bonferoni factor for lofreq_snpcaller.py
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
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301 USA.



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
__version__ = "0.0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = ""
__license__ = ""
__credits__ = [""]
__status__ = ""


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')



def read_exclude_pos_file(fexclude):
    """Parse file containing ranges of positions to exclude and return
    positions as list.
    """

    excl_pos = []
    fhandle = open(fexclude, 'r')
    for line in fhandle:
        if line.startswith('#'):
            continue
        if len(line.strip()) == 0:
            continue

        start = int(line.split()[0])
        end = int(line.split()[1])
        # ignore the rest

        assert start < end, ("Invalid position found in %s" % fexclude)

        excl_pos.extend(range(start, end))
    fhandle.close()

    return excl_pos

    

def determine_bonf_factor(fbam, chrom, num_excl_pos=0):
    """Determine Bonferroni correction factor from sequence length in
    given BAM files
    """

    bam_header = pileup.header(fbam)
    if bam_header == False:
        LOG.critical("samtools header parsing failed test")
        raise ValueError
    sq = pileup.sq_from_header(bam_header)
    assert chrom in sq, (
        "Couldn't find chromosome '%s' in BAM file '%s'" % (
            chrom, fbam))
    sq_len = pileup.len_for_sq(bam_header, chrom)
    # one for each nonconsensus nucleotide:
    bonf_factor = (sq_len-num_excl_pos) * 3

    return bonf_factor




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

    parser.add_option("-i", "--input",
                      dest="fbam", # type="string|int|float"
                      help="BAM input.")
    parser.add_option("-e", "--exclude",
                      dest="fexclude", # type="string|int|float"
                      help="Optional: Exclude positions listed in this file"
                      " format is: start end [comment ...]"
                      " , with zero-based, half-open coordinates")
    parser.add_option("-c", "--chrom",
                      dest="chrom", # type="string|int|float"
                      help="Chromsome to use from BAM.")
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
        pileup.LOG.setLevel(logging.INFO)

    if opts.debug:
        LOG.setLevel(logging.DEBUG)
        pileup.LOG.setLevel(logging.DEBUG)

    if not opts.fbam:
        parser.error("BAM input file argument missing.")
        sys.exit(1)
    if not os.path.exists(opts.fbam):
        LOG.fatal("%s does not exist'.\n" % (opts.fbam))
        sys.exit(1)

    if not opts.chrom:
        parser.error("Chromosome argument missing.")
        sys.exit(1)


    # exclude positions
    #
    excl_pos = []
    if opts.fexclude:
        excl_pos = read_exclude_pos_file(opts.fexclude)
        LOG.info("Parsed %d positions from %s" % (
            len(excl_pos), opts.fexclude))

    print determine_bonf_factor(opts.fbam, opts.chrom, len(excl_pos))




if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
