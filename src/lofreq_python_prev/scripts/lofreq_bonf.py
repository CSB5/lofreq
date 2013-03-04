#!/usr/bin/env python
"""Determine Bonferoni factor for lofreq_snpcaller.py

Bonferroni factor is either determined from:
- A bed file listing regions to analyse
or
- the header of a BAM file
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
from optparse import OptionParser, SUPPRESS_HELP

#--- third-party imports
#
# /

#--- project specific imports
#
from lofreq import sam

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011, 2012 Genome Institute of Singapore"
__license__ = "GPL2"


#global logger
# http://docs.python.org/library/logging.html
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

    parser.add_option("", "--bam",
                      dest="fbam", # type="string|int|float"
                      help="BAM input.")
    parser.add_option("", "--exclude",
                      dest="fexclude", # type="string|int|float"
                      help=SUPPRESS_HELP)
    #help="Optional/Deprecated: Exclude positions listed in this file"
    #                  " format is: start end [comment ...]"
    #                  " , with zero-based, half-open coordinates"
    #                  " (clashes with --bed)")
    parser.add_option("", "--chrom",
                      dest="chrom", # type="string|int|float"
                      help=SUPPRESS_HELP)
    #help="Optional: Chromsome to use from BAM"
    #                  " (needed for --exclude)")
    parser.add_option("", "--bed",
                      dest="fbed", # type="string|int|float"
                      help="Optional: Bed file listing positions used"
                      " in later pileup")
    #                    (clashes with --exclude and bam)")
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
        sam.LOG.setLevel(logging.INFO)

    if opts.debug:
        LOG.setLevel(logging.DEBUG)
        sam.LOG.setLevel(logging.DEBUG)

    if not opts.fbam and not opts.fbed:
        parser.error("BAM input file and bed input file argument missing.")
        sys.exit(1)
    if opts.fbam and not os.path.exists(opts.fbam):
        LOG.fatal("%s does not exist'.\n" % (opts.fbam))
        sys.exit(1)


    if opts.fexclude and not opts.chrom:
        parser.error("Chromosome argument missing." 
                     " Needed when using exclude file")
        sys.exit(1)

    if opts.fbed and opts.fexclude:
        parser.error("Can only use either exclude file or"
                     " bed file ('include')")
        sys.exit(1)

    if opts.fbed and opts.chrom:
        parser.error("chrom argument should only be used "
                     " with exclude file, not with bed file ('include')")
        sys.exit(1)
        
    
    bonf_factor = sam.auto_bonf_factor(
        opts.fbam, opts.fbed, opts.fexclude, opts.chrom)
    print bonf_factor


if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
