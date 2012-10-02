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


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')



def read_bed_coords(fbed):
    """Reads coordinates from bed file and returns them in a dict with
    chromsomes as key. Ranges are a tuple of start and end pos
    (0-based; half-open just like bed and python ranges)
    """

    bed_coords = dict()
    
    fh = open(fbed, 'r')
    for line in fh:
        if line.startswith('#'):
            continue
        if len(line.strip()) == 0:
            continue
        try:
            (chrom, start, end) = line.split("\t")[0:3]
        except IndexError:
            sys.stderr.write(
                "FATAL: Failed to parse bed line: %s\n" % line)
            raise ValueError
        start = int(start)
        end = int(end)
        assert end >= start, (
            "Start value (%d) not lower equal end value (%d)."
            " Parsed from file %s" % (
                start, end, fbed))
        if not bed_coords.has_key(chrom):
            bed_coords[chrom] = []
        bed_coords[chrom].append((start, end))
    fh.close()
    return bed_coords

    

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
        assert start < end, (
            "Invalid position found in %s" % fexclude)        
        excl_pos.extend(range(start, end))
        LOG.debug("got excl pos: %d-%d" % (start, end))
    fhandle.close()


    # remove duplicate ones (F and R)
    excl_pos = set(excl_pos)

    return excl_pos

    

def sum_chrom_len(fbam, chrom_list=None):
    """
    Return length of all chromsomes. If chrom_list is not empty then
    only consider those. Length is extracted from BAM file (fbam)
    header
    """
    
    sum_sq_len = 0
    bam_header = pileup.header(fbam)
    if bam_header == False:
        LOG.critical("samtools header parsing failed test")
        raise ValueError
    sq = pileup.sq_from_header(bam_header)

    # use all if not set
    if not chrom_list:
        chrom_list = sq
        
    for chrom in chrom_list:
        assert chrom in sq, (
        "Couldn't find chromosome '%s' in BAM file '%s'" % (
            chrom, fbam))
        
        sq_len = pileup.len_for_sq(bam_header, chrom)
        sum_sq_len += sq_len
        
    return sum_sq_len

    



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
                      " , with zero-based, half-open coordinates"
                      " (clashes with --bed)")
    parser.add_option("-c", "--chrom",
                      dest="chrom", # type="string|int|float"
                      help="Optional: Chromsome to use from BAM"
                      " (needed for --exclude)")
    parser.add_option("-b", "--bed",
                      dest="fbed", # type="string|int|float"
                      help="Optional: Bed file listing positions used"
                      " in later pileup (clashes with --exclude)")
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
        

    
    # exclude positions
    #
    if opts.fexclude:
        excl_pos = []
        excl_pos = read_exclude_pos_file(opts.fexclude)
        LOG.info("Parsed %d positions from %s" % (
            len(excl_pos), opts.fexclude))
        
        sum_sq_len = sum_chrom_len(opts.fbam, [opts.chrom])
        sum_sq_len -= len(excl_pos)

    elif opts.fbed:

        bed_coords = read_bed_coords(opts.fbed)
        sum_sq_len = 0
        for (chrom, ranges) in bed_coords.iteritems():
            for r in ranges:
                LOG.debug("bed coord range for %s: %d-%d" % (
                    chrom, r[0], r[1]))
                diff = r[1]-r[0]
                assert diff > 0
                sum_sq_len += diff
    else:
        # look at all
        sum_sq_len = sum_chrom_len(opts.fbam)

    bonf_factor = sum_sq_len * 3
    print bonf_factor


if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
