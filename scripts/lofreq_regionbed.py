#!/usr/bin/env python
"""Creates a bed file listing all chromsomes and there lengths which
can be used as an input to lofreq_snpcaller.py -l (i.e. samtools
mpileup -l)
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

MYNAME = os.path.basename(sys.argv[0])




def chrom_and_len_from_bam(bam):
    """Return a list of chromosome/reference sequences listed in BAM
    file and their length
    """

    ret = []
    header = sam.sam_header(bam)
    if header == False:
        LOG.critical("parsing samtools header from %s failed" % bam)
        raise ValueError
    sq_list = sam.sq_list_from_header(header)
    for sq in sq_list:
        sq_len = sam.len_for_sq(header, sq)
        ret.append((sq, sq_len))
    return ret

    
def bam_header_to_bed_region(bam, fh_bed):
    """Parses chromsome and length from BAM file and convert to bed
    file
    """ 
    
    for (sq, length) in chrom_and_len_from_bam(bam):
        fh_bed.write("%s\t%d\t%d\n" % (sq, 0, length))
    
   
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
    parser.add_option("-i", "--bam",
                      dest="bam", # type="string|int|float",
                      help="BAM input file")
    parser.add_option("-o", "--bed",
                      dest="bed", # type="string|int|float",
                      default="-",
                      help="bed output file ('-' for stdout, which is default)")
    
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


    # checks for mandatory files
    for (filename, descr, direction) in [(opts.bam, "BAM file", 'in'),
                                         (opts.bed, "bed file", 'out')]:
        if not filename:
            parser.error("%s argument missing." % descr)
            sys.exit(1)
        if filename == '-':
            continue
        if direction == 'in' and not os.path.exists(filename):
            sys.stderr.write(
                "file '%s' does not exist.\n" % filename)
            sys.exit(1)
        if direction == 'out' and os.path.exists(filename):
            sys.stderr.write(
                "Refusing to overwrite existing file '%s'.\n" % filename)
            sys.exit(1)
            
    
    if opts.bed == '-':
        fh = sys.stdout
    else:
        fh = open(opts.bed, 'w')

    bam_header_to_bed_region(opts.bam, fh)

    if fh != sys.stdout:
        fh.close()


        
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")


    
