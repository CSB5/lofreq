#!/usr/bin/env python
"""Generic utils
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
from math import log10
import sys

#--- third-party imports
#
# /


#--- project specific imports
#
# nothing should go here by definition

    
__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011, 2012 Genome Institute of Singapore"
__license__ = "GPL2"




#def mean_and_stdv(x):
#    """
#    Calculate mean and standard deviation of data x[]:
#    mean = {\sum_i x_i \over n}
#    std = sqrt(\sum_i (x_i - mean)^2 \over n-1)
# 
#    Based on
#    http://www.physics.rutgers.edu/~masud/computing/WPark_recipes_in_python.html
#    """
# 
#    num = len(x)
#    assert num != 0
#    if num == 1:
#        return (x[0], 0.0)
#        
#    mean = sum(x)/float(num)
#    std = sum([(a-mean)**2 for a in x])
#    std = sqrt(std / float(num-1))
# 
#    return mean, std



def prob_to_phredqual(prob):
    """
    Turns an error probability into a phred value
    
    >>> prob_to_phredqual(0.01)
    20
    
    """

    assert prob >= 0.0, (
        "Probability can't be smaller than 0 but got %f" % prob)
    try:
        return int(round(-10.0 * log10(prob)))
    except ValueError:
        # prob is zero
        return sys.maxint


def phredqual_to_prob(phredqual):
    """
    Turns a phred quality into an error probability

    >>> '%.2f' % phredqual_to_prob(20)
    '0.01'

    """

    assert phredqual >= 0, ("Phred-quality must be >= 0, but is %s" % phredqual)
    return 10**(-phredqual/10.0)


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
        # http://genome.ucsc.edu/FAQ/FAQformat.html#format1
        # 4: name, score, strand...
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
    
    Orientation agnostic!
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
        # ignore  rest of line        
        assert start < end, (
            "Invalid position found in %s" % fexclude)
        excl_pos.extend(range(start, end))
    fhandle.close()
    # remove duplicates (on F and R)
    return set(excl_pos)



    
if __name__ == '__main__':
    import doctest
    doctest.testmod()        
            
