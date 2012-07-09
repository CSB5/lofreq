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
# /

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
    """

    assert phredqual >= 0, ("Phred-quality must be >= 0, but is %s" % phredqual)
    return 10**(-phredqual/10.0)



#def count_bases(bases, allowed_bases = "ACGT"):
#    """
#    Counts bases and return counts as dict with base as key and count
#    as value. Also returns consensus base. Only counts allowed_bases.
#    No case conversion will be done!
#
#    Alternative collections.counter is only available in Python 2.7
#
#    doctest:
#    >>> count_bases("AAAA")[1]
#    'A'
#    >>> count_bases("AAAC")[1]
#    'A'
#    >>> count_bases("AACC")[1]
#    'N'
#    >>> count_bases("QQQQ")[1]
#    '-'
#    """
#
#    basecounts = dict()
#    for b in allowed_bases:
#        basecounts[b] = bases.count(b)
#
#    # sorted (ascending) list of key/value tuples. take last one
#    # [-1] and only its key [0]
#    #
#    if sum(basecounts.values()) == 0:
#        # empty? return gap
#        cons_base = "-"
#    else:
#        # ascending order
#        sorted_counts = sorted(basecounts.items(),
#                               key=lambda x: x[1])
#        if sorted_counts[-1][1] == sorted_counts[-2][1]:
#            # tie? return N
#            cons_base = 'N'
#        else:
#            # return "true" consensus
#            cons_base = sorted_counts[-1][0]
#            
#    return (basecounts, cons_base)
#

if __name__ == '__main__':
    import doctest
    doctest.testmod()        
            
