#!/usr/bin/env python
"""Generic utils
"""



#--- standard library imports
#
from __future__ import division
from math import sqrt, log10


#--- third-party imports
#
# /


#--- project specific imports
#
# /


__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = ""
__license__ = ""
__credits__ = [""]
__status__ = ""




def mean_and_stdv(x):
    """
    Calculate mean and standard deviation of data x[]:
    mean = {\sum_i x_i \over n}
    std = sqrt(\sum_i (x_i - mean)^2 \over n-1)

    Based on http://www.physics.rutgers.edu/~masud/computing/WPark_recipes_in_python.html
    """

    num = len(x)
    assert num != 0
    if num == 1:
        return (x[0], 0.0)
        
    mean = sum(x)/float(num)
    std = sum([(a-mean)**2 for a in x])
    std = sqrt(std / float(num-1))

    return mean, std



def prob_to_phredqual(prob):
    """
    Turns an error probability into a phred value
    """

    return int(round(-10.0 * log10(prob)))



def count_bases(bases, allowed_bases = "ACGT"):
    """
    Counts bases and return counts as dict with base as key and count
    as value. Also returns consensus base. Only counts allowed_bases.
    No case conversion will be done!
    """

    basecounts = dict()
    for b in allowed_bases:
        basecounts[b] = bases.count(b)

    # sorted (ascending) list of key/value tuples. take last one
    # [-1] and only its key [0]
    if sum(basecounts.values()) != 0:
        cons_base = sorted(basecounts.items(),
                           key=lambda x: x[1])[-1][0]
    else:
        cons_base = "-"
        
    return (basecounts, cons_base)



