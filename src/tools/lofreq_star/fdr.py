"""FDR routines
"""

#--- standard library imports
#

import sys
try:
    from itertools import izip
except ImportError:
    def izip(a, b):
        return zip(a, b)


#--- third-party imports
#
# /

#--- project specific imports
#
# /


__author__ = "Grace Hui Ting Yeo"
__email__ = "yeohtg@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "The MIT License"


if sys.version_info >= (3, 0):
    def xrange(*args, **kwargs):
        return iter(range(*args, **kwargs))


def fdr(pvals, a=0.05, n=None):
    """
    Implementation of the Benjamini-Hochberg procedure.
    Takes a list of p-values and returns a list of the indices of those p-values that pass.
    Does not adjust p-values.
    See http://sas-and-r.blogspot.sg/2012/05/example-931-exploring-multiple-testing.html
    for pseudocode.

    Test data from : http://udel.edu/~mcdonald/statmultcomp.html
    >>> import random
    >>> pvals = [0.6, 0.07, 0.49, 0.2, 0.48, 0.74, 0.68, 0.01, 0.97, 0.38, 0.032, 0.07]
    >>> random.shuffle(pvals)
    >>> sorted([pvals[i] for i in fdr(pvals, a=0.20)])
    [0.01, 0.032]
    >>> fdr([])
    []
    >>> fdr([1])
    []
    """

    if n != None:
        assert n >= len(pvals)
    else:
<<<<<<< HEAD
        n=len(pvals)

    sorted_pvals_indices = sorted(range(len(pvals)), key=lambda k:pvals[k])
    t = next((rank for rank, spi in zip(range(len(pvals), 0, -1),
                                         reversed(sorted_pvals_indices))
=======
        n = len(pvals)

    sorted_pvals_indices = sorted(xrange(len(pvals)), key=lambda k: pvals[k])
    t = next((rank for rank, spi in izip(xrange(len(pvals), 0, -1),
                                         reversed(sorted_pvals_indices))
>>>>>>> py3
              if pvals[spi] < rank*a/n), None)
    if t:
        return sorted_pvals_indices[:t]
    return []
