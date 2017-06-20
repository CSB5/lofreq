#!/usr/bin/env python
"""Commonly used multiple correction routines

Original source: multiple_testing.py from goatools (see below).
https://github.com/tanghaibao/goatools
f75455067a7f7aad66f5b229ab514977b70c34d9

AW:
- Modified to get rid of numpy dependence.
- Added n argument (for input of clipped pvalues)

Original Authors:
- Haibao Tang (tanghaibao),
- Brent Pedersen (brentp),
- Aurelien Naldi (aurelien-naldi)
Email: tanghaibao@gmail.com
License: BSD
"""

__author__ = "Haibao Tang, Brent Pedersen, Aurelien Naldi"
__email__ = "tanghaibao@gmail.com"
#__copyright__ = ""
__license__ = "BSD"

from itertools import groupby


class AbstractCorrection(object):
    
    def __init__(self, pvals, a=.05, n=None):
        self.pvals = self.corrected_pvals = list(pvals)

        # number of multiple tests
        if n:
            assert n>len(pvals)
            self.n = n 
        else:
            self.n = len(self.pvals)
        # type-1 error cutoff for each test   
        self.a = a                  

        self.set_correction()

    def set_correction(self):
        # the purpose of multiple correction is to lower the alpha
        # instead of the canonical value (like .05)
        pass


    
class Bonferroni(AbstractCorrection):
    """http://en.wikipedia.org/wiki/Bonferroni_correction
    >>> ["%.4f" % v for v in Bonferroni([0.01, 0.01, 0.03, 0.05, 0.005], a=0.05).corrected_pvals]
    ['0.0500', '0.0500', '0.1500', '0.2500', '0.0250']
    """
    
    def set_correction(self):
        self.corrected_pvals = [pv * self.n
                                for pv in self.corrected_pvals]


        
class Sidak(AbstractCorrection):
    """http://en.wikipedia.org/wiki/Bonferroni_correction
    >>> ["%.8f" % v for v in Sidak([0.01, 0.01, 0.03, 0.05, 0.005], a=0.05).corrected_pvals]
    ['0.04898974', '0.04898974', '0.14696923', '0.24494871', '0.02449487']
    """
    def set_correction(self):
        if self.n != 0:
            correction = self.a * 1. / (1 - (1 - self.a) ** (1. / self.n))
        else:
            correction = 1
        self.corrected_pvals = [pv * correction
                                for pv in self.corrected_pvals]


        
class HolmBonferroni(AbstractCorrection):
    """http://en.wikipedia.org/wiki/Holm-Bonferroni_method
    given a list of pvals, perform the Holm-Bonferroni correction
    and return the indexes from original list that are significant.
    (cant use p-value as that may be repeated.)
    >>> ["%.4f" % v for v in HolmBonferroni([0.01, 0.01, 0.03, 0.05, 0.005], a=0.05).corrected_pvals]
    ['0.0400', '0.0400', '0.0600', '0.0500', '0.0250']
    """

    def set_correction(self):
        if len(self.pvals):
            for (i, c) in self.generate_significant():
                self.corrected_pvals[i] *= c
        
    def generate_significant(self):
        pvals = self.pvals
        pvals_idxs = zip(pvals, range(len(pvals)))
        pvals_idxs = sorted(pvals_idxs)

        #lp = len(self.pvals)
        lp = self.n

        for pval, idxs in groupby(pvals_idxs, lambda x: x[0]):
            idxs = list(idxs)
            for p, i in idxs:
                if p * 1. / lp < self.a:
                    yield (i, lp)
            lp -= len(idxs) 

# also in the original file, but removed here:
#class FDR
#def calc_qval

if __name__ == '__main__':
    import doctest
    doctest.testmod()        
