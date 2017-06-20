#!/usr/bin/env python
"""Generic utils for LoFreq
"""


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011 Genome Institute of Singapore"
__license__ = "The MIT License"



#--- standard library imports
#
from math import log10, log
import sys
from time import strftime
import string

MAX_INT = 2147483647
# instead of sys.maxint

#--- third-party imports
#
from scipy.stats import chi2


#--- project specific imports
#
# nothing should go here by definition




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



def now():
    return strftime("%Y-%m-%d %H:%M:%S")



def fisher_comb(pv1, pv2):
    """
    Fisher's method for combining p-values
    
    See for example
    http://en.wikipedia.org/wiki/Fisher's_method
    and
    breseq-0.18b:polymorphism_statistics.r
    """
    
    if pv1 == 0 or pv2 == 0:
        # not sure if this is correct.
        # see also http://stats.stackexchange.com/questions/58537/fishers-method-when-p-value-0
        return 0.0
    
    comb_log = -2.0 * (log(pv1) + log(pv2))
    # http://stackoverflow.com/questions/11725115/p-value-from-chi-sq-test-statistic-in-python
    comb_pv = 1.0 - chi2.cdf(comb_log, 4)    
    return comb_pv



def complement(strand, na_type='DNA'):
    """return complement of nucleic acid seqeunce

    original source http://stackoverflow.com/questions/1738633/more-pythonic-way-to-find-a-complementary-dna-strand
    Nadia Alramli

    Added DNA/RNA handling

    >>> complement("UTAGC")
    'AATCG'
    >>> complement("utagc")
    'aatcg'
    >>> complement("UTAGC", na_type="RNA")
    'AAUCG'
    >>> complement("utagc", na_type="RNA")
    'aaucg'
    """

    if na_type == 'DNA':
        if sys.version_info[0] > 2:
            tr = bytes.maketrans(b'UTAGCutagc', b'AATCGaatcg')
        else:
            tr = string.maketrans('UTAGCutagc', 'AATCGaatcg')
    elif na_type == 'RNA':
        if sys.version_info[0] > 2:
            tr = bytes.maketrans(b'UTAGCutagc', b'AAUCGaaucg')
        else:
            tr = string.maketrans('UTAGCutagc', 'AAUCGaaucg')
    else:
        raise ValueError("Unknown NA type %s" % na_type)
    return strand.translate(tr)



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
        #return sys.maxint
        return MAX_INT


    
def phredqual_to_prob(phredqual):
    """
    Turns a phred quality into an error probability

    >>> '%.2f' % phredqual_to_prob(20)
    '0.01'

    """

    assert isinstance(phredqual, int)
    #assert phredqual >= 0, ("Phred-quality must be >= 0, but is %s" % phredqual)
    # also works for phredqual=0
    return 10**(-phredqual/10.0)

    
if __name__ == '__main__':
    import doctest
    doctest.testmod()        
            
