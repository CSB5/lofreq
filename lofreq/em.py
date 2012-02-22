#!/usr/bin/env python
"""Prediction of low frequency variants by Expectation-Maximization.
"""


#--- standard library imports
#
from __future__ import division
import logging
import os

#--- third-party imports
#
# default is to use our ownbinomial extension instead of introducing a
# scipy dependency.
#
USE_SCIPY = False
if USE_SCIPY:
    from scipy import stats
    #from scipy.stats.distributions import binom


#--- project specific imports
#
# /
from lofreq import snp
from lofreq.utils import count_bases
if not USE_SCIPY:
    from lofreq_ext import binom_sf


__author__ = ["Grace Yeo", "Andreas Wilm"]
__version__ = "0.1.20111004"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = ""
__license__ = ""
__credits__ = [""]
__status__ = ""


# initial error probabilities. usually converge real fast.
# format: [from_base][snp_base]
DEFAULT_ERROR_PROBS = {
    'A': {'C':0.001, 'G':0.001, 'T':0.001},
    'C': {'A':0.001, 'G':0.001, 'T':0.001},
    'G': {'A':0.001, 'C':0.001, 'T':0.001},
    'T': {'A':0.001, 'C':0.001, 'G':0.001},
    
    # N is for the 4-parameter model
    'N': {'A':0.001, 'C':0.001, 'G':0.001, 'T':0.001}
    }


# 4 or 12 parameter model
DEFAULT_NUM_PARAM = 12

# maximum number of em iterations (usually converges much faster)
MAX_ITER = 100

# covergence if new probs don't differ by more than this threshold
DEFAULT_CONVERGENCE_EPSILON = 1e-9

# bonferroni factor
DEFAULT_BONF_FACTOR = 1.0

# significance threshold
DEFAULT_SIG_THRESH = 0.01


LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

if USE_SCIPY:
    LOG.warn("Using scipy functions instead of internal ones")


def expectation(base_counts, ref_seq,
                snp_base, snp_base_error_prob,
                bonf_factor, sig_thresh,
                from_base='N'):
    """Returns SNP calls for seeing given snp-base at each position in
    region. SNP-calls are tuples of pvalues and positions. P-Values
    are Bonferroni corrected and filtered according to significance
    threshold. Allowed SNP calls: from_base to N

    Calculated pvalue is the tail of the binomial distribution,
    which is defined by the coverage at that position, and the
    assumed error probability for the snp_base.
    """

    assert len(base_counts) == len(ref_seq)
    
    snp_calls = []
    for (pos, base_counts) in enumerate(base_counts):
        ref_base = ref_seq[pos]

        if snp_base == ref_base:
            continue    

        if from_base != 'N':
            if ref_base != from_base:
                continue

        coverage = sum(base_counts.values())
        if coverage == 0:
            continue
        
        if USE_SCIPY:
            # old (less accurate): complementary cumulative
            # distribution function
            # pvalue = (1 - bdist.cdf(base_counts[snp_base] - 1))

            # binom.sf
            #bdist = binom(coverage, snp_base_error_prob)
            # complementary cumulative distribution function
            #pvalue = bdist.sf(base_counts[snp_base]) * bonf_factor

            pvalue = stats.binom.sf(base_counts[snp_base]-1,
                coverage, snp_base_error_prob) * bonf_factor

        else:
            pvalue  = binom_sf(coverage, base_counts[snp_base]-1,
                               snp_base_error_prob) * bonf_factor
            
        if pvalue < sig_thresh:
            snp_calls.append((pos, pvalue))
        #LOG.debug("pos %d coverage=%d base_counts[snp_base]=%d pvalue %f" % (pos, coverage, base_counts[snp_base], pvalue))
    return snp_calls



def maximization(base_counts, ref_seq,
                 snp_base, snp_cols,
                 from_base='N'):
    """Recalculate the error probability for given (SNP) base by
    returning the relative frequency of the SNP base (changing from
    from_base if given). Ignore already called SNP
    positions.

    base_counts are the base counts for all alignment
    positions in ref_seq.
    
    snp_cols stores columns for which this snp_base was already
    predicted to be a SNP.
    """

    assert len(base_counts) == len(ref_seq)
    if len(snp_cols):
        assert max(snp_cols) < len(ref_seq)

    # total coverage (base count sum) in region
    cov_sum = 0
    # snp base count
    snp_count_sum = 0 
    
    # loop over all columns/positions, ignoring the ones where the snp
    # base to test is actually the consensus and the ones where we
    # know it's a snp. restrict this to a certain src base if
    # wanted.
    #
    for col_num in xrange(len(ref_seq)):
        col_base_counts = base_counts[col_num]
        ref_base = ref_seq[col_num]

        # NOTE: in Grace's original code the snp_col test was there,
        # but broken. remove or col_num in snp_cols and you get her
        # (broken) version
        #
        if snp_base == ref_base or col_num in snp_cols:
            continue
        if from_base != 'N':
            if ref_base != from_base:
                continue
        
        cov_sum += sum(col_base_counts.values())
        snp_count_sum += col_base_counts[snp_base]

    if cov_sum > 0.0:
        return snp_count_sum/float(cov_sum)
    else:
        return 0.0
    


class EmBasedSNPCaller(object):
    """
    EM based SNP caller.
    """

    def __init__(self,
                 num_param = DEFAULT_NUM_PARAM,
                 bonf_factor = DEFAULT_BONF_FACTOR,
                 sig_thresh = DEFAULT_SIG_THRESH,
                 convergence_epsilon = DEFAULT_CONVERGENCE_EPSILON):
        """
        init function
        """

        assert convergence_epsilon > 0.0 and convergence_epsilon < 1.0, (
            "Invalid convergence_epsilon")
        assert num_param in [4, 12], (
            "Unsupported number of parameters chosen")

        self.bonf_factor = bonf_factor
        self.sig_thresh = sig_thresh
        self.num_param = num_param
        self.convergence_epsilon = convergence_epsilon
        
        self.error_probs = dict()



    def set_default_error_probs(self):
        """Initializes error probabilities to their default
        """

        self.error_probs = dict()
        for src_base in DEFAULT_ERROR_PROBS.keys():
            self.error_probs[src_base] = dict()
            for dst_base in DEFAULT_ERROR_PROBS[src_base].keys():
                prob = DEFAULT_ERROR_PROBS[src_base][dst_base]
                self.error_probs[src_base][dst_base] = prob



    def load_error_probs(self, error_prob_file):
        """Loads error probabilties from a file
        """

        self.set_default_error_probs()

        fhandle = open(error_prob_file, 'r')
        for line in fhandle:
            line = line.rstrip(os.linesep)
            if len(line.strip()) == 0 or line.startswith("#"):
                continue

            line_split = line.split()
            src_base = line_split[0]
            assert src_base in 'ACGTN'
            #assert len(line_split[1:])%2 == 0
            for i in range(1, len(line_split), 2):
                dst_base = line_split[i]
                assert dst_base != src_base
                assert dst_base in 'ACGTN'
                prob = float(line_split[i+1])
                self.error_probs[src_base][dst_base] = prob
                
        fhandle.close()
        

    def save_error_probs(self, error_prob_file):
        """Saves error probabilties from a file
        """
       
        assert len(self.error_probs), (
            "Error probabilities haven't"
            " been initialized. Please call em_training() or"
            " load_error_probs() first.")

        assert not os.path.exists(error_prob_file), (
            "Cowardly refusing to overwrite existing file %s" % (
                error_prob_file))
        fhandle = open(error_prob_file, 'w')
        for src_base in sorted(self.error_probs.keys()):
            fhandle.write("%s" % src_base)
            for dst_base in sorted(self.error_probs[src_base].keys()):
                prob = self.error_probs[src_base][dst_base]
                fhandle.write(" %s %f" % (dst_base, prob))
            fhandle.write("\n")
        fhandle.close()
        

        
    def em_training(self, base_counts, ref_seq):
        """Runs EM iteratively to compute error probabilities (snps
        are returned).
     
        Error probabilities are a dict of dicts, where the first key is
        the consensus base and the second is the snp. SNPs are returned as
        list of tuples of pos, variant and pvalue.
     
        Can be run on a small subset and returned error probs can than
        be used for new and final expectation call to get all SNPs by
        calling expectation again.
        """
     
        assert len(base_counts) == len(ref_seq)

        final_snp_calls = []

        self.set_default_error_probs()

        if self.num_param == 12:
            from_bases = self.error_probs.keys()
            from_bases.remove('N')
        elif self.num_param == 4:
            from_bases = 'N'
        else:
            raise ValueError, ("Unsupported number of parameters chosen")

        for from_base in from_bases:
            for snp_base in self.error_probs[from_base].keys():
                error_prob = self.error_probs[from_base][snp_base]
     
                num_iter = 0
                while True:
     
                    # E step
                    snp_calls = expectation(base_counts, ref_seq,
                                            snp_base, error_prob,
                                            self.bonf_factor, self.sig_thresh,
                                            from_base)
                    # M step
                    snp_pos = [pos for (pos, pvalue) in snp_calls]
                    new_error_prob = maximization(base_counts, ref_seq,
                                                  snp_base, snp_pos,
                                                  from_base)     
         
                    # break if error probabilities converge or if maximum
                    # number of iterations reached
                    #
                    if abs(new_error_prob - error_prob) < self.convergence_epsilon \
                       or num_iter == MAX_ITER:
         
                        if num_iter == MAX_ITER:
                            reason = "maximum number of iterations reached"
                        else:
                            reason = "error probabilities converged"
                        LOG.info("EM iteration for %s>%s done (%s)."
                                 " final error prob for %s>%s = %f" % (
                                     from_base, snp_base,
                                     reason, from_base,
                                     snp_base, new_error_prob))
                        
                        self.error_probs[from_base][snp_base] = new_error_prob
                        final_snp_calls.extend([(pos, snp_base, pvalue)
                                           for (pos, pvalue) in snp_calls])
     
                        break # while True
                    
                    # update the base error and run EM again
                    #
                    num_iter += 1
                    LOG.info("EM iteration no %d: updating error_prob for"
                             " %s>%s from %f to %f..." % (
                                 num_iter, from_base, snp_base,
                                 error_prob, new_error_prob))
                    error_prob = new_error_prob

        return (final_snp_calls)
     



    def call_snp_in_column(self, col_coord, col_bases, ref_base=None):
        """Call SNP in a column given read bases. If ref_base is not
        given a consensus base will be derived.
        """

        assert len(self.error_probs), (
            "Error probabilities haven't"
            " been initialized. Please call em_training() or"
            " load_error_probs() first.")

        snp_calls = []

        (col_base_counts, cons_base_est) = count_bases(col_bases)
        if not ref_base:
            ref_base = cons_base_est
        if col_base_counts.has_key('N'):
            del col_base_counts['N']

        coverage = sum(col_base_counts.values())
        if coverage == 0:
            return []

        if self.num_param == 12:
            from_bases = self.error_probs.keys()
            from_bases.remove('N')
        elif self.num_param == 4:
            from_bases = 'N'
        else:
            raise ValueError, ("Unsupported number of parameters chosen")

        for from_base in from_bases:
            for (snp_base, error_prob) in self.error_probs[from_base].iteritems():
                for (rel_pos, pvalue) in expectation([col_base_counts],
                                                     ref_base, snp_base, error_prob,
                                                     self.bonf_factor, self.sig_thresh,
                                                     from_base):
                    info_dict = dict()
                    info_dict['pvalue'] = pvalue
                    info_dict['em-error-prob-%s-to-%s' % (
                        from_base, snp_base)] = error_prob
                    info_dict['coverage'] = coverage
                    for (k, v) in col_base_counts.iteritems():
                        info_dict["basecount-%s" % k] = v
                    snpcall = snp.ExtSNP(col_coord,
                                         ref_base,
                                         snp_base,
                                         col_base_counts[snp_base]/coverage,
                                         info_dict)
                    LOG.info("LofreqNQ SNP %s." % (snpcall))
                    snp_calls.append(snpcall)

        return snp_calls



if __name__ == "__main__":
    # FIXME Add tests
    pass
