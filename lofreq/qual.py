#!/usr/bin/env python
"""Quality aware SNP caller.

"""


#--- standard library imports
#
from __future__ import division
import logging

#--- third-party imports
#
# /

#--- project specific imports
#
# /
from lofreq_ext import snpcaller_qual
from lofreq import snp
from lofreq import conf

__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = ""
__license__ = ""
__credits__ = [""]
__status__ = ""


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')



class QualBasedSNPCaller(object):
    """
    Quality/probabilty based SNP caller.
    """

    def __init__(self,
                 noncons_default_qual = conf.NONCONS_DEFAULT_QUAL,
                 noncons_filter_qual = conf.NONCONS_FILTER_QUAL,
                 ign_bases_below_q = conf.DEFAULT_IGN_BASES_BELOW_Q,
                 bonf_factor = 1,
                 sig_thresh = conf.DEFAULT_SIG_THRESH):
        """
        init function

        noncons_default_qual is the assumed quality of non-consensus bases,
        for which it is to be determined if they are an error or a SNP

        noncons_filter_qual is a threshold. non-consensus basepairs
        below this quality will be completely ignored.
        """

        assert isinstance(noncons_default_qual, int)
        assert isinstance(noncons_filter_qual, int)
        assert isinstance(bonf_factor, int)

        self.replace_noncons_quals = True
        self.noncons_default_qual = noncons_default_qual
        self.noncons_filter_qual = noncons_filter_qual
        self.ign_bases_below_q = ign_bases_below_q
        self.bonf_factor = bonf_factor
        self.sig_thresh = sig_thresh
        LOG.debug("New QualBasedSNPCaller: noncons default qual = %d"
                  ", noncons filter qual = %d"
                  ", ign_bases_below_q = %d"
                  ", bonferroni factor = %d"
                  ", sign. level = %f"
                  " and replace_noncons_quals = %s" % (
                self.noncons_default_qual,
                self.noncons_filter_qual,
                self.ign_bases_below_q,
                self.bonf_factor,
                self.sig_thresh,
                self.replace_noncons_quals))



    def call_snp_in_column(self, col_coord, base_qual_hist, ref_base):
        """Call SNP in one pileup colum given it's bases and qualities.
        """

        ret_snp_calls = []

        assert len(base_qual_hist.keys())==4, (
            "Expected exactly four bases as keys for base_qual_hist")
        for base in base_qual_hist.keys():
            assert base in 'ACGT', (
                "Only allowed bases/keys are A, C, G or T, but not %s" % base)
        assert ref_base in 'ACGT', (
            "consensus base must be one of A, C, G or T, but not %s" % ref_base)


        # count non-consensus bases and remove if their quality is
        # below noncons_filter_qual or ign_bases_below_q
        #
        noncons_counts_dict = dict()
        noncons_filterq = max(self.noncons_filter_qual, self.ign_bases_below_q)
        for b in base_qual_hist.keys():
            if b == ref_base:
                continue
            noncons_counts_dict[b]  = sum(c for (q, c) in base_qual_hist[b].iteritems() 
                                     if q >= noncons_filterq)

        # return if no non-consbases left after filtering
        if sum(noncons_counts_dict.values()) == 0:
            LOG.debug("Consensus bases only. Early exit...")
            return []

        # get list of consensus qualities but remove if quality is
        # below ign_bases_below_q
        #
        cons_quals = []
        for (q, c) in base_qual_hist[ref_base].iteritems():
            if q < self.ign_bases_below_q:
                continue
            cons_quals.extend([q]*c)
        cons_count = len(cons_quals)

        # add a default value for each non consensus base
        # (WARNING: reusing cons_quals!)
        base_quals = cons_quals
        base_quals.extend([self.noncons_default_qual] * sum(noncons_counts_dict.values()))
      
        
        # need noncons counts and bases in order since only counts are
        # handed down to snpcaller_qual
        noncons_bases = []
        noncons_counts = []
        for (b, c) in noncons_counts_dict.iteritems():
            noncons_bases.append(b)
            noncons_counts.append(c)
        noncons_bases = tuple(noncons_bases)
        noncons_counts = tuple(noncons_counts)

        # using sorted base_quals might in theory be better for speedup and
        # numerical stability. in practice I don't see a difference
        pvalues = snpcaller_qual(sorted(base_quals), noncons_counts,
                                 self.bonf_factor, self.sig_thresh)
        # reported pvalues are already bonferroni corrected

        # setup info dictionary shared between different alleles for
        # this position
        #
        info_dict = dict()

        # note: coverage after filtering!
        coverage = sum(noncons_counts_dict.values()) + cons_count
        info_dict['coverage'] = coverage

        for (k, v) in noncons_counts_dict.iteritems():
            info_dict["basecount-%s" % k] = v
        info_dict["basecount-%s" % ref_base] = cons_count

        # check pvalues for all possible mutations and report if significant
        #
        for (base, count, pvalue) in zip(noncons_bases, noncons_counts, pvalues):
            if pvalue < self.sig_thresh:
                info_dict['pvalue'] = pvalue
                snpcall = snp.ExtSNP(col_coord, ref_base, base,
                                     count/float(coverage), info_dict)
                ret_snp_calls.append(snpcall)


        return ret_snp_calls


        
if __name__ == "__main__":
    pass
