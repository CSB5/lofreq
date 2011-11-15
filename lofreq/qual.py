#!/usr/bin/env python
"""Quality/probabilty based SNP caller.

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
from lofreq.utils import count_bases


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


DEFAULT_SIG_THRESH = 0.01
NONCONS_DEFAULT_QUAL = 20
NONCONS_FILTER_QUAL = 20
DEFAULT_IGN_BASES_BELOW_Q = 3



class QualBasedSNPCaller(object):
    """
    Quality/probabilty based SNP caller.
    """

    def __init__(self,
                 ign_bases_below_q = DEFAULT_IGN_BASES_BELOW_Q,
                 noncons_default_qual = NONCONS_DEFAULT_QUAL,
                 noncons_filter_qual = NONCONS_FILTER_QUAL,
                 bonf_factor = 1.0,
                 sig_thresh = DEFAULT_SIG_THRESH):
        """
        init function

        noncons_default_qual is the assumed quality of non-consensus bases,
        for which it is to be determined if they are an error or a SNP

        noncons_filter_qual is a threshold. non-consensus basepairs
        below this quality will be completely ignored.
        """

        self.replace_noncons_quals = True
        self.noncons_default_qual = noncons_default_qual
        self.noncons_filter_qual = noncons_filter_qual
        self.bonf_factor = bonf_factor
        self.sig_thresh = sig_thresh
        # Illumina indicates errors with Q2 (or lower). Those bases
        # should not be used for downstream analysis. Therefore we
        # could use 3 as cutoff. GATK does by default not correct Q<5,
        self.ign_bases_below_q = ign_bases_below_q

        LOG.debug("New QualBasedSNPCaller: ign_bases_below_q = %d"
                 ", noncons default qual = %d"
                 ", noncons filter qual = %d"
                 ", bonferroni factor = %f"
                 ", sign. level = %f"
                 " and replace_noncons_quals = %s" % (
                     self.ign_bases_below_q,
                     self.noncons_default_qual,
                     self.noncons_filter_qual,
                     self.bonf_factor,
                     self.sig_thresh,
                     self.replace_noncons_quals))



    def call_snp_in_column(self, col_coord, col_bases, col_base_quals,
                           cons_base=None):
        """
        Call SNP in one pileup colum given it's bases and qualities.
        """

        #debug_cols = [3857]
        debug_cols = []
        # for testing:
        #col_base_quals = [random.randrange(5, 30)
        #                  for x in range(len(col_bases))]

        ret_snp_calls = []

        assert len(col_bases) == len(col_base_quals), (
            "Need a quality for each base")

        (raw_base_counts, dummy) = count_bases(col_bases.upper())
        if col_coord+1 in debug_cols:
            LOG.critical("Raw base counts at %d are=%s" % (
                col_coord+1, raw_base_counts))


        # make a cleaned copy of input bases and qualities. remove N's
        # and their qualities and uppercase bases and ignore anything
        # with quality below ign_bases_below_q
        #
        bases_and_quals = [(b, q)
                           for (b, q) in zip(col_bases.upper(), col_base_quals)
                           if q > self.ign_bases_below_q]
        ign_bases_below_q_count = len(col_base_quals)-len(bases_and_quals)

        bases_and_quals = [(b, q)
                           for (b, q) in bases_and_quals
                           if b != "N"]


        #LOG.critical("Using fixed values (see mqual_to_em4/fixed-em-probs)")
        #bases_and_quals = [(b, 30)
        #                  for (b, q) in bases_and_quals]

        #LOG.critical("Using fixed values (see mqual_to_em4/fixed-em-probs)")
        #bases_and_quals_rpl = [(b, 30)
        #                   for (b, q) in bases_and_quals
        #                   if b =='A']
        #bases_and_quals_rpl.extend([(b, 31)
        #                   for (b, q) in bases_and_quals
        #                   if b =='C'])
        #bases_and_quals_rpl.extend([(b, 32)
        #                   for (b, q) in bases_and_quals
        #                   if b =='T'])
        #bases_and_quals_rpl.extend([(b, 33)
        #                   for (b, q) in bases_and_quals
        #                   if b =='G'])
        #bases_and_quals = bases_and_quals_rpl


        # remove non-consensus bases (incl. their qualities) if their
        # quality is below filter-threshold
        #
        if not cons_base:
            (dummy, cons_base_est) = count_bases(
                [b for (b, q) in bases_and_quals])
            cons_base = cons_base_est
        bases_and_quals = [(b, q) for (b, q) in bases_and_quals
                           if b == cons_base or q >= self.noncons_filter_qual]


        # construct a list of consensus quality values and
        # non-consensus quality values replaced by self.noncons_qual
        # (or something equivalent, eg the cons average)
        #
        if self.replace_noncons_quals:
            cons_base_quals = [q for (b, q) in bases_and_quals if b==cons_base]
            num_noncons = len(bases_and_quals) - len(cons_base_quals)

            # used to be just: base_quals.extend(num_noncons *
            # [self.noncons_default_qual])
            #
            noncons_base_quals = num_noncons * [self.noncons_default_qual]
            # Use average consensus quality instead?
            if False:
                avg_cons_qual = sum(cons_base_quals)/len(cons_base_quals)
                noncons_base_quals = num_noncons * [avg_cons_qual]

            base_quals = cons_base_quals
            base_quals.extend(noncons_base_quals)

        else:
            base_quals = [q for (b, q) in bases_and_quals]
            num_noncons = len([b for (b, q) in bases_and_quals
                               if b != cons_base])

        # return if no non-consbases left after filtering
        if num_noncons <= 0:
            LOG.debug("Consensus bases only. early exit...")
            return []


        # get non consensus counts and bases
        #
        (base_counts, dummy) = count_bases(
            [b for (b, q) in bases_and_quals])
        LOG.debug("testing col %d with base_count %s" % (
            col_coord, base_counts))

        coverage = sum(base_counts.values())
        # the following might look ugly, but we need a tuple of counts in the
        # same order as noncons_bases for calling the snp caller
        noncons_bases = [b for b in base_counts.keys() if b != cons_base]
        noncons_counts = tuple([base_counts[b] for b in noncons_bases])

        # print quality histogram
        if False:
            qual_list = [q for (b, q) in bases_and_quals]
            for q in range(min(qual_list), max(qual_list)+1):
                LOG.debug("quality histogram at %d: %d x %d" % (
                    col_coord+1, qual_list.count(q), q))


        # using sorted base_quals might in theory be better for speedup and
        # numerical stability. in practice I don't see a difference
        pvalues = snpcaller_qual(sorted(base_quals), noncons_counts,
                                 self.bonf_factor, self.sig_thresh)
        # reported pvalues are already correct with bonferroni factor

        # setup info dictionary shared between different alleles for
        # this position
        #
        info_dict = dict()
        for k, v in raw_base_counts.iteritems():
            info_dict["raw-basecount-%s" % k] = v
        info_dict['raw-coverage'] = sum(raw_base_counts.values())

        for k, v in base_counts.iteritems():
            info_dict["filtered-basecount-%s" % k] = v
        info_dict['filtered-coverage'] = sum(base_counts.values())

        info_dict['ign-bases-below-qual-%d' % self.ign_bases_below_q] = ign_bases_below_q_count

        # check pvalues for all possible mutations and report if significant
        #
        for (base, count, pvalue) in zip(noncons_bases, noncons_counts, pvalues):
            if pvalue < self.sig_thresh:
                info_dict['pvalue'] = pvalue

                snpcall = snp.ExtSNP(col_coord,
                                     cons_base, base,
                                     count/coverage, info_dict)
                ret_snp_calls.append(snpcall)
                LOG.info("Qual Detected SNP %s." % (snpcall))


        #debug_cols = []
        if col_coord+1 in debug_cols:
            (tmp_base_counts, tmp_cons_base_est) = count_bases(
                [b for (b, q) in bases_and_quals])
            LOG.critical("Filtered base counts at %d are=%s" % (
                col_coord+1, tmp_base_counts))
            for base in set([b for (b, q) in bases_and_quals]):
                quals = [q for (b, q) in bases_and_quals if b == base]
                LOG.critical("Average qualities for %s: %f" % (
                    base, sum(quals)/len(quals)))
            LOG.critical("Debugging based on hardcoded SNP positions:"
                         " zero-offset col_coord=%d" % col_coord)
            import pdb
            pdb.set_trace()

        return ret_snp_calls


        
if __name__ == "__main__":
    # FIXME add tests
    pass
