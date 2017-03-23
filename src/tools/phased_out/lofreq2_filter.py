#!/usr/bin/env python
"""Apply number of filters to given list of SNVs.

Each filter is applied to all SNVs, i.e. not just the previously
PASSED ones!
"""
from builtins import zip
from builtins import range


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2013 Genome Institute of Singapore"
__license__ = "GPL2"


#--- standard library imports
#
import sys
import logging
import os
# optparse deprecated from Python 2.7 on. need optparse here to mess
# with the default options if needed.
from optparse import OptionParser#, SUPPRESS_HELP
import gzip


#--- third-party imports
#
#/

#--- project specific imports
#
try:
    import lofreq2_local
except ImportError:
    pass

try:
    from lofreq_star import vcf
except ImportError:
    sys.stderr.write("FATAL(%s): Couldn't find LoFreq's vcf module."
                     " Are you sure your PYTHONPATH is set correctly (= %s)?\n" % (
                         (sys.argv[0], os.environ['PYTHONPATH'])))
    sys.exit(1)
from lofreq_star import multiple_testing
from lofreq_star import fdr
from lofreq_star.utils import prob_to_phredqual, phredqual_to_prob, MAX_INT


# invocation of ipython on exceptions
#import sys, pdb
#from IPython.core import ultratb
#sys.excepthook = ultratb.FormattedTB(mode='Verbose',
#                                     color_scheme='Linux', call_pdb=1)


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')



def win_filter(snvs_on_cur_chrom, win_size, vcf_info_id):
    """Makes snv INFO[vcf_info_id] with 0 if there is a neigbouring
    snv within win_size, otherwise 1
    """

    # make sure snvs are sorted
    # snvs_on_cur_chrom = sorted(snvs_on_cur_chrom, 
    #                            key = lambda x: x.POS)
    # disabled because possible overkill: lofreq produces ordered
    # lists by default since BAM file is sorted

    for (ci, cur_snv) in enumerate(snvs_on_cur_chrom):
        
        # prev_snv: snv at < pos on same chrom
        prev_snv = None
        for pi in reversed(range(ci)):
            tmp = snvs_on_cur_chrom[pi]
            assert tmp.POS <= cur_snv.POS
            assert tmp.CHROM == cur_snv.CHROM
            if tmp.POS != cur_snv.POS:
                prev_snv = tmp
                break
            
        # next_snv: snv at > pos on same chrom
        next_snv = None
        for ni in range(ci+1, len(snvs_on_cur_chrom)):
            tmp = snvs_on_cur_chrom[ni]
            assert tmp.POS >= cur_snv.POS
            assert tmp.CHROM == cur_snv.CHROM
            if tmp.POS != cur_snv.POS:
                next_snv = tmp
                break
            
        LOG.debug("prev_snv=%d cur_snv=%d next_snv=%d" % (
            prev_snv.POS if prev_snv else -1,
            cur_snv.POS,
            next_snv.POS if next_snv else -1))

        cur_snv.INFO[vcf_info_id] = 1 # pass by default
        if prev_snv != None:
            if cur_snv.POS-prev_snv.POS <= win_size:
                cur_snv.INFO[vcf_info_id] = 0
        if next_snv != None:
            if next_snv.POS-cur_snv.POS <= win_size:
                cur_snv.INFO[vcf_info_id] = 0
            
                
    
def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog [Options]\n" \
      + "\n" + __doc__
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true",
                      dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true",
                      dest="debug",
                      help="enable debugging")
    parser.add_option("-i", "--vcf_in",
                      dest="vcf_in",
                      help="Input vcf file (gzip supported; - for stdin).")
    default = "-"
    parser.add_option("-o", "--outfile",
                      dest="vcf_out",
                      default=default,
                      help="Output vcf file (gzip supported; - for stdout)."
                      " Default = %s)" % default)

    parser.add_option("-p", "--pass-only",
                      action="store_true",
                      dest="pass_only",
                      help="Only print PASSed variants")

    default = "holm-bonf"
    parser.add_option("", "--strandbias",
                      default=default,
                      help="Filter variants with strandbias."
                      " Valid values are 'bonf' (Bonferroni),"
                      " 'holm-bonf' (Holm-Bonferroni), an integer value"
                      " or 'off'. If 'bonf' or 'holm-bonf', variants with"
                      " accordingly corrected strand-bias pvalue"
                      " < strandbias-alpha will be filtered. If an int"
                      " was given, variants with strand-bias phred-scores"
                      " larger than this value will be filtered."
                      " (default: %s)" % default)
    default = 0.05
    parser.add_option("", "--strandbias-alpha",
                      default=default,
                      type="float",
                      help="Alpha/significance-level for strandbias testing."
                      " (applies only to 'bonf' and 'holm-bonf'; "
                      " default: %s)" % default)
    
    default = 10
    parser.add_option("", "--min-cov",
                      dest="min_cov",
                      type='int',
                      default=default,
                      help="Filter variant if coverage is"
                      " below this value (int; default = %d)" % default)
    parser.add_option("", "--max-cov",
                      dest="max_cov",
                      type='int',
                      help="Filter variant if coverage is"
                      " above this cap (int)")
    parser.add_option("", "--min-af",
                      dest="min_af",
                      type="float",
                      help="Filter if (allele) freq is"
                      " below this threshold (float)")

    parser.add_option("", "--snv-qual",
                      help="Filter variants based on quality. Valid values"
                      " are 'fdr', 'bonf', 'holmbonf' or an integer value."
                      " If FDR Benjamini-Hochberg correction will be used."
                      " If 'bonf' Bonferroni- and if 'holm-bonf'"
                      " Holm-Bonferroni-correction will be used."
                      " If an int was given, variants with a phred-score"
                      " below this value will be filtered")
    parser.add_option("", "--snv-qual-alpha",
                      type="float",
                      help="Alpha/significance threshold for multiple testing"
                      " correction routines during SNV quality filtering."
                      " Only applies to 'bonf', 'holm-bonf' and 'fdr'")
    parser.add_option("", "--snv-qual-numtests",
                      type="int",
                      help="Set number of tests for multiple testing"
                      " correction routines during SNV quality filtering."
                      " Only applies to 'fdr', 'bonf' and 'holm-bonf'."
                      " Defaults to number of pvalues.")

    parser.add_option("", "--window-size",
                      dest="window_size",
                      type='int',
                      help='Ignore variants, if another'
                      ' variant is present within a window of this size'
                      ' (ignoring multi-allelic vars at same pos).')

    #parser.add_option("--force",
    #                  #help=SUPPRESS_HELP,
    #                  dest="force_overwrite", action="store_true")

    return parser



def main():
    """main function
    """

    tmp_vcf_markup = []

    parser = cmdline_parser()

    # WARNING: undocumented arg to remove all defaults (and the reason
    # why we have to use OptParse)
    if '--no-defaults' in sys.argv:
        for (k, v) in list(parser.defaults.items()):
            parser.defaults[k] = None
        sys.argv = [x for x in sys.argv if x != "--no-defaults"]

    (opts, args) = parser.parse_args()

    if len(args):
        parser.error("Unrecognized arguments found: %s." % (
            ' '.join(args)))
        sys.exit(1)


    if opts.verbose:
        LOG.setLevel(logging.INFO)
    if opts.debug:
        LOG.setLevel(logging.DEBUG)

    for (in_file, descr) in [(opts.vcf_in, "VCF")]:
        if not in_file:
            parser.error("%s input file argument missing." % descr)
            sys.exit(1)
        if not os.path.exists(in_file) and in_file != "-":
            sys.stderr.write(
                "file '%s' does not exist.\n" % in_file)
            sys.exit(1)

    for (out_file, descr) in [(opts.vcf_out, "VCF output file")]:
        if not out_file:
            parser.error("%s output file argument missing." % descr)
            sys.exit(1)
        if os.path.exists(out_file) and out_file!="-":
            sys.stderr.write("Cowardly refusing to overwrite existing"
                             " output file '%s'.\n" % out_file)
            sys.exit(1)


    if opts.vcf_in == '-':
        vcf_reader = vcf.VCFReader(sys.stdin)
    else:
        if opts.vcf_in[-3:] == '.gz':
            vcf_reader = vcf.VCFReader(gzip.open(opts.vcf_in,'r'))
        else:
            vcf_reader = vcf.VCFReader(open(opts.vcf_in,'r'))
    snvs = [r for r in vcf_reader]
    LOG.info("Parsed %d SNVs from %s" % (len(snvs), opts.vcf_in))


    
    # list of tuples: first element is a filter func, which takes a
    # snv and a filter-id as input. second is the filter id. variant
    # will be marked as filtered if func returns True
    filters = []

    
    if opts.min_af != None:
        vcf_filter = vcf._Filter(
            id=("minaf%f" % opts.min_af).rstrip('0'),
            desc="Minimum allele frequency")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

        filters.append((
            lambda s, f_id: f_id if s.INFO['AF'] < opts.min_af else None,
            vcf_filter.id
            ))


    if opts.max_cov != None:
        if not all(['DP' in s.INFO for s in snvs]):
            LOG.error("At least one SNV was not annotated with depth info (DP)"
                      " (was this file produced with LoFreq?).")
            sys.exit(1)

        vcf_filter = vcf._Filter(
            id="maxcov%d" % opts.max_cov,
            desc="Maximum coverage")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

        filters.append((
            lambda s, f_id: f_id if s.INFO['DP'] > opts.max_cov else None,
            vcf_filter.id
            ))


    if opts.min_cov != None:
        if not all(['DP' in s.INFO for s in snvs]):
            LOG.error("At least one SNV was not annotated with depth info (DP)"
                      " (was this file produced with LoFreq?).")
            sys.exit(1)

        vcf_filter = vcf._Filter(
            id="mincov%d" % opts.min_cov,
            desc="Minimum coverage")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

        filters.append((
            lambda s, f_id: f_id if s.INFO['DP'] < opts.min_cov else None,
            vcf_filter.id
            ))

    # structured as opts.snv_qual filtering, but keeps corrected
    # values.
    if opts.strandbias != None:

        if opts.strandbias in ['bonf', 'holm-bonf']:
            if not opts.strandbias_alpha:
                LOG.fatal("Need alpha/significance threshold for strandbias"
                          " multiple testing correction")
                sys.exit(1)

            vcf_filter = vcf._Filter(
                id="strandbias%s" % opts.strandbias.replace("-", ""),
                desc="Strand-bias filter (%s corrected < %g)" % (
                    opts.strandbias, opts.strandbias_alpha))
            vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

            if opts.strandbias == 'bonf':
                vcf_info_id = "SBBC"
            elif opts.strandbias == 'holm-bonf':
                vcf_info_id = "SBHBC"
            else:
                raise ValueError
            vcf_info = vcf._Info(
                id=vcf_info_id, num=1, type='Integer',
                desc="Strand-bias %s corrected" % opts.strandbias)
            vcf_reader.infos[vcf_info.id] = vcf_info

            try:
                pvals = (phredqual_to_prob(s.INFO['SB']) for s in snvs)
            except (KeyError, AssertionError) as e:
                LOG.error("At least one SNV was not annotated properly with"
                          " strandbias info (SB)"
                          " (was this file produced with LoFreq?)"
                          " You will need to switch strandbias filtering off")
                sys.exit(1)

            if opts.strandbias == 'bonf':
                corr_pvals = multiple_testing.Bonferroni(
                    pvals).corrected_pvals
            elif opts.strandbias == 'holm-bonf':
                corr_pvals = multiple_testing.HolmBonferroni(
                    pvals).corrected_pvals
            else:
                raise ValueError
            for (cp, s) in zip(corr_pvals, snvs):
                s.INFO[vcf_info.id] = prob_to_phredqual(cp)
                if s.INFO[vcf_info.id] > MAX_INT:
                    s.INFO[vcf_info.id] = MAX_INT

            filters.append((
                lambda s, f_id: f_id if s.INFO[vcf_info.id] > prob_to_phredqual(opts.strandbias_alpha) else None,
                vcf_filter.id
                ))

        # int
        elif opts.strandbias != 'off':
            try:
                max_strandbias_phred = int(opts.strandbias)
                assert max_strandbias_phred >= 0
            except (ValueError, AssertionError) as e:
                LOG.fatal("Invalid strandbias argument: %s" % (opts.strandbias))
                sys.exit(1)

            vcf_filter = vcf._Filter(
                max_strandbias_phred = int(
                id="sbp%d" % opts.max_strandbias_phred,
                desc="Phred-based strand-bias filter (max)"))
            vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

            filters.append((
                lambda s, f_id: f_id if float(s.INFO['SB']) > opts.max_strandbias_phred else None,
                vcf_filter.id
                ))
            

    # structured as opts.strandbias filtering, but doesn't keep
    # corrected values.
    if opts.snv_qual != None:

        if opts.snv_qual in ['bonf', 'holm-bonf', 'fdr']:
            if not opts.snv_qual_alpha:
                LOG.fatal("Need alpha/significance threshold for snv quality"
                          " multiple testing correction")
                sys.exit(1)

            vcf_filter = vcf._Filter(
                id="snvqual%s" % opts.snv_qual.replace("-", ""),
                desc="SNV quality filter (%s corrected < %g)" % (
                    opts.snv_qual, opts.snv_qual_alpha))
            vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

            vcf_info_id = "SNVQUALPASS" # tmp markup
            tmp_vcf_markup.append(vcf_info_id)

            pvals = []
            pidx = []
            for (i, s) in enumerate(snvs):
                # if qual is not NA, convert to pvalue, else don't
                # use filter (set filter to NA)
                if s.QUAL != '.':
                    pvals.append(phredqual_to_prob(s.QUAL))
                    pidx.append(i)
                    s.INFO[vcf_info_id] = 0
                else:
                    s.INFO[vcf_info_id] = '.'

            if opts.snv_qual == 'bonf':
                for (i, p) in enumerate(
                        multiple_testing.Bonferroni(
                            pvals, n=opts.snv_qual_numtests).corrected_pvals):
                    if p <= opts.snv_qual_alpha:
                        snvs[pidx[i]].INFO[vcf_info_id] = 1

            elif opts.snv_qual == 'holm-bonf':
                for (i, p) in enumerate(
                        multiple_testing.HolmBonferroni(
                            pvals, n=opts.snv_qual_numtests).corrected_pvals):
                    if p <= opts.snv_qual_alpha:
                        snvs[pidx[i]].INFO[vcf_info_id] = 1
 
            elif opts.snv_qual == 'fdr':
                for i in fdr.fdr(pvals, a=opts.snv_qual_alpha, 
                                 n=opts.snv_qual_numtests):
                    snvs[pidx[i]].INFO[vcf_info_id] = 1

            else:
                raise ValueError

            filters.append((
                lambda s, f_id: f_id if s.INFO[vcf_info_id] != '.' and s.INFO[vcf_info_id] == 0 else None,
                vcf_filter.id
                ))

        elif opts.snv_qual != 'off':
            try:
                min_qual = int(opts.snv_qual)
                assert min_qual >= 0
            except (ValueError, AssertionError) as e:
                LOG.fatal("Invalid snv quality argument: %s" % (opts.snv_qual))
                sys.exit(1)

            vcf_filter = vcf._Filter(
                id="minqual%d" % min_qual,
                desc="Minimum SNV quality")
            vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

            filters.append((
                lambda s, f_id: f_id if s.QUAL != '.' and s.QUAL < min_qual else None,
                vcf_filter.id
                ))


    if opts.window_size != None:
        vcf_filter = vcf._Filter(
            id="snvwin%d" % opts.window_size,
            desc="SNV window filter (SNVs within %d bp distance)" % (
                opts.window_size))
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

        vcf_info_id = "SNVWINPASS" # tmp markup
        tmp_vcf_markup.append(vcf_info_id)

        snvs_on_cur_chrom = []
        last_chrom = None
        seen_chroms = []
        for (i, cur_snv) in enumerate(snvs): # assumes snvs are sorted by chrom
            if i == 0:
                last_chrom = cur_snv.CHROM
                
            if cur_snv.CHROM != last_chrom:
                assert cur_snv.CHROM not in seen_chroms, (
                    "SNV input not ordered by chromosome."
                    " Sure this file was procuced by LoFreq?")
                win_filter(snvs_on_cur_chrom, opts.window_size, vcf_info_id)
                seen_chroms.append(last_chrom)
                last_chrom = cur_snv.CHROM
                snvs_on_cur_chrom = [cur_snv]
                
            else:
                snvs_on_cur_chrom.append(cur_snv)

        # don't forget last chrom
        win_filter(snvs_on_cur_chrom, opts.window_size, vcf_info_id)

        
        filters.append((
            lambda s, f_id: f_id if s.INFO[vcf_info_id] != '.' and s.INFO[vcf_info_id] == 0 else None,
            vcf_filter.id
            ))
            

    # The actual filtering: if filter function returns 1 the
    # corresponding snv has to be filtered
    #
    # FIXME can't this be done easier with map()?
    #
    if len(filters) == 0:
        LOG.error("No filters activated.")
        sys.exit(1)

    #import pdb; pdb.set_trace()
    for (filter_func, filter_id) in filters:
        for (i, s) in enumerate(snvs):
            f = filter_func(s, filter_id)
            if f:
                # just s = s.__replace() can't work
                if s.FILTER == '.' or s.FILTER == 'PASS':
                    snvs[i] = s._replace(FILTER=f)
                else:
                    snvs[i] = s._replace(FILTER="%s;%s" % (s.FILTER, f))

                        
    
    # should all also work if we get already PASSed input

    n_passed = 0
    for (i, s) in enumerate(snvs):
        if s.FILTER == '.':
            snvs[i] = s._replace(FILTER="PASS")
            n_passed += 1
    LOG.info("%d SNVs passed all filters." % n_passed)

    # remove temporary markup
    for tmpkey in tmp_vcf_markup:
        for s in snvs:
            if tmpkey in s.INFO:
                del s.INFO[tmpkey]

    if opts.pass_only:
        snvs = (s for s in snvs if s.FILTER == 'PASS')

    if opts.vcf_out == '-':
        fh_out = sys.stdout
    else:
        if opts.vcf_out[-3:] == '.gz':
            fh_out = gzip.open(opts.vcf_out, 'w')
        else:
            fh_out = open(opts.vcf_out, 'w')

    vcf_writer = vcf.VCFWriter(fh_out)
    vcf_writer.meta_from_reader(vcf_reader)
    vcf_writer.write(snvs)

    if fh_out != sys.stdout:
        fh_out.close()


if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
