#!/usr/bin/env python
"""Apply number of filters to given list of SNVs.

Note, that each filter is applied to all SNVs, i.e. not just the
previously PASSED ones!
"""


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2013 Genome Institute of Singapore"
__license__ = "Free for non-commercial use"


#--- standard library imports
#
import sys
import logging
import os
# optparse deprecated from Python 2.7 on
from optparse import OptionParser, SUPPRESS_HELP
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
from lofreq_star.utils import prob_to_phredqual, phredqual_to_prob

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
                      help="Output vcf file (gzip supported; - for stdout). Default = %s)" % default)

    parser.add_option("-p", "--pass-only",
                      action="store_true",
                      dest="pass_only",
                      help="Only print PASSed variants")

    default = "holm-bonf"
    parser.add_option("", "--strandbias",
                      default="holm-bonf",
                      help="Filter variant with strandbias."
                      " Valid values are 'bonf' (Bonferroni),"
                      " 'holm-bonf' (Holm-Bonferroni), an integer or off."  
                      " If 'bonf' or 'holm-bonf', variants with accordingly"
                      " corrected strand-bias pvalue <0.05 will be filtered."
                      " Otherwise, a variant will be filtered if the strand-bias"
                      " phred-score is larger than the given int."
                      " (default: %s)" % default)
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
                      help="Optional: Filter variant if coverage is"
                      " above this cap (int)")
    parser.add_option("", "--min-af",
                      dest="min_af", 
                      type="float",
                      help="Optional: Filter if (allele) freq is"
                      " below this threshold (float)")
    parser.add_option("", "--snv-phred", 
                      dest="min_snv_phred",
                      type='int',
                      help="Optional: Filter variant if its phred-score"
                      " is below this value (int)")
    parser.add_option("", "--snv-fdr", 
                      dest="snv_fdr",
                      type='float',
                      help="Optional: Set alpha for FDR (Benjamini-Hochberg) of SNV pvalues")
    #parser.add_option("", "--window-size",
    #                  dest="window_size",
    #                  type='int',
    #                  help="Optional: Filter variants if at least"
    #                  " one or more were called within this window size")

    parser.add_option("--force", help=SUPPRESS_HELP,
                      dest="force_overwrite", action="store_true") 

    return parser



def main():


    parser = cmdline_parser()            

    # warning: undocumented arg to remove all defaults
    if '--no-defaults' in sys.argv:
        for (k, v) in parser.defaults.items(): 
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

    for (in_file, descr) in [(opts.vcf_in, "VCF input file")]:
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
            sys.stderr.write(
                "Cowardly refusing to overwrite existing output file '%s'.\n" % out_file)
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

    if opts.strandbias and opts.strandbias == 'bonf':            
        vcf_filter = vcf._Filter(
            id="sbb", 
            desc="Strand-bias filter on Bonferroni corrected p-values")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer
        
        vcf_info = vcf._Info(
            id="SBBC", num=1, type='Integer',
            desc="Strand-bias Bonferroni corrected")        
        vcf_reader.infos[vcf_info.id] = vcf_info
        
        pvals = [phredqual_to_prob(s.INFO['SB']) for s in snvs]
        corr_pvals = multiple_testing.Bonferroni(pvals).corrected_pvals
        for (cp, s) in zip(corr_pvals, snvs):
            s.INFO[vcf_info.id] = prob_to_phredqual(cp)
            if s.INFO[vcf_info.id] > sys.maxint:
                s.INFO[vcf_info.id] = sys.maxint
                
        filters.append((
            lambda s, f_id: f_id if s.INFO["SBBC"] > 13 else None,
            vcf_filter.id
            ))

        
    elif opts.strandbias and opts.strandbias == 'holm-bonf':
        vcf_filter = vcf._Filter(
            id="sbh", 
            desc="Strand-bias filter on Holm-Bonferroni corrected p-values")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer
        
        vcf_info = vcf._Info(
            id="SBHC", num=1, type='Integer',
            desc="Strand-bias Holm-Bonferroni corrected")        
        vcf_reader.infos[vcf_info.id] = vcf_info
        
        pvals = [phredqual_to_prob(s.INFO['SB']) for s in snvs]
        corr_pvals = multiple_testing.HolmBonferroni(pvals).corrected_pvals
        #import pdb; pdb.set_trace()
        for (cp, s) in zip(corr_pvals, snvs):
            s.INFO[vcf_info.id] = prob_to_phredqual(cp)
            if s.INFO[vcf_info.id] > sys.maxint:
                s.INFO[vcf_info.id] = sys.maxint
            
        filters.append((
            lambda s, f_id: f_id if s.INFO["SBHC"] > 13 else None,
            vcf_filter.id
            ))                

    elif opts.strandbias and opts.strandbias != 'off':
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
        vcf_filter = vcf._Filter(
            id="maxcov%d" % opts.max_cov, 
            desc="Maximum coverage")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

        filters.append((
            lambda s, f_id: f_id if s.INFO['DP'] > opts.max_cov else None,
            vcf_filter.id
            ))

        
    if opts.min_cov != None:  
        vcf_filter = vcf._Filter(
            id="mincov%d" % opts.min_cov, 
            desc="Minimum coverage")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

        filters.append((
            lambda s, f_id: f_id if s.INFO['DP'] < opts.min_cov else None,
            vcf_filter.id
            ))

    if opts.min_snv_phred != None:  
        vcf_filter = vcf._Filter(
            id="minqual%d" % opts.min_snv_phred, 
            desc="Minimum SNV quality")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

        filters.append((
            lambda s, f_id: f_id if s.QUAL != '.' and s.QUAL < opts.min_snv_phred else None,
            vcf_filter.id
            ))

    if opts.snv_fdr:
        vcf_filter = vcf._Filter(
            id=("fdr%f" % opts.snv_fdr).rstrip('0'),
            desc="SNV Phred FDR")
        vcf_reader.filters[vcf_filter.id] = vcf_filter# reader serves as template for writer

        # we make this only temporary part of the vcf records
        # and delete it later, so no need to add to header
        
        #vcf_info = vcf._Info(
        #    id="FDR", num=1, type='Integer',
        #    desc="Passed FDR test")        
        #vcf_reader.infos[vcf_info.id] = vcf_info
        vcf_info_id = "FDRPASS"
        pvals = []
        pidx = []
        for (i, s) in enumerate(snvs):
            if s.QUAL != '.':
                pvals.append(phredqual_to_prob(s.QUAL))
                pidx.append(i)
                s.INFO[vcf_info_id] = 0
            else:
                s.INFO[vcf_info_id] = '.'

        for i in fdr.fdr(pvals, opts.snv_fdr):
            snvs[pidx[i]].INFO[vcf_info_id] = 1
                
        filters.append((
            lambda s, f_id: f_id if s.INFO["FDRPASS"] != '.' and s.INFO["FDRPASS"] == 0 else None,
            vcf_filter.id
            ))                
        
        
    if len(filters) == 0:
        LOG.error("No filters activated. Will exit now.")
        sys.exit(1)

    # The actual filtering:
    #
    # FIXME can't this be done easier with map()?
    #
    for (filter_func, filter_id) in filters:
        for (i, s) in enumerate(snvs):
            f = filter_func(s, filter_id)
            if f:
                # just s = s.__replace() can't work
                if s.FILTER == '.' or s.FILTER == 'PASS':
                    snvs[i] = s._replace(FILTER=f)
                else:
                    snvs[i] = s._replace(FILTER="%s,%s" % (s.FILTER, f))

    #if opts.window_size != None:  
    #    raise NotImplementedError # FIXME


    # should all also work if we get already PASSed input  
    
    n_passed = 0
    for (i, s) in enumerate(snvs):
        if s.FILTER == '.':
            snvs[i] = s._replace(FILTER="PASS")
            n_passed += 1
    LOG.info("%d SNVs passed all filters." % n_passed)

    # remove temporary fdr markup
    for tmpkey in ['FDRPASS']:
        for s in snvs:
            if s.INFO.has_key(tmpkey):
                del s.INFO[tmpkey]

    if opts.pass_only:
        snvs = [s for s in snvs if s.FILTER == 'PASS']
        
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
    # FIXME need test cases
    LOG.info("Successful program exit")
