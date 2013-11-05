#!/usr/bin/env python
"""Cluster SNVs based on SNV freqs confidence interval
"""


__author__ = "Andreas Wilm, Niranjan Nagarajan"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2013 Genome Institute of Singapore"
__license__ = "GPL2"



# --- standard library imports
#
import sys
import logging
import os
import argparse
from math import sqrt


#--- third-party imports
#
# /

#--- project specific imports
#
# legacy snp format
HAVE_SNP_MODULE = False
try:
    from lofreq import snp
    HAVE_SNP_MODULE = True
except ImportError:
    pass
# vcf format
HAVE_VCF_MODULE = False
try:
    import lofreq2_local
    from lofreq_star import vcf
    HAVE_VCF_MODULE = True
except ImportError:
    pass    
if HAVE_SNP_MODULE == False and HAVE_VCF_MODULE == False:
    sys.stderr.write("Couldn't import any of LoFreq SNP format modules\n")
    sys.exit(1)
SUPPORTED_FORMATS = []
if HAVE_SNP_MODULE:
    SUPPORTED_FORMATS.append('snp')
if HAVE_VCF_MODULE:
    SUPPORTED_FORMATS.append('vcf')
    

    
__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2013 Genome Institute of Singapore"
__license__ = "GPL2"



#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')



# invocation of ipython on exceptions
#import sys, pdb
#from IPython.core import ultratb
#sys.excepthook = ultratb.FormattedTB(mode='Verbose',
#                                     color_scheme='Linux', call_pdb=1)

class MetaVar(object):
    """Wrapper for SNV and VCF format variants to look the same
    """
    
    def __init__(self, vcf_var=None, snp_var=None):
        """
        """

        if vcf_var and snp_var:
            raise ValueError, ("Can only take one: vcf- or snp-var")
        self.repr = None
        self.coverage = None
        self.freq = None
        self.var_count = None
        self.max_ci = None
        self.min_ci = None
        if vcf_var:
            self.add_vcf_var(vcf_var)
        elif snp_var:
            self.add_snp_var(snp_var)


    def add_vcf_var(self, vcf_var):
        """Add variant in vcf format
        """
        self.repr = "%s %d %c>%c %f" % (vcf_var.CHROM, vcf_var.POS,
                                        vcf_var.REF, ','.join(vcf_var.ALT),
                                        vcf_var.INFO['AF'])
        self.coverage = vcf_var.INFO['DP']
        self.freq = vcf_var.INFO['AF']
        self.var_count = int(self.coverage * self.freq)
        (self.min_ci, self.max_ci) = self.compute_ci(
            self.coverage, self.var_count)
        LOG.info('CI for %s: %f--%f' % (
            self.repr, self.min_ci, self.max_ci))

        
    def add_snp_var(self, snp_var):
        """Add variant in legacy SNP format
        """
        self.repr = "%s %d %c>%c %f" % (snp_var.chrom, snp_var.pos+1,
                                        snp_var.wildtype, snp_var.variant,
                                        snp_var.freq)
        self.coverage = int(snp_var.info['coverage'])
        self.freq = snp_var.freq
        self.var_count = int(self.coverage * self.freq)
        (self.min_ci, self.max_ci) = self.compute_ci(
            self.coverage, self.var_count)
        LOG.info('CI for %s: %f--%f' % (
            self.repr, self.min_ci, self.max_ci))

        
    @staticmethod
    def compute_ci(coverage, var_count):
        """Compute confidnce interval:
        
        Agresti-Coull Interval at the 0.05 level
        http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Agresti-Coull_Interval
        
        n~ = n + 4
        p~ = 1/n~ * (X + 4/2)
        ci: p~ +- 2*sqrt(1/n~ * p~ * (1-p~)
        """
        n_t = float(coverage + 4)
        p_t = (var_count + 2) / n_t
        ci = 2 * sqrt(p_t * (1-p_t) / n_t)
        min_ci = p_t - 3*ci
        if min_ci < 0.0:
            min_ci = 0.0
        max_ci = p_t + 3*ci

        return (min_ci, max_ci)
    

        

def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/dev/howto/argparse.html
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--verbose",
                      action="store_true", 
                      dest="verbose",
                      help="be verbose")
    parser.add_argument("--debug",
                      action="store_true", 
                      dest="debug",
                      help="enable debugging")
    parser.add_argument("-i", "--variants",
                      dest="var_file",
                      help="variant input file (supported formats: %s)" % (
                          ', '.join(SUPPORTED_FORMATS)))
    parser.add_argument("-o", "--out",
                      dest="cluster_file",
                      default="-",
                      help="Cluster output file (- for stdout = default)")

    return parser


    
def main():
    """The main function
    """

    parser = cmdline_parser()
    args = parser.parse_args()

    # FIXME catch unrecognized args (not just (len(args)

    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)


    for (in_file, descr) in [(args.var_file, "variant file")]:
        if not in_file:
            parser.error("%s input file argument missing." % descr)
            sys.exit(1)
        if not os.path.exists(in_file) and in_file != "-":
            sys.stderr.write(
                "file '%s' does not exist.\n" % in_file)
            sys.exit(1)
            
    for (out_file, descr) in [(args.cluster_file, "cluster output file")]:
        if not out_file:
            parser.error("%s output file argument missing." % descr)
            sys.exit(1)
        if os.path.exists(out_file) and out_file!="-":
            sys.stderr.write(
                "Cowardly refusing to overwrite existing"
                " output file '%s'.\n" % out_file)
            sys.exit(1)


    # A lot of code for supporting legacy SNP format. 
    #
    # FIXME this and MetaVar() should use vcf by default and just
    # convert snp to vcf
    #
    is_vcf = False
    if HAVE_VCF_MODULE:
        if args.var_file == '-':
            vcf_fh = sys.stdin
        else:
            vcf_fh = open(args.var_file)
            # FIXME gzip support
        vcf_reader = vcf.VCFReader(vcf_fh)
        try:
            var_list = [MetaVar(vcf_var=r)
                        for r in vcf_reader]
            is_vcf = True
        except:
            pass
        if vcf_fh != sys.stdin:
            vcf_fh.close()

    is_snp = False
    if not is_vcf and HAVE_SNP_MODULE:
        try:
            var_list = [MetaVar(snp_var=s) 
                        for s in snp.parse_snp_file(args.var_file)]
            is_snp = True
        except IndexError:
            pass

    if not is_snp and not is_vcf:
        LOG.error("Can't parse %s. Tried the following formats: %s" % (
            args.var_file, ', '.join(SUPPORTED_FORMATS)))
        sys.exit(1)

    
    LOG.info("Parsed %d SNPs from %s" % (len(var_list), args.var_file))

    
    var_list =  sorted(var_list, key=lambda x: x.freq, reverse=True)

    if args.cluster_file == '-':
        fh_out = sys.stdout
    else:
        fh_out = open(args.cluster_file, 'w')

        
    if len(var_list)==0:
        fh_out.write("No SNPs <-> no clusters!\n")
        if fh_out != sys.stdout:
            print "No SNPs <-> no clusters!"
            fh_out.close()
        sys.exit(0)

        
    cluster = dict()
    clu_no = 0
    seed = var_list[0]
    #cluster[clu_no,'members'] = ["%s %f" % (seed.repr, seed.freq)]
    cluster[clu_no,'members'] = ["%s" % (seed.repr)]
    cluster[clu_no,'min'] = seed.min_ci
    cluster[clu_no,'max'] = seed.max_ci

    for var in var_list[1:]:
        LOG.debug("checking %s %f: max_ci %f vvar. clu_min %f" % (
            var.repr, var.freq, var.max_ci, cluster[clu_no,'min']))
        if var.max_ci > cluster[clu_no,'min']:
            #cluster[clu_no,'members'].append("%s %f" % (var.repr, var.freq))
            cluster[clu_no,'members'].append("%s" % (var.repr))
        else:
            clu_no += 1
            seed = var
            #cluster[clu_no,'members'] = ["%s %f" % (seed.repr, seed.freq)]
            cluster[clu_no,'members'] = ["%s" % (seed.repr)]
            cluster[clu_no,'min'] = seed.min_ci
            cluster[clu_no,'max'] = seed.max_ci

        
    for i in range(clu_no+1):
        fh_out.write("cluster %d (freq. range: %f - %f): %s\n" % (
            i+1, cluster[i,'min'], cluster[i,'max'], 
            ', '.join(cluster[i,'members'])))
        
    if fh_out != sys.stdout:
        fh_out.close()
    print "%d clusters found (written to %s)" % (clu_no+1, fh_out.name)
 

if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
        
