#!/usr/bin/env python
"""Cluster SNVs based on SNV freqs confidence interval
"""


__author__ = "Andreas Wilm, Niranjan Nagarajan"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2013,2014 Genome Institute of Singapore"
__license__ = "The MIT License"



# --- standard library imports
#
import sys
import logging
import os
import argparse
from math import sqrt
from collections import namedtuple
from itertools import groupby

#--- third-party imports
#
# /

#--- project specific imports
#
# James Casbon's pyvcf
import vcf
    


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


CI = namedtuple('CI', ['min', 'max'])

# invocation of ipython on exceptions
#import sys, pdb
#from IPython.core import ultratb
#sys.excepthook = ultratb.FormattedTB(mode='Verbose',
#                                     color_scheme='Linux', call_pdb=1)

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

    return CI._make([min_ci, max_ci])

    
def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence

    Brent Pedersen: https://www.biostars.org/p/710/
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        #header = header.next()[1:].strip()
        header = header.next()[1:].strip().split(" ")[0]
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq
        

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
                      help="variant input file (vcf format)")
    parser.add_argument("-r", "--ref",
                      dest="reffa",
                        help="Reference fasta file (for reconstructing cluster sequence)")
    parser.add_argument("-o", "--out",
                      dest="cluster_file",
                      default="-",
                      help="Cluster output file (- for stdout = default)")

    return parser


def vcf_var_to_str(v):
    return "%s %d %s>%s %f" % (
        v.CHROM, v.POS, v.REF, ','.join(["%s" % x for x in v.ALT]), v.INFO['AF'])

    
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


    if args.cluster_file == '-':
        fh_out = sys.stdout
    else:
        fh_out = open(args.cluster_file, 'w')


    if args.reffa:
        refno = 0
        for refname, refseq in fasta_iter(args.reffa):
            if refno > 0:
                sys.stderr.write("Only supporting one sequence\n")
                sys.exit(1)
            refno += 1
    else:
        refseq = ""


    if args.var_file == '-':
        vcf_fh = sys.stdin
    else:
        vcf_fh = vcf.VCFReader(filename=args.var_file)
    var_list = [v for v in vcf_fh]

    if any([not v.is_snp for v in var_list]):
        sys.stderr.write("WARNING: Only supporting SNPs! Automatically removing others\n")
        var_list = [v for v in var_list if v.is_snp]
         
    LOG.info("Parsed %d SNPs from %s" % (len(var_list), args.var_file))
    
    

    assert all([x.INFO.has_key('AF') and x.INFO.has_key('DP')
                for x in var_list])
    var_list = sorted(var_list, key=lambda x: x.INFO['AF'], reverse=True)
    ci_list = [compute_ci(v.INFO['DP'], int(v.INFO['AF'] * v.INFO['DP'])) 
               for v in var_list]
    var_ci_list = zip(var_list, ci_list)
    del var_list, ci_list# paranoia


        
    if len(var_ci_list)==0:
        fh_out.write("No variants <-> no clusters!\n")
        if fh_out != sys.stdout:
            fh_out.close()
        sys.exit(0)

        
    cluster = dict()
    clu_no = 0
    seed_var, seed_ci = var_ci_list[0]
    #cluster[clu_no,'members'] = ["%s %f" % (seed.repr, seed.freq)]
    cluster[clu_no,'members'] = [seed_var]
    cluster[clu_no,'min'] = seed_ci.min
    cluster[clu_no,'max'] = seed_ci.max

    for var, ci in var_ci_list[1:]:
        LOG.debug("checking %s %f: max_ci %f vvar. clu_min %f" % (
            var, var.INFO['AF'], ci.max, cluster[clu_no,'min']))
        if ci.max > cluster[clu_no,'min']:
            #cluster[clu_no,'members'].append("%s %f" % (var.repr, var.freq))
            cluster[clu_no,'members'].append(var)
        else:
            clu_no += 1
            seed = var
            #cluster[clu_no,'members'] = ["%s %f" % (seed.repr, seed.freq)]
            cluster[clu_no,'members'] = [seed]
            cluster[clu_no,'min'] = ci.min
            cluster[clu_no,'max'] = ci.max

        
    for i in range(clu_no+1):
        fh_out.write("# cluster %d (freq. range: %f - %f): %s\n" % (
            i+1, cluster[i,'min'], cluster[i,'max'], 
            ', '.join([vcf_var_to_str(x) for x in cluster[i,'members']])))

        # write sequence as well if we have a reference
        if refseq:
            haplotype = refseq
            for v in sorted(cluster[i,'members'], key = lambda v: v.POS):
                # FIXME random order for multi-allelic psositions
                assert v.CHROM == refname
                assert refseq[v.POS-1] == v.REF# use refseq to not break for multi-allelic positions
                assert len(v.ALT)==1, ("Support for 1 base alt only")
                alt = str(v.ALT[0])
                idx = v.POS-1
                haplotype = haplotype[:idx] + alt + haplotype[idx + 1:]
            fh_out.write(">haplotype-cluster-{}\n{}\n".format(i+1, haplotype))

        
    if fh_out != sys.stdout:
        fh_out.close()
    print "%d clusters found (written to %s)" % (clu_no+1, fh_out.name)
 

if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
        
