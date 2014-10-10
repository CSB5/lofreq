#!/usr/bin/env python
"""If you know about false negative somatic calls, find where they were lost along the way
"""

import sys
import argparse

import vcf



__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "GPL2"



def cmdline_parser():
    """
    creates an OptionParser instance
    """
    
    parser = argparse.ArgumentParser(description=__doc__)
    
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        dest="verbose",
                        help="be verbose")
    parser.add_argument("--fn",
                        required=True,
                        dest="vcf_fn",
                        help="FN vcf file")
    parser.add_argument("--n-rlx",
                        required=True,
                        dest="vcf_nrlx",
                        help="Normal relaxed vcf file")
    parser.add_argument("--n-str",
                        required=True,
                        dest="vcf_nstr",
                        help="Normal stringent vcf file")
    parser.add_argument("--t-rlx",
                        required=True,
                        dest="vcf_trlx",
                        help="Tumor relaxed vcf file")
    parser.add_argument("--t-str",
                        required=True,
                        dest="vcf_tstr",
                        help="Tumor stringent vcf file")
    parser.add_argument("--s-raw",
                        required=True,
                        dest="vcf_sraw",
                        help="Somatic raw vcf file")
    parser.add_argument("--s-final",
                        required=True,
                        dest="vcf_sfinal",
                        help="Somatic final vcf file")
    parser.add_argument("--s-final-wo-dbsnp",
                        required=True,
                        dest="vcf_sfinal_wo_dbsnp",
                        help="Somatic final vcf file without dbSNP")
    return parser



def main():
    """main function
    """
    
    vcf_fh = dict()
    #vcf_files = dict()

    parser = cmdline_parser()
    args = parser.parse_args()
        
    for (k, v) in [
            ('FN', args.vcf_fn),
            ('normal_rlx', args.vcf_nrlx),
            ('normal_str', args.vcf_nstr),
            ('tumor_rlx', args.vcf_trlx),
            ('tumor_str', args.vcf_tstr),
            ('somatic_raw', args.vcf_sraw),
            ('somatic_final', args.vcf_sfinal),
            ('somatic_final_minus_dbsnp', args.vcf_sfinal_wo_dbsnp)]:
        #vcf_files[k] = v
        try:
            vcf_fh[k] = vcf.VCFReader(filename=v)
        except:
            sys.stderr.write("Reading %s failed\n" % v)
            raise
    
    sys.stderr.write("Analyzing FN %s and friends\n" % vcf_fh['FN'].filename)
    
    ORDER = ['normal_rlx', 'normal_str', 'tumor_rlx', 'tumor_str', 'somatic_raw', 'somatic_final', 'somatic_final_minus_dbsnp']
    
    
    print "#CHROM\tPOS\tREF\tALT\t%s" % ('\t'.join(ORDER))
    for fn in vcf_fh['FN']:
        present_in = dict()
        for k in ORDER:
            present_in[k] = 0
            for t in vcf_fh[k].fetch(fn.CHROM, fn.POS-1, fn.POS):
                assert len(fn.REF) == len(t.REF)
                assert len(fn.ALT)==1
                assert len(t.ALT)==1            
                if t.ALT[0] == fn.ALT[0]:
                    if t.QUAL:
                        q = t.QUAL
                    else:
                        q = "."
                    try:
                        present_in[k] = "Q=%s;SB=%s;DP=%d;AF=%f" % (q, t.INFO['SB'], t.INFO['DP'], t.INFO['AF'])
                    except KeyError:
                        sys.stderr.write("Key Error. Dropping to debugger\n")
                        import pdb; pdb.set_trace()
                    break
        print "%s\t%s\t%s\t%s\t%s" % (
            fn.CHROM, fn.POS, fn.REF, fn.ALT[0], '\t'.join(["%s" % present_in[k] for k in ORDER]))
    
    
if __name__ == "__main__":
    main()    
