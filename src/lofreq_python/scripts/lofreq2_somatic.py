#!/usr/bin/env python
"""LoFreq* Somatic SNV Caller.
"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2013 Genome Institute of Singapore"
__license__ = "Free for non-commercial use"
#
# FIXME:update-copyright
#


#--- standard library imports
#
import sys
import logging
import os
import argparse
import subprocess
import datetime

#--- third-party imports
#
#/

#--- project specific imports
#
# sets PATH so that local scripts/binary is used if present, i.e.
# stuff can be run without installing it
try:
    import lofreq2_local
except:
    pass    

    
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

VCF_NORMAL_EXT = "normal.vcf"
VCF_TUMOR_EXT = "tumor.vcf"
VCF_RAW_EXT = "lofreq_somatic_raw.vcf"
VCF_FILTERED_EXT = "lofreq_somatic_filtered.vcf"
VCF_FINAL_EXT = "lofreq_somatic_final.vcf"



def timestamp():
    """Generate a timestamp string
    """

    return datetime.datetime.now().strftime("%Y%m%d-%H%M%S")


def cmdline_parser():
    """Returns an argparse instance
    """

    # http://docs.python.org/dev/howto/argparse.html
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-v", "--verbose", 
                        action="store_true",
                        help="be verbose")
    parser.add_argument("--debug", 
                        action="store_true",
                        help="enable debugging")
    parser.add_argument("-n", "--normal", 
                        required=True,
                        help="Normal BAM file")
    parser.add_argument("-t", "--tumor", 
                        required=True,
                        help="Tumor BAM file")
    parser.add_argument("-o", "--outprefix", 
                        help="Prefix for output files."
                        " Final somatic SNV calls will be stored in"
                        " PREFIX+%s. If empty"
                        " will be set to tumor BAM file." % VCF_FINAL_EXT)
    parser.add_argument("-f", "--ref", 
                        required=True,
                        help="Reference fasta file")
    parser.add_argument("-l", "--bed", 
                        help="BED file listing regions to restrict analysis to")

    default = 0.001
    parser.add_argument("-s", "--normal-sig", 
                        type=float,
                        default=default,
                        help="Significance threshold / evalue for"
                        " SNV prediction on normal sample"
                        " (default: %f)" % default)
    default = 10
    parser.add_argument("-S", "--tumor-sig", 
                        type=float,
                        default=default,
                        help="Significance threshold / evalue for"
                        " SNV prediction on tumor sample"
                        " (default: %f)" % default)
    default = 20
    parser.add_argument("-m,", "--mq-filter", 
                        type=int,
                        default=default,
                        help="mapping quality filter for"
                        " tumor sample (default=%d)" % default)
    parser.add_argument("--reuse-normal-vcf", 
                        help="Expert only: reuse already compute normal"
                        " prediction (PREFIX+%s)" % VCF_NORMAL_EXT)

    return parser



def somatic(bam_n, bam_t, ref, out_prefix, bed=None,
            sig_n=0.001, sig_t=10, mq_filter_n=20,
            reuse_normal=None):
    """Core of the somatic SNV callers, which calls all the necessary
    parts and glues them together
    """

    infiles = [bam_n, bam_t, ref]
    if bed:
        infiles.append(bed)
    for inf in infiles:
        assert os.path.exists(inf), ("File %s does not exist" % inf)

    vcf_n = out_prefix + VCF_NORMAL_EXT
    if reuse_normal:
        assert os.path.exists(reuse_normal), (
            "%s does not exist" % reuse_normal)
        vcf_n = reuse_normal
        
    vcf_t = out_prefix + VCF_TUMOR_EXT
    vcf_som_raw = out_prefix + VCF_RAW_EXT
    vcf_som_filtered = out_prefix + VCF_FILTERED_EXT
    vcf_som_final = out_prefix + VCF_FINAL_EXT
    
    outfiles = [vcf_t, vcf_som_raw, vcf_som_filtered, vcf_som_final]
    if not reuse_normal:
        outfiles.append(vcf_n)
    for outf in outfiles:
        assert not os.path.exists(outf), (
            "Cowardly refusing to overwrite already existing file %s" % outf)
            
    somatic_commands = []

    if not reuse_normal:
        cmd = ['lofreq', 'call', '-f', ref]
        if bed:
            cmd.extend(['-l', bed])
        cmd.extend(['-b', "%d" % 1, '-s', "%f" % sig_n, '-o', vcf_n, bam_n])
        somatic_commands.append(cmd)
    
    cmd = ['lofreq', 'call', '-f', ref,]
    if bed:
        cmd.extend(['-l', bed])
    cmd.extend(['-b', 'dynamic', '-s', "%f" % sig_t, 
                '-m', "%d" % mq_filter_n, '-o', vcf_t, bam_t])
    somatic_commands.append(cmd)

    # FIXME both of the above can run in theory be simultaneously
    
    cmd = ['lofreq', 'vcfset', '-1', vcf_t, '-2', vcf_n, 
           '-a', 'complement', '-o', vcf_som_raw]
    somatic_commands.append(cmd)
    
    cmd = ['lofreq', 'filter', '-i', vcf_som_raw, 
           '--min-cov', "%d" % 10, '--strandbias-holmbonf', 
           '-p', '-o', vcf_som_filtered]
    somatic_commands.append(cmd)

    # filter and uniq could be combined into one
    
    if False:
        cmd = ['lofreq', 'uniq', '-v', vcf_som_filtered,
               '-o', vcf_som_final, bam_n]
    else:
        LOG.warn("'uniq' is temporarily disabled")
        cmd = ['cp', vcf_som_filtered, vcf_som_final]
    somatic_commands.append(cmd)
    
    # FIXME gzip in and output should be supported internally 
    cmd = ['gzip', vcf_n, vcf_t, vcf_som_raw, vcf_som_filtered]
    somatic_commands.append(cmd)

    for (i, cmd) in enumerate(somatic_commands):
        LOG.info("Running cmd %d out of %d: %s" % (
            i+1, len(somatic_commands), ' '.join(cmd)))
        #LOG.critical("DEBUG continue before executing %s" % ' '.join(cmd)); continue

        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError as e:
            LOG.fatal("The following command failed: %s (%s)" % (
                ' '.join(cmd), str(e)))
            raise
        except OSError as e:
            LOG.fatal("The following command failed: %s (%s)" % (
                ' '.join(cmd), str(e)))
            LOG.fatal("Looks like the lofreq binary is not in your PATH (which is: %s)" % (os.environ["PATH"]))
            sys.exit(1)



def main():
    """The main function
    """
    
    parser = cmdline_parser()
    args = parser.parse_args()
    
    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)


    for (in_file, descr) in [(args.normal, "BAM file for normal tissue"),
                             (args.tumor, "BAM file for tumor tissue")]:
        if not in_file:
            LOG.error("%s input file argument missing." % descr)
            #parser.print_help()
            sys.exit(1)
        if not os.path.exists(in_file): # and in_file != "-":
            LOG.error("file '%s' does not exist.\n" % in_file)
            #parser.print_help()
            sys.exit(1)
            
    if args.outprefix:
        outprefix = args.outprefix
    else:
        outprefix = args.tumor.replace(".bam", "")

    if args.reuse_normal_vcf:
        if not os.path.exists(args.reuse_normal_vcf):
            LOG.error("file '%s' does not exist.\n" % (args.reuse_normal_vcf))
            sys.exit(1)

    LOG.debug("args = %s" % args)

    somatic(args.normal, args.tumor, args.ref, outprefix, args.bed,
            args.normal_sig, args.tumor_sig, args.mq_filter, 
            args.reuse_normal_vcf)


    
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
