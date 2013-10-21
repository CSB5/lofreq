#!/usr/bin/env python
"""LoFreq* Somatic SNV Caller.
"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2013 Genome Institute of Singapore"
__license__ = "GPL2"


#--- standard library imports
#
import sys
import logging
import os
import argparse
import subprocess
import tempfile

#--- third-party imports
#
#/

#--- project specific imports
#
# sets PATH so that local scripts/binary is used if present, i.e.
# stuff can be run without installing it
try:
    import lofreq2_local
except ImportError:
    pass



#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')



class SomaticSNVCaller(object):
    """Somatic SNV caller using LoFreq

    FIXME add gzip support
    """

    VCF_NORMAL_EXT = "normal.vcf"
    VCF_TUMOR_PREFILTER_EXT = "tumor-prefilter.vcf"
    VCF_TUMOR_EXT = "tumor.vcf"
    VCF_RAW_EXT = "somatic_raw.vcf"
    VCF_FILTERED_EXT = "somatic_filtered.vcf"
    VCF_FINAL_EXT = "somatic_final.vcf"

    LOFREQ = 'lofreq'

    DEFAULT_ALPHA_N = 0.001
    DEFAULT_ALPHA_T = 0.000001
    DEFAULT_CORR_T = 'fdr'
    DEFAULT_MQ_FILTER_T = 13
    DEFAULT_BAQ_OFF = False
    DEFAULT_SRC_QUAL_ON = False
    DEFAULT_SRC_QUAL_IGN_VCF = None

    MIN_COV = 10

    def __init__(self, bam_n=None, bam_t=None,
                 ref=None, outprefix=None,
                 bed=None, reuse_normal_vcf=None):
        """init function
        """

        assert all([bam_n, bam_t, ref, outprefix]), (
            "Missing mandatory arguments")

        # make sure infiles exist and save them
        #
        infiles = [bam_n, bam_t, ref]
        if bed:
            infiles.append(bed)
        for f in infiles:
            assert os.path.exists(f), (
                "File %s does not exist" % f)

        self.bam_n = bam_n
        self.bam_t = bam_t
        self.ref = ref
        self.bed = bed
        self.outprefix = outprefix

        self.outfiles = []

        self.vcf_n = None
        self.vcf_t_prefilter = None
        self.vcf_t = None
        self.vcf_som_raw = None
        self.vcf_som_filtered = None
        self.vcf_som_final = None

        self.set_output(reuse_normal_vcf)


        self.alpha_n = None
        self.alpha_t = None
        self.corr_t = None
        self.mq_filter_t = None
        self.baq_off = None
        self.src_qual_on = None
        self.src_qual_ign_vcf = None

        self.set_default_params()


    def set_output(self, reuse_normal_vcf):
        """setup output file names and make sure they don't exist
        """

        assert self.outprefix

        if reuse_normal_vcf:
            assert os.path.exists(reuse_normal_vcf)
            self.vcf_n = reuse_normal_vcf
        else:
            self.vcf_n = self.outprefix + self.VCF_NORMAL_EXT + ".gz"
            assert not os.path.exists(self.vcf_n), (
                "Cowardly refusing to overwrite already existing file %s" % (
                    self.vcf_n))

        self.vcf_t = self.outprefix + self.VCF_TUMOR_EXT + ".gz"
        self.vcf_t_prefilter = self.outprefix + self.VCF_TUMOR_PREFILTER_EXT + ".gz"
        self.vcf_som_raw = self.outprefix + self.VCF_RAW_EXT + ".gz"
        self.vcf_som_filtered = self.outprefix + self.VCF_FILTERED_EXT + ".gz"
        self.vcf_som_final = self.outprefix + self.VCF_FINAL_EXT

        self.outfiles = [self.vcf_t_prefilter, self.vcf_t, self.vcf_som_raw,
                         self.vcf_som_filtered, self.vcf_som_final]
        # vcf_n already checked
        for f in self.outfiles:
            assert not os.path.exists(f), (
                "Cowardly refusing to overwrite already existing file %s" % f)


    def set_default_params(self):
        """Set up default parameters
        """

        self.alpha_n = self.DEFAULT_ALPHA_N
        self.alpha_t = self.DEFAULT_ALPHA_T
        self.corr_t = self.DEFAULT_CORR_T
        self.mq_filter_t = self.DEFAULT_MQ_FILTER_T
        self.baq_off = self.DEFAULT_BAQ_OFF
        self.src_qual_on = self.DEFAULT_SRC_QUAL_ON



    @staticmethod
    def subprocess_wrapper(cmd, close_tmp=True):
        """Wrapper for subprocess.check_call

        Returns (rewound) fh for cmd stdout and stderr if close_tmp is
        False. Caller will then have to closer upon which the files
        will be deleted automaitcally.

        """

        assert isinstance(cmd, list)
        fh_stdout = tempfile.TemporaryFile()
        fh_stderr = tempfile.TemporaryFile()

        try:
            LOG.info("Executing %s", ' '.join(cmd))
            subprocess.check_call(cmd, stdout=fh_stdout, stderr=fh_stderr)
        except subprocess.CalledProcessError as e:
            LOG.fatal("The following command failed: %s (%s)" % (
                ' '.join(cmd), str(e)))
            LOG.fatal("An error message indicating the source of"
                      " this error should have bee printed above")
            raise
        except OSError as e:
            LOG.fatal("The following command failed: %s (%s)" % (
                ' '.join(cmd), str(e)))
            LOG.fatal("An error message indicating the source of"
                      " this error should have bee printed above")
            LOG.fatal("Looks like the lofreq binary is not in your PATH")
            raise

        if close_tmp:
            fh_stdout.close()
            fh_stderr.close()
            return (None, None)
        else:
            # will be destroyed upon closing, i.e. caller has to close!
            fh_stdout.seek(0)
            fh_stderr.seek(0)
            return (fh_stdout, fh_stderr)


    def call_normal(self):
        """Relaxed call of variants on normal sample
        """

        if os.path.exists(self.vcf_n):
            LOG.info('Reusing %s' % self.vcf_n)
            return

        cmd = [self.LOFREQ, 'call']
        cmd.extend(['-f', self.ref])
        #if self.baq_off:
        #    cmd.append('-B')
        cmd.append('--verbose')
        if self.bed:
            cmd.extend(['-l', self.bed])
        cmd.append('--no-default-filter')# no filtering wanted
        cmd.extend(['-b', "%d" % 1, '-s', "%f" % self.alpha_n])
        cmd.extend(['-m', "%d" % self.mq_filter_t])
        cmd.extend(['-o', self.vcf_n])
        cmd.append(self.bam_n)

        # cmd = ['valgrind', '--tool=memcheck', '--leak-check=full'] + cmd

        (o, e) = self.subprocess_wrapper(cmd, close_tmp=False)

        olines = o.readlines()
        elines = e.readlines()
        #for l in elines:
        #    print l
        o.close()
        e.close()


    def call_tumor(self):
        """Variant call on tumor sample
        """

        cmd = [self.LOFREQ, 'call']
        # coverage is filtered later anyway, but ignoring it during call
        # makes things faster and avoids trouble if user forgots to give
        # bed-file etc
        cmd.extend(['-C', "%d" % self.MIN_COV])
        cmd.extend(['-f', self.ref])
        if self.baq_off:
            cmd.append('-B')
        if self.src_qual_on:
            cmd.append('-S')
        if self.src_qual_ign_vcf:
            cmd.extend(['-V', self.src_qual_ign_vcf])
            
        cmd.append('--verbose')
        if self.bed:
            cmd.extend(['-l', self.bed])
        cmd.append('--no-default-filter')# filtering explicitely
        cmd.extend(['-b', "%d" % 1, '-s', "%f" % 0.001]) # should probably roughly match alpha_n
        cmd.extend(['-m', "%d" % self.mq_filter_t])
        cmd.extend(['-o', self.vcf_t_prefilter])
        cmd.append(self.bam_t)

        #cmd = ['valgrind', '--tool=memcheck', '--leak-check=full'] + cmd

        (o, e) = self.subprocess_wrapper(cmd, close_tmp=False)

        olines = o.readlines()
        elines = e.readlines()
        o.close()
        e.close()
        num_tests = -1
        for l in elines:
            #print l
            if l.startswith('Number of tests performed'):
                num_tests = int(l.split(':')[1])
                break
        if num_tests == -1:
            LOG.error("Couldn't parse number of tests from lofreq call output"
                      " (which was: %s)" % (elines))
            raise ValueError

        cmd = [self.LOFREQ, 'filter', '-i', self.vcf_t_prefilter,
               '--snv-qual', "%s" % self.corr_t,
               '--snv-qual-alpha', '%f' % self.alpha_t,
               '--snv-qual-numtests', '%d' % num_tests,
               '-o', self.vcf_t]

        self.subprocess_wrapper(cmd)


    def complement(self):
        """Produce complement of tumor and normal variants and filter
        them
        """

        cmd = [self.LOFREQ, 'vcfset', '-1', self.vcf_t,
               '-2', self.vcf_n,
               '-a', 'complement',
               '-o', self.vcf_som_raw]
        self.subprocess_wrapper(cmd)

        # apply filter to complement
        #
        cmd = [self.LOFREQ, 'filter', '-i', self.vcf_som_raw,
               '--min-cov', "%d" % self.MIN_COV,
               '--strandbias', 'holm-bonf',
               '-p', '-o', self.vcf_som_filtered]
        self.subprocess_wrapper(cmd)


    def uniq(self):
        """Run LoFreq uniq as final check on somatic variants
        """

        cmd = [self.LOFREQ, 'uniq',
               '--uni-freq', "0.5",
               '-v', self.vcf_som_filtered,
               '-o', self.vcf_som_final,
               self.bam_n]
        self.subprocess_wrapper(cmd)


    def run(self):
        """Run the whole somatic SNV calling pipeline
        """

        for (k, v) in [(x, self.__getattribute__(x)) for x in dir(self)
                       if not x.startswith('_')]:
            if callable(v):
                continue
            LOG.debug("%s %s" % (k, v))
        #import pdb; pdb.set_trace()

        if self.src_qual_ign_vcf and not self.src_qual_on:
            LOG.fatal("ign-vcf file was provided, but src-qual is off")
            sys.exit(1)
        self.call_normal()
        self.call_tumor()
        self.complement()
        self.uniq()

        # FIXME replace source line in final output with sys.argv?

        #cmd = ['gzip']
        #gzip_files = list(self.outfiles)
        #gzip_files.append(self.vcf_n) # not part of outfiles by default
        #gzip_files.remove(self.vcf_som_final)
        #for f in set(gzip_files):
        #    cmd.append(f)
        #self.subprocess_wrapper(cmd)




def cmdline_parser():
    """Returns an argparse instance
    """

    # http://docs.python.org/dev/howto/argparse.html
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Be verbose")
    parser.add_argument("--debug",
                        action="store_true",
                        help="Enable debugging")
    parser.add_argument("-n", "--normal",
                        required=True,
                        help="Normal BAM file")
    parser.add_argument("-t", "--tumor",
                        required=True,
                        help="Tumor BAM file")
    parser.add_argument("-o", "--outprefix",
                        help="Prefix for output files."
                        " Final somatic SNV calls will be stored in"
                        " PREFIX+%s." % (SomaticSNVCaller.VCF_FINAL_EXT))
    parser.add_argument("-f", "--ref",
                        required=True,
                        help="Reference fasta file")
    parser.add_argument("-l", "--bed",
                        help="BED file listing regions to restrict analysis to")

    default = 'fdr'
    choices = ['bonf', 'holm-bonf', 'fdr']
    parser.add_argument("--tumor-mtc",
                        #required=True,
                        default=default,
                        choices = choices,
                        help="Type of multiple testing correction on tumor"
                        " (default: %s)" % default)

    default = 0.01
    parser.add_argument("--tumor-alpha",
                        #required=True,
                        default=default,
                        type=float,
                        help="Significance threshold (alpha) for SNV pvalues"
                        " in tumor sample (default: %f)" % default)

    default = 0.001
    parser.add_argument("--normal-alpha",
                        #required=True,
                        default=default,
                        type=float,
                        help="Significance threshold (alpha) for SNV pvalues"
                        "  in normal sample (non-Bonferroni corrected)"
                        " (default: %f)" % default)

    default = 13
    parser.add_argument("-m,", "--mq-filter",
                        type=int,
                        default=default,
                        help="Mapping quality filter (default=%d)" % default)

    parser.add_argument("-B", "--baq-off",
                        action="store_true",
                        help="Disable BAQ computation")

    parser.add_argument("-S", "--src-qual",
                        action="store_true",
                        help="Enable use of source quality")

    parser.add_argument("-V", "--ign-vcf",
                        help="Ignore variants in this vcf file for source quality computation")

    parser.add_argument("--reuse-normal-vcf",
                        help="Expert only: reuse already computed"
                        " normal prediction")

    return parser



def main():
    """The main function
    """

    parser = cmdline_parser()
    args = parser.parse_args()

    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)
        import sys, pdb
        from IPython.core import ultratb
        sys.excepthook = ultratb.FormattedTB(mode='Verbose',
                                             color_scheme='Linux', call_pdb=1)


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

    if args.reuse_normal_vcf:
        if not os.path.exists(args.reuse_normal_vcf):
            LOG.error("file '%s' does not exist.\n" % (
                args.reuse_normal_vcf))
            sys.exit(1)

    LOG.debug("args = %s" % args)

    somatic_snv_caller = SomaticSNVCaller(
        bam_n = args.normal,
        bam_t = args.tumor,
        ref = args.ref,
        outprefix = args.outprefix,
        bed = args.bed,
        reuse_normal_vcf = args.reuse_normal_vcf)

    somatic_snv_caller.alpha_n = args.normal_alpha
    somatic_snv_caller.alpha_t = args.tumor_alpha
    somatic_snv_caller.corr_t = args.tumor_mtc
    if args.mq_filter:
        somatic_snv_caller.mq_filter_t = args.mq_filter
    if args.baq_off:
        somatic_snv_caller.baq_off = True
    if args.src_qual:
        somatic_snv_caller.src_qual_on = True
    if args.ign_vcf:
        somatic_snv_caller.src_qual_ign_vcf = args.ign_vcf

    somatic_snv_caller.run()


if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
