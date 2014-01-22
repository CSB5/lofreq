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



def bam_index_exists(bam):
    """check if an index for an BAM file exists
    """
    for f in [bam + ".bai", os.path.splitext(bam)[0] + ".bai"]:
        if os.path.exists(f):
            return True
    return False



class SomaticSNVCaller(object):
    """Somatic SNV caller using LoFreq
    """

    VCF_NORMAL_RLX_EXT = "normal_relaxed.vcf"
    VCF_NORMAL_RLX_LOG_EXT = "normal_relaxed.log"
    VCF_TUMOR_MAPERRPROF_EXT = "maperrprof.txt"
    VCF_TUMOR_RLX_EXT = "tumor_relaxed.vcf"
    VCF_TUMOR_RLX_LOG_EXT = "tumor_relaxed.log"
    VCF_TUMOR_STR_EXT = "tumor_stringent.vcf"
    VCF_SOMATIC_RAW_EXT = "somatic_raw.vcf"
    VCF_SOMATIC_FINAL_EXT = "somatic_final.vcf"

    VCF_GERMLINE_EXT = "germline.vcf"

    LOFREQ = 'lofreq'

    DEFAULT_ALPHA_N = 0.01;# i.e. qual
    DEFAULT_ALPHA_T = 0.01;# i.e. qual
    DEFAULT_MTC_T = 'bonf'
    DEFAULT_MTC_ALPHA_T = 1
    DEFAULT_MQ_FILTER_T = 13
    DEFAULT_BAQ_OFF = False
    DEFAULT_MQ_OFF = False
    DEFAULT_SRC_QUAL_ON = False
    DEFAULT_SRC_QUAL_IGN_VCF = None
    DEFAULT_ALN_ERR_PROF_ON = False
    DEFAULT_MIN_COV = 10
    DEFAULT_USE_ORPHAN = False
    DEFAULT_NUM_THREADS = 1

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

        # setup output files
        #
        self.outfiles = []
        if reuse_normal_vcf:
            assert os.path.exists(reuse_normal_vcf)
            self.vcf_n_rlx = reuse_normal_vcf
        else:
            self.vcf_n_rlx = self.outprefix + self.VCF_NORMAL_RLX_EXT + ".gz"
            self.vcf_n_rlx_log = self.outprefix + self.VCF_NORMAL_RLX_LOG_EXT
            for f in [self.vcf_n_rlx, self.vcf_n_rlx_log]:
                assert not os.path.exists(f), (
                "Cowardly refusing to overwrite already existing file %s" % (f))

        self.vcf_t_maperrprof = self.outprefix + self.VCF_TUMOR_MAPERRPROF_EXT
        self.vcf_t_rlx = self.outprefix + self.VCF_TUMOR_RLX_EXT + ".gz"
        self.vcf_t_rlx_log = self.outprefix + self.VCF_TUMOR_RLX_LOG_EXT
        self.vcf_t_str = self.outprefix + self.VCF_TUMOR_STR_EXT + ".gz"
        self.vcf_som_raw = self.outprefix + self.VCF_SOMATIC_RAW_EXT + ".gz"
        self.vcf_som_fin = self.outprefix + self.VCF_SOMATIC_FINAL_EXT
        self.vcf_germl = self.outprefix + self.VCF_GERMLINE_EXT + ".gz"

        self.outfiles = [self.vcf_t_maperrprof, self.vcf_t_rlx,
                         self.vcf_t_str, self.vcf_som_raw,
                         self.vcf_som_fin, self.vcf_germl]
        # vcf_n already checked
        for f in self.outfiles:
            assert not os.path.exists(f), (
                "Cowardly refusing to overwrite already existing file %s" % f)

        # other params
        self.alpha_n = self.DEFAULT_ALPHA_N
        self.alpha_t = self.DEFAULT_ALPHA_T
        self.mtc_t = self.DEFAULT_MTC_T
        self.mtc_alpha_t = self.DEFAULT_MTC_ALPHA_T
        self.mq_filter_t = self.DEFAULT_MQ_FILTER_T
        self.baq_off = self.DEFAULT_BAQ_OFF
        self.mq_off = self.DEFAULT_MQ_OFF
        self.src_qual_on = self.DEFAULT_SRC_QUAL_ON
        self.aln_err_prof_on = self.DEFAULT_ALN_ERR_PROF_ON
        self.src_qual_ign_vcf = self.DEFAULT_SRC_QUAL_IGN_VCF
        self.min_cov = self.DEFAULT_MIN_COV
        self.use_orphan = self.DEFAULT_USE_ORPHAN
        self.num_threads = self.DEFAULT_NUM_THREADS



    @staticmethod
    def subprocess_wrapper(cmd, close_tmp=True):
        """Wrapper for subprocess.check_call

        Returns (rewound) fh for cmd stdout and stderr if close_tmp is
        False. Caller will then have to closer upon which the files
        will be deleted automatically.
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

        if os.path.exists(self.vcf_n_rlx):
            LOG.info('Reusing %s' % self.vcf_n_rlx)
            return


        # one could argue that BAQ (and MQ) should always be off for
        # normal, as it only every reduces the chance of calls. tests
        # on melanoma CHH930 50-50 g10ac045 orphan show that the ones
        # called extra, do not overlap with any of the somatic_raw
        # files. therefore we switched the behaviour back to the
        # previous default which was to use the same as in tumor

        # FIXME aln_err_prof?

        if self.num_threads < 2:
            cmd = [self.LOFREQ, 'call']
        else:
            cmd = [self.LOFREQ, 'call-parallel', '--pp-threads', "%d" % self.num_threads]
        cmd.extend(['-f', self.ref])
        if self.baq_off:
            cmd.append('-B')
        #if self.mq_off:
        cmd.append('-J')
        if self.use_orphan:
            cmd.append('--use-orphan')

        cmd.append('--verbose')
        if self.bed:
            cmd.extend(['-l', self.bed])
        cmd.append('--no-default-filter')# no filtering wanted
        cmd.extend(['-b', "%d" % 1, '-s', "%f" % self.alpha_n])
        cmd.extend(['-m', "%d" % self.mq_filter_t])
        cmd.extend(['-o', self.vcf_n_rlx])
        cmd.append(self.bam_n)

        # cmd = ['valgrind', '--tool=memcheck', '--leak-check=full'] + cmd

        (o, e) = self.subprocess_wrapper(cmd, close_tmp=False)
        fh = open(self.vcf_n_rlx_log, 'w')
        fh.write('# %s\n' % ' '.join(cmd))
        olines = o.readlines()
        elines = e.readlines()
        for l in elines:
            fh.write("stderr: %s" % l)
            LOG.info("cmd stderr: %s" % l.rstrip())
        for l in olines:
            fh.write("stdout: %s" % l)
        fh.close()
        o.close()
        e.close()


    def call_tumor(self):
        """Variant call on tumor sample
        """

        if self.aln_err_prof_on:
            cmd = [self.LOFREQ, 'bamstats']
            cmd.extend(['-f', self.ref])
            if self.bed:
                cmd.extend(['-l', self.bed])
            # -q ?
            cmd.extend(['-m', "%d" % self.mq_filter_t])
            cmd.extend(['-o' , self.vcf_t_maperrprof])
            cmd.append(self.bam_t)
            self.subprocess_wrapper(cmd)

        if self.num_threads < 2:
            cmd = [self.LOFREQ, 'call']
        else:
            cmd = [self.LOFREQ, 'call-parallel', '--pp-threads', "%d" % self.num_threads]
        cmd.extend(['-f', self.ref])
        if self.baq_off:
            cmd.append('-B')
        if self.mq_off:
            cmd.append('-J')
        if self.use_orphan:
            cmd.append('--use-orphan')

        cmd.append('--verbose')
        if self.bed:
            cmd.extend(['-l', self.bed])
        cmd.append('--no-default-filter')# filtering explicitely
        cmd.extend(['-b', "%d" % 1, '-s', "%f" % self.alpha_t])
        cmd.extend(['-m', "%d" % self.mq_filter_t])
        cmd.extend(['-o', self.vcf_t_rlx])

        # coverage is filtered later anyway, but ignoring it during call
        # makes things faster and avoids trouble if user forgots to give
        # bed-file etc
        if self.aln_err_prof_on:
            cmd.extend(['-A' , self.vcf_t_maperrprof])
        cmd.extend(['-C', "%d" % self.min_cov])
        if self.src_qual_on:
            cmd.append('-S')
        if self.src_qual_ign_vcf:
            cmd.extend(['-V', self.src_qual_ign_vcf])
        cmd.append(self.bam_t)

        #cmd = ['valgrind', '--tool=memcheck', '--leak-check=full'] + cmd

        (o, e) = self.subprocess_wrapper(cmd, close_tmp=False)
        fh = open(self.vcf_t_rlx_log, 'w')
        fh.write('# %s\n' % ' '.join(cmd))
        olines = o.readlines()
        elines = e.readlines()
        for l in elines:
            fh.write("stderr: %s" % l)
            LOG.info("cmd stderr: %s" % l.rstrip())
        for l in olines:
            fh.write("stdout: %s" % l)
        fh.close()
        o.close()
        e.close()

        num_tests = -1
        for l in elines:
            if l.startswith('Number of substitution tests performed'):
                num_tests = int(l.split(':')[1])
                break
        if num_tests == -1:
            LOG.error("Couldn't parse number of tests from lofreq call output"
                      " (which was: %s)" % (elines))
            raise ValueError

        cmd = [self.LOFREQ, 'filter', '-i', self.vcf_t_rlx,
               '--snv-qual', "%s" % self.mtc_t,
               '--snv-qual-alpha', '%f' % self.mtc_alpha_t,
               '--snv-qual-numtests', '%d' % num_tests,
               '--pass-only', '-o', self.vcf_t_str]

        self.subprocess_wrapper(cmd)


    def call_germline(self):
        """Call germline variants by taking the intersection between
        the stringent tumor and relaxed normal calls

        FIXME this is ad-hoc
        """

        cmd = [self.LOFREQ, 'vcfset', '-1', self.vcf_n_rlx,
               '-2', self.vcf_t_str,
               '-a', 'intersect',
               '-o', self.vcf_germl]
        # FIXME no further filtering and using vcf_n_rlx entries
        self.subprocess_wrapper(cmd)


    def complement(self):
        """Produce complement of tumor and normal variants and filter
        them

        Also adds SOMATIC tag
        """

        assert not os.path.exists(self.vcf_som_raw)
        vcf_som_raw_fh =  open(self.vcf_som_raw, 'w')

        cmd = [self.LOFREQ, 'vcfset', '-1', self.vcf_t_str,
               '-2', self.vcf_n_rlx,
               '-a', 'complement', '-o', '-']

        (cmd_o, cmd_e) = self.subprocess_wrapper(cmd, close_tmp=False)
        for l in cmd_e:
            LOG.warn("cmd stderr: %s" % l.rstrip())
        cmd_e.close()

        # quick and dirty vcf agnostic way of tagging them as SOMATIC
        for l in cmd_o:
            if not l.startswith('#'):
                l = l.rstrip() + ";SOMATIC"  + "\n"
            vcf_som_raw_fh.write("%s" % l)
        cmd_o.close()
        vcf_som_raw_fh.close()


    def uniq(self):
        """Run LoFreq uniq as final check on somatic variants
        """

        cmd = [self.LOFREQ, 'uniq',
               '--uni-freq', "0.5",
               '-v', self.vcf_som_raw,
               '-o', self.vcf_som_fin]
        if self.use_orphan:
            cmd.append('--use-orphan')
        cmd.append(self.bam_n)

        (o, e) = self.subprocess_wrapper(cmd, close_tmp=False)
        for l in e.readlines():
            LOG.warn("uniq stderr: %s" % l)
        o.close()
        e.close()


    def run(self):
        """Run the whole somatic SNV calling pipeline
        """

        if not bam_index_exists(self.bam_n):
            LOG.fatal("Normal BAM file is not indexed."
                      " Please create the index first"
                      " with e.g. samtools index %s" % (self.bam_n))
            return False


        for (k, v) in [(x, self.__getattribute__(x)) for x in dir(self)
                       if not x.startswith('_')]:
            if callable(v):
                continue
            LOG.debug("%s %s" % (k, v))
        #import pdb; pdb.set_trace()

        if self.src_qual_ign_vcf and not self.src_qual_on:
            LOG.fatal("ign-vcf file was provided, but src-qual is off")
            sys.exit(1)
        try:
            self.call_normal()
            self.call_tumor()
            self.complement()
            self.uniq()
            self.call_germline()
        except:
            return False

        # FIXME replace source line in final output with sys.argv?
        return True


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
                        help="Prefix for output files. Final somatic SNV"
                        " calls will be stored in PREFIX+%s." % (
                            SomaticSNVCaller.VCF_SOMATIC_FINAL_EXT))
    parser.add_argument("-f", "--ref",
                        required=True,
                        help="Reference fasta file")
    parser.add_argument("-l", "--bed",
                        help="BED file listing regions to restrict analysis to")


    default = 0.01
    parser.add_argument("--tumor-alpha",
                        #required=True,
                        default=default,
                        type=float,
                        help="Advanced: Significance threshold (alpha)"
                        " for SNV pvalues in (relaxed) tumor vcf"
                        " (default: %f)" % default)

    default = 0.01
    parser.add_argument("--normal-alpha",
                        #required=True,
                        default=default,
                        type=float,
                        help="Advanced: Significance threshold (alpha) for SNV pvalues"
                        "  in (relaxed) normal vcf (default: %f)" % default)

    default = 'bonf'
    choices = ['bonf', 'holm-bonf', 'fdr']
    parser.add_argument("--tumor-mtc",
                        #required=True,
                        default=default,
                        choices = choices,
                        help="Type of multiple testing correction for tumor"
                        " (default: %s)" % default)

    default = 10
    parser.add_argument("--tumor-mtc-alpha",
                        #required=True,
                        default=default,
                        type=float,
                        help="Multiple testing correction alpha for tumor"
                        " (default: %f)" % default)

    default = 13
    parser.add_argument("-m,", "--mq-filter",
                        type=int,
                        default=default,
                        help="Ignore reads with mapping quality below this value (default=%d)" % default)

    parser.add_argument("-B", "--baq-off",
                        action="store_true",
                        help="Disable BAQ computation in tumor")

    parser.add_argument("-J", "--mq-off",
                        action="store_true",
                        help="Disable use of mapping quality in LoFreq's model")

    parser.add_argument("-A", "--aln-err-prof",
                        action="store_true",
                        help="Advanced: Use alignment error profile")

    parser.add_argument("-S", "--src-qual",
                        action="store_true",
                        help="Advanced: Enable use of source quality in tumor (see also -V)")
    default = "normal"
    parser.add_argument("-V", "--ign-vcf",
                        default = default,
                        help="Ignore variants in this vcf-file for"
                        " source quality computation in tumor (only used if -S was given)."
                        " Use 'None' for NA and 'normal' for predictions"
                        " in normal sample (default = %s)" % default)

    parser.add_argument("--reuse-normal-vcf",
                        help="Advanced: reuse already computed"
                        " vcf for normal sample")

    parser.add_argument("--use-orphan",
                        action="store_true",
                        help="Advanced: use orphaned/anomalous reads from read pairs")

    parser.add_argument("--threads",
                        type=int,
                        dest="num_threads",
                        help="Advanced: use this many threads for each call")

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
        import pdb
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

    # check if outdir exists
    outdir = os.path.dirname(args.outprefix)
    if outdir != "" and not os.path.exists(outdir):
        LOG.error("The directory part of the given output prefix points"
                  " to a non-existing directory: '%s').\n" % (outdir))
        sys.exit(1)



    somatic_snv_caller = SomaticSNVCaller(
        bam_n = args.normal,
        bam_t = args.tumor,
        ref = args.ref,
        outprefix = args.outprefix,
        bed = args.bed,
        reuse_normal_vcf = args.reuse_normal_vcf)

    somatic_snv_caller.alpha_n = args.normal_alpha
    somatic_snv_caller.alpha_t = args.tumor_alpha
    somatic_snv_caller.mtc_t = args.tumor_mtc
    somatic_snv_caller.mtc_alpha_t = args.tumor_mtc_alpha
    somatic_snv_caller.num_threads = args.num_threads
    if args.mq_filter:
        somatic_snv_caller.mq_filter_t = args.mq_filter
    if args.baq_off:
        somatic_snv_caller.baq_off = True
        
    if args.src_qual:
        somatic_snv_caller.src_qual_on = True
        if args.ign_vcf:
            if args.ign_vcf.lower() not in ["none", "na"]:
                if args.ign_vcf == "normal":
                    somatic_snv_caller.src_qual_ign_vcf = somatic_snv_caller.vcf_n_rlx
                else:
                    somatic_snv_caller.src_qual_ign_vcf = args.ign_vcf
        
    if args.aln_err_prof:
        somatic_snv_caller.aln_err_prof_on = True
    if args.use_orphan:
        somatic_snv_caller.use_orphan = True

    if not somatic_snv_caller.run():
        LOG.fatal("Somatic SNV caller failed. Exiting")
        sys.exit(1)


if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
