#!/usr/bin/env python
"""LoFreq* Somatic SNV Caller: Predict somatic variants from a paired
normal/disease sample.

The script will produce several output files using the prefix specified.
"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2013,2014 Genome Institute of Singapore"
__license__ = "The MIT License"


#--- standard library imports
#
import sys
import logging
import os
import argparse
import subprocess
import tempfile
from socket import gethostname

#--- third-party imports
#
#/

#--- project specific imports
#
# sets PATH so that local scripts/binary is used if present, i.e.
# stuff can be run without installing it
#
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

    VCF_NORMAL_RLX_EXT = "normal_relaxed.vcf.gz"
    VCF_NORMAL_RLX_LOG_EXT = "normal_relaxed.log"
    VCF_NORMAL_STR_EXT = "normal_stringent.snvs.vcf.gz"
    VCF_INDELS_NORMAL_STR_EXT = "normal_stringent.indels.vcf.gz"
    #
    VCF_TUMOR_RLX_EXT = "tumor_relaxed.vcf.gz"
    VCF_TUMOR_RLX_LOG_EXT = "tumor_relaxed.log"
    VCF_TUMOR_STR_EXT = "tumor_stringent.snvs.vcf.gz"
    VCF_INDELS_TUMOR_STR_EXT = "tumor_stringent.indels.vcf.gz"
    #
    VCF_SOMATIC_RAW_EXT = "somatic_raw.snvs.vcf.gz"
    VCF_INDELS_SOMATIC_RAW_EXT = "somatic_raw.indels.vcf.gz"
    VCF_SOMATIC_FINAL_EXT = "somatic_final.snvs.vcf.gz"
    VCF_INDELS_SOMATIC_FINAL_EXT = "somatic_final.indels.vcf.gz"
    VCF_SOMATIC_FINAL_WO_DBSNP_EXT = "somatic_final_minus-dbsnp.snvs.vcf.gz"
    VCF_INDELS_SOMATIC_FINAL_WO_DBSNP_EXT = "somatic_final_minus-dbsnp.indels.vcf.gz"
    #
    VCF_GERMLINE_EXT = "germline.snvs.vcf.gz"
    VCF_GERMLINE_INDELS_EXT = "germline.indels.vcf.gz"

    LOFREQ = 'lofreq'

    # call parameters for relaxed calls in normal and tumor
    DEFAULT_ALPHA_N = 0.10
    DEFAULT_ALPHA_T = 0.01
    # tumor only
    DEFAULT_SRC_QUAL_ON = True
    DEFAULT_SRC_QUAL_IGN_VCF = None
    DEFAULT_MIN_COV = 7
    DEFAULT_USE_ORPHAN = False# always on for normal
    DEFAULT_BAQ_OFF = False# always off for normal

    # stringent parameters for tumor
    DEFAULT_MTC_T = 'bonf'
    DEFAULT_MTC_ALPHA_T = 1
    DEFAULT_INDEL_MTC_T = 'bonf'
    DEFAULT_INDEL_MTC_ALPHA_T = 0.01# conservative value reduces dep on dbsnp

    # stringent parameters for normal (only used for sq)
    DEFAULT_MTC_N = 'fdr'
    DEFAULT_MTC_ALPHA_N = 0.01

    # uniq parameters
    DEFAULT_SNV_UNIQ_MTC = 'fdr'
    DEFAULT_SNV_UNIQ_MTC_ALPHA = 0.001
    DEFAULT_INDEL_UNIQ_MTC = 'fdr'
    DEFAULT_INDEL_UNIQ_MTC_ALPHA = 0.0001

    # misc
    DEFAULT_CALL_INDELS = False
    DEFAULT_NUM_THREADS = 1
    DEFAULT_DO_GERMLINE = False
    DEFAULT_SB_MTC_ALPHA = 0.001
    DEFAULT_MAX_COV = 100000

    def __init__(self, bam_n, bam_t, ref, outprefix,
                 bed=None, dbsnp=None, continue_interrupted=False):
        """init function
        """

        assert all([bam_n, bam_t, ref, outprefix]), (
            "Missing mandatory arguments")

        # make sure infiles exist
        #
        infiles = [bam_n, bam_t, ref]
        if bed:
            infiles.append(bed)
        if dbsnp:
            infiles.append(dbsnp)
        for f in infiles:
            assert os.path.exists(f), (
                "File %s does not exist" % f)

        self.bam_n = bam_n
        self.bam_t = bam_t
        self.ref = ref
        self.bed = bed
        self.dbsnp = dbsnp
        self.outprefix = outprefix

        # continue interrupted program. use with caution. existing
        # files will be treated as coming from a successfully
        # completed process run with the same options as here
        self.continue_interrupted = continue_interrupted

        # setup output files
        #
        self.vcf_n_rlx = self.outprefix + self.VCF_NORMAL_RLX_EXT
        self.vcf_n_rlx_log = self.outprefix + self.VCF_NORMAL_RLX_LOG_EXT
        self.vcf_n_str = self.outprefix + self.VCF_NORMAL_STR_EXT
        self.vcf_indels_n_str = self.outprefix + self.VCF_INDELS_NORMAL_STR_EXT
        #
        self.vcf_t_rlx = self.outprefix + self.VCF_TUMOR_RLX_EXT
        self.vcf_t_rlx_log = self.outprefix + self.VCF_TUMOR_RLX_LOG_EXT
        self.vcf_t_str = self.outprefix + self.VCF_TUMOR_STR_EXT
        self.vcf_indels_t_str = self.outprefix + self.VCF_INDELS_TUMOR_STR_EXT
        #
        self.vcf_som_raw = self.outprefix + self.VCF_SOMATIC_RAW_EXT
        self.vcf_indels_som_raw = self.outprefix + self.VCF_INDELS_SOMATIC_RAW_EXT
        self.vcf_som_fin = self.outprefix + self.VCF_SOMATIC_FINAL_EXT
        self.vcf_indels_som_fin = self.outprefix + self.VCF_INDELS_SOMATIC_FINAL_EXT
        self.vcf_som_fin_wo_dbsnp = self.outprefix + self.VCF_SOMATIC_FINAL_WO_DBSNP_EXT
        self.vcf_indels_som_fin_wo_dbsnp = self.outprefix + self.VCF_INDELS_SOMATIC_FINAL_WO_DBSNP_EXT
        #
        self.vcf_germl = self.outprefix + self.VCF_GERMLINE_EXT
        self.vcf_germl_indels = self.outprefix + self.VCF_GERMLINE_INDELS_EXT

        self.call_rlx_extra_args = None

        # make sure output files don't exist if we are not in
        # 'continue' mode
        #
        self.outfiles = []
        self.outfiles = [self.vcf_n_rlx, self.vcf_n_rlx_log, self.vcf_n_str, self.vcf_indels_n_str,
                         self.vcf_t_rlx, self.vcf_t_rlx_log, self.vcf_t_str, self.vcf_indels_t_str,
                         self.vcf_som_raw, self.vcf_som_fin,
                         self.vcf_indels_som_raw, self.vcf_indels_som_fin,
                         self.vcf_som_fin_wo_dbsnp, self.vcf_indels_som_fin_wo_dbsnp,
                         self.vcf_germl, self.vcf_germl_indels]
        if not self.continue_interrupted:
            for f in self.outfiles:
                assert not os.path.exists(f), (
                    "Cowardly refusing to overwrite already existing file %s" % f)

        # other params
        self.alpha_n = self.DEFAULT_ALPHA_N
        self.alpha_t = self.DEFAULT_ALPHA_T

        self.mtc_t = self.DEFAULT_MTC_T
        self.mtc_alpha_t = self.DEFAULT_MTC_ALPHA_T
        self.indel_mtc_t = self.DEFAULT_MTC_T
        self.indel_mtc_alpha_t = self.DEFAULT_MTC_ALPHA_T

        # stringent normal (SQ ign. only)
        self.mtc_n = self.DEFAULT_MTC_N
        self.mtc_alpha_n = self.DEFAULT_MTC_ALPHA_N

        self.snv_uniq_mtc = self.DEFAULT_SNV_UNIQ_MTC
        self.snv_uniq_mtc_alpha = self.DEFAULT_SNV_UNIQ_MTC_ALPHA
        self.indel_uniq_mtc = self.DEFAULT_INDEL_UNIQ_MTC
        self.indel_uniq_mtc_alpha = self.DEFAULT_INDEL_UNIQ_MTC_ALPHA

        self.src_qual_on = self.DEFAULT_SRC_QUAL_ON
        self.src_qual_ign_vcf = self.DEFAULT_SRC_QUAL_IGN_VCF
        self.min_cov = self.DEFAULT_MIN_COV
        self.use_orphan = self.DEFAULT_USE_ORPHAN
        self.baq_off = self.DEFAULT_BAQ_OFF
        self.num_threads = self.DEFAULT_NUM_THREADS
        self.call_indels = self.DEFAULT_CALL_INDELS
        self.do_germline = self.DEFAULT_DO_GERMLINE
        self.sb_mtc_alpha = self.DEFAULT_SB_MTC_ALPHA
        self.max_cov = self.DEFAULT_MAX_COV


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
            LOG.fatal("The following command failed with code %d: %s" % (
                e.returncode, ' '.join(cmd)))
            try:
                fh_stderr.seek(0)
                LOG.fatal("Received the following on stderr:")
                for line in fh_stderr:
                    sys.stderr.write(line + "\n")
            except:
                pass
            raise
        except OSError as e:
            LOG.fatal("The following command failed: %s (%s)" % (
                ' '.join(cmd), str(e)))
            LOG.fatal("Maybe the lofreq binary is not in your PATH")
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


    @staticmethod
    def num_tests_from_log(stream):
        """Extract number of performed SNV and indel tests from log file"""

        num_subst_tests = -1
        num_indel_tests = -1
        for l in stream:
            if l.startswith('Number of substitution tests performed'):
                num_subst_tests = int(l.split(':')[1])
            elif l.startswith('Number of indel tests performed'):
                num_indel_tests = int(l.split(':')[1])
            if num_subst_tests != -1 and num_indel_tests != -1:
                break
        if num_subst_tests == -1 and num_indel_tests == -1:
            LOG.error("Couldn't parse number of tests from reused log")
            raise ValueError
        return (num_subst_tests, num_indel_tests)


    def call_rlx(self, sample_type):
        """Relaxed calling of variants in normal or tumor. Calls indels and
        substitutions! Can be prevented by setting call_rlx_extra_args
        accordingly.

        """

        assert sample_type in ['normal', 'tumor']

        # shared arguments for both sample types
        #
        if self.num_threads < 2:
            cmd = [self.LOFREQ, 'call']
        else:
            cmd = [self.LOFREQ, 'call-parallel',
                   '--pp-threads', "%d" % self.num_threads]
        cmd.extend(['-d', "%d" % int(self.max_cov*1.01)])
        cmd.extend(['-f', self.ref])
        cmd.append('--verbose')
        cmd.append('--no-default-filter')# we filter later explicitely
        cmd.extend(['-b', "%d" % 1])# bonferroni factor 1
        if self.bed:
            cmd.extend(['-l', self.bed])

        if self.call_indels:
            cmd.append("--call-indels")
        if self.call_rlx_extra_args:
            cmd.extend(self.call_rlx_extra_args)

        # sample type specific arguments
        #
        if sample_type == "normal":
            cmd.extend(['-a', "%f" % self.alpha_n])
            cmd.append('--use-orphan')
            cmd.append('-B')# BAQ off
            cmd.append('-N')# MQ off
            cmd.append('-A')# IDAQ off

            out_vcf = self.vcf_n_rlx
            out_log = self.vcf_n_rlx_log
            cmd.append(self.bam_n)

        elif sample_type == "tumor":
            cmd.extend(['-a', "%f" % self.alpha_t])
            if self.use_orphan:
                cmd.append('--use-orphan')
            if self.baq_off:
                cmd.append('-B')
            cmd.extend(['-C', "%d" % self.min_cov])
            if self.src_qual_on:
                cmd.append('-s')
            if self.src_qual_ign_vcf:
                cmd.extend(['-S', self.src_qual_ign_vcf])
            out_vcf = self.vcf_t_rlx
            out_log = self.vcf_t_rlx_log
            cmd.append(self.bam_t)

        else:
            raise ValueError(sample_type)

        # out_vcf now set
        cmd.extend(['-o', out_vcf])


        # before we actually do anything check existance of output
        # files and whether we should reuse them
        #
        if self.continue_interrupted:
            if os.path.exists(out_vcf):
                assert os.path.exists(out_log), (
                    "%s exists but %s is missing." % (out_vcf, out_log))
                LOG.info("Skipping rlx call on %s" % sample_type)

                LOG.info("Parsing number of tests from log file %s" % out_log)
                fh = open(out_log, 'r')
                elines = [l.decode().replace("stderr: ", "") for l in fh.readlines()]
                fh.close()
                return self.num_tests_from_log(elines)
            else:
                assert not os.path.exists(out_log)

        (o, e) = self.subprocess_wrapper(cmd, close_tmp=False)
        fh = open(out_log, 'w')
        fh.write('# %s\n' % ' '.join(cmd))
        olines = [l.decode() for l in o.readlines()]
        elines = [l.decode() for l in e.readlines()]
        for l in elines:
            fh.write("stderr: %s" % l)
            LOG.info("cmd stderr: %s" % l.rstrip())
        for l in olines:
            fh.write("stdout: %s" % l)
        fh.close()
        o.close()
        e.close()

        return self.num_tests_from_log(elines)


    def rlx_to_str(self, sample_type, num_tests):
        """Using tumor filtering settings to create stringent calls
        from relaxed calls
        """

        num_snv_tests, num_indel_tests = num_tests 
        assert sample_type in ['normal', 'tumor']

        # filtering stringently using tumor stringent settings
        if sample_type == "normal":
            vcf_rlx = self.vcf_n_rlx
            vcf_str = self.vcf_n_str
            vcf_indels_str = self.vcf_indels_n_str

            mtc = self.mtc_n
            mtc_alpha = self.mtc_alpha_n
            indel_mtc = mtc
            indel_mtc_alpha = mtc_alpha

        elif sample_type == "tumor":
            vcf_rlx = self.vcf_t_rlx
            vcf_str = self.vcf_t_str
            vcf_indels_str = self.vcf_indels_t_str

            mtc = self.mtc_t
            mtc_alpha = self.mtc_alpha_t
            indel_mtc = self.indel_mtc_t
            indel_mtc_alpha = self.indel_mtc_alpha_t
        else:
            raise ValueError(sample_type)

        # filter indels and snvs separately
        #
        filter_base_cmd = [
            self.LOFREQ, 'filter', '-i', vcf_rlx,
            '--sb-mtc', 'fdr', '--sb-alpha', '%f' % self.sb_mtc_alpha,
            '--cov-max', "%d" % self.max_cov,
            '--cov-min', '%d' % self.min_cov]
        filter_snv_cmd = filter_base_cmd + [
            '--only-snvs',
            '--snvqual-mtc', "%s" % mtc,
            '--snvqual-alpha', '%f' % mtc_alpha,
            '--snvqual-ntests', '%d' % num_snv_tests]
        filter_indel_cmd = filter_base_cmd + [
            '--only-indels',
            '--indelqual-mtc', "%s" % indel_mtc,
            '--indelqual-alpha', '%f' % indel_mtc_alpha,
            '--indelqual-ntests', '%d' % num_indel_tests]

        for (vcf_out, cmd) in [(vcf_str, filter_snv_cmd),
                               (vcf_indels_str, filter_indel_cmd)]:
            if self.continue_interrupted and os.path.exists(vcf_out):
                LOG.info('Reusing %s' % (vcf_out))
            else:
                cmd = cmd + ["-o", vcf_out]
                self.subprocess_wrapper(cmd)


    def call_germline(self):
        """Call germline variants by taking the intersection between
        the stringent tumor and relaxed normal calls

        WARNING this is ad-hoc. There is no further downstream
        filtering and we're using the meta-info from the vcf_n_rlx
        entries.
        """

        cmd = [self.LOFREQ, 'vcfset',
               '-a', 'intersect',
               '-1', self.vcf_n_rlx, '-2', self.vcf_t_str,
               '-o', self.vcf_germl]
        cmd = [self.LOFREQ, 'vcfset',
               '-a', 'intersect',
               '-1', self.vcf_n_rlx, '-2', self.vcf_indels_t_str,
               '-o', self.vcf_germl_indels]
        self.subprocess_wrapper(cmd)


    def remove_normal(self):
        """Produce complement of tumor and normal variants and add SOMATIC tag
        """

        vcfset_base_cmd = [self.LOFREQ, 'vcfset', '-a', 'complement',
                            '-2', self.vcf_n_rlx, '--add-info', 'SOMATIC']
        vcfset_snv_cmd = vcfset_base_cmd + [
            '--only-snvs', '-1', self.vcf_t_str, '-o', self.vcf_som_raw]
        vcfset_indels_cmd = vcfset_base_cmd + [
            '--only-indels', '--only-pos', '-1', self.vcf_indels_t_str,
            '-o', self.vcf_indels_som_raw]

        for (vcf_out, cmd) in [(self.vcf_som_raw, vcfset_snv_cmd),
                               (self.vcf_indels_som_raw, vcfset_indels_cmd)]:
            if self.continue_interrupted and os.path.exists(vcf_out):
                LOG.info('Reusing %s' % self.vcf_som_raw)
                continue
            else:
                assert not os.path.exists(vcf_out), (
                    "%s already exists. Please remove or run me with"
                    " --continue if you want to reuse this file" % vcf_out)
            self.subprocess_wrapper(cmd)


    def uniq(self):
        """Run LoFreq uniq as final check on somatic variants
        """

        uniq_base_cmd = [self.LOFREQ, 'uniq', '--uni-freq', "0.5", "--is-somatic"]


        uniq_snv_cmd = uniq_base_cmd + [
            "-v", self.vcf_som_raw,
            '--uniq-mtc', self.snv_uniq_mtc,
            '--uniq-alpha', "%s" % self.snv_uniq_mtc_alpha]
        uniq_indels_cmd = uniq_base_cmd + [
            "-v", self.vcf_indels_som_raw,
            '--uniq-mtc', self.indel_uniq_mtc,
            '--uniq-alpha', "%s" % self.indel_uniq_mtc_alpha]

        for (vcf_out, cmd) in [(self.vcf_som_fin, uniq_snv_cmd),
                               (self.vcf_indels_som_fin, uniq_indels_cmd)]:

            if self.continue_interrupted and os.path.exists(vcf_out):
                LOG.info('Reusing %s' % vcf_out)
                continue
            else:
                assert not os.path.exists(vcf_out), (
                    "%s already exists. Please remove or run me with"
                    " --continue if you want to reuse this file" % vcf_out)

            cmd.extend(['-o', vcf_out])
            cmd.append(self.bam_n)

            (o, e) = self.subprocess_wrapper(cmd, close_tmp=False)
            for l in e.readlines():
                LOG.warn("uniq stderr: %s" % l.decode())
            o.close()
            e.close()


    def remove_dbsnp(self):
        """Remove dbSNP from 'final' somatic calls
        """

        complement_base_cmd = [self.LOFREQ, 'vcfset',
                               '-a', 'complement',
                               '-2', self.dbsnp]
        complement_snv_cmd = complement_base_cmd + [
            '-1', self.vcf_som_fin, '--only-snvs']
        complement_indels_cmd = complement_base_cmd + [
            '-1', self.vcf_indels_som_fin, '--only-pos', '--only-indels']

        for (vcf_out, cmd) in [(self.vcf_som_fin_wo_dbsnp, complement_snv_cmd),
                               (self.vcf_indels_som_fin_wo_dbsnp, complement_indels_cmd)]:

            if self.continue_interrupted and os.path.exists(vcf_out):
                LOG.info('Reusing %s' % vcf_out)
                return
            else:
                assert not os.path.exists(vcf_out), (
                    "%s already exists. Please remove or"
                    " run me with --continue if you want to reuse this file" % vcf_out)

            cmd.extend(["-o", vcf_out])
            self.subprocess_wrapper(cmd)


    def run(self):
        """Run the whole somatic SNV calling pipeline

        Will raise an exception on error
        """

        LOG.info("Running on %s" % gethostname())

        # sanity checks
        #
        for b in [self.bam_n, self.bam_t]:
            if not bam_index_exists(b):
                LOG.fatal("BAM file %s is not indexed."
                          " Please create the index first"
                          " with e.g. samtools index (or use lofreq)" % (b))
                return ValueError

        if self.src_qual_ign_vcf and not self.src_qual_on:
            LOG.fatal("ign-vcf file was provided, but src-qual is off")
            return ValueError


        for (k, v) in [(x, self.__getattribute__(x)) for x in dir(self)
                       if not x.startswith('_')]:
            if callable(v):
                continue
            LOG.debug("%s %s" % (k, v))
        #import pdb; pdb.set_trace()


        try:
            (num_subst_tests, num_indel_tests) = self.call_rlx("normal")
            self.rlx_to_str("normal", (num_subst_tests, num_indel_tests))

            (num_subst_tests, num_indel_tests) = self.call_rlx("tumor")
            self.rlx_to_str("tumor", (num_subst_tests, num_indel_tests))
        except:
            #return False
            raise

        self.remove_normal()
        self.uniq()
        if self.dbsnp:
            self.remove_dbsnp()

        if self.do_germline:
            self.call_germline()

        # FIXME add source line (sys.argv) in final outputs



def cmdline_parser():
    """Returns an argparse instance
    """

    # http://docs.python.org/dev/howto/argparse.html
    parser = argparse.ArgumentParser(prog="lofreq somatic",
                                     description=__doc__)

    basic = parser.add_argument_group('Basic Options')

    basic.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Be verbose")
    basic.add_argument("-n", "--normal",
                        required=True,
                        help="Normal BAM file")
    basic.add_argument("-t", "--tumor",
                        required=True,
                        help="Tumor BAM file")
    basic.add_argument("-o", "--outprefix",
                        required=True,
                        help="Prefix for output files")
    basic.add_argument("-f", "--ref",
                        required=True,
                        help="Reference fasta file")
    basic.add_argument("-l", "--bed",
                        help="BED file listing regions to restrict analysis to")
    basic.add_argument("-d", "--dbsnp",
                        help="vcf-file (bgzipped and index with tabix)"
                       " containing known germline variants (e.g. dbsnp for human")

    default = SomaticSNVCaller.DEFAULT_NUM_THREADS
    basic.add_argument("--threads",
                        type=int,
                        default=default,
                        dest="num_threads",
                        help="Use this many threads for each call")

    ###

    advanced = parser.add_argument_group('Advanced Options (PLEASE read the documentation before changing any of these)')

    default = SomaticSNVCaller.DEFAULT_MTC_T
    choices = ['bonf', 'holm-bonf', 'fdr']
    advanced.add_argument("--tumor-mtc",
                        #required=True,
                        default=default,
                        choices=choices,
                        help="Type of multiple testing correction for tumor"
                        " (default: %s)" % default)

    default = SomaticSNVCaller.DEFAULT_MTC_ALPHA_T
    advanced.add_argument("--tumor-mtc-alpha",
                        #required=True,
                        default=default,
                        type=float,
                        help="Multiple testing correction alpha for tumor"
                        " (default: %f)" % default)

    default = SomaticSNVCaller.DEFAULT_INDEL_MTC_T
    choices = ['bonf', 'holm-bonf', 'fdr']
    advanced.add_argument("--indel-tumor-mtc",
                        #required=True,
                        default=default,
                        choices=choices,
                        help="Type of multiple testing correction for tumor"
                        " (default: %s)" % default)

    default = SomaticSNVCaller.DEFAULT_INDEL_MTC_ALPHA_T
    advanced.add_argument("--indel-tumor-mtc-alpha",
                        #required=True,
                        default=default,
                        type=float,
                        help="Multiple testing correction alpha for tumor"
                        " (default: %f)" % default)

    advanced.add_argument("--call-indels",
                        action="store_true",
                        help="Also call indels (see documentation  on how to preprocess your BAM files)")


    default = SomaticSNVCaller.DEFAULT_MIN_COV
    advanced.add_argument("--min-cov",
                        type=int,
                        default=default,
                        help="Minimum coverage for somatic calls"
                        " (default: %d)" % default)

    advanced.add_argument("--germline",
                        action="store_true",
                        help="Also list germline calls in separate file")


    ###

    experts = parser.add_argument_group('Experts (PLEASE do not use/change these, unless you know exactly what you are doing and'
                                        ' if you change them nevertheless, light a candle first)')

    default = SomaticSNVCaller.DEFAULT_ALPHA_N
    experts.add_argument("--normal-alpha",
                        #required=True,
                        default=default,
                        type=float,
                        help=argparse.SUPPRESS,
                        #help="Significance threshold (alpha) for SNV pvalues"
                        #"  in (relaxed) normal vcf"
                        #" (default: %f)" % default
                     )
    default = SomaticSNVCaller.DEFAULT_ALPHA_T
    experts.add_argument("--tumor-alpha",
                         #required=True,
                         default=default,
                         type=float,
                         help=argparse.SUPPRESS,
                         #help="Significance threshold (alpha) for SNV pvalues"
                         #"  in (relaxed) tumor vcf"
                         #" (default: %f)" % default
                     )
    default = "normal"
    experts.add_argument("-S", "--ign-vcf",
                        default=default,
                        help="Ignore variants in this vcf-file for source"
                        " quality computation in tumor (collides with "
                        " --no-src-qual). Default is to use (stringently"
                         " filtered) predictions in normal sample")

    experts.add_argument("--use-orphan",
                              action="store_true",
                              help="Use orphaned/anomalous reads from pairs"
                              " in all samples")
    experts.add_argument("--baq-off",
                              action="store_true",
                              help="Switch use of BAQ off in all samples")
    experts.add_argument("--call-rlx-extra-args",
                              dest="call_rlx_extra_args",
                              help="Extra arguments to call_rlx (replace dashes with @)")

    experts.add_argument("--no-src-qual",
                        action="store_true",
                        help="Disable use of source quality in tumor (see also -V)")
    experts.add_argument("--debug",
                          action="store_true",
                          help="Enable debugging")
    experts.add_argument("--continue",
                              dest="continue_interrupted",
                              action="store_true",
                              help="continue interrupted run. Will reuse"
                              " existing files, assuming they are complete"
                              " and created with identical options!")

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

    LOG.debug("args = %s" % args)

    # check if outdir exists
    outdir = os.path.dirname(args.outprefix)
    if outdir != "" and not os.path.exists(outdir):
        LOG.error("The directory part of the given output prefix points"
                  " to a non-existing directory: '%s').\n" % (outdir))
        sys.exit(1)

    if not args.dbsnp:
        LOG.warn("No dbsnp file given. Using dbsnp is highly recommended"
                 " when dealing with human data.")
    elif not os.path.exists(args.dbsnp + ".tbi"):
        LOG.warn("Looks like dbsnp was not indexed. Please run bgzip and tabix"
                 " on your dbsnp vcf if 'lofreq somatic' fails and rerun with"
                 " --continue")
    try:
        somatic_snv_caller = SomaticSNVCaller(
            bam_n=args.normal, bam_t=args.tumor, ref=args.ref,
            outprefix=args.outprefix, bed=args.bed, dbsnp=args.dbsnp,
            continue_interrupted=args.continue_interrupted)
    except AssertionError as e:
        LOG.fatal("%s" % str(e))
        sys.exit(1)

    somatic_snv_caller.alpha_n = args.normal_alpha
    somatic_snv_caller.alpha_t = args.tumor_alpha
    somatic_snv_caller.mtc_t = args.tumor_mtc
    somatic_snv_caller.mtc_alpha_t = args.tumor_mtc_alpha
    somatic_snv_caller.indel_mtc_t = args.indel_tumor_mtc
    somatic_snv_caller.indel_mtc_alpha_t = args.indel_tumor_mtc_alpha
    somatic_snv_caller.num_threads = args.num_threads
    somatic_snv_caller.min_cov = args.min_cov
    if args.baq_off:
        somatic_snv_caller.baq_off = True
    else:
        somatic_snv_caller.baq_off = False
    if args.use_orphan:
        somatic_snv_caller.use_orphan = True
    else:
        somatic_snv_caller.use_orphan = False
    if args.call_indels:
        somatic_snv_caller.call_indels = True
    if args.call_rlx_extra_args:
        extra_args = args.call_rlx_extra_args.replace('@', '-').split(" ")
        somatic_snv_caller.call_rlx_extra_args = extra_args

    if args.no_src_qual:
        somatic_snv_caller.src_qual_on = False
    else:
        somatic_snv_caller.src_qual_on = True
        if args.ign_vcf:
            if args.ign_vcf == "normal":
                somatic_snv_caller.src_qual_ign_vcf = ",".join([
                    somatic_snv_caller.vcf_n_str,
                    somatic_snv_caller.vcf_indels_n_str
                    ])
            else:
                somatic_snv_caller.src_qual_ign_vcf = args.ign_vcf

    somatic_snv_caller.do_germline = args.germline

    try:
        somatic_snv_caller.run()
    except:
        LOG.fatal("Somatic SNV caller failed. Exiting")
        #raise
        sys.exit(1)


if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
