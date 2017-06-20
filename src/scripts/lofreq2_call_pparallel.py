#!/usr/bin/env python
"""Parallel wrapper for 'lofreq call': Runs one thread per seq/chrom
listed in header (used as region to make use of indexing feature) and
bed file (if given) and combines results at the end.
"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2013, 2014 Genome Institute of Singapore"
__license__ = "The MIT License"


#--- standard library imports
#
import sys
import logging
import subprocess
import multiprocessing
import tempfile
import shutil
import os
import gzip
from collections import namedtuple

#--- third-party imports
#
#/

#--- project specific imports
#
# sets PATH so that scripts/binary presentin src dir are used first if
# present, i.e. stuff can be run without installing it
try:
    import lofreq2_local
except ImportError:
    pass

# py2to3
MAXINT = 9223372036854775807

Region = namedtuple('Region', ['chrom', 'start', 'end'])
# coordinates in Python-slice / bed format, i.e. zero-based half-open


# global logger
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

BIN_PER_THREAD = 2


def prob_to_phredqual(prob):
    """WARNING: near-identical copy from utils.py. copied here to make
    script independent.

    Turns an error probability into a phred value

    >>> prob_to_phredqual(0.01)
    20
    """

    from math import log10

    assert prob >= 0.0 and prob <= 1.0, (
        "Probability must be >=0 and <=1, but got %f" % prob)
    try:
        return int(round(-10.0 * log10(prob)))
    except ValueError:
        # prob is zero
        return MAXINT

    
def split_region_(start, end):
    """split region (given in zero-based half-open start and end
    coordinates) in two halves
    """
    l = end - start
    assert l > 1, ("Region too small to be split: %d--%d" % (start, end))
    m = l//2# explicit integer divison
    return ((start, start+m), (start+m, end))


def region_length(reg):
    return reg.end-reg.start


def split_region(reg):
    """split region (given in zero-based half-open start and end
    coordinates) in two halves
    """
    return [Region(reg.chrom, x[0], x[1])
            for x in split_region_(reg.start, reg.end)]


def read_bed_coords(fbed):
    """Fault-resistant reading of coordinates from bed file. Yields
    regions as chrom, start, end tuple with zero-based half-open
    coordinates. Based on the implementation in LoFreq 0.6.0
    """

    with open(fbed, 'r') as fh:
        for line in fh:
            if line.startswith('#') or len(line.strip()) == 0:
                continue
            # bed should use tab as delimiter. use whitespace if tab fails.
            if len(line.split('\t')) >= 3:
                (chrom, start, end) = line.split("\t")[0:3]
            elif len(line.split()) >= 3:
                (chrom, start, end) = line.split()[0:3]
            else:
                start = end = "NAN"# caught later
            try:
                # float conversion for support of scientific notation
                (start, end) = [int(float(x)) for x in [start, end]]
            except ValueError:
                if line.startswith('browser') or line.startswith('track'):
                    continue
                else:
                    #import pdb; pdb.set_trace()
                    raise ValueError("Couldn't parse the following line"
                                     " from bed-file %s: %s" % (fbed, line))

            if end <= start or end < 0 or start < 0:
                LOG.fatal("Invalid coordinates start=%d end=%d read from %s" % (
                    start, end, fbed))
                raise ValueError

            yield (chrom, start, end)


def total_num_tests_from_logs(log_files):
    """Extract number of performed tests from all log files and
    returns their sum (for multiple testing correction)
    """

    total_num_snv_tests = 0
    total_num_indel_tests = 0
    for f in log_files:
        fh = open(f, 'r')
        num_snv_tests_found = False
        num_indel_tests_found = False

        for line in fh:
            if line.startswith('Number of substitution tests performed'):
                num_snv_tests = int(line.split(':')[1])
                total_num_snv_tests += num_snv_tests
                num_snv_tests_found = True
            if line.startswith('Number of indel tests performed'):
                num_indel_tests = int(line.split(':')[1])
                total_num_indel_tests += num_indel_tests
                num_indel_tests_found = True
        if not num_snv_tests_found:
            LOG.fatal("Didn't find number of snv tests in log-file %s" % (f))
            return (-1, -1)
        if not num_indel_tests_found:
            LOG.fatal("Didn't find number of indel tests in log-file %s" % (f))
            return (-1, -1)

        fh.close()

    return (total_num_snv_tests, total_num_indel_tests)


def concat_vcf_files(vcf_files, vcf_out, source=None):
    """FIXME source unused
    """

    assert not os.path.exists(vcf_out)

    cmd = ['lofreq', 'vcfset', '-a', 'concat', '-o', vcf_out, '-1']
    cmd.extend(vcf_files)
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        LOG.fatal("The following command failed with return code %d: %s" % (
            e.returncode, ' '.join(cmd)))
        sys.exit(1)


def sq_list_from_bam_samtools(bam):
    """Extract SQs listed in BAM head using samtools

    Elements of returned list is a 2-tuple with sq, length.
    """

    assert os.path.exists(bam), ("BAM file %s does not exist" % bam)
    cmd = 'samtools view -H %s' % (bam)
    LOG.debug("cmd=%s" % cmd)
    process = subprocess.Popen(cmd.split(),
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = process.communicate()

    retcode = process.returncode
    if retcode != 0:
        LOG.fatal("%s exited with error code '%d'." \
                  " Command was '%s'. stderr was: '%s'" % (
                      cmd.split()[0], retcode, cmd, stderrdata))
        raise OSError
    if sys.version_info[0] > 2:
        stdoutdata = stdoutdata.decode()
    stdout_lines = str.splitlines(stdoutdata)
        
    sq_list = []

    for line in stdout_lines:
        if not line.startswith("@SQ"):
            continue
        line_split = line.rstrip().split()
        sn_field = [x for x in line_split if x.startswith("SN:")][0]
        sq = sn_field[3:]
        ln_field = [x for x in line_split if x.startswith("LN:")][0]
        ln = int(ln_field[3:])
        sq_list.append((sq, ln))

    if len(sq_list) == 0:
        LOG.error("No mapping reads in index for %s found."
                  " Reindexing should solve this."
                  " Trying samtools instead" % (bam))
        sys.exit(1)

    return sq_list


def sq_list_from_bam(bam):
    """Extract SQs listed in BAM. Elements of returned list is a
    3-tuple with sq, length and number of mapped reads.
    """

    assert os.path.exists(bam), ("BAM file %s does not exist" % bam)
    cmd = 'lofreq idxstats %s' % (bam)
    LOG.debug("cmd=%s" % cmd)
    process = subprocess.Popen(cmd.split(),
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = process.communicate()

    retcode = process.returncode
    if retcode != 0:
        LOG.fatal("%s exited with error code '%d'." \
                  " Command was '%s'. stderr was: '%s'" % (
                      cmd.split()[0], retcode, cmd, stderrdata))
        raise OSError
    if sys.version_info[0] > 2:
        stdoutdata = stdoutdata.decode()
    stdout_lines = str.splitlines(stdoutdata)

    # orig src: sq_list_from_header()
    sq_list = []
    for line in stdout_lines:
        # chrom len #mapped #unmapped
        (sq, sqlen, n_mapped, n_unmapped) = line.rstrip().split()
        (sqlen, n_mapped, n_unmapped) = [int(x) for x in [
            sqlen, n_mapped, n_unmapped]]
        if sq != '*':
            sq_list.append((sq, sqlen, n_mapped))

    return sq_list


def bins_from_bamheader(bam):
    """Returns regions/bins determine by chromosomes listed in bam
    file. Will skip chromosomes with no reads mapped.
    """

    sq_list = sq_list_from_bam(bam)

    # get list of chromosome and their length. if supported also get
    # number of mapped reads to remove chromosome with no reads mapped
    #
    # have three elements and 3rd is #reads
    if len(sq_list[0]) == 3:
        # remove those with no reads mapped
        sq_list = [x for x in sq_list if x[2] > 0]
        if len(sq_list) == 0:
            LOG.warning("Looks like the index for %s is a bit old"
                        " (idxstats reports no reads mapped). Reindexing"
                        " should solve this. Will continue by calling samtools,"
                        " so no need to worry for now though." % (bam))
            sq_list = sq_list_from_bam_samtools(bam)
            if len(sq_list) == 0:
                LOG.fatal("Sorry, samtools failed as well :(")
                sys.exit(1)

    if len(sq_list) == 0:
        LOG.fatal("Oops. Found no chromsomes in header of %s"
                  " that have any reads mapped!?" % bam)
        sys.exit(1)

    return [(x[0], 0, x[1]) for x in sq_list]


def lofreq_cmd_per_bin(lofreq_call_args, bins, tmp_dir):
    """Returns argument for one lofreq call per bins (Regions()).
    Order is by length byt file naming is according to input order
    """

    # longest bins first, but keep input order as index so that we can
    # use this as file name and only need to concatenate later and
    # output will be sorted by input order

    enum_bins = sorted(enumerate(bins),
                       key=lambda eb: region_length(eb[1]), reverse=True)

    for (i, b) in enum_bins:
        LOG.debug("length sorted bin keeping input index #%d: %s" % (i, b))
        # maintain region order by using index
        reg_str = "%s:%d-%d" % (b.chrom, b.start+1, b.end)
        cmd = ' '.join(lofreq_call_args)
        cmd += ' --no-default-filter'# needed here whether user-arg or not
        cmd += ' -r "%s" -o %s/%d.vcf.gz > %s/%d.log 2>&1' % (
            reg_str, tmp_dir, i, tmp_dir, i)
        #LOG.warn("DEBUG: yielding %s" % cmd)
        yield cmd


def work(cmd):
    """Command caller wrapper for multiprocessing"""

    #print "DEBUG", os.environ["PATH"]
    #from subprocess import Popen, PIPE
    #call(["which", "lofreq"])
    #which = Popen(['which', 'lofreq'], stdout=PIPE).stdout.read()
    #which = Popen("which lofreq", shell=True, stdout=PIPE).stdout.read()
    #LOG.warn("Executing (lofreq=%s): %s" % (which, cmd))
    # res = subprocess.call("lofreq version", shell=True)
    # cmd = 'valgrind --tool=memcheck ' + cmd
    res = subprocess.call(cmd, shell=True)
    if res:
        LOG.fatal("Following command failed with status %d: %s" % (res, cmd))
        # can't exit here: sys.exit(1)
    return res


def main():
    """The main function
    """

    orig_argv = list(sys.argv[1:])

    #
    # 1. parse pparallel specific args: get and remove from list
    #

    # poor man's usage
    #
    if '-h' in orig_argv:
        sys.stderr.write(__doc__ + "\n")
        sys.stderr.write("All arguments except --pp-threads (mandatory),"
                         " --pp-debug, --pp-verbose\nand --pp-dryrun will"
                         " be passed down to 'lofreq call'. Make sure that"
                         " the\nremaining args are valid 'lofreq call'"
                         " args as no syntax check will be\nperformed.\n")
        sys.exit(1)

    verbose = True
    try:
        idx = orig_argv.index('--pp-verbose')
        orig_argv = orig_argv[0:idx] +  orig_argv[idx+1:]
        verbose = True
    except (IndexError, ValueError):
        pass
    if verbose:
        LOG.setLevel(logging.INFO)

    debug = False
    try:
        idx = orig_argv.index('--pp-debug')
        orig_argv = orig_argv[0:idx] +  orig_argv[idx+1:]
        debug = True
    except (IndexError, ValueError):
        pass
    if debug:
        LOG.setLevel(logging.DEBUG)

    dryrun = False
    try:
        idx = orig_argv.index('--pp-dryrun')
        orig_argv = orig_argv[0:idx] +  orig_argv[idx+1:]
        dryrun = True
    except (IndexError, ValueError):
        pass

    # number of threads
    #
    num_threads = -1
    try:
        idx = orig_argv.index('--pp-threads')
        num_threads = int(orig_argv[idx+1])
        orig_argv = orig_argv[0:idx] +  orig_argv[idx+2:]
    except (ValueError, IndexError) as e:
        LOG.fatal("Parallel wrapper requires --pp-threads"
                  " [int] as arg (number of threads)")
        sys.exit(1)
    if num_threads > multiprocessing.cpu_count():
        LOG.fatal("Requested number of threads higher than number"
                 " of CPUs. Will reduce value to match number of CPUs")
        #num_threads = multiprocessing.cpu_count()
        sys.exit(1)


    # reference. if not indexed, multiple threads might want to index it
    # before actually running call, which creates a race condition
    idx = -1
    for arg in ['-f', '--ref']:
        if arg in orig_argv:
            idx = orig_argv.index(arg)
            break
    ref = orig_argv[idx+1]
    if not os.path.exists(ref + ".fai"):
        LOG.fatal("Index for reference %s missing. Use samtools or lofreq faidx %s" % (ref, ref))
        sys.exit(1)


    # Doh!
    if 'call' in orig_argv:
        LOG.fatal("argument 'call' not needed")
        sys.exit(1)

    lofreq_call_args = list(orig_argv)
    #LOG.warn("lofreq_call_args = %s" % ' '.join(lofreq_call_args))


    #
    # 2. check for disallowed args
    #

    # using region ourselves
    #
    # FIXME (re-) use of region could easily be merged into main logic
    # by turning it into a region and intersecting with the rest
    #
    for disallowed_arg in ['--plp-summary-only', '-r', '--region']:
        if disallowed_arg in lofreq_call_args:
            LOG.fatal("%s not allowed in pparallel mode" % disallowed_arg)
            sys.exit(1)


    #
    # 3. modify args that we use
    #

    # get final/original output file name and remove arg
    #
    final_vcf_out = "-"# default
    idx = -1
    for arg in ['-o', '--out']:
        if arg in lofreq_call_args:
            idx = lofreq_call_args.index(arg)
            break
    if idx != -1:
        final_vcf_out = lofreq_call_args[idx+1]
        if os.path.exists(final_vcf_out):
            LOG.fatal("Cowardly refusing to overwrite already existing"
                      " VCF output file %s" % final_vcf_out)
            sys.exit(1)
        lofreq_call_args = lofreq_call_args[0:idx] +  lofreq_call_args[idx+2:]

    # bed-file
    #
    bed_file = None
    idx = -1
    for arg in ['-l', '--bed']:
        if arg in lofreq_call_args:
            idx = lofreq_call_args.index(arg)
            break
    if idx != -1:
        bed_file = lofreq_call_args[idx+1]
        if not os.path.exists(bed_file):
            LOG.fatal("Bed-file %s does not exist" % bed_file)
            sys.exit(1)
        lofreq_call_args = lofreq_call_args[0:idx] +  lofreq_call_args[idx+2:]


    # parse significance and bonf factor. needed for final filtering
    # if bonf is computed dynamically.
    #
    # determine bonf option
    #
    bonf_opt = 'dynamic'# NOTE: needs to be default is in lofreq call
    idx = -1
    for arg in ['-b', '--bonf']:
        if arg in lofreq_call_args:
            idx = lofreq_call_args.index(arg)
            break
    if idx != -1:
        bonf_opt = lofreq_call_args[idx+1]

    if bonf_opt == 'auto':
        raise NotImplementedError(
            'FIXME bonf "auto" handling not implemented')

    sig_opt = 0.01# WARN: needs to be default is in lofreq call
    idx = -1
    for arg in ['-a', '--sig']:
        if arg in lofreq_call_args:
            idx = lofreq_call_args.index(arg)
            break
    if idx != -1:
        sig_opt = float(lofreq_call_args[idx+1])

    # determine whether no-default-filter was given
    #
    no_default_filter = False
    if '--no-default-filter' in lofreq_call_args:
        no_default_filter = True


    # prepend actual lofreq command
    #
    # lofreq2_local makes automatically sure we call the correct binary
    lofreq_call_args.insert(0, 'lofreq')
    lofreq_call_args.insert(1, 'call')

    bam = None
    for arg in lofreq_call_args:
        ext = os.path.splitext(arg)[1].lower()
        if ext in [".bam", ".sam"] and os.path.exists(arg):
            bam = arg
            break
    if not bam:
        LOG.fatal("Couldn't determine BAM file from argument list"
                  " or file doesn't exist")
        sys.exit(1)

    # now use one thread per region. output is numerated per thread
    # (%d.log and %d.vcf.gz) and goes into tmp_dir
    #
    tmp_dir = tempfile.mkdtemp(prefix='lofreq2_call_parallel')
    LOG.debug("tmp_dir = %s" % tmp_dir)
    LOG.debug("bonf_opt = %s" % bonf_opt)
    LOG.debug("sig_opt = %s" % sig_opt)
    LOG.debug("final_vcf_out = %s" % final_vcf_out)
    LOG.debug("num_threads = %s" % num_threads)
    LOG.debug("no_default_filter = %s" % (no_default_filter))
    LOG.debug("lofreq_call_args = %s" % (lofreq_call_args))
    #import pdb; pdb.set_trace()

    LOG.info("Using %d threads with following basic args: %s\n" % (
            num_threads, ' '.join(lofreq_call_args)))


    # At this stage all basic args are known. In theory there are
    # three major variables that determine the splitting logic:
    #
    # - bed file
    # - region arg
    # - sq+length from bam header
    #
    # We should in theory use intersection of all three, bed, region and sq
    # (e.g. using pybedtools or a lightweight alternative)
    # but for now we disallow regions (need it for ourselves; see above)

    bam_bins = [Region._make(x) for x in bins_from_bamheader(bam)]
    if bed_file:
        bed_bins = [Region._make(x) for x in read_bed_coords(bed_file)]

        # if the number of regions is huge and they are scattered
        # across all major chromosomes/sequences, then it's much
        # faster to use the bam_bins as region and bed as extra
        # argument
        bam_sqs = set([b[0] for b in bam_bins])
        bed_sqs = set([b[0] for b in bed_bins])
        if len(bed_bins) > 100*len(bam_bins) and len(bed_sqs) > len(bam_sqs)/10.0:
            bed_sqs = set([b[0] for b in bed_bins])
            bins = [b for b in bam_bins if b[0] in bed_sqs]
            lofreq_call_args.extend(['-l', bed_file])
        else:
            bins = bed_bins
    else:
        bins = bam_bins

    for (i, b) in enumerate(bins):
        LOG.debug("initial bins: #%d %s %d %d len %d" % (
            i, b.chrom, b.start, b.end, region_length(b)))

    # split greedily into bins such that nregions ~ 2*threads:
    # keep more bins than threads to make up for differences in regions
    # even after split
    #
    total_length = sum([region_length(b) for b in bins])
    while True:
        #  inefficient but doesn't matter in practice: should split
        #  max and insert new elements
        # intelligently to avoid sorting whole list.
        bins = sorted(bins, key=lambda b: region_length(b))
        biggest = bins[-1]
        biggest_length = region_length(biggest)

        LOG.debug("biggest_length=%d total_length/(%d*num_threads)=%f" % (
            biggest_length, BIN_PER_THREAD, total_length/(BIN_PER_THREAD*num_threads)))
        if biggest_length < total_length/(BIN_PER_THREAD*num_threads):
            break
        elif biggest_length < 100:
            LOG.warn("Regions getting too small to be efficiently processed")
            break

        biggest = bins.pop()
        (b1, b2) = split_region(biggest)
        bins.extend([b1, b2])

    for (i, b) in enumerate(bins):
        LOG.debug("bins after splitting: #%d %s %d %d len %d" % (
            i, b.chrom, b.start, b.end, region_length(b)))


    # need to make sure bins are order as chromosome order in BAM
    # header (might not be the case in bed-file either which samtools
    # parses the whole BAM when given a bed), otherwise output will
    # not be sorted since it just concatenates
    #
    # sort first by start position and then by predefined chromsome
    # order
    bins = sorted(bins, key=lambda b: b.start)
    sq_list = sq_list_from_bam(bam)
    LOG.debug("sq_list  %s" % sq_list)
    sdict = dict()
    for (i, sq) in enumerate(sq_list):
        sdict[sq[0]] = i
    bins = sorted(bins, key=lambda b: sdict[b.chrom])

    for (i, b) in enumerate(bins):
        LOG.debug("bins after chrom ordering: #%d %s %d %d len %d" % (
            i, b.chrom, b.start, b.end, region_length(b)))

    #bins = [Region('chr22', 0, 50000000)]# TMPDEBUG
    cmd_list = list(lofreq_cmd_per_bin(lofreq_call_args, bins, tmp_dir))
    #FIXME assert len(cmd_list) > 1, (
    #    "Oops...did get %d instead of multiple commands to run on BAM: %s" % (len(cmd_list), bam))
    LOG.info("Adding %d commands to mp-pool" % len(cmd_list))

    #import pdb; pdb.set_trace()
    LOG.debug("cmd_list = %s" % cmd_list)
    if dryrun:
        for cmd in cmd_list:
            print("%s" % cmd)
        LOG.critical("dryrun ending here")
        sys.exit(1)

    #def mycallback(x):
    #    if x:
    #        LOG.warn("multiprocessing.Pool call result: %s" % x)
    #        pool.terminate()

    pool = multiprocessing.Pool(processes=num_threads)
    results = pool.map(work, cmd_list, chunksize=1)
    #results = pool.map_async(work, cmd_list, chunksize=1, callback=mycallback)
    pool.close()# not adding any more
    pool.join()# wait until all done

    if any(results):
        # rerrors printed in work()
        LOG.fatal("Some commands in pool failed. Can't continue")
        sys.exit(1)

    # concat the output by number
    #
    vcf_concat = os.path.join(tmp_dir, "concat.vcf.gz")
    # maintain order
    vcf_files = [os.path.join(tmp_dir, "%d.vcf.gz" % no)
                 for no in range(len(cmd_list))]
    if not all([os.path.exists(f) for f in vcf_files]):
        LOG.fatal("Missing some vcf output from threads")
        sys.exit(1)
    concat_vcf_files(vcf_files, vcf_concat,
                     "##source=%s" % ' '.join(sys.argv))

    # filtering
    #
    log_files = [os.path.join(tmp_dir, "%d.log" % no)
                 for no in range(len(cmd_list))]
    num_snv_tests, num_indel_tests = total_num_tests_from_logs(log_files)
    if num_snv_tests == -1 or num_indel_tests == -1:
        sys.exit(1)
    # same as in lofreq_call.c and used by lofreq2_somatic.py
    sys.stderr.write("Number of substitution tests performed: %d\n" % num_snv_tests)
    sys.stderr.write("Number of indel tests performed: %d\n" % num_indel_tests)

    cmd = ['lofreq', 'filter', '-i', vcf_concat, '-o', final_vcf_out]
    if no_default_filter:
        cmd.append('--no-defaults')

    if bonf_opt == 'dynamic':
        # if bonf was computed dynamically, use bonf sum
        sub_bonf = num_snv_tests
        indel_bonf = num_indel_tests
        if sub_bonf == 0:
            sub_bonf = 1
        if indel_bonf == 0:
            indel_bonf = 1
        sub_phredqual = prob_to_phredqual(sig_opt/float(sub_bonf))
        indel_phredqual = prob_to_phredqual(sig_opt/float(indel_bonf))
        cmd.extend(['--snvqual-thresh', "%s" % sub_phredqual])
        cmd.extend(['--indelqual-thresh', "%s" % indel_phredqual])

    elif bonf_opt == 'auto':
        raise NotImplementedError

    if bonf_opt not in ['auto', 'dynamic'] and no_default_filter:
        # if bonf_opt was a fixed int, then it was already used properly and
        # there's no need to filter against snv qual. if furthermore,
        # --no-defaults is given we then don't filter at all
        #
        LOG.info("Copying concatenated vcf file to final destination")
        LOG.debug("vcf_concat=%s final_vcf_out=%s" % (vcf_concat, final_vcf_out))

        if final_vcf_out == "-":
            fh_in = gzip.open(vcf_concat, 'r')
            fh_out = sys.stdout
            shutil.copyfileobj(fh_in, fh_out)
            fh_in.close()
        else:
            # check again if final output doesn't exist, just to be sure
            if os.path.exists(final_vcf_out):
                LOG.fatal("Cowardly refusing to overwrite %s with %s" % (
                    final_vcf_out, vcf_concat))
                sys.exit(1)

            shutil.copy(vcf_concat, final_vcf_out)

            # try to copy index as well if exists
            tidx_src = vcf_concat + ".tbi"
            if os.path.exists(tidx_src):
                tidx_dst = final_vcf_out + ".tbi"
                shutil.copy(tidx_src, tidx_dst)

    else:
        cmd = ' '.join(cmd)# subprocess.call takes string
        LOG.info("Executing %s\n" % (cmd))
        if subprocess.call(cmd, shell=True):
            LOG.fatal("Final filtering command failed."
            " Commmand was %s" % (cmd))
            sys.exit(1)

    # remove temp files/dir
    if False:
        LOG.warn("Not deleting tmp dir %s" % tmp_dir)
    else:
        shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    main()
