#!/usr/bin/env python
"""Parallel wrapper for LoFreq* SNV Caller: Runs one thread per
seq/chrom listed in header (used as region to make use of indexing
feature) and bed file (if given) and combines results at the end.
"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2013 Genome Institute of Singapore"
__license__ = "Free for non-commercial use"


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




def prob_to_phredqual(prob):
    """WARNING: copy from utils.py. copied here to make script independent
    
    Turns an error probability into a phred value
    
    >>> prob_to_phredqual(0.01)
    20
    """

    from math import log10
    MAX_INT = 2147483647
    # instead of sys.maxint
    
    assert prob >= 0.0, (
        "Probability can't be smaller than 0 but got %f" % prob)
    try:
        return int(round(-10.0 * log10(prob)))
    except ValueError:
        # prob is zero
        #return sys.maxint
        return MAX_INT

                                
def total_num_tests_from_logs(log_files):
    """FIXME:add-doc
    """

    total_num_tests = 0
    for f in log_files:
        fh = open(f, 'r')
        num_tests_found = False

        for line in fh:
            if line.startswith('Number of substitution tests performed'):
                num_tests = int(line.split(':')[1])
                total_num_tests += num_tests
                num_tests_found = True
                break
        if not num_tests_found:
            LOG.fatal("Didn't find number of tests in log-file %s" % (f))
            return -1

        fh.close()

    return total_num_tests



def concat_vcf_files(vcf_files, vcf_concat, source=None):
    """Keeps only head of first vcf file (with '##source=lofreq call'
    replaced by source + \n if given) and write content of all to
    vcf_concat

    """

    # FIXME add gzip support to concat_vcf_files() or make vcfset subcommand"
    assert not os.path.exists(vcf_concat)
    fh_out = open(vcf_concat, 'w')

    for (i, f) in enumerate(vcf_files):
        fh = open(f, 'r')
        for line in fh:
            if i > 0 and line.startswith('#'):
                continue
            if source and line.startswith('##source=lofreq call'):
                fh_out.write(source + "\n")
                continue
            fh_out.write(line)
        fh.close()
    fh_out.close()



#def read_bed_coords(fbed):
#    """Reads coordinates from bed file and returns them as list of
#    tuples (chrom, start, end). Start and end pos are 0-based;
#    half-open just like bed and python ranges
#
#    Based on lofreq 0.6.0
#    """
#
#    with open(fbed, 'r') as fh:
#        for line in fh:
#            if line.startswith('#'):
#                continue
#            if len(line.strip()) == 0:
#                continue
#            try:
#                (chrom, start, end) = line.split("\t")[0:3]
#            except IndexError:
#                LOG.fatal(
#                "FATAL: Failed to parse bed line: %s\n" % (line))
#                raise ValueError
#
#            # http://genome.ucsc.edu/FAQ/FAQformat.html#format1
#            # 4: name, score, strand...
#            start = int(float(start))
#            end = int(float(end))
#            # stupid float is necessary for scientific notation, e.g. 1.25e+08
#            # the only alternative is to use Decimal
#            if end <= start :
#                LOG.fatal("Start value (%d) not lower end value (%d)."
#                " Parsed from file %s" % (
#                    start, end, fbed))
#                raise ValueError
#            yield (chrom, start, end)
#


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
    (stdoutdata, stderrdata) =  process.communicate()

    retcode = process.returncode
    if retcode != 0:
        LOG.fatal("%s exited with error code '%d'." \
                  " Command was '%s'. stderr was: '%s'" % (
                      cmd.split()[0], retcode, cmd, stderrdata))
        raise OSError
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
    (stdoutdata, stderrdata) =  process.communicate()

    retcode = process.returncode
    if retcode != 0:
        LOG.fatal("%s exited with error code '%d'." \
                  " Command was '%s'. stderr was: '%s'" % (
                      cmd.split()[0], retcode, cmd, stderrdata))
        raise OSError
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


def lofreq_cmd_per_sq(bam, lofreq_call_args, tmp_dir):
    """Returns argument for one lofreq call per chromosome. Commands
    will be inversely sorted by number of mapped reads, if index
    supports idxstats otherwise by length. Output vcf file names are
    in order of chromosome.
    """

    # FIXME: this can be improved a lot. When asked for x threads we
    # should simply bin the BAM file such that we get exactly x
    # threads. The number of mapped reads per sq should be taken into
    # account for the binning process and to make it even more
    # complicated the regions in the optional bed file could be taken
    # into account as well. 

    sq_list = sq_list_from_bam(bam)
    if len(sq_list) < 2:
        LOG.warn("Not more than one SQ found in BAM header."
                 " No point in running in parallel mode.")
        sys.exit(1)

    # A quick&dirty way to make sure all processors are always busy is
    # to sort chromosomes by coverage/nreads/length (but keep original
    # order for output prefixing)
    #
    # sq_list contains extra info len and, where possible, number of
    # mapped reads as two or three-tuple. if three elems then third is
    # n_mapped. use that for sorting otherwise fall back to chromosome
    # length
    #
    if len(sq_list[0])==3:
        # have three elements and 3rd is #reads
        sort_idx = 2
        # remove those with no reads mapped
        sq_list = [x for x in sq_list if x[2] > 0]

        if len(sq_list) == 0:
            LOG.warning("Looks like the index for %s is a bit old"
                        " (idxstats reports no reads mapped). Reindexing"
                        " should solve this. Can continue without problem,"
                        " so no need to worry for now though." % (bam))
            sq_list = sq_list_from_bam_samtools(bam)
            if len(sq_list) == 0:
                LOG.fatal("Sorry, fallback solution failed as well :(")
                sys.exit(1)
            
    if len(sq_list[0])!=3:
        # only have two elements and 2nd is chrom length or fallback
        sort_idx = 1
        
    # sort list by sort_idx and prepend original index, which is used
    # for vcf file naming to create the same order as in the bam
    enum_sq_list = sorted(enumerate(sq_list), key=lambda x: x[1][sort_idx], reverse=True)
    # keep only original index and sq name
    enum_sq_list = [(x[0], x[1][0]) for x in enum_sq_list]
    #import pdb; pdb.set_trace()
    #from IPython import embed; embed()
    # if all the above is too compliated just use enumerate(sq_list) below

    #import pdb; pdb.set_trace()
    
    for (i, sq) in enum_sq_list:
        # maintain region order by using index
        reg_str = "%s" % sq
        # which surprisingly works without pos:pos
        cmd = ' '.join(lofreq_call_args)
        cmd += ' --no-default-filter'# needed here whether user-arg or not
        cmd += " -r %s -o %s/%d.vcf > %s/%d.log 2>&1" % (
            reg_str, tmp_dir, i, tmp_dir, i)
        #LOG.warn("DEBUG: yielding %s" % cmd)
        yield cmd


def work(cmd):
    """Command caller wrapper for multiprocessing"""

    #http://stackoverflow.com/questions/884650/how-to-spawn-parallel-child-processes-on-a-multi-processor-system
    #return subprocess.check_output(['which', 'lofreq'])
    #return "Would execute: %s" % (cmd)
    # FIXME any way to capture stderr?
    return subprocess.call(cmd, shell=True)



def main():
    """The main function
    """

    orig_argv = list(sys.argv[1:])

    # parse pparallel specific args: get and remove from list
    #

    # poor man's usage
    #
    if '-h' in orig_argv:
        sys.stderr.write(__doc__ + "\n")
        sys.stderr.write("All arguments except --pp-threads (mandatory),"
                  " --pp-debug, --pp-verbose\nand --pp-dryrun"
                  " will be passed down to 'lofreq call'."
                  "Make sure that the\nremaining args are valid 'lofreq call'"
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

    # Doh!
    if 'call' in orig_argv:
        LOG.fatal("argument 'call' not needed")
        sys.exit(1)

    lofreq_call_args = list(orig_argv)
    #LOG.warn("lofreq_call_args = %s" % ' '.join(lofreq_call_args))


    # check for disallowed args
    #
    # use region ourselves
    #
    for disallowed_arg in ['--plp-summary-only', '-r', '--region']:
        if disallowed_arg in lofreq_call_args:
            LOG.fatal("regions argument -r not allowed in pparallel mode")
            sys.exit(1)


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
        # Just need derivation here and no more filtering at end

    sig_opt = 0.05# NOTE: needs to be default is in lofreq call
    idx = -1
    for arg in ['-s', '--sig']:
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


    # now use one thread per region. output is numerated per thread
    # (%d.log and %d.vcf) and goes into tmp_dir
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

    bam = None
    for arg in lofreq_call_args:
        ext = os.path.splitext(arg)[1].lower()
        if ext in [".bam", ".sam"] and os.path.exists(arg):
            bam = arg
            break
    if not bam:
        LOG.fatal("Could determine BAM file from argument list"
                  " or file doesn't exist")
        sys.exit(1)

    cmd_list = list(lofreq_cmd_per_sq(bam, lofreq_call_args, tmp_dir))
    assert len(cmd_list)>1, (
        "Oops...did get %d instead of multiple commands to run on BAM: %s" % (len(cmd_list), bam))
    LOG.info("Adding %d commands to mp-pool" % len(cmd_list))
    LOG.debug("cmd_list = %s" % cmd_list)
    if dryrun:
        for cmd in cmd_list:
            print "%s" % (cmd)
        LOG.critical("dryrun ending here")
        sys.exit(1)

    pool = multiprocessing.Pool(processes=num_threads)
    results = pool.map(work, cmd_list, chunksize=1)
    # report failures and exit if found
    if any(results):
        #for (res) in results:
        #    LOG.fatal("At least one process reported: %s", res)
        for (res, cmd) in zip(results, cmd_list):
            if res != 0:
                LOG.fatal("The following command failed"
                          " with code %d: %s" %  (res, cmd))
        LOG.fatal("Can't continue")
        sys.exit(1)


    # concat the output
    #
    vcf_concat = os.path.join(tmp_dir, "concat.vcf")
    # maintain order
    vcf_files = [os.path.join(tmp_dir, "%d.vcf" % no)
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
    num_tests = total_num_tests_from_logs(log_files)
    if num_tests == -1:
        sys.exit(1)
    # same as in lofreq_snpcaller.c and used by lofreq2_somatic.py
    sys.stderr.write("Number of substitution tests performed: %d\n" % num_tests)

    cmd = ['lofreq', 'filter', '--only-passed', '-i', vcf_concat, '-o', final_vcf_out]
    if no_default_filter:
        cmd.append('--no-defaults')

    if bonf_opt == 'dynamic':
        # if bonf was computed dynamically, use bonf sum
        bonf = num_tests
        phredqual = prob_to_phredqual(sig_opt/float(bonf))
        cmd.extend(['--snvqual-thresh', "%s" % phredqual])

    elif bonf_opt == 'auto':
        raise NotImplementedError
    
    if bonf_opt not in ['auto', 'dynamic'] and no_default_filter:
        # if bonf_opt was a fixed int, then it was already used properly and
        # there's no need to filter against snv qual. if furthermore,
        # --no-defaults is given we then don't filter at all
        LOG.info("Copying concatenated vcf file to final destination")
        LOG.debug("vcf_concat=%s final_vcf_out=%s" % (vcf_concat, final_vcf_out))
        fh_in = open(vcf_concat, 'r')
        if final_vcf_out == "-":
            fh_out = sys.stdout
        else:
            # check again to make race conditions less likely
            if os.path.exists(final_vcf_out):
                LOG.fatal("Cowardly refusing to overwrite %s with %s" % (
                    final_vcf_out, vcf_concat))
                sys.exit(1)
                                
            if final_vcf_out[-3:] == '.gz':
                fh_out = gzip.open(final_vcf_out, 'w')
            else:
                fh_out = open(final_vcf_out, 'w')
        shutil.copyfileobj(fh_in, fh_out)
        fh_in.close()
        if fh_out != sys.stdout:
            fh_out.close
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
