#!/usr/bin/env python
"""Parallel wrapper for LoFreq* SNV Caller.

The idea is to run one thread per seq/chrom listed in header (used as
region to make use of indexing feature) and bed file (if given).
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
from lofreq_star.utils import prob_to_phredqual

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




def total_num_tests_from_logs(log_files):
    """FIXME:add-doc
    """

    total_num_tests = 0
    for f in log_files:
        fh = open(f, 'r')
        num_tests_found = False

        for line in fh:
            if line.startswith('Number of tests performed'):
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

    LOG.warn("FIXME add gzip support to concat_vcf_files()"
             " or make vcfset subcommand")
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



def read_bed_coords(fbed):
    """Reads coordinates from bed file and returns them as list of
    tuples (chrom, start, end). Start and end pos are 0-based;
    half-open just like bed and python ranges

    Based on lofreq 0.6.0
    """

    bed_coords = []

    with open(fbed, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            if len(line.strip()) == 0:
                continue
            try:
                (chrom, start, end) = line.split("\t")[0:3]
            except IndexError:
                LOG.fatal(
                "FATAL: Failed to parse bed line: %s\n" % (line))
                raise ValueError

            # http://genome.ucsc.edu/FAQ/FAQformat.html#format1
            # 4: name, score, strand...
            start = int(float(start))
            end = int(float(end))
            # stupid float is necessary for scientific notation, e.g. 1.25e+08
            # the only alternative is to use Decimal
            if end <= start :
                LOG.fatal("Start value (%d) not lower end value (%d)."
                " Parsed from file %s" % (
                    start, end, fbed))
                raise ValueError
            yield (chrom, start, end)


def sq_list_from_bam(bam):
    """Extract SQs listed in BAM to which at least one read maps
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
        (n_mapped, n_unmapped) = map(int, [n_mapped, n_unmapped])
        if sq != '*' and n_mapped > 0:
            sq_list.append(sq)
    return sq_list


def lofreq_cmd_per_sq(sq_list, lofreq_call_args, tmp_dir):
    """Order in sq_list should be the same as in BAM file1"""

    for (i, sq) in enumerate(sq_list):
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

        
    # poor man's usage
    #
    if '-h' in orig_argv:
        LOG.fatal("Calls 'lofreq call' per chrom/seq in bam"
                         " and combines result at the end.\n"
                         "All arguments except '--pp-threads",
                         "'--pp-debug', '--pp-verbose' and --pp-dryrun"
                         "will be passed down to 'lofreq call'.")
        sys.exit(1)
        
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
        raise NotImplementedError, (
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
        LOG.fatal("Could determine BAM file from argument list")
        sys.exit(1)
        
    sq_list = sq_list_from_bam(bam)
    if len(sq_list) == 1:
        LOG.warn("Only one SQ found in BAM header. No point in running in parallel mode.")
        sys.exit(1)
    #import pdb; pdb.set_trace()
        
    cmd_list = list(lofreq_cmd_per_sq(sq_list, lofreq_call_args, tmp_dir))
    if dryrun:
        for cmd in cmd_list:
            print "%s" % (cmd)
        LOG.critical("dryrun ending here")
        sys.exit(1)
      
    results = []
    pool = multiprocessing.Pool(processes=num_threads)
    p = pool.map_async(work, cmd_list, callback=results.extend)
    p.wait()
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


    # if bonf was computed dynamically, parse bonf from each run's log
    # and filter using their sum
    #
    if bonf_opt == 'dynamic':
        log_files = [os.path.join(tmp_dir, "%d.log" % no)
                     for no in range(len(cmd_list))]
        bonf = total_num_tests_from_logs(log_files)
        if bonf == -1:
            sys.exit(1)
        phredqual = prob_to_phredqual(sig_opt/float(bonf))
        cmd = ['lofreq2_filter.py', '-p', '-i', vcf_concat, '-o', final_vcf_out]
        cmd.extend(['--snv-qual', "%s" % phredqual])
        if no_default_filter:
            cmd.append('--no-defaults')
        cmd = ' '.join(cmd)# subprocess.call takes string
        LOG.info("Executing %s\n" % (cmd))
        if subprocess.call(cmd, shell=True):
            LOG.fatal("Final filtering command failed."
                      " Commmand was %s" % (cmd))
            sys.exit(1)

    elif bonf_opt == 'auto':
        assert NotImplementedError
        
    # otherwise bonf_opt was a fixed int and was already used properly
    else:
        LOG.info("Copying concatenated vcf file to final destination\n")
        fh_in = open(vcf_concat, 'r')
        if final_vcf_out == "-":
            fh_out = sys.stdout
        else:
            fh_out = open(final_vcf_out, 'w')
        shutil.copyfileobj(fh_in, fh_out)
        fh_in.close()
        if fh_out != sys.stdout:
            fh_out.close

    # remove temp files/dir
    if False:
        LOG.warn("Not deleting tmp dir %s" % tmp_dir)
    else:
        shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    sys.stderr.write("NOTE: Running this only makes sense,"
                     " if you have multiple SQs and -b is fixed or coverage is low.\n")
    # otherwise runtime optimization through dyn. bonf. kicks in
    
    LOG.warn("Largely untested!")
    main()
