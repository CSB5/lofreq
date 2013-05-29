#!/usr/bin/env python
"""Parallel wrapper for LoFreq* SNV Caller.
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
import subprocess
import multiprocessing
import tempfile
import shutil
import os
#import itertools

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


def bonf_sum_from_dynamic(log_files):
    """FIXME:add-doc
    """
 
    bonf_sum = 0
    for f in log_files:
        fh = open(f, 'r')
        bonf_found = False
        for line in fh:
            if line.startswith('Dynamic Bonferroni factor:'):
                bonf = int(line.rstrip().split()[-1])
                bonf_sum += bonf
                bonf_found = True
                break
        fh.close()
        assert bonf_found, ("Didn't find Bonferroni factor in %s" % f)

    return bonf_sum


        
def concat_vcf_files(vcf_files, vcf_concat, source=None):
    """Keeps only head of first vcf file (with '##source=lofreq call'
    replaced by source + \n if given) and write content of all to
    vcf_concat"""

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
    
    fh = open(fbed, 'r')
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
        bed_coords.append((chrom, start, end))
    fh.close()
    
    return bed_coords



def work(cmd):
    """Command caller wrapper for multiprocessing"""
    
    #http://stackoverflow.com/questions/884650/how-to-spawn-parallel-child-processes-on-a-multi-processor-system
    #return subprocess.check_output(['which', 'lofreq'])
    #return "Would execute: %s" % (cmd)
    return subprocess.call(cmd, shell=True)



def main():
    """The main function
    """
    
    debug = True
    verbose = True


    lofreq_call_args = list(sys.argv[1:])
    #LOG.warn("original args = %s" % ' '.join(lofreq_call_args))

    
    # FIXME usage
    #
    if '-h' in lofreq_call_args:
        sys.stderr.write("Calls 'lofreq call' per region and combines result at the end.\n"
                         "All arguments except '-n threads' will be passed down to 'lofreq call'.\n"
                         "Requires a bed-file as input\n")
        sys.exit(1)


    # debugging is switched on for lofreq call switch it on here as
    # well
    if '--debug' in lofreq_call_args:
        debug = True

        
    # disallow regions. we need them
    #
    if '-r' in lofreq_call_args:
        LOG.fatal("Can't work on regions and bed file at the same time")
        sys.exit(1)

        
    # get num threads and remove from arg list
    #
    try:
        idx = lofreq_call_args.index('-n')
        num_threads = int(lofreq_call_args[idx+1])
        lofreq_call_args = lofreq_call_args[0:idx] +  lofreq_call_args[idx+2:]
    except (ValueError, IndexError) as e:
        LOG.fatal("Parallel wrapper requires -n"
                     " [int] as arg (number of threads)")
        sys.exit(1)
    if num_threads > multiprocessing.cpu_count():
        LOG.warn("Requested number of threads higher than number"
                 " of CPUs. Will reduce value to match number of CPUs")
        num_threads = multiprocessing.cpu_count()
        
        
    # get bed file and parse regions from it
    #
    try:
        idx = lofreq_call_args.index('-l')
        bed_file = lofreq_call_args[idx+1]
        lofreq_call_args = lofreq_call_args[0:idx] +  lofreq_call_args[idx+2:]
    except (ValueError, IndexError) as e:
        LOG.fatal("Parallel wrapper requires a bed-file as argument")
        sys.exit(1)
    try:
        regions = read_bed_coords(bed_file)
    except ValueError:
        LOG.fatal("Parsing of %s failed" % bed_file)
        sys.exit(1)

        
    # parse final/original output file name
    #
    final_vcf_out = "-" # default
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

        
    # parse sig. threshold (needed for final filtering ig bonfis
    # computed dynamically)
    #
    sig_opt = 0.05 # WARN: needs to be same default as in 'lofreq call'
    idx = -1
    for arg in ['-s', '--sig']:
        if arg in lofreq_call_args:
            idx = lofreq_call_args.index(arg)
            break
    if idx != -1:
        sig_opt = float(lofreq_call_args[idx+1])

        
    # determine bonf option
    #
    try:
        bonf_opt = lofreq_call_args[lofreq_call_args.index('-b')+1]
    except (IndexError, ValueError) as e:
        LOG.fatal("Need explicitely set Bonferroni option")
        sys.exit(1)

        
    # need verbose to see dynamic Bonferroni factor in log (if set)
    #
    if '--verbose' not in lofreq_call_args:
        lofreq_call_args.append("--verbose")

        
    # prepend actual lofreq command
    #
    lofreq_call_args.insert(0, 'lofreq')
    lofreq_call_args.insert(1, 'call')
    
    
    # now use one thread per region. output is numerated per thread
    # (%d.log and %d.vcf) and goes into tmp_dir
    # 
    tmp_dir = tempfile.mkdtemp(prefix='lofreq2_call_parallel')
    if debug:
        sys.stderr.write("DEBUG: bonf_opt = %s\n" % bonf_opt)
        sys.stderr.write("DEBUG: sig_opt = %s\n" % sig_opt)
        sys.stderr.write("DEBUG: final_vcf_out = %s\n" % final_vcf_out)
        sys.stderr.write("DEBUG: num_threads = %s\n" % num_threads)
        sys.stderr.write("DEBUG: tmp_dir = %s\n" % tmp_dir)
    #import pdb; pdb.set_trace()
    if verbose:
        sys.stderr.write("Using %d threads with basic args: %s\n" % (
            num_threads, ' '.join(lofreq_call_args)))
    pool = multiprocessing.Pool(processes=num_threads)
    thread_cmds = []
    for (i, reg) in enumerate(regions):
        reg_str = "%s:%d-%d" % reg
        this_cmd = ' '.join(lofreq_call_args)
        this_cmd += " -r %s -o %s/%d.vcf > %s/%d.log 2>&1" % (
            reg_str, tmp_dir, i, tmp_dir, i)
        thread_cmds.append(this_cmd)
    results = []
    p = pool.map_async(work, thread_cmds,
                        callback=results.extend)
    p.wait()
    # report failures and exit if found
    if any(results):    
        for (res, cmd) in zip(results, thread_cmds):
            if res != 0:
                LOG.fatal("The following command failed: %s" %  cmd)
        LOG.fatal("Can't continue")
        sys.exit(1)


    # concat the output
    #
    vcf_concat = os.path.join(tmp_dir, "concat.vcf")
    # don't use glob. this way we get a check for all files for free
    vcf_files = [os.path.join(tmp_dir, "%d.vcf" % no) 
                 for no in range(len(regions))]
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
                     for no in range(len(regions))]
        bonf_sum = bonf_sum_from_dynamic(log_files)
        phredqual = prob_to_phredqual(sig_opt/float(bonf_sum))
        cmd = "lofreq2_filter.py -p -i %s -o %s --snv-phred %d" % (
            vcf_concat, final_vcf_out, phredqual)
        if verbose:
            sys.stderr.write("Executing %s\n" % ' '.join(cmd))
        if subprocess.call(cmd, shell=True):
            LOG.fatal("Final filtering command failed."
                      " Commmand was %s" % (cmd))
            sys.exit(1)
    else:
        # just copy file
        if verbose:
            sys.stderr.write("Copying concatenated vcf file to final destination\n")
        shutil.copy(vcf_concat, final_vcf_out)

        
    # remove temp files/dir
    if False:
        LOG.warn("Not deleting tmp dir %s" % tmp_dir)
    else:
        shutil.rmtree(tmp_dir)

    
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
