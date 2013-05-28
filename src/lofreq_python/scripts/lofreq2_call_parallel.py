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



# originally from lofreq 0.6.0
def read_bed_coords(fbed):
    """Reads coordinates from bed file and returns them as list of
    tuples (chrom, start, end). Start and end pos are 0-based;
    half-open just like bed and python ranges
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
    #http://stackoverflow.com/questions/884650/how-to-spawn-parallel-child-processes-on-a-multi-processor-system
    #return subprocess.check_output(['which', 'lofreq'])
    #return "Would execute: %s" % (cmd)
    return subprocess.call(cmd, shell=True)



def main():
    """The main function
    """

    lofreq_call_args = sys.argv[1:]
    #LOG.warn("original args = %s" % ' '.join(lofreq_call_args))

    # FIXME usage
    #
    if '-h' in lofreq_call_args:
        sys.stderr.write("Calls 'lofreq call' per region and combines result at the end.\n"
                         "All arguments except '-n threads' will be passed down to 'lofreq call'.\n"
                         "Requires a bed-file as input\n")
        sys.exit(1)

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
        LOG.critical("Parallel wrapper requires -n [int] as arg (number of threads)")
        sys.exit(1)

    # get bed file and parse regions from it
    #
    try:
        idx = lofreq_call_args.index('-l')
        bed_file = lofreq_call_args[idx+1]
        lofreq_call_args = lofreq_call_args[0:idx] +  lofreq_call_args[idx+2:]
    except (ValueError, IndexError) as e:
        LOG.critical("Parallel wrapper requires a bed-file as argument")
        sys.exit(1)
    try:
        regions = read_bed_coords(bed_file)
    except ValueError:
        LOG.fatal("Parsing of %s failed" % bed_file)
        sys.exit(1)
        
    # parse final/original output file name
    #
    orig_out = "-" # default
    idx = -1
    try:
        idx = lofreq_call_args.index('-o')
    except (IndexError, ValueError) as e:
        try:
            idx = lofreq_call_args.index('--out')
        except (IndexError, ValueError) as e:
            pass
    if idx != -1:
        orig_out = lofreq_call_args[idx+1]
        lofreq_call_args = lofreq_call_args[0:idx] +  lofreq_call_args[idx+2:]

    # determine bonf option
    #
    bonf_opt = None
    try:
        bonf_opt = lofreq_call_args[lofreq_call_args.index('-b')+1]
    except (IndexError, ValueError) as e:
        LOG.fatal("Need explicitely set Bonferroni option")
        sys.exit(1)

    # need verbose to see dynamic Bonferroni factor in log (if set)
    #
    if '--verbose' not in lofreq_call_args:
        lofreq_call_args.append("--verbose")

    lofreq_call_args.insert(0, 'lofreq')
    lofreq_call_args.insert(1, 'call')
    

    # now use one thread per region
    # 
    LOG.critical("bonf_opt = %s" % bonf_opt)
    LOG.critical("orig_out = %s" % orig_out)
    LOG.critical("num_threads = %s" % num_threads)
    #import pdb; pdb.set_trace()
    pool = multiprocessing.Pool(processes=num_threads)
    thread_cmds = []
    tmp_dir = tempfile.mkdtemp(prefix='lofreq2_call_parallel')
    LOG.critical("tmp_dir = %s" % tmp_dir)

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

    if any(results):    
        for (res, cmd) in zip(results, thread_cmds):
            if res != 0:
                LOG.fatal("The following command failed: %s" %  cmd)
        LOG.fatal("Can't continue")
        sys.exit(1)

    LOG.fatal("FIXME: unfinished")
    import pdb; pdb.set_trace()
    """
    concat all vcf-files (keep only one header)
    save to tmp file
    if bonf = dyn: 
      parse bonf from all logs ("Dynamic Bonferroni factor: %lld") and filter with sum to orig_out
    else
      cp tmp file to orig_out
    """          

    shutil.rmtree(tmp_dir)

    
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
