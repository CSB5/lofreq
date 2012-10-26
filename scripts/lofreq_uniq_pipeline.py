#!/usr/bin/env python
"""Create a diff of two SNV files, i.e. extract SNVs uniquely
predicted in only one or the other, or alternatively common to both.
This is a pre-processing step for 'somatic' calls, which should then
be completed by lofreq_uniq.py
"""

# Copyright (C) 2011, 2012 Genome Institute of Singapore
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.




#--- standard library imports
#
from __future__ import division
import sys
import logging
import os
import subprocess
# optparse deprecated from Python 2.7 on
from optparse import OptionParser
import shutil

#--- third-party imports
#
# /

#--- project specific imports
#
# /

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011, 2012 Genome Institute of Singapore"
__license__ = "GPL2"


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

MYNAME = os.path.basename(sys.argv[0])

def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog [Options]\n" \
            + "\n" + __doc__
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true", dest="debug",
                      help="enable debugging")
    parser.add_option("", "--bam1",
                      dest="bam1", # type="string|int|float"
                      help="First BAM file (e.g. normal tissue)")
    parser.add_option("", "--bam2",
                      dest="bam2", # type="string|int|float"
                      help="Second BAM file (e.g. cancer tissue)")
    parser.add_option("", "--ref",
                      dest="reffa", # type="string|int|float",
                      help="Reference fasta file")
    parser.add_option("", "--bed",
                      dest="regionbed", # type="string|int|float",
                      help="Region bed file to restrict"
                      " analysis to (see lofreq_regionbed.py)")
    parser.add_option("", "--sb-filter",
                      dest="sbfilter",
                      action="store_true",
                      help="Optional: Apply strand-bias filtering")
    MIN_COV = 10
    parser.add_option("", "--min-cov",
                      dest="mincov",
                      type="int",
                      default=MIN_COV,
                      help="Coverage filter."
                      " Ignore SNVs with coverage below this value"
                      " (default is %d)" % (MIN_COV))
    parser.add_option("-o", "--outdir",
                      dest="outdir", # type="string|int|float",
                      help="Output directory")
    
    return parser



def bam_to_result_filename(bam, what, outdir=None):
    """Derive result filename for step 'what' from BAM filename.
    """
    
    basename = os.path.splitext(os.path.basename(bam))[0]
    if outdir:
        basename = os.path.join(outdir, basename)
    
    if what == 'raw-snv':
        return basename + "-raw.snp"
    elif what == 'filtered-snv':
        return basename + ".snp"
    elif what == 'diff-snv':
        return basename + "-diff.snp"
    elif what == 'uniq-snv':
        return basename + "-uniq.snp"
    else:
        raise ValueError, (
            "Don't how to produce file name for %s" % (what))


def get_bonf(regionbed):
    """Determine Bonferroni factor
    """

    bonf = None
    cmd_list = ['lofreq_bonf.py', '--bed', regionbed]
    try:
        p = subprocess.Popen(cmd_list, 
                             stdout=subprocess.PIPE, 
                             stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
    except OSError:
        LOG.error("Couldn't execute '%s'" % (stderr))
        raise
    
    if p.returncode != 0:
        LOG.fatal("%s failed and produced the following on stderr: %s" % (
            cmd_list[0], stderr))
        LOG.fatal("Could your bed file be malformatted,"
                  " e.g. use space instead of tabs as delimiter?")
        sys.exit(1)
    try:
        bonf = int(stdout)
    except:
        LOG.error("Couldn't parse result from '%s'" % stdout)
        raise
    return bonf



def snv_call_wrapper(bam, snv_raw, reffa, regionbed):
    """Wrapper for lofreq_snpvaller.py
    """

    if os.path.exists(snv_raw):
        LOG.warn("Reusing %s." % snv_raw)
        return snv_raw

    bonf = get_bonf(regionbed)
    LOG.info("Will use Bonferroni factor %d"
             " for SNV quality filtering of %s" % (bonf, bam))

    cmd_list = ['lofreq_snpcaller.py', '--bonf', "%d" % bonf, 
                '-f', reffa, '-l', regionbed,
                '-b', bam, '-o', snv_raw]
    try:
        LOG.info("Executing: %s" % ' '.join(cmd_list))
        subprocess.check_call(cmd_list)
    except subprocess.CalledProcessError:
        LOG.fatal("The following command failed: %s" % (
            ' '.join(cmd_list)))
        raise
    return snv_raw



def snv_filter_wrapper(snv_raw, snv_filt, mincov=None, sbfilter=None):
    """Wrapper for lofreq_filter.py
    """

    if os.path.exists(snv_filt):
        LOG.warn("Reusing %s" % snv_filt)
        return snv_filt

    if not mincov and not sbfilter:
        LOG.info("No filtering requested. Will just copy files")
        shutil.copy(snv_raw, snv_filt)
        return snv_filt

    cmd_list = ['lofreq_filter.py']
    if sbfilter:
        cmd_list.append('--strandbias-holmbonf')
    if mincov:
        cmd_list.extend(['--min-cov', "%d" % mincov])
        cmd_list.extend(['-i', snv_raw, '-o', snv_filt])
    try:
        LOG.info("Executing: %s" % ' '.join(cmd_list))                
        subprocess.check_call(cmd_list)
    except subprocess.CalledProcessError:
        LOG.fatal("The following command failed: %s" % (
            ' '.join(cmd_list)))
        raise
    return snv_filt


def snv_diff_wrapper(snv_to_diff, snv_ref, snv_diff_out):
    """Wrapper for lofreq_diff.py
    """

    if os.path.exists(snv_diff_out):
        LOG.warn("Reusing %s" % snv_diff_out)
        return snv_diff_out
            
    cmd_list = ['lofreq_diff.py', '-s', snv_to_diff, '-t', snv_ref, 
                '-m', 'uniq_to_1', '-o', snv_diff_out]
    try:
        LOG.info("Executing: %s" % ' '.join(cmd_list))
        subprocess.check_call(cmd_list)
    except subprocess.CalledProcessError:
        LOG.fatal("The following command failed: %s" % (
            ' '.join(cmd_list)))
        raise
    return snv_diff_out


def snv_uniq_wrapper(snv_diff, snv_uniq, reffa, ref_bam):
    """Wraper for lofreq_uniq.py
    """

    if os.path.exists(snv_uniq):
        raise ValueError, (
            "Cowardly refusing to overwrite existing file %s" % snv_uniq)

            
    cmd_list = ['lofreq_uniq.py', '-d', snv_diff, '-r', reffa, 
                '-b', ref_bam, '-o', snv_uniq]
    try:
        LOG.info("Executing: %s" % ' '.join(cmd_list))
        subprocess.check_call(cmd_list)
    except subprocess.CalledProcessError:
        LOG.fatal("The following command failed: %s" % (
            ' '.join(cmd_list)))
        raise
    return snv_uniq


def main():
    """
    The main function
    """
    
    parser = cmdline_parser()
    (opts, args) = parser.parse_args()

    if len(args):
        parser.error("Unrecognized arguments found: %s." % (
            ' '.join(args)))
        sys.exit(1)

    if opts.verbose:
        LOG.setLevel(logging.INFO)

    if opts.debug:
        LOG.setLevel(logging.DEBUG)


    # checks for mandatory files
    for (filename, descr, direction) in [
            (opts.bam1, "First BAM file", 'in'),
            (opts.bam2, "Second BAM file", 'in'),
            (opts.reffa, "Reference fasta file", 'in'),
            (opts.regionbed, "Region bed file", 'in')
            ]:
        if not filename:
            parser.error("%s argument missing." % descr)
            sys.exit(1)
        if filename == '-':
            continue
        if direction == 'in' and not os.path.exists(filename):
            sys.stderr.write(
                "file '%s' does not exist.\n" % filename)
            sys.exit(1)
        if direction == 'out' and os.path.exists(filename):
            sys.stderr.write(
                "Refusing to overwrite existing file '%s'.\n" % filename)
            sys.exit(1)
            
    if not opts.outdir:
        parser.error("Missing output directory argument\n")
        sys.exit(1)

        
    if not os.path.isdir(opts.outdir):
        os.mkdir(opts.outdir)


    # snv calls
    #
    #
    for bam in [opts.bam1, opts.bam2]:
        snv_raw = bam_to_result_filename(bam, 'raw-snv', opts.outdir)
        snv_call_wrapper(bam, snv_raw, opts.reffa, opts.regionbed)

        
    # snv filtering
    #
    #
    for bam in [opts.bam1, opts.bam2]:
        snv_raw_in = bam_to_result_filename(bam, 'raw-snv', opts.outdir)
        snv_filt_out = bam_to_result_filename(bam, 'filtered-snv', opts.outdir)
        snv_filter_wrapper(snv_raw_in, snv_filt_out, opts.mincov, opts.sbfilter)


    # diff & uniq
    #
    #
    for (bam_to_diff, bam_ref) in [(opts.bam1, opts.bam2),
                                   (opts.bam2, opts.bam1)]:

        snv_to_diff = bam_to_result_filename(
            bam_to_diff, 'filtered-snv', opts.outdir)
        snv_ref = bam_to_result_filename(
            bam_ref, 'filtered-snv', opts.outdir)
        snv_diff_out = bam_to_result_filename(
            bam_to_diff, 'diff-snv', opts.outdir)
        snv_diff_wrapper(snv_to_diff, snv_ref, snv_diff_out)

        snv_diff = snv_diff_out        
        snv_uniq_out = bam_to_result_filename(
            bam_to_diff, 'uniq-snv', opts.outdir)
        snv_uniq_wrapper(snv_diff, snv_uniq_out, opts.reffa, bam_ref)

    print "Done. %s will now contain the result files, which are" % (opts.outdir)
    print " *-raw.snp: Unfiltered SNV predictions"
    print " *.snp: Filtered SNV predictions"
    print " *-diff.snp: SNVs only predicted in one file (possibly due to coverage etc.)"
    print " *-uniq.snp: SNVs uniquely predicted in one file (the main result)"

if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")


    
