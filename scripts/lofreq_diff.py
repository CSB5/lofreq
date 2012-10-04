#!/usr/bin/env python
"""Create a diff of two SNV files, i.e. extract SNVs uniquely
predicted in only one or the other, or alternatively common to both.
This is a pre-processing step for 'somatic' calls, which should then
be completed by lofreq_unique.py
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
# optparse deprecated from Python 2.7 on
from optparse import OptionParser

#--- third-party imports
#
# /

#--- project specific imports
#
from lofreq import snp


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



def snvs_to_dict(snv_file):
    """
    Load SNVs from file and convert to a dictionary with identifiers
    as keys and the complete entry as value
    """
    
    snvs  = snp.parse_snp_file(snv_file)
    LOG.info("Parsed %d SNVs from %s" % (len(snvs), snv_file))
    snvs_dict = dict()
    for s in snvs:
        snvs_dict[s.identifier()] = s
    return snvs_dict



def run_diff(snv_file_1, snv_file_2, mode):
    """FIXME:add-missing-doc-str
    """

    snvs_1_dict = snvs_to_dict(snv_file_1)
    snvs_2_dict = snvs_to_dict(snv_file_2)

    snvs_1_ids = set(snvs_1_dict.keys())
    snvs_2_ids = set(snvs_2_dict.keys())
    

    if mode == 'uniq_to_1':
        out_snvs_ids = snvs_1_ids.difference(snvs_2_ids)
        out_snvs = [snvs_1_dict[s] for s in out_snvs_ids]
    elif mode == 'uniq_to_2':
        out_snvs_ids = snvs_2_ids.difference(snvs_1_ids)
        out_snvs = [snvs_2_dict[s] for s in out_snvs_ids]
    elif mode == 'common':
        out_snvs_ids = snvs_1_ids.intersection(snvs_2_ids)
        out_snvs = []
        # can't reuse the original SNV entries: need to compute avg freqs etc.
        for s in out_snvs_ids:
            chrom = snvs_1_dict[s].chrom
            pos = snvs_1_dict[s].pos
            wt = snvs_1_dict[s].wildtype
            var = snvs_1_dict[s].variant
            freq = (snvs_1_dict[s].freq + snvs_2_dict[s].freq)/2.0
            d = dict()
            d['common-snv-src-1'] = snv_file_1
            d['common-snv-src-2'] = snv_file_2
            out_snvs.append(snp.ExtSNP(pos, wt, var, freq, d, chrom))
    else:
        raise ValueError, "Unknown mode of action %s" % mode

    return out_snvs


    
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
    parser.add_option("-o", "--output",
                      dest="out", # type="string|int|float"
                      default='-',
                      help="Output file name ('-' for stdout; default)")
    parser.add_option("-s", "--snv1",
                      dest="snv_file_1", # type="string|int|float",
                      help="SNV file 1")
    parser.add_option("-t", "--snv2",
                      dest="snv_file_2", # type="string|int|float",
                      help="SNV file 2")
    choices = ['uniq_to_1', 'uniq_to_2', 'common']
    parser.add_option("-m", "--mode",
                      choices=choices,
                      dest="mode",
                      help="Mode of action. One of %s" % ', '.join(choices))
    
    return parser



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
    for (filename, descr, direction) in [(opts.snv_file_1, "First SNV file", 'in'),
                                         (opts.snv_file_2, "Second SNV file", 'in'),
                                         (opts.out, "Output file", 'out')]:
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
            
    if not opts.mode:
        parser.error("No mode of action chosen\n")
        sys.exit(1)


    out_snvs = run_diff(opts.snv_file_1, opts.snv_file_2, opts.mode)

    if opts.out == '-':
        fh = sys.stdout
    else:
        fh = open(opts.out, 'w')

    snp.write_snp_file(fh, out_snvs)

    if fh != sys.stdout:
        fh.close()


        
if __name__ == "__main__":
    LOG.debug("FIXME: add support for vcf in and out")
    main()
    LOG.info("Successful program exit")


    
