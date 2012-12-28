#!/usr/bin/env python
"""Convert SNP positions to alignment positions or unaligned positions
of a target sequence"""


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
import os
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser, SUPPRESS_HELP
import difflib

#--- third-party imports
#
#import Bio
from Bio import SeqIO

#--- project specific imports
#
from lofreq import snp
from lofreq import posmap

                                                        
__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__license__ = "GPL2"



# global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')
 
    
def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog: " + __doc__ + "\n" \
            "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("", "--verbose",
                      dest="verbose",
                      action="store_true",
                      help=SUPPRESS_HELP) #"be verbose")
    parser.add_option("", "--debug",
                      dest="debug",
                      action="store_true", 
                      help=SUPPRESS_HELP) #"debugging")
    parser.add_option("-a", "--msa",
                      dest="msa_file",
                      help="Multiple sequence alignment")
    parser.add_option("-i", "--snp-in",
                      dest="snp_in",
                      help="SNP file")
    parser.add_option("-o", "--snp-out",
                      dest="snp_out",
                      default='-',
                      help="SNP output file ('- for stdout = default)'")
    parser.add_option("-s", "--seqid",
                      dest="seq_id",
                      help="Sequence id of interest (for which SNP"
                      " file was produced). If empty aligned positions are assumed!")
    parser.add_option("-m", "--map-to-id",
                      dest="map_to_id",
                      help="Convert SNP positions to unaligned"
                        " positions of this seq instead of alignment"
                        " coordinates")
    parser.add_option("", "--force-overwrite",
                      action="store_true",
                      dest="force_overwrite",
                      help="Force overwriting of output file")
    default = "fasta"
    parser.add_option("", "--fmt",
                      dest="aln_fmt",
                      default=default,
                      help="Alignment format (must be supported by"
                      " Biopython; default is '%s')" % default)
    return parser



def map_pos_of_snps(snp_list_in, seq_id, map_to_id, pos_map): 
    """The actual main function which will convert snp positions of
    snps in snp_list_in coming from seq_id to map_to_id (or alignment
    if None) based on PosMap pos_map. If seq_id is empty positions are
    assumed to be aligned
    """

    assert seq_id or map_to_id, ("Need src or dst id")
    
    conv_pos_map = pos_map.convert(seq_id, map_to_id)
            
    snp_list_offset = []
    for s in snp_list_in:
        try:
            LOG.debug("Changing pos from %d to %d for %s" % (
                s.pos, conv_pos_map[s.pos], s))
            s.info['orig-pos'] = s.pos+1
            s.pos = conv_pos_map[s.pos]
        except KeyError:
            if map_to_id:
                LOG.fatal("Position %d in %s has no equivalent in %s."
                          % (s.pos+1, seq_id, map_to_id))
            else:
                LOG.fatal("INTERNAL ERROR: Position %d in %s seems invalid."
                          % (s.pos+1, seq_id))
            raise
        #FIXME
        if s.pos == -1:
            LOG.fatal("No match for this SNP. This needs to be fixed by the developers")
            sys.exit(1)
            
        if map_to_id:
            s.info['offset'] = map_to_id
        else: 
            s.info['offset'] = "aligned"
        snp_list_offset.append(s)

        #print "old", s
        #print "new", s
        #print

    return snp_list_offset

    

def main():
    """
    The main function
    """


    parser = cmdline_parser()
    (opts, args) = parser.parse_args()

    if opts.verbose:
        LOG.setLevel(logging.INFO)
    if opts.debug:
        LOG.setLevel(logging.DEBUG)


    if len(args) != 0:
        parser.error("Unrecognized args found")
        sys.exit(1)

    # file check
    for (filename, descr, direction, mandatory) in [
            (opts.snp_in, "SNP input file", 'in', True),
            (opts.snp_out, "SNP output file", 'out', True),
            (opts.msa_file, "MSA file", 'in', True),
            ]:

        if not mandatory and not filename:
            continue
            
        if not filename:
            parser.error("%s argument missing." % descr)
            sys.exit(1)
        if filename == '-':
            continue
        
        if direction == 'in' and not os.path.exists(filename):
            LOG.fatal(
                "file '%s' does not exist.\n" % filename)
            sys.exit(1)
        if direction == 'out' and os.path.exists(filename) \
          and not opts.force_overwrite:
            LOG.fatal(
                "Refusing to overwrite existing file '%s'.\n" % filename)
            sys.exit(1)

    fh_in = open(opts.msa_file, 'r')
    msa_seqs = SeqIO.to_dict(SeqIO.parse(fh_in, opts.aln_fmt))
    fh_in.close()
    if len(msa_seqs)<2:
        LOG.fatal("%s contains only one sequence\n" % (opts.msa_file))
        sys.exit(1)
        
        

    seq_ids_to_check = []
    if opts.seq_id:
        seq_ids_to_check.append(opts.seq_id)
    if opts.map_to_id:
        seq_ids_to_check.append(opts.map_to_id)
    for seq_id in seq_ids_to_check:
        if not seq_id in msa_seqs.keys():
            bestmatch = None
            matches = difflib.get_close_matches(
                seq_id, msa_seqs.keys(), 1, 0.5)
            if matches:
                bestmatch = matches[0]
            LOG.fatal("Couldn't find seq '%s' in MSA file"
            " (best match was: %s)" % (
                seq_id, bestmatch))
            sys.exit(1)
        
    pos_map = posmap.PosMap(msa_seqs.values())
    snp_list = snp.parse_snp_file(opts.snp_in)
    if opts.debug:
        pos_map.output()
    snp_list_offset = map_pos_of_snps(
        snp_list, opts.seq_id, opts.map_to_id, pos_map)

    if opts.snp_out == '-':
        fh_out = sys.stdout
    else:
        fh_out = open(opts.snp_out, 'w')
    snp.write_snp_file(fh_out, snp_list_offset)
    if fh_out != sys.stdout:
        fh_out.close()

    
if __name__ == "__main__":
    main()
    LOG.debug("Successful exit")
    # FIXME in all other scripts LOG.info works but here it's always
    # printed. Couldn't find out why
    #LOG.critical("v=%s d=%s" % (opts.verbose, opts.debug))
    #LOG.critical("LOG.isEnabledFor(logging.INFO)=%s" % LOG.isEnabledFor(logging.INFO))
    #LOG.critical("LOG.getEffectiveLevel()=%s" % LOG.getEffectiveLevel())
