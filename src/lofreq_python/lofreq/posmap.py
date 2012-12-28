#!/usr/bin/env python
"""Class for converting alignmed positions or unaligned positions
between source to target sequence"""


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


#--- third-party imports
#
# /

#--- project specific imports
#
# /

                                                        
__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011, 2012 Genome Institute of Singapore"
__license__ = "GPL2"



# global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

 

class PosMap(object):
    """Position map class

    All zero offset. Gapped pos translate into -1 (instead of None to
    avoid mixup with 0)

    >>> from lofreq import posmap
    >>> from Bio.Seq import Seq
    >>> from Bio.SeqRecord import SeqRecord
    >>> sr1 = SeqRecord(Seq("GA-TC-G"))
    >>> sr1.id = "s1"
    >>> sr2 = SeqRecord(Seq("--GATCG"))
    >>> sr2.id = "s2"
    >>> pm = posmap.PosMap([sr1, sr2])
    >>> 
    >>> aln_to_s1 = pm.convert(None, "s1")
    >>> aln_to_s1
    {0: 0, 1: 1, 2: -1, 3: 2, 4: 3, 5: -1, 6: 4}
    >>>
    >>> aln_to_s2 = pm.convert(None, "s2")
    >>> aln_to_s2
    {0: -1, 1: -1, 2: 0, 3: 1, 4: 2, 5: 3, 6: 4}
    >>>
    >>> s1_to_s2 = pm.convert("s1", "s2")
    >>> s1_to_s2
    {0: -1, 1: -1, 2: 1, 3: 2, 4: 4}
    >>> sorted(s1_to_s2.keys()) == range(len(sr1.seq.ungap("-")))
    True
    >>>
    >>> s2_to_s1 = pm.convert("s2", "s1")
    >>> s2_to_s1 
    {0: -1, 1: 2, 2: 3, 3: -1, 4: 4}
    >>> sorted(s2_to_s1.keys()) == range(len(sr2.seq.ungap("-")))
    True
    """

    def __init__(self, seqrecs=None):
        """
        """

        self.seq_ids = []
        self.pos_map = dict()

        if seqrecs:
            self.generate(seqrecs)

            
    @staticmethod
    def isgap(res, gap_chars = "-"):
        """Return true if given residue is a gap character

        Don't use . or ~ here. . is often used in literature so denote
        consensus and ~ is often used by us internally to signify
        'dunno', which is different from N.
        """
        return (res in gap_chars)
    

    def generate(self, seqrecs):
        """Computes a position map, which is a dict with aligned
        positions as main key. Sequence ids are 2nd dim key and their
        corresponding unaligned position is the value. Gapped
        positions get a -1 assigned.
        
        NOTE: the format is terribly inefficient. Should rather be
        spit out as blocks/ranges for liftover chains (troublesome for
        >2 though)")
        """
       
        self.pos_map = dict()
        self.seq_ids = [s.id for s in seqrecs]
        
        aln_len = len(seqrecs[0].seq)
        for s in seqrecs:
            assert len(s.seq) == aln_len, (
                "Looks like your seqs are not aligned")
        # we use the ids as identifiers so they have to be unique
        assert len(set([s.id for s in seqrecs])) == len(seqrecs)
            
        for apos in xrange(aln_len):
            self.pos_map[apos] = dict()
        for sr in seqrecs:
            # running count of gaps for this seqrec
            gap_count = 0
            for apos in xrange(aln_len):
                res = sr.seq[apos]  
                if not self.isgap(res):
                    self.pos_map[apos][sr.id] = apos - gap_count
                else:
                    self.pos_map[apos][sr.id] = -1
                    gap_count += 1
                    
    
    def output(self, handle=sys.stdout):
        """Print position map
        """

        print "pos map: all pos zero-offset"
        print "aln-pos\t%s" % ('\t'.join(self.seq_ids))
        for aln_pos in sorted(self.pos_map.keys()):
            line = "%d" % aln_pos
            for s in self.seq_ids:
                p = self.pos_map[aln_pos][s]
                line += " %d" % (p)
            handle.write("%s\n" % (line))
            

            
    def parse(self, pos_map_file):
        """Parse position map from file

        NOTE: offset is untouched, i.e. as in file
        """

        LOG.critical("Untested function")
        handle = open(pos_map_file, 'r')
        line = handle.readline()
        header = line.rstrip().split('\t')
        assert header[0] == 'aln-pos', (
            "Was expecting first field to be aln-pos, but it's '%s'" % (
                header[0]))
    
        self.pos_map = dict()
        for line in handle:
            positions = [int(x) for x in line.rstrip().split('\t')]
            assert len(positions) == len(header)
    
            aln_pos = positions[0]
            self.pos_map[aln_pos] = dict()
            for (pos, sid) in zip(positions[1:], header[1:]):
                self.pos_map[aln_pos][sid] = pos
        handle.close()
    
    

    def convert(self, src=None, query=None):
        """Mangles input pos_map and returns a dict containing
        unaligned position matching between src and query ids.

        If query is None, then aligned positions for src are returned.
        Likewise if src is None, then unaligned positions for aligned
        pos are returned
        """
                
        if query and src:
            pos_map = dict([(v[src], v[query]) 
                      for (k, v) in self.pos_map.iteritems()])
        elif src:
            pos_map = dict([(v[src], k) 
                      for (k, v) in self.pos_map.iteritems()])
            #for (k, v) in d.iteritems():
            #    assert k<=v
        elif query:            
            pos_map = dict([(k, v[query]) 
                         for (k, v) in self.pos_map.iteritems()])
        else:
            raise ValueError

        # remove gapped positions
        if pos_map.has_key(-1):
            del pos_map[-1]
            
        return pos_map

    
if __name__ == "__main__":
    import doctest
    doctest.testmod()        


    
