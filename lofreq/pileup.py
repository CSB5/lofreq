#!/usr/bin/env python
"""
Helper functions for samtools' m/pileup

Should be replaced with PySam in the future once mpileup and all its options are supported properly
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
import subprocess
import logging
import re
import copy
from itertools import chain

#--- third-party imports
#
# /

#--- project specific imports
#
from lofreq import utils


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011, 2012 Genome Institute of Singapore"
__license__ = "GPL2"


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


VALID_BASES = ['A', 'C', 'G', 'T', 'N']
        
    
class PileupColumn():
    """
    Pileup column class. Parses samtools m/pileup output.

    See http://samtools.sourceforge.net/pileup.shtml

    Using namedtuples seemed too inflexible:
    namedtuple('PileupColumn', 'chrom coord ref_base coverage read_bases base_quals')
    pileup_column = PileupColumn(*(line.split('\t')))
    pileup_column = PileupColumn._make((line.split('\t')))
    """


    def __init__(self, line=None):
        """
        """
        # chromosome name
        self.chrom = None
        
        # 0-based coordinate (in: 1-based coordinate)
        self.coord = None
        
        # reference base
        self.ref_base = None

        # locally determined consensus base
        self.cons_base = None
        
        # the number of reads covering the site
        self.coverage = None

        self._bases_and_quals = dict()
        for b in VALID_BASES:
            self._bases_and_quals[b.upper()] = dict()
            self._bases_and_quals[b.lower()] = dict()

        self.num_ins_events = 0
        self.num_del_events = 0
        self.num_read_starts = 0
        self.num_read_ends = 0

        if line:
            self.parse_line(line)

            
#    def determine_cons(self):
#        """
#        Quality aware
#        """
#        
#        #cons_dict = dict()
#        #for base in self._bases_and_quals.keys():
#        #    for q in self._bases_and_quals[base]:
#        #        cons_dict[b.upper()] = cons_dict[b.upper()].get(q, 0) + (1.0 - phredqual_to_prob(self._bases_and_quals[b][q]))
#        #import pdb; pdb.set_trace()
#        # FIXME untested, unfinished

        
    def parse_line(self, line):
        """
        Split a line of pileup output and set values accordingly

        From http://samtools.sourceforge.net/pileup.shtml:
        
        At the read base column, a dot stands for a match to the
        reference base on the forward strand, a comma for a match on
        the reverse strand, `ACGTN' for a mismatch on the forward
        strand and `acgtn' for a mismatch on the reverse strand. A
        pattern `\+[0-9]+[ACGTNacgtn]+' indicates there is an
        insertion between this reference position and the next
        reference position. The length of the insertion is given by
        the integer in the pattern, followed by the inserted sequence.
        Here is an example of 2bp insertions on three reads:

        seq2 156 A 11  .$......+2AG.+2AG.+2AGGG    <975;:<<<<<


        Similarly, a pattern `-[0-9]+[ACGTNacgtn]+' represents a
        deletion from the reference. Here is an exmaple of a 4bp
        deletions from the reference, supported by two reads:

        seq3 200 A 20 ,,,,,..,.-4CACC.-4CACC....,.,,.^~. ==<<<<<<<<<<<::<;2<<

        Also at the read base column, a symbol `^' marks the start of
        a read segment which is a contiguous subsequence on the read
        separated by `N/S/H' CIGAR operations. The ASCII of the
        character following `^' minus 33 gives the mapping quality. A
        symbol `$' marks the end of a read segment. Start and end
        markers of a read are largely inspired by Phil Green's CALF
        format. These markers make it possible to reconstruct the read
        sequences from pileup. SAMtools can optionally append mapping
        qualities to each line of the output. This makes the output
        much larger, but is necessary when a subset of sites are
        selected.
        """

        assert self.coord == None, (
            "Seems like I already read some values")

        line_split = line.split('\t')
        assert len(line_split) == 6, (
            "Couldn't parse pileup line: '%s'" % line)

        self.chrom = line_split[0]
        # in: 1-based coordinate
        self.coord = int(line_split[1]) - 1
        self.ref_base = line_split[2].upper() # paranoia upper()
        self.coverage = int(line_split[3])

        bases = line_split[4]

        # convert quals immediately to phred scale
        quals = [ord(c)-33 for c in line_split[5]]

        # convert special reference markup to actual reference
        bases = bases.replace(".", self.ref_base.upper())
        bases = bases.replace(",", self.ref_base.lower())

        # NOTE: we are not using start/end info, so delete it to avoid
        # confusion. we are not using indel info, so delete it to avoid
        # confusion. deletion on reference ('*') have qualities which
        # will be deleted as well.
        bases = self.rem_startend_markup(bases)
        (bases, quals) = self.rem_indel_markup(bases, quals)
        assert len(bases) == len(quals), (
            "Mismatch between number of parsed bases and quality values at %s:%d\n"
            % (self.chrom, self.coord+1))

        if len(bases) != self.coverage-self.num_del_events:
            LOG.warn("Mismatch between number of bases (= %d) and samtools coverage value (= %d)."
                     " Ins/del events: %d/%d. Cleaned base_str is '%s'. Line was '%s'" % (
                    len(bases), self.coverage, self.num_ins_events, self.num_del_events, bases, line))

        for (i, b) in enumerate(bases):
            if b.upper() not in VALID_BASES: # paranoia. once encountered gaps in pileup
                continue
            q = quals[i]
            self._bases_and_quals[b][q] = self._bases_and_quals[b].get(q, 0) + 1

        # FIXME: this doesn't take qualities into account
        # implement function here and get rid of util import
        (base_counts, cons_base) = utils.count_bases(bases.upper())
        # use pileup refbase on tie and ambigiouity 
        if cons_base == '-' or cons_base == 'N':
            cons_base = self.ref_base
            
        self.cons_base = cons_base

        
            
    def rem_startend_markup(self, bases_str):
        """
        Remove end and start (incl mapping) markup from read bases string

        From http://samtools.sourceforge.net/pileup.shtml:

        ...at the read base column, a symbol `^' marks the start of a
        read segment which is a contiguous subsequence on the read
        separated by `N/S/H' CIGAR operations. The ASCII of the
        character following `^' minus 33 gives the mapping quality. A
        symbol `$' marks the end of a read segment.
        """
        
        org_len = len(bases_str)
        bases_str = bases_str.replace('$', '')
        self.num_read_ends = org_len-len(bases_str)

        org_len = len(bases_str)
        bases_str = re.sub('\^.', '', bases_str)
        self.num_read_starts = org_len-len(bases_str)

        return bases_str

    
    def rem_indel_markup(self, bases_str, quals):
        """
        Remove indel markup from read bases string

        From http://samtools.sourceforge.net/pileup.shtml:

        A pattern '\+[0-9]+[ACGTNacgtn]+' indicates there is an
        insertion between this reference position and the next
        reference position. Similarly, a pattern
        '-[0-9]+[ACGTNacgtn]+' represents a deletion from the
        reference. The deleted bases will be presented as '*' in the
        following lines.
        """
       
        # First the initial +- markup for which no quality value
        # exists: find out how many insertion/deletions happened
        # first, so that you can then delete the right amount of
        # nucleotides afterwards.
        #
        while True:
            match = re.search('[-+][0-9]+', bases_str)
            if not match:
                break

            if bases_str[match.start()] == '+':
                self.num_ins_events += 1
            else:
                assert bases_str[match.start()] == '-'

            num = int(bases_str[match.start()+1:match.end()])
            left = bases_str[:match.start()]
            right = bases_str[match.end()+num:]
            bases_str = left + right


        # now delete the deletion on the reference marked as stars
        # (which have quality values; see also
        # http://seqanswers.com/forums/showthread.php?t=3388)
        # and return

        self.num_del_events = bases_str.count('*')

        quals = [(q) for (b, q) in zip(bases_str, quals)
                      if b != '*']
        bases_str = ''.join([b for b in bases_str if b != '*'])

        return (bases_str, quals)


    def get_counts_for_base(self, base, min_qual=3, keep_strand_info=True):
        """Count base (summarise histograms) and return as fw, rv
        count dict. If keep_strand_info is false, then counts are
        returned as sum of fw and rv
        """

        fw_count = 0
        b = base.upper()
        fw_count += sum([c for (q, c) in self._bases_and_quals[b].iteritems()
                      if q >= min_qual])

        rv_count = 0
        b = base.lower()
        rv_count += sum([c for (q, c) in self._bases_and_quals[b].iteritems()
                      if q >= min_qual])
        if keep_strand_info:
            return (fw_count, rv_count)
        else:
            return sum([fw_count, rv_count])


    def get_all_base_counts(self, min_qual=3, keep_strand_info=True):
        """Frontend to get_count_for_base: Count bases (summarise
        histograms) and return as dict with (uppercase) bases as keys.
        Values will be an int (sum of fw and rv) unless
        keep_strand_info is False (returns sum of both)
        """

        base_counts = dict()
        for base in VALID_BASES:
            base_counts[base] = self.get_counts_for_base(base, min_qual, keep_strand_info)
    
        return base_counts


    def get_base_and_qual_hist(self, keep_strand_info=True):
        """Return a copy of base/quality histograms. If
        keep_strand_info is False, then only uppercase bases will be
        used as keys and values are counts summarised for fw and rv
        strand
        """
        
        if keep_strand_info:
            return copy.deepcopy(self._bases_and_quals)

        # a bit more tricky...
        base_and_qual_hists = dict()
        for base in self._bases_and_quals:
            # don't merge twice
            if base.islower():
                continue
                
            # like dict.update() but add instead of replace
            fw_dict = self._bases_and_quals[base.upper()]
            rv_dict = self._bases_and_quals[base.lower()]
            qual_union = set(fw_dict.keys() + rv_dict.keys())
            base_and_qual_hists[base] = dict([(q, fw_dict.get(q, 0) + rv_dict.get(q, 0))
                                              for q in qual_union])

        return base_and_qual_hists

            
def sq_from_header(header):
    """
    Parse sequence name/s from header. Will return a list. Not sure if
    several names are allowed.
    """

    sq_list = []
    for line in header:
        line_split = line.split()
        try:
            if line_split[0] != "@SQ":
                continue
            if line_split[1].startswith("SN:"):
                sq_list.append(line_split[1][3:])
        except IndexError:
            continue
    return sq_list



def len_for_sq(header, sq):
    """
    Parse sequence length from header.
    """

    for line in header:
        line_split = line.split('\t')
        try:
            if line_split[0] != "@SQ":
                continue
            if not line_split[1].startswith("SN:"):
                continue
            if not line_split[1][3:] == sq:
                continue

            # right line, but which is the right field?
            for field in line_split[2:]:
                if field.startswith("LN:"):
                    return int(field[3:])
        except IndexError:
            continue
    return None

        
            
    
def header(fbam, samtools_binary="samtools"):
    """
    Calls 'samtools -H view', parse output and return
    
    Arguments:
    - fbam:
    is the bam file to parse
    - samtools_binary:
    samtools binary name/path
    
    Results:
    Returns the raw header as list of lines
    
    NOTE:
    Results might change or exeuction fail depending on used samtools version.
    Tested on Version: 0.1.13 (r926:134)
    """
    
    cmd = '%s view -H %s' % (
        samtools_binary, fbam)

    # http://samtools.sourceforge.net/pileup.shtml
    LOG.debug("calling: %s" % (cmd))
    # WARNING: "The data read is buffered in memory, so do not use
    # this method if the data size is large or unlimited." Only other
    # option I see is to write to file.
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
                       
    for line in stderrdata.split("\n"):
        if not len(line):
            continue
        if line == "[mpileup] 1 samples in 1 input files":
            continue
        if line == "[fai_load] build FASTA index.":
            continue
        else:
            LOG.warn("Unhandled line on stderr detected: %s" % (line))

    return stdoutdata.split("\n")
        
