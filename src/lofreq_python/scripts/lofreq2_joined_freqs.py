#!/usr/bin/env python
"""Report joined nucleotide frequencies for given positions. Only
positions on the same read or read pair will be counted. Given
positions can be annotated with SNV info in the form of ref-alt bases,
in which case different base-call quality filtering mechanisms can be
used and any bases not in the ref-alt combo will be ignored.
"""
 

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2013 Genome Institute of Singapore"
__license__ = "Free for non-commercial use"


#--- standard library imports
#
import sys
import os
import logging
import argparse
from collections import namedtuple

#--- third-party imports
#
# FIXME get rid of pysam dep
try:
    import pysam
except ImportError:
    sys.stderr.write("FATAL(%s): This lofreq utility script relies"
                     " on Pysam (http://code.google.com/p/pysam/)"
                     " which seems to be missing.")
    sys.exit(1)    


#--- project specific imports
#
#try:
#    import lofreq2_local
#except ImportError:
#    pass    
#
#try:
#    from lofreq_star import vcf
#except ImportError:
#    sys.stderr.write(
#        "FATAL(%s): Couldn't find LoFreq's vcf module."
#        " Are you sure your PYTHONPATH is set correctly (= %s)?\n" % (
#            (sys.argv[0], os.environ['PYTHONPATH'])))
#    sys.exit(1)
#

# global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN, 
                    format='%(levelname)s [%(asctime)s]: %(message)s')

SnvPos = namedtuple('SnvPos', ['pos', 'ref', 'alt'])


class ParseSnvPos(argparse.Action):
    def __init__(self,
                 option_strings,
                 dest,
                 nargs=None,
                 const=None,
                 default=None,
                 type=None,
                 choices=None,
                 required=False,
                 help=None,
                 metavar=None):
        argparse.Action.__init__(self,
                                 option_strings=option_strings,
                                 dest=dest,
                                 nargs=nargs,
                                 const=const,
                                 default=default,
                                 type=type,
                                 choices=choices,
                                 required=required,
                                 help=help,
                                 metavar=metavar,
                                 )
        #print 'Initializing CustomAction'
        for (name, value) in sorted(locals().items()):
            if name == 'self' or value is None:
                continue
            #print '  %s = %r' % (name, value)
        return

    def __call__(self, parser, namespace, values, option_string=None):
        #print 'Processing CustomAction for "%s"' % self.dest        
        assert isinstance(values, list)
        snv_pos = []
        for v in values:
            if not ':' in v:
                snv_pos.append(SnvPos(int(v)-1, 'N', 'N'))
            elif '-' in v:
                (pos, bases) = v.split(':')
                pos = int(pos) - 1 
                (refbase, altbase) = bases.upper().split('-')
                snv_pos.append(SnvPos(pos, refbase, altbase))
            else:
                raise ValueError, (
                    "Was expecting a snv_position in the form of"
                    " pos:ref-alt where ref and alt can be ommited")
            #print "Saving snv_pos %s" % (str(snv_pos[-1]))
        setattr(namespace, self.dest, snv_pos)

        
def cmdline_parser():
    """
    creates an argparse instance
    """

    # http://docs.python.org/library/optparse.html
    parser = argparse.ArgumentParser(
        description=__doc__)

    parser.add_argument("--verbose",
                      dest="verbose",
                      action="store_true",
                      help="Enable verbose output")
    
    parser.add_argument("--debug",
                      dest="debug",
                      action="store_true", 
                      help=argparse.SUPPRESS) #"debugging")
    
    parser.add_argument("-b", "--bam",
                      dest="bam",
                      required=True,
                      help="Mapping input file (BAM)")
    
    parser.add_argument("-c", "--chrom",
                        dest="chrom",
                        required=True,
                        help="Positions are on this chromosome/sequence")
    
    parser.add_argument("-p", "--snv-positions", 
                        dest="snv_positions",
                        metavar='NN',
                        required=True,
                        action=ParseSnvPos,
                        nargs='+',
                        help='List of positions in the form of'
                        ' pos:alt-ref or simply pos (which is the'
                        ' same as pos:N-N')
    
    parser.add_argument("-f", "--fasta",
                        dest="ref_fa",
                        help="Will print bases at given positions in"
                        " reference fasta file")
    
    default = 1
    parser.add_argument("-M", "--min-mq",
                        dest="min_mq",
                        type=int,
                        default=default,
                        help="Ignore reads with MQ smaller than this value"
                        " (default=%d)" % default)
    
    default = 3
    parser.add_argument("-B", "--min-bq",
                        dest="min_bq",
                        type=int,
                        default=default,
                        help="Ignore any base with BQ smaller than this value"
                        " (default=%d)" % default)

    default = 20
    parser.add_argument("-A", "--min-altbq",
                        dest="min_altbq",
                        type=int,
                        default=default,
                        help="Ignore any alternate base with BQ smaller"
                        " than this value. Works only if ref-alt"
                        " information was given")# (default=%d)" % default)
    
    return parser



def joined_counts(sam, chrom, snv_positions, min_mq=2, min_bq=3, min_altbq=20):
    """
    - sam: samfile object (pysam)
    - chrom: chromsome/sequence of interest
    - positions: list of SnvPos. Use with alt == 'N' if no altbq specific behaviour is needed
    - min_mq: filter reads with mapping qualiy below this value
    - min_bq: filter bases with call quality below this value
    
    Note, will only report counts that overlap *all* positions.
    """

    min_pos =  min([sp.pos for sp in snv_positions])
    max_pos =  max([sp.pos for sp in snv_positions])

    skip_stats = dict()
    for x in ['dups', 'anomalous', 'qcfail', 'secondary', 'below_mq_min']:
        skip_stats[x] = 0
    
    assert len(snv_positions)>=2 and min(snv_positions)>=0
    assert chrom in sam.references
    assert max_pos < sam.lengths[sam.references.index(chrom)]
    
    pos_overlap = dict()

    # NOTE: we actually don't really need to fetch all pos. one after
    # max should be enough. but what if circular? best to be on the
    # safe side. For full control over BQ, BAQ etc we would have to
    # use mpileup() instead of fetch().
    
    for alnread in sam.fetch(chrom, min_pos, max_pos+1):
        # FIXME +1 necessary?
        assert not alnread.is_unmapped # paranoia

        if alnread.is_duplicate:
            skip_stats['dups'] += 1
            continue
        
        if alnread.is_paired and not alnread.is_proper_pair:
            # check both as is_proper_pair might contain nonsense
            # value if not paired
            skip_stats['anomalous'] += 1
            continue

        if alnread.is_qcfail:
            skip_stats['qcfail'] += 1
            continue

        if alnread.is_secondary:
            skip_stats['secondary'] += 1
            continue

        if alnread.mapq < min_mq:
            skip_stats['below_mq_min'] += 1
            continue
            

        # create a map of ref position (key; long int(!)) and the
        # corresponding query (clipped read) nucleotides (value)
        #
        # pos_nt_map = dict([(rpos, alnread.query[qpos])
        #                   for (qpos, rpos) in alnread.aligned_pairs 
        #                  if qpos!=None and rpos!=None])
        #                   # can't just use if qpos and rpos since
        #                   # they might be 0
        aln_pairs = [(qpos, rpos) for (qpos, rpos) in alnread.aligned_pairs
                     if qpos!=None and rpos!=None]
                     # can't just use if qpos and rpos since
                     # they might be 0
        #pos_nt_map = dict([(rpos, alnread.query[qpos])
        #                   for (qpos, rpos) in aln_pairs
        #                   if ord(alnread.qqual[qpos])-33 >= min_bq])
        #                   # MQs come as ints in pysam, but BQs are
        #                   # ASCII-encoded?

        # create a dictionary with ref pos as key and read base and
        # it's phred qual as value (tuple)
        pos_nt_map = dict([(rpos, 
                            (alnread.query[qpos], ord(alnread.qqual[qpos])-33))
                            for (qpos, rpos) in aln_pairs])

        #import pdb; pdb.set_trace()

        # remove positions where bq is below min_bq
        rem_below_bq = [k for (k, v) in pos_nt_map.items() 
                        if v[1] < min_bq]
        for k in rem_below_bq:
            del pos_nt_map[k]

        # remove positions where bq of an alt base is below min_altbq
        for snv_pos in snv_positions:
            if snv_pos.alt == 'N':
                continue
            if pos_nt_map.has_key(snv_pos.pos):
                base = pos_nt_map[snv_pos.pos][0]
                qual = pos_nt_map[snv_pos.pos][1]
                if base == snv_pos.alt and qual < min_altbq:
                    #LOG.critical("Removing %d because %c is"
                    #             " alt with q %d" % (
                    #                 snv_pos.pos, base, qual))
                    #import pdb; pdb.set_trace()
                    del pos_nt_map[snv_pos.pos]
                    
        # also remove positions where base is neither alt nor ref base
        for snv_pos in snv_positions:
            if snv_pos.ref == 'N' and snv_pos.alt == 'N':
                continue
            if pos_nt_map.has_key(snv_pos.pos):
                base = pos_nt_map[snv_pos.pos][0]
                if base not in [snv_pos.ref, snv_pos.alt]:
                    del pos_nt_map[snv_pos.pos]
                            
        # NOTE: the above filters effectively mean that if a base is
        # removed and overlaps a position of interest, the read and
        # it's partner will never be counted
        
        # create a read-id which is identical for PE reads. NOTE: is
        # there another effective way to identify pairs? Using the
        # name only might lead to silent errors
        if alnread.qname[-2] in [".", "/", "#"]:
            read_id_base = alnread.qname[:-2]
        else:
            read_id_base = alnread.qname

        # record read-id (same for PE reads) and nucleotide for each
        # given position
        for pos in [sv.pos for sv in snv_positions]:
            if pos_nt_map.has_key(pos):
                if not pos_overlap.has_key(pos):
                    pos_overlap[pos] = dict()
                pos_overlap[pos][read_id_base] = pos_nt_map[pos][0]

                
    # create a list of read-id's overlapping all given pos
    overlapping_all = frozenset(pos_overlap[snv_positions[0].pos].keys())
    for snv_pos in snv_positions[1:]:
        overlapping_all = overlapping_all.intersection(
            pos_overlap[snv_pos.pos].keys())
    overlapping_all = list(overlapping_all)

    counts = dict()
    for read_id in overlapping_all:
        key = ''.join([pos_overlap[sp.pos][read_id] for sp in snv_positions])
        counts[key] = counts.get(key, 0) + 1
    
    print "# ignored reads: %s" % (
        ', '.join(["%s: %d" % (k ,v) for (k, v) in skip_stats.items()]))

    counts_sum = sum(counts.values())
    print "# %d reads overlapped with given positions %s"  % (
        counts_sum, ', '.join([str(sp.pos+1) for sp in snv_positions]))

    return counts


    
def main():
    """
    The main function
    """

    parser = cmdline_parser()
    args = parser.parse_args()
    
    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)
        
    # file check
    if not os.path.exists(args.bam):
        LOG.fatal("file '%s' does not exist.\n" % args.bam)
        sys.exit(1)
    if args.ref_fa and not os.path.exists(args.ref_fa):
        LOG.fatal("file '%s' does not exist.\n" % args.ref_fa)
        sys.exit(1)

    if len(args.snv_positions) < 2:
        LOG.fatal("Need at least two positions of interest")
        parser.print_usage(sys.stderr)
        sys.exit(1)

    if args.ref_fa:
        ref_bases = ""
        fastafile = pysam.Fastafile(args.ref_fa)
        for pos in [sp.pos for sp in args.snv_positions]:
            # region = "%s:%d-%d" % (args.chrom, pos+1, pos+1)
            # region needs +1 again which is not intuitive
            # refbase = fastafile.fetch(region=region)
            b = fastafile.fetch(args.chrom, pos, pos+1)
            if not b:
                LOG.fatal("Couldn't fetch region from fastafile %s."
                          " Possible chromsome/sequence name"
                          " mismatch" % (args.ref_fa))
                ref_bases += '-'
            ref_bases += b
        print "# ref %s" % (ref_bases)
            
    sam = pysam.Samfile(args.bam, "rb")
    if args.chrom not in sam.references:
        LOG.fatal("Chromosome/Sequence %s not found in %s" % (
            args.chrom, args.bam))
    counts = joined_counts(sam, args.chrom, args.snv_positions, 
                           args.min_mq, args.min_bq)
    counts_sum = sum(counts.values())
    print "# bases counts freq"
    #for k in sorted(counts, key=counts.get):
    for k in sorted(counts, key=lambda x: counts[x]):
        if counts[k]:
            print "%s %d %.4f" % (k, counts[k], counts[k]/float(counts_sum))

    
if __name__ == "__main__":
    main()
    LOG.critical("TESTS TESTS TEST:"
                 " counts, MQ filter, BQ filter, position overlap")
    LOG.info("Successful exit")
