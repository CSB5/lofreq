#!/usr/bin/env python
"""
Create summary of Base calls per pileup colum given by stdin filtered
at different quality levels
"""

import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser


from lofreq import pileup



#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

BASES = ['A', 'C', 'G', 'T', 'N']


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

    parser.add_option("-i", "--input",
                      dest="pileup", # type="string|int|float"
                      default="-",
                      help="Pileup (- for stdin)")
    parser.add_option("-q", "--qual",
                      dest="qual",
                      type="int",
                      help="Quality threshold")
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

    if not opts.qual and opts.qual != 0:
        parser.error("Missing quality threshold argument.")
        sys.exit(1)

    if opts.verbose:
        LOG.setLevel(logging.INFO)
        pileup.LOG.setLevel(logging.INFO)

    if opts.debug:
        LOG.setLevel(logging.DEBUG)
        pileup.LOG.setLevel(logging.DEBUG)


    if opts.pileup == "-":
        LOG.info("Reading from stdin")
        fh = sys.stdin
    else:
        LOG.info("Reading from %s" % opts.pileup)
        fh = open(opts.pileup, 'r')

    qual = int(opts.qual)


    print "coord\tbases (Q>=%d)\ttotal\tfw\trv" % qual
    for line in fh:
        if len(line.strip())==0:
            continue
        
        pcol = pileup.PileupColumn()
        pcol.parse_line(line)
        
        base_counts = pcol.get_all_base_counts(qual)
        for base in sorted(base_counts.keys()):
            print "%d\t%s\t\t%d\t%d\t%d" % (pcol.coord+1, base,
                                            sum(base_counts[base]),
                                            base_counts[base][0],
                                            base_counts[base][1])

    if fh != sys.stdin:
        fh.close()
    



if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
