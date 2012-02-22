#!/usr/bin/env python
"""
Create summary of Base calls per pileup colum given by stdin filtered
at different quality levels
"""

import sys
from lofreq import pileup

BASES = ['A', 'C', 'G', 'T', 'N']

fh = sys.stdin
try:
    qual = int(sys.argv[1])
except IndexError:
    sys.stderr.write("FATAL: Need quality cutoff argument\n")
    sys.exit(1)
print "coord\tbase\tQ>%d\tcov\tfw\trv" % qual
for line in fh:
    if len(line.strip())==0:
        continue
    
    pcol = pileup.PileupColumn(line)
    pcol.parse_line(line, keep_strand_info=True, delete_raw_values=False)
    
    bq = [(b, q) for (b, q) in zip(pcol.read_bases, pcol.base_quals)]
    res = dict()
    for base in BASES:
        cov = len([(b, q) for (b, q) in bq if b.upper()==base.upper() and q>qual])
        fw = len([(b, q) for (b, q) in bq if b==base.upper() and q>qual])
        rv = len([(b, q) for (b, q) in bq if b==base.lower() and q>qual])
        print "%d\t%s\t\t%d\t%d\t%d" % (pcol.coord+1, base, cov, fw, rv)
        
if fh != sys.stdin:
    fh.close()
