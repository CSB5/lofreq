#!/usr/bin/env python
"""Experimental implementation of various quality bias checks
"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "WTFPL http://www.wtfpl.net/"


#--- standard library imports
#
import sys
from collections import namedtuple
import gzip

#--- third-party imports
#
#/


VCFEntry = namedtuple('VCFEntry', ['chrom', 'pos', 'dbsnpid', 'ref', 'alt', 'qual', 'filter', 'info'])


def write_var(var, fh=sys.stdout):    
    var = var._replace(pos=str(var.pos))
    fh.write("%s\n" % '\t'.join(var))

    
def vcf_line_to_var(line):
    fields = line.rstrip().split('\t')[:8]
    e = VCFEntry._make(fields)
    return e._replace(pos=int(e.pos))


def var_len(var):
    return len(var.alt)-len(var.ref)    


def af_from_var(var):
    for f in var.info.split(';'):
        if f.startswith('AF='):
            return float(f[3:]) 
    return None


def overlap(e1, e2):
    assert e1.pos <= e2.pos, ("unsorted input: %s > %s" % (e1, e2))
    if e1.pos + var_len(e1) >= e2.pos or e2.pos + var_len(e2) <= e1.pos:
        #print "overlap between %d %s %s and %d %s %s" % (e1.pos, e1.ref, e1.alt, e2.pos, e2.ref, e2.alt)
        return True
    else:
        return False


if len(sys.argv) != 2:
    sys.stderr.write("FATAL: Need (one) vcf file as only argument\n")
    sys.exit(1)
    
vcf = sys.argv[1]
if vcf == "-":
    fh = sys.stdin
elif vcf.endswith(".gz"):
    fh = gzip.open(vcf)
else:
    fh = open(vcf)

prev_vars = []
for line in fh:
    line = line.rstrip()
    if line.startswith('#'):
        print line
        continue
    
    cur_var = vcf_line_to_var(line)
    #if cur_var.pos==2114100:
    #    import pdb; pdb.set_trace()
        
    #print "len(prev_vars)=%d" % (len(prev_vars))
    if len(prev_vars):
        if cur_var.chrom != prev_vars[-1].chrom or not overlap(prev_vars[-1], cur_var):
            #if len(prev_vars)>1:
                #import pdb; pdb.set_trace()
            # pick highest qual/af from stack and empty stack
            picked_var = sorted(prev_vars, key=lambda e: af_from_var(e), reverse=True)[0]
            write_var(picked_var)
            prev_vars = []
    prev_vars.append(cur_var)

# don't forget remaining ones
picked_var = sorted(prev_vars, key=lambda e: af_from_var(e), reverse=True)[0]
write_var(picked_var)

    
if fh != sys.stdout:
    fh.close()
    
#print "%d prev_vars left" % (len(prev_vars))

