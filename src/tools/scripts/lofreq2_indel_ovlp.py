#!/usr/bin/env python
"""Removes overlapping indels
"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "The MIT License"


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


#def var_len(var):
#    return abs(len(var.alt)-len(var.ref))


def af_from_var(var):
    for f in var.info.split(';'):
        if f.startswith('AF='):
            return float(f[3:]) 
    return None


def qual_from_var(var):
    """takes care of missing values, int conversion and ties in comparisons
    """
    if var.qual==".":
        if sys.version_info >= (3, 0):
            return sys.maxsize
        else:
            return sys.maxint    
    else:
        # add AF to deal with ties
        return int(var.qual)+af_from_var(var)


def overlap(v1, v2):
    """determine whether affected positions of two variants overlap
    """

    #if v1.pos==4589049:
    #    import pdb; pdb.set_trace()
    pos1 = set([v1.pos+i for i in range(max([len(v1.ref), len(v1.alt)]))])
    pos2 = set([v2.pos+i for i in range(max([len(v2.ref), len(v2.alt)]))])
    return len(pos1.intersection(pos2))>0

def main():
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
    
    #pic_best_func = af_from_var
    pick_best_func = qual_from_var

    prev_vars = []
    for line in fh:
        line = line.rstrip()
        if line.startswith('#'):
            print(line)
            continue
        
        cur_var = vcf_line_to_var(line)
        if False:
            sys.stderr.write("INFO: looking at %d:%s>%s\n" % (cur_var.pos, cur_var.ref, cur_var.alt))
            sys.stderr.write("INFO: on stack: %s\n" % (', '.join(["%d:%s>%s" % (v.pos, v.ref, v.alt) for v in prev_vars])))
        if len(prev_vars):
            if cur_var.chrom != prev_vars[-1].chrom or not overlap(prev_vars[-1], cur_var):
                # pick highest qual/af from stack and empty stack
                picked_var = sorted(prev_vars, key=lambda e: pick_best_func(e), reverse=True)[0]
                #if len(prev_vars)>1:
                #    print "picked %s from %s" % (picked_var, prev_vars)
                write_var(picked_var)
                prev_vars = []
        prev_vars.append(cur_var)
    
    # don't forget remaining ones
    picked_var = sorted(prev_vars, key=lambda e: pick_best_func(e), reverse=True)[0]
    write_var(picked_var)
    
        
    if fh != sys.stdout:
        fh.close()
        
    #print "%d prev_vars left" % (len(prev_vars))

if __name__ == "__main__":
    main()
