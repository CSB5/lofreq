#!/bin/bash

source lib.sh || exit 1


# testing depth functions

bam="../../lofreq-test-data/denv2-multiplex-replicates/ATCACG_1.bam"
reg="consensus:1-100"

# Python API versus samtools and awk
py_depth=$(python -c "from lofreq_ext import depth_stats; print '%.2f %d' % (depth_stats(\"$bam\", region=\"$reg\"))")
st_depth=$(samtools depth -r $reg $bam | awk '{s+=$NF; n+=1} END {printf "%.2f %d", s/n, n}')
#echo "DEBUG: py_depth=$py_depth st_depth=$st_depth" 1>&2
if [ "$py_depth" != "$st_depth" ]; then
    echoerror "Samtools and Python depth_stats results differ"
else
    echook "Samtools and Python depth_stats returned identical result"
fi

# auto bonferroni using Python API versus samtools sys call
autobonfdepth_py=$(python -c 'from lofreq import sam; print sam.auto_bonf_factor_from_depth("../../lofreq-test-data/denv2-multiplex-replicates/ATCACG_1.bam")')
autobonfdepth_st=$(python -c 'from lofreq import sam; print sam.__auto_bonf_factor_from_depth("../../lofreq-test-data/denv2-multiplex-replicates/ATCACG_1.bam", samtools="samtools")')
#echo "DEBUG: autobonfdepth_py=$autobonfdepth_py autobonfdepth_st=$autobonfdepth_st" 1>&2
if [ "$autobonfdepth_py" != "$autobonfdepth_st" ]; then
    echoerror "Samtools and Python auto_bonf_factor_from_depth() implementations differ"
else
    echook "Samtools and Python auto_bonf_factor_from_depth() give identical result"
fi
