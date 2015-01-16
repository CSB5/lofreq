#!/bin/bash

source lib.sh || exit 1


python -c 'import sys; sys.path.insert(0, "../src/scripts/");import lofreq2_call_pparallel; print "\n".join([str(x) for x in lofreq2_call_pparallel.read_bed_coords("data/reg.bed")])' > /dev/null
if [ $? -eq 0 ]; then
	echook "bed reading function works"
else
	echoerror "bed reading function works"
	exit 1
fi

