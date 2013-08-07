#!/bin/bash

# Test whether we can detect if the wrong reference was given

source lib.sh || exit 1


bam=data/denv2-pseudoclonal/denv2-pseudoclonal.bam
reffa=data/denv2-simulation/denv2-refseq.fa


cmd="$LOFREQ call -f $reffa $bed $bam"
if eval $cmd 2>/dev/null; then
    echoerror "LoFreq should have failed but didn't. Command was $cmd"
    exit 1
else
    echook "LoFreq detected use of wrong reference"
fi

