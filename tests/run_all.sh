#!/bin/bash

# FIXME missing option to run on cluster
for f in $(ls *sh | grep -v run_all.sh | grep -v lib.sh); do
	echo "*** Running $f";
	./$f || echo "FAILED: $f" ;
	echo;
done
