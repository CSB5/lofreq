#!/bin/bash

for f in $(ls *sh | grep -v run_all.sh | grep -v lib.sh); do
	echo "*** Running $f";
	./$f || echo "FAILED: $f" ;
	echo;
done
