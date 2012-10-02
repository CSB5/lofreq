#!/bin/bash

echo "*** Testing if scripts are in path"
for s in lofreq_snpcaller.py lofreq_filter.py lofreq_bonf.py; do
	which $s >/dev/null || exit 1
	echo "$s: OK"	
done

echo
echo "*** Running internal sensitivity test"
lofreq_snpcaller.py --test-sens || exit 1



 