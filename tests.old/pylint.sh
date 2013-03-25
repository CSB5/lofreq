#!/bin/bash

source lib.sh || exit 1

#files_to_test=$(grep 'scripts/' ../src/lofreq_python/setup.py | tr -d ",'" | tr -d '[\t ]' | sed -e 's,^,../,')
files_to_test=$(grep '^[^#].*\.py'  ../src/lofreq_python/Makefile.am | grep -v PYTHON | cut -f 2 -d = | tr -d '\' | tr -d '[\t ]')


for pylint in pylint pylint-2.6 pylint-2.7; do
   which $pylint >/dev/null 2>&1 || continue
   echoinfo "Using $(which $pylint)"
   log=$(mktemp -t ${pylint}.XXXXX)
   for f in $files_to_test; do
       echoinfo "Testing $f"
       $pylint -E --rcfile pylint.rc $f >> $log
   done 
   if [ -s $log ]; then
       echoerror "pylint produced errors:"
       cat $log
       exit 1
   else
       echook "pylint produced no errors"
   fi
   rm $log
done


