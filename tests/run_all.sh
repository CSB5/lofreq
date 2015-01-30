#!/bin/bash

source lib.sh >/dev/null || exit 1

if hostname | grep -q aquila; then
	on_cluster=1
	ln -sf /mnt/projects/wilma/lofreq/testing/data .	
else
	on_cluster=0
	ln -sf /mnt/pnsg10_projects/wilma/lofreq/testing/data .
fi

#mail="-m bes -M wilma@gis.a-star.edu.sg"
mail=""

log=$(basename $0 .sh).$(date +%Y%m%d-%H%M).log

for f in $(ls *sh | grep -v run_all*sh | grep -v lib.sh); do
	if [ $on_cluster -eq 1 ]; then
		echo "*** Scheduling $f"
		name="lofreq-testing-$(basename $f)"
#cat<<EOF
		qsub -N $name -pe OpenMP $threads $mail -l h_vmem=8G -l h_rt=3:00:00 -o $log -j y -V -b y -cwd "bash $f"
#EOF
	else
		echo "*** Running $f";	
#cat<<EOF
		./$f || echo "FAILED: $f" ;
#EOF
	fi
	echo
done

if [ $on_cluster -eq 1 ]; then
	echo "After all jobs completed check $log"
fi
