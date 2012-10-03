for f in $(ls *sh | grep -v wrapper.sh); do
	echo "*** Running $f";
	./$f || echo "FAILED: $f" ;
	echo;
done
