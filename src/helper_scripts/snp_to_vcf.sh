#!/bin/bash

snv=$1
test -n "$snv" || exit 1
test -s "$snv" || exit 1

echo "##fileformat=VCFv4.0"
#echo "##fileDate=20130305"
echo "##source=\"$(basename $0)\""
printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
awk '/^[^#]/ {split($3, b, ">"); phred="."; printf "%s\t%d\t.\t%c\t%c\t%s\tPASS\tAF=%f\n", $1, $2, b[1], b[2], phred, $4}' $snv
