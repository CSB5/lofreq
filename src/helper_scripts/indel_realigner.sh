#!/bin/bash

# defaults. have to go before usage()
JAVA_EXTRA_ARGS='-Xmx8g'
known_vcf=""

GATK_DIR_DEFAULT=/mnt/software/stow/GenomeAnalysisTK-2.7-2-g6bda569/
test -z "$GATK_DIR" && export GATK_DIR=$GATK_DIR_DEFAULT
gatk_jar=$GATK_DIR/GenomeAnalysisTK.jar


usage() {
    # keep in sync with arg parsing below
cat <<EOF
$(basename $0): local realignment of reads, misaligned due to indels

Wrapper for GATK's indel realignment functionality.

Needs GATK version 2! You will have to set the GATK_DIR environment
variable to point to your GATK installation.

  Options:
        -i | --bam_in   : input BAM file (indexed)
        -r | --ref_fa   : reference input fasta
        -o | --bam_out  : output BAM file
        -s | --vcf      : already existing vcf file, e.g. for dbSNP
        -h | --help )   : display this help
EOF
}

# defaults
while [ "$1" != "" ]; do
    case $1 in
        -i | --bam_in )
            shift
            bam_in=$1
            ;;
        -o | --bam_out )
            shift
            bam_out=$1
            ;;
        -r | --ref_fa )
            shift
            ref_fa=$1
            ;;
        -s | --vcf )
            shift
            known_vcf=$1
            ;;
        -h | --help ) 
            usage
            exit 0
            ;;
        * ) 
            echo "FATAL: unknown argument \"$1\""
            usage
            exit 1
    esac
    shift
done


# make sure all necessary args where given and files exist
#
if [ ! -e $bam_in ] || [ -z $bam_in ]; then
    echo "FATAL: bam input file \"$bam_in\" missing" 1>&2
    echo
    usage
    exit 1
fi
if [ ! -e $ref_fa ] || [ -z $ref_fa ]; then
    echo "FATAL: reference fasta file \"$ref_fa\" missing" 1>&2
    echo
    usage
    exit 1
fi
if [ -z $bam_out ]; then
    echo "FATAL: missing output BAM file argument" 1>&2;
    usage
    exit 1
fi
if [ -e $bam_out ]; then
    echo "FATAL: refusing to overwrite existing BAM file \"$bam_out\"" 1>&2
    usage
    exit 1
fi


# check for required jar files
#
if [ ! -s $gatk_jar ]; then
    echo "FATAL: couldn't find GATK's $(basename $gatk_jar). Please set the GATK_DIR environemnt variable to your GATK installation" 1>&2
    exit 1
fi


# check for executables
#
for dep in samtools; do
    if ! which $dep >/dev/null 2>&1; then
        echo "FATAL: couldn't find $dep. Please make sure it's in your path" 1>&2
        exit 1
    fi
done

test -e ${ref_fa}.fai || \
    samtools faidx $ref_fa || exit 1

known_vcf_arg=""
if [ ! -z $known_vcf ]; then
    if [ ! -s $known_vcf ]; then
        echo "FATAL: Non-existing vcf file: $known_vcf" 1>&2
        exit 1
    else
         known_vcf_arg="--known $known_vcf"
    fi
fi


# Expecting version 2 of GATK
#
# 1: RealignerTargetCreator
# http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_indels_RealignerTargetCreator.html
# 2: IndelRealigner
# http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_indels_IndelRealigner.html


intervals=$(dirname $bam_out)/$(basename $bam_in .bam).intervals
log=${intervals}.log
if ! java $JAVA_EXTRA_ARGS -jar $gatk_jar \
    -I $bam_in \
    -R $ref_fa \
    -T RealignerTargetCreator \
    -o $intervals \
    $known_vcf_arg > $log 2>&1; then
    echo "ERROR: GATK's RealignerTargetCreator failed. See $log" 1>&2
    exit 1
fi

log=${bam_out}.log
if ! java $JAVA_EXTRA_ARGS -jar $gatk_jar \
   -I $bam_in \
   -R $ref_fa \
   -T IndelRealigner \
   -targetIntervals $intervals \
   -o $bam_out \
    $known_vcf_arg > $log 2>&1; then
    echo "ERROR: GATK's RealignerTargetCreator failed. See $log" 1>&2
    exit 1
fi

samtools index $bam_out 

gzip -f $intervals

#EOF

echo "Successful exit."
