#!/bin/bash


# defaults. have to go before usage()
PICARDDIR_DEFAULT="/mnt/software/stow/picard-1.74/bin/"
threads=4
illumina=0
keep_temp=0
force=0
debug=0
usage() {
    # keep in sync with arg parsing below
cat <<EOF
$(basename $0): wrapper for unique BWA mapping

Performs a mapping of given SR/PE reads with BWA, keeping only
uniquely mapping reads and adding read groups.

Intermediate (sai) files will be kept in a temp dir within /tmp (or
TMPDIR if set). Need to be able to construct index in reference fasta
directory (if missing).

You can set the PICARDDIR environment variable to point to
something else than the default, which is $PICARDIR_DEFAULT

  Mandatory options:
    -f | --fastq1   : Input fastq[.gz] file
    -r | --ref      : Reference fasta file
    -o | --outbam   : Output BAM file
  Optional:
    -h | --help     : Display this help
    -g | --fastq2   : Fastq[.gz], second in pair (optional)
    -k | --keep     : Keep temp. directory
    -t | --threads  : Number of threads to use (default=$threads)
         --illumina : Phred qualities are ASCII64, ie. Illumina 1.3-1.7 (check with FastQC)
         --force    : Force overwriting of files
EOF
}

while [ "$1" != "" ]; do
    case $1 in
        -f | --fastq1 )
            shift
            fastq1=$1
            ;;
        -g | --fastq2 )
            shift
            fastq2=$1
            ;;
        -h | --help )
            usage
            exit
            ;;
            --illumina )
            illumina=1
            ;;
        -k | --keep )
            keep_temp=1
            ;;
        -r | --ref )           
            shift
            reffa=$1
            ;;
        -t | --threads )
            shift
            threads=$1
            ;;
        -o | --outbam )
            shift
            outbam=$1
            ;;
            --force )
            force=1
            ;;
        --debug )
            debug=1
            ;;
        * ) 
            echo "FATAL: unknown argument \"$1\""
            usage
            exit 1
    esac
    shift
done

test -z "$PICARDDIR" && export PICARDDIR=$PICARDDIR_DEFAULT

# make sure all necessary args where given and files exist
#
if [ ! -e $fastq1 ] || [ -z $fastq1 ]; then
    echo "FATAL: fastq file \"$fastq1\" missing" 1>&2
    echo
    usage
    exit 1
fi
if [ ! -z $fastq2 ] && [ ! -s $fastq2 ]; then
    echo "FATAL: fastq file \"$fastq2\" missing" 1>&2
    usage
    exit 1
fi
if [ ! -e $reffa ] || [ -z $reffa ]; then
    echo "FATAL: reference fasta file \"$reffa\" missing" 1>&2
    usage
    exit 1
fi
if [ -z $outbam ]; then
    echo "FATAL: missing output BAM file argument" 1>&2;
    usage
    exit 1
fi
if [ -e $outbam ]  && [ $force -ne 1 ]; then
    echo "FATAL: refusing to overwrite existing BAM file \"$outbam\"" 1>&2
    usage
    exit 1
fi




# check for binaries
for bin in bwa java samtools; do
    if ! which $bin >/dev/null 2>&1; then
        echo "FATAL: couldn't find $bin. make sure it's in your path" 1>&2
        exit 1
    fi
done

# check for picard needed for readgroup adding (needed for GATK)
picard_readgroup_jar=${PICARDDIR}/AddOrReplaceReadGroups.jar
if [ ! -s $picard_readgroup_jar ]; then
    echo "FATAL: couldn't find Picard's $(basename picard_readgroup_jar). Please set PICARDDIR to your Picard installation" 1>&2
    exit 1
fi



# ----

tmpdir=$(mktemp --tmpdir -d "$(basename $0).XXXXXX")

bwa_aln_extra_args="-t $threads -q 3"
# -t: number of threads
# -q: quality threshold for read trimming down to 35bp
if [ $illumina -eq 1 ]; then
    bwa_aln_extra_args="$bwa_aln_extra_args -I"
    # -I: if phred qualities are ascii64, ie. Illumina 1.3-1.7
fi

bwa_samse_extra_args=""
bwa_sampe_extra_args="-s"
# -s disable Smith-Waterman for the unmapped mate

if [ $debug -eq 1 ]; then
cat <<EOF
DEBUG
 threads=$threads
 reffa=$reffa
 fastq1=$fastq1
 fastq2=$fastq2
 bwa_aln_extra_args=$bwa_aln_extra_args
 bwa_samse_extra_args=$bwa_samse_extra_args
 bwa_sampe_extra_args=$bwa_sampe_extra_args
 illumina=$illumina
 tmpdir=$tmpdir
EOF
fi

# create reference index if missing
test -s ${reffa}.pac || bwa index $reffa || exit 1


# align (create sai's)
sais=""
for fastq in $fastq1 $fastq2; do
    sai=$tmpdir/$(basename $fastq | sed -e 's,.gz$,,' | sed -e 's,.fastq$,,' | sed -e 's,.txt$,,').sai
    sais="$sais $sai"
#cat <<EOF
    test $debug -eq 1 && echo "DEBUG: bwa aln $bwa_aln_extra_args -f $sai $reffa $fastq"
    bwa aln $bwa_aln_extra_args -f $sai $reffa $fastq || exit 1    
#EOF
done

# SR/PE specific options to bwa sam and samtools
if [ -z $fastq2 ]; then
    args="samse $bwa_samse_extra_args $reffa $sais $fastq1"
    # remove unmapped reads from single end mapping
else
    args="sampe $bwa_sampe_extra_args $reffa $sais $fastq1 $fastq2"
    # keep only reads mapped in proper pair
fi

# samtools filtering: used to use -f 0x2 for paired-end reads and -F
# 0x4 for single-end. This however is not correct. Noticed this when
# picard failed with """Exception in thread "main"
# net.sf.samtools.SAMFormatException: SAM validation error: ERROR:
# ..., MAPQ should be 0 for unmapped read.""" which results in a
# truncated BAM file (and samtools sort won't fail). The reason there
# was that bwa sometimes happily maps a read spannning genome start
# and end (Burrows Wheeler!), even uniquely, but then sets it to
# unmapped later. So you can have a uniquely mapped read with a MAQP>0
# that's unmapped! Either filter -F 12 (or 0xC) which works with single-end
# or paired-end or Picard VALIDATION_STRINGENCY=LENIENT (or both) to
# prevent it.

#cat <<EOF
test $debug -eq 1 && echo "DEBUG: bwa $args | grep..."
bwa $args | \
    grep '\(^@\|XT:A:U\)' | \
    samtools view -b -S -F 12 - | \
    java -jar $picard_readgroup_jar \
        INPUT=/dev/stdin OUTPUT=/dev/stdout \
        VALIDATION_STRINGENCY=LENIENT \
        RGCN="GIS" RGPL="illumina" RGLB="library-placeholder" \
        RGPU="platform-unit-placeholder" RGSM="sample-name-placeholder" | \
    samtools sort - ${outbam%.bam} || exit 1
#EOF

# user might want to move the file, so don't: samtools index $outbam
if [ $keep_temp -ne 1 ]; then
    test -d $tmpdir && rm -rf $tmpdir
else
    echo "Keeping $tmpdir"
fi

echo "Successful exit. Mapping written to see $outbam"
