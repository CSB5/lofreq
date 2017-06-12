#!/bin/bash
make clean >/dev/null

wget -q -O - http://downloads.sourceforge.net/project/samtools/samtools/1.1/samtools-1.1.tar.bz2 | tar -xj
cd samtools-1.1 || exit
echo -e "Making samtools"
make -j 4 &> log_make_samtools.log
cd .. || exit
if [ "$(uname)" == 'Darwin' ]; then
    glibtoolize
elif [ "$(uname)" == 'Linux' ]; then
    libtoolize
fi

echo -e "Bootstrap and configure"
./bootstrap && ./configure SAMTOOLS="${PWD}/samtools-1.1/" HTSLIB="${PWD}/samtools-1.1/htslib-1.1/" &> log_configure.log

if [ "$(uname)" == 'Darwin' ]; then
    echo -e "Making objects for MacOsX"
    make -j 4 &> log_make_lofreq.log
    pushd src/lofreq
    rm lofreq
    # build lofreq linking static libraries (libbam.a, libhts.a, ...)
    cc -std=gnu99 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall \
    -o lofreq \
    -lz -lpthread -lm ../../samtools-1.1/libbam.a ../../samtools-1.1/htslib-1.1/libhts.a \
    /opt/macports_local/lib/libz.a ../cdflib90/libcdf.a \
    bam_index.o bam_md_ext.o bedidx.o binom.o fet.o kprobaln_ext.o lofreq_alnqual.o \
    lofreq_call.o lofreq_checkref.o lofreq_filter.o lofreq_indelqual.o lofreq_index.o \
    lofreq_main.o lofreq_uniq.o lofreq_vcfset.o lofreq_viterbi.o log.o multtest.o \
    plp.o samutils.o snpcaller.o utils.o vcf.o viterbi.o || exit 1
    # check with otool -L lofreq
    popd
elif [ "$(uname)" == 'Linux' ]; then
    # assuming GIS setup with locally statically compiled libz
    echo -e "Making objects for linux"
    make -j 4 &> log_make_lofreq.log
else
    echo "Unknown system: $(uname)" 1>&2
    exit 1
fi

echo -e "Objects done, now packaging"
# all systems:
version=$(grep '^AC_INIT' configure.ac | cut -f 2 -d , | tr -d '[\[\]] ')
reldir=lofreq_star-$version
test -d "$reldir" && exit 1
mkdir "$reldir"
# use tar to preserve directory structure as if untouched src tree
# but removing unwanted stuff
tar c src/cdflib90.README src/lofreq/lofreq "$(find src/scripts/ -name \*py -or -name \*py.README | grep -v '/build/' | grep -v setup)" | tar x --strip-components 1 -C "$reldir" || exit 1
cp README.md LICENSE "$reldir"
tar cvzf "$(basename $reldir)".tgz "$reldir" && rm -rf "${reldir:?}/"* && rmdir "$reldir"
echo "release in $(basename $reldir).tgz. rename to architecture. unpack and test" 1>&2
