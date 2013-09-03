autoreconf
make clean >/dev/null


if [ $(uname) == 'Darwin' ]; then
    # assuming MacOsX with MacPorts:
    ./configure CFLAGS='-I/opt/local/include' LDFLAGS='-L/opt/local/lib' >/dev/null
    make -j 2 >/dev/null
    pushd src/lofreq
    rm lofreq
    # repeat last line from make but replace -lz with MacPorts static one
    cc -std=gnu99 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall -I../libbam/ -I../lofreq_core/ -I/opt/local/include -D_THREAD_SAFE -pthread -D_THREAD_SAFE -pthread -o lofreq lofreq_main.o lofreq_uniq.o lofreq_snpcaller.o /opt/local/lib/libz.a  -L/opt/local/lib ../lofreq_core/liblofreq_core.a ../libbam/libbam.a ../cdflib90/libcdf.a  -lm -pthread || exit 1
    # check with otool -L lofreq
    popd
elif [ $(uname) == 'Linux' ]; then
    # assuming GIS setup with locally statically compiled libz
    ./configure CFLAGS='-I/home/wilma/local/include' LDFLAGS='-L/home/wilma/local/lib' 
    make -j 2
else
    echo "Unknown system: $(uname)" 1>&2
    exit 1
fi

# all systems:
version=$(grep '^AC_INIT' configure.ac | cut -f 2 -d , | tr -d '[\[\]] ')
reldir=lofreq_star-$version
test -d $reldir && exit 1
mkdir $reldir
# use tar to preserve directory structure as if untouched src tree
# but removing unwanted stuff
tar c src/zlib.LICENSE src/libbam/AUTHORS src/libbam/COPYING src/cdflib90.README src/lofreq/lofreq $(find src/lofreq_python/ -name \*py -or -name \*py.README | grep -v '/build/' | grep -v setup) | tar x --strip-components 1 -C $reldir || exit 1
cp README LICENSE $reldir
tar cvzf $(basename $reldir).tgz $reldir && rm -rf $reldir/* && rmdir $reldir

echo "release in $(basename $reldir).tgz. rename to architecture. unpack and test" 1>&2
