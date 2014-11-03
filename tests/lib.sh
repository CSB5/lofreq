echoerror() {
    echo "ERROR: $@" 1>&2
}
echook() {
    echo "OK: $@" 1>&2
}
echowarn() {
    echo "WARN: $@" 1>&2
}
echoinfo() {
    echo "INFO: $@" 1>&2
}
echodebug() {
    echo "DEBUG: $@" 1>&2
}

# md5sum is md5 on mac
md5=$(which md5sum 2>/dev/null || which md5)

seq=$(which seq 2>/dev/null || which gseq)

ncpus=$(sysctl -2 hw.ncpu 2>/dev/null || grep -c ^processor /proc/cpuinfo 2>/dev/null || echo 1)
# use 1/8 of available cpus at max but 4 min for parallel tasks
threads=$(echo $ncpus | awk '{n=$1/8; if (n<4) {n=4}; print n}')

LOFREQ=../src/lofreq/lofreq
#LOFREQ=../lofreq_star-2.0.0-beta/lofreq/lofreq


