echoerror() {
    echo "ERROR: $@" 1>&2
}
echook() {
    echo "OK: $@" 1>&2
}
echowarn() {
    echo "WARN: $@" 1>&2
}

# md5sum is md5 on mac
md5=$(which md5sum 2>/dev/null || which md5)
