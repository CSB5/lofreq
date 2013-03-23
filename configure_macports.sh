test -n "$PREFIX" && prefixarg="--prefix $PREFIX"
./configure CFLAGS='-I/opt/local/include'  LDFLAGS='-L/opt/local/lib' $prefixarg
