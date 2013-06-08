AM_CFLAGS = -Wall -I../cdflib90/ -I../libbam/ -fPIC @AM_CFLAGS@
AM_LDFLAGS = -shared -fPIC @AM_LDFLAGS@
noinst_LIBRARIES = liblofreq_core.a
liblofreq_core_a_SOURCES = bed.c binom.c fet.c log.c plp.c snpcaller.c utils.c vcf.c \
	bed.h binom.h fet.h log.h plp.h snpcaller.h utils.h vcf.h
liblofreq_core_a_LIBADD = ../cdflib90/libcdf.a ../libbam/libbam.a
EXTRA_DIST = fet.c.AUTHORS fet.c.COPYING

#install-exec-hook:
#	 cd ../lofreq_python && LDFLAGS='-L$(prefix)/lib' $(PYTHON) setup.py install --prefix '$(prefix)'
#
#clean-local:
#	 cd ../lofreq_python && $(PYTHON) setup.py clean
#	 cd ../lofreq_python && rm -rf ./build

