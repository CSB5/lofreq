AM_CFLAGS = -Wall -I../cdflib90/ -I../libbam/ -fPIC
AM_LDFLAGS = -shared -fPIC
noinst_LIBRARIES = liblofreq_core.a
liblofreq_core_a_SOURCES = bam2depth.c fet.c log.c sam_view.c snpcaller.c utils.c vcf.c \
	bam2depth.h binom.h fet.h log.h sam_view.h snpcaller.h utils.h vcf.h
liblofreq_core_a_LIBADD = ../libbam/libbam.a
EXTRA_DIST = fet.c.AUTHORS fet.c.COPYING

#install-exec-hook:
#	 cd ../lofreq_python && LDFLAGS='-L$(prefix)/lib' $(PYTHON) setup.py install --prefix '$(prefix)'
#
#clean-local:
#	 cd ../lofreq_python && $(PYTHON) setup.py clean
#	 cd ../lofreq_python && rm -rf ./build

