# libbcf convenience library (doesn't get installed: noinst)
#
AM_CPPFLAGS = -I.. -D_FILE_OFFSET_BITS=64 -D_USE_KNETFILE -Wall
noinst_LTLIBRARIES = libbcf.la
libbcf_la_SOURCES =  bcf.c vcf.c bcfutils.c prob1.c em.c kfunc.c \
	kmin.c index.c fet.c mut.c bcf2qcall.c \
	bcf.h kmin.h prob1.h

# header inclusion:
# for f in libbcf_la_SOURCES; do grep '^#include' src/lofreq_samtools/bcftools/$f; done | sort -u
