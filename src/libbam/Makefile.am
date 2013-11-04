AM_CFLAGS = -Wall -fPIC -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE @AM_CFLAGS@
# -D_USE_KNETFILE defined in configure.ac to be passed down to all
AM_LDFLAGS = -shared -fPIC @AM_LDFLAGS@
# -D_USE_KNETFILE defined by toplevel configure to make it active across sub-sources
KNETFILE = knetfile.c
noinst_LIBRARIES = libbam.a

# use
# cfiles=$(grep -A2 libbam_la_SOURCES Makefile.am | sed -e 's,.*=,,' | tr -d '[\n\\]')
# cfiles="$cfiles knetfile.c"
# grep '^#include "' $cfiles | cut -f 2 -d ' ' | tr -d '"' | sort -u 
# to get
BAM_H = bam.h bam_endian.h bgzf.h faidx.h kaln.h khash.h knetfile.h \
	kprobaln.h kseq.h ksort.h kstring.h razf.h sam.h sam_header.h

# original LJOBS s/.o/.c and KNETFILE instead of KNETFILE_O
libbam_a_SOURCES = bgzf.c kstring.c bam_aux.c bam.c bam_import.c \
	sam.c bam_index.c bam_pileup.c bam_lpileup.c bam_md.c razf.c \
	faidx.c bedidx.c $(KNETFILE) sam_header.c \
	kprobaln.c $(BAM_H)
	

