#SUBDIRS = bcftools
AM_CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall
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
	kprobaln.c

#libbam_la_LIBADD = bcftools/libbcf.la
#libbam_la_LDFLAGS = -static
#library_includedir=$(includedir)/clustalo/
#library_include_HEADERS = clustal-omega-config.h clustal-omega.h


#bin_PROGRAMS = lofreq_samtools
#lofreq_samtools_SOURCES = bam_plcmd.c sam_view.c \
#	bam_rmdup.c bam_rmdupse.c bam_mate.c bam_stat.c bam_color.c \
#	bamtk.c kaln.c bam2bcf.c bam2bcf_indel.c errmod.c sample.c \
#	bam2depth.c
#lofreq_samtools_LDADD = bcftools/libbcf.la libbam.la

