AM_CFLAGS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall -I../cdflib90/ -I../uthash $(HTSLIB_CPPFLAGS) @AM_CFLAGS@
AM_LDFLAGS = $(LDFLAGS_for_htslib) @AM_LDFLAGS@
bin_PROGRAMS = lofreq
lofreq_SOURCES = bam_md_ext.c bam_md_ext.h \
bedidx.c bam_index.c \
binom.c binom.h \
defaults.h \
fet.c fet.h \
kprobaln_ext.c kprobaln_ext.h \
log.c log.h \
lofreq_alnqual.c lofreq_alnqual.h \
lofreq_index.c lofreq_index.h \
lofreq_uniq.h lofreq_uniq.c \
lofreq_checkref.h lofreq_checkref.c \
lofreq_indelqual.h lofreq_indelqual.c \
lofreq_main.c \
lofreq_viterbi.c lofreq_viterbi.h \
lofreq_vcfset.c lofreq_vcfset.h \
lofreq_filter.c lofreq_filter.h  \
lofreq_call.c lofreq_call.h \
multtest.c multtest.h \
plp.c plp.h \
samutils.h samutils.c \
snpcaller.h snpcaller.c \
utils.c utils.h \
vcf.c vcf.h \
viterbi.c viterbi.h
#lofreq_bamstats.h lofreq_bamstats.c

if HTSLIB_IS_SRC

include $(HTSDIR)/htslib_static.mk
LDFLAGS_for_htslib = $(HTSLIB_static_LDFLAGS)
LIBS_for_htslib = $(HTSDIR)/libhts.a $(HTSLIB_static_LIBS)

else

LDFLAGS_for_htslib = $(HTSLIB_LDFLAGS)
LIBS_for_htslib = -lhts

endif

# note: order matters
lofreq_LDADD = $(LIBS_for_htslib) ../cdflib90/libcdf.a
