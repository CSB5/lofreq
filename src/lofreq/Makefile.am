AM_CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall \
	-I../cdflib90/ -I../libbam/ -I../uthash @AM_CFLAGS@
bin_PROGRAMS = lofreq
lofreq_SOURCES = binom.c binom.h \
defaults.h \
fet.c fet.h \
log.c log.h \
lofreq_index.c lofreq_index.h \
lofreq_uniq.h lofreq_uniq.c \
lofreq_bamstats.h lofreq_bamstats.c \
lofreq_checkref.h lofreq_checkref.c \
lofreq_indelqual.h lofreq_indelqual.c \
lofreq_main.c \
lofreq_viterbi.c lofreq_viterbi.h \
lofreq_vcfset.c lofreq_vcfset.h \
lofreq_filter.c lofreq_filter.h  \
lofreq_snpcaller.c lofreq_snpcaller.h \
multtest.c multtest.h \
plp.c plp.h \
samutils.h samutils.c \
snpcaller.h snpcaller.c \
utils.c utils.h \
vcf.c vcf.h \
viterbi.c viterbi.h


# note: order matters
lofreq_LDADD = ../libbam/libbam.a ../cdflib90/libcdf.a

