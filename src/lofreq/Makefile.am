AM_CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall \
	-I../cdflib90/ -I../libbam/ -I../uthash @AM_CFLAGS@
bin_PROGRAMS = lofreq
lofreq_SOURCES = bed.h bed.c \
fet.c fet.h \
log.c log.h	\
lofreq_index.c lofreq_index.h \
lofreq_uniq.h  lofreq_uniq.c \
samutils.h samutils.c \
utils.c utils.h \
binom.c binom.h \
lofreq_bamstats.h lofreq_bamstats.c \
lofreq_checkref.h lofreq_checkref.c \
lofreq_indel_quality.h lofreq_indel_quality.c \
multtest.c multtest.h \
vcf.c vcf.h \
viterbi.c viterbi.h \
lofreq_viterbi.c lofreq_viterbi.h \
lofreq_vcfset.c lofreq_vcfset.h \
snpcaller.h snpcaller.c \
lofreq_filter.c lofreq_filter.h  \
lofreq_snpcaller.c lofreq_snpcaller.h \
plp.c plp.h \
lofreq_main.c defaults.h


# note: order matters
lofreq_LDADD = ../libbam/libbam.a ../cdflib90/libcdf.a

