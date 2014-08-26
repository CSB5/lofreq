AM_CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall \
	-I../cdflib90/ -I../libbam/ -I../uthash @AM_CFLAGS@
bin_PROGRAMS = lofreq
lofreq_SOURCES = fet.c lofreq_filter.h lofreq_snpcaller.h log.c plp.h	\
utils.c fet.h lofreq_index.c lofreq_uniq.c log.h samutils.c utils.h		\
binom.c lofreq_bamstats.c lofreq_checkref.c lofreq_checkref.h lofreq_index.h lofreq_uniq.h multtest.c		\
samutils.h vcf.c binom.h lofreq_bamstats.h lofreq_main.c				\
lofreq_vcfset.c multtest.h snpcaller.c vcf.h defaults.h					\
lofreq_filter.c lofreq_snpcaller.c lofreq_vcfset.h plp.c snpcaller.h

# note: order matters
lofreq_LDADD = ../libbam/libbam.a ../cdflib90/libcdf.a

