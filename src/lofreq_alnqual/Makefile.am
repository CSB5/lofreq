AM_CFLAGS=-Wall -I../libbam -I../lofreq -I@HTSLIB@ -I@SAMTOOLS@ @AM_CFLAGS@
bin_PROGRAMS = lofreq2_alnqual
lofreq2_alnqual_SOURCES = bam_md_ext.c bam_md_ext.h kprobaln_ext.c kprobaln_ext.h lofreq2_alnqual.c
lofreq2_alnqual_LDADD =  @HTSLIB@/libhts.a @SAMTOOLS@/libbam.a

