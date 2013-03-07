# Top-level Makefile.am for LoFreq
#
# based on Makefile.am for Clustal Omega

# we do not use GNU-style files: NEWS README AUTHORS ChangeLog
# AUTOMAKE_OPTIONS = foreign

#SUBDIRS = src/cdflib90 src/libbam src/helper_scripts src/lofreq_core src/lofreq_samtools src/lofreq_snpcaller  src/lofreq_python
SUBDIRS = src/libbam src/lofreq_core src/lofreq src/lofreq_python
bug-tests: all
	cd tests && $(SHELL) run_all.sh;

