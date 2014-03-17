# Top-level Makefile.am for LoFreq
#
# based on Makefile.am for Clustal Omega

# we do not use GNU-style files: NEWS README AUTHORS ChangeLog
# AUTOMAKE_OPTIONS = foreign

if ENABLE_HELPER
HELPER = src/helper_scripts
endif

if ENABLE_TOOLS
TOOLS = src/tools
endif

SUBDIRS = src/cdflib90 $(HELPER) src/libbam src/lofreq src/scripts $(TOOLS) src/uthash
bug-tests: all
	cd tests && $(SHELL) run_all.sh;

