---
layout: post
title: LoFreq version 0.6.0 released
---
Changes:
- New option --cons-as-ref for cases where you want to call a consensus
      base per position and then call SNVs against it instead of calling SNVs
      against
      the given reference
- Default Bonferroni factor reverted to auto instead of auto-ign-zero-cov
- Removed Biopython dependency
- Added CONSVAR INFO field to vcf for denoting majority changes with regard
      to the reference base
- User supplied CFLAGS and LDFLAGS are passed down to Python's extension
      build as well
