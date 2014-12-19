---
layout: post
title: Release of 2.1.1 (bug fix release)
---

A last minute bug sneaked into version 2.1.0 (just released yesterday). It was
 triggered when calling SNVs and indels simultaneously with `lofreq
 call` (not `somatic`) and would have resulted in a segfault in most
 cases (at the automatically run filtering stage).

The actual reason was broken indexing during the FDR correction for
 strand-bias filtering when indels and SNVs where present in the same
 vcf file. The somatic pipeline was not affected by this because it
 keeps separate files for SNVs and indels.

Andreas
