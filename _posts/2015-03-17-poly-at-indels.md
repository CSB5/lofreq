---
layout: post
title: Indels in Poly-AT repeats
---

In Illumina data one can often observe low allele frequency indels in
poly-AT regions, which are likely false positives. GATK's BQSR rarely
sets the corresponding indel qualities low enough though and so these
indels often get predicted by LoFreq and they are not automatically
removed. Newer versions of LoFreq (>=2.1.2) will indicate the length of
homopolymer runs in which an indel was predicted (see extra VCF field
called HRUN, akin to IonTorrent's HRUN). Users might want to remove
indels  associated with high values, especially when the predicted indel
ends in an A or T.

Andreas



