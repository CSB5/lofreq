---
layout: post
title: Where are the FORMAT and SAMPLE fields?
---

We recently got asked a lot why LoFreq's VCF output has no FORMAT and
SAMPLE columns. The reason for their absence is that they represent
genotyping information and current LoFreq versions don't call
genotypes. This shouldn't stop you from using LoFreq, as genotype
information is often not even needed (depending on your
analysis). Some downstream tools might required it (e.g. vcf_melt),
even though these columns are optional according to the
[VCF specification (see e.g. section 1.3)](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

As a workaround, you can just add fake columns in cases where you know
that the information is actually not required and for somatic
samples you can use
([lofreq2\_add\_sample.py](https://github.com/CSB5/lofreq/blob/master/src/tools/scripts/lofreq2_add_sample.py)
, which comes with LoFreq (note the pysam dependency). The next versions of LoFreq (2.2) will be
able to call genotypes.
