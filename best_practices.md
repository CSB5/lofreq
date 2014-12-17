---
layout: page
title: 'Best Practices'
---

# Best Practices for LoFreq-Star

#### Aligning your reads

LoFreq works with any mapper that produces BAM. However, we highly recommend to use BWA-MEM. 

If your aligner does not produce mapping qualities or you don't trust them, then switch the use of mapping quality in LoFreq off (option `-J` for `lofreq call`). Two examples where this is likely appropriate are (1) BWASW, which seems to underestimate mapping qualities and (2) Bowtie, which reports only either 0 or 255, under some settings.

Some mappers have trouble placing the second read in paired-end sequencing. The result are so called orphaned reads. LoFreq ignores these by default (so does `samtools`) but you can include them if needed (option `--use-orphan` for `lofreq call`). Note, we have observed high number of orphaned reads in alignments created with Mosaik depending on the used parameters.


Note, we do not anymore recommend to remove non-uniquely mapped reads (i.e. those with mapping quality zero) as we did for previous version.
  
#### Pre-processing your BAM file:

We suggest to pre-process your BAM file by following [GATK's best practices protocol](http://www.broadinstitute.org/gatk/guide/best-practices). In a nutshell, this requires running the following steps: MarkDuplicates (Picard), indel realignment (GATK) and base quality recalibration (GATK).

If you are working on any other organism than human, GATK might complain about a missing vcf-file containing known variants (dbSNP). You can either use an empty vcf-file (only containing a header; even though that's cheating) or a vcf-file produced by running LoFreq on the raw BAM file (likely better alternative). Once done, run LoFreq again on the recalibrated file.


#### Under Construction: Suggested settings for certain scenarios

- BAQ, bonf, -S/-V, -l sens vs spec etc
