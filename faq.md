---
layout: page
title: FAQ
---



##### Do I need to filter LoFreq predictions?

You usually don't. Predicted variants are already filtered using
default parameters (which include coverage, strand-bias, snv-quality
etc). If you however, need extra conservative settings you can use
`lofreq filter` to add some more filter criteria.

##### What happens if i didn't recalibrate base-qualities?

You can run LoFreq on non-recalibrated data, even though we recommend
to follow GATK's best practices for pre-processing of your BAM file
for Illumina data. Non-recalibrated Illumina data might lead to a
slight increase in false positive calls but we've rarely seen bad
examples.

##### How about recalibrating non-human data?

We recommend GATK base-recalibration even for non-human data (or
targeted sequencing), even though GATK requires the input of known
variant sites (a circular problem!?), which are not known for many
organisms (use dbsnp for human data). One option is to use initial
LoFreq predictions as vcf input to GATK and then run LoFreq again or
simply use an empty vcf file.

##### Regions for Targeted Resequencing, Exome etc.

If you want to run LoFreq only on certain regions use the appropriate
bed-file as input with `-l regions.bed`. Not doing so will negatively
affect LoFreq's sensitivity and it might call variants outside of the
desired regions.

##### PCR amplified data (amplicons)

If your data was heavily PCR-amplified LoFreq will likely call SNVs in
primer regions, due to ambiguous primer positions, primer impurities
etc. The best way to deal with this is to create a bed-file that only
lists non-primer regions and use that as input to LoFreq with `-l
regions.bed`.

Do not use Picard's Markduplicate on very high coverage data. It will
mark actual reads as duplicates.

Another problem with heavily PCRed input is that PCR artifacts will
show up as low-frequency SNVs. Their allele frequency will usually be
low (unless a misamplification happened in early cycles).
Computational tools will be unable to distinguish these from true
low-frequency variants, since the mis-incorporated bases look real to
the sequencing machine. To get rid of these you will either have to
run your samples in duplicates (before amplification) or remove SNVs
below the expected PCR error frequency.

##### Choice of mapper: I get few SNVs when running LoFreq on a BAM file with Bowtie/BWA-SW

BWA-SW assigns very low mapping qualities to mapped reads. Consider
disabling the use of mapping quality for `lofreq call` with `-N`.

The same is true for Bowtie, depending on the parameters you used.
Bowtie will sometimes produce alignments with reads that have either
only mapping quality 0 (non unique) or 255 (not available). It is
unclear to us why that happens, or what exactly that means (see also
[this thread](http://seqanswers.com/forums/showthread.php?t=3142)).
LoFreq incorporates mapping quality by default into its variant
calling model, which doesn't help in this scenario. We would suggest
to switch off the use of mapping quality in such cases as well (as
above)

##### `call-parallel` fails with a cryptic error message 

This is easiest to debug by running the single-threaded version `call`
first with the same parameters, check the error message and fix the
command accordingly. Then execute `call-parallel` again with fixed
parameters.
