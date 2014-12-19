---
layout: page
title: FAQ
---



### Do I need to filter LoFreq predictions?

You usually don't. Predicted variants are already filtered using
default parameters (which include coverage, strand-bias, snv-quality
etc). If you however need extra conservative settings you can use
`lofreq filter` to add some more filter criteria.

### What if I didn't/can't recalibrate base-qualities?

You can run LoFreq on non-recalibrated data, even though we recommend
to follow GATK's best practices for pre-processing of your BAM file
for Illumina data. Non-recalibrated Illumina data might lead to a
slight increase in false positive calls but we've rarely seen really bad
examples.

### How about recalibrating non-human data?

We recommend GATK base-recalibration even for non-human data (or
targeted sequencing), even though GATK requires the input of known
variant sites (a circular problem actually), which are not known for
many organisms (use dbsnp for human data). One option is to run LoFreq
first and use its predictions as "known" variant input to GATK and
then run LoFreq again. The other alternative is to simply use an empty
vcf file.

### Targeted Resequencing, Exome etc.

If you want to run LoFreq only on certain regions use the appropriate
bed-file as input with `-l regions.bed`. Not doing so will negatively
affect LoFreq's sensitivity (because the automatically applied
Bonferroni correction will be too harsh) and it also might call variants
outside of the desired regions.

### PCR amplified data (amplicons)

If your data was heavily PCR-amplified LoFreq will likely call SNVs in
primer regions, due to ambiguous primer positions, primer impurities
etc. The best way to deal with this is to create a bed-file that only
lists non-primer regions and use that as input to LoFreq with `-l
regions.bed`.

Do not use Picard's Markduplicate on very high coverage data. It will
mark actual reads as duplicates and greatly reduce the effective
sequencing depth.

Another problem with heavily PCRed input is that PCR artifacts will
show up as low-frequency SNVs. Their allele frequency will usually be
low (unless a misamplification happened in early cycles).
Computational tools will be unable to distinguish these from true
low-frequency variants, since the mis-incorporated bases look real to
the sequencing machine. To get rid of these you will either have to
run your samples in duplicates (before amplification) or remove SNVs
below the expected PCR error frequency.

### Choice of mapper: Bowtie and BWA-SW

BWA-SW assigns very low mapping qualities to mapped reads. LoFreq
incorporates mapping quality by default into its variant calling
model, which doesn't help in this scenario and will reduce its
sensitivity. Consider disabling the use of mapping quality for `lofreq
call` with `-N`.

The mapping quality issue also applies to Bowtie, depending on the
parameters you used. Bowtie will sometimes produce alignments with
reads that have either only mapping quality 0 (non unique) or 255 (not
available). It is unclear to us why that happens, or what exactly that
means (see also
[this thread](http://seqanswers.com/forums/showthread.php?t=3142)).
We suggest to switch off the use of mapping quality in such
cases as well (see above).

### Cryptic error messages: `call-parallel` fails

`call-parallel` is really just a wrapper to `lofreq call`.
This is easiest to debug it is to run `lofreq call`
first with the same parameters, check the error message and fix the
parameter choice. Then execute `call-parallel` again with fixed
parameters.

### Why the Star/Asterisk in LoFreq*?


The star is used in programming as a wildcard character matching
any pattern. In LoFreq* it's supposed to match any sequencing assay
([*Seq](http://liorpachter.wordpress.com/2013/08/19/genesis-of-seq/)
as Lior Pachter would call it), because it is generic enough to
be applied to a number of different sequencing assays without change
in parameters.

### The multitude of parameters are confusing

Don't worry: there is rarely need to change the default parameters.
Only change them if following the instructions and recipes given here
or if you know **exactly** what you are doing.

### What's the minimum coverage?

In theory LoFreq doesn't have a minimum (or maximum coverage). The
maximum is restricted by runtime (it gets slow for a coverage of a
million and above). However, at very low coverages sampling biases can
happen which is why we set a minimum coverage of 10 by default.
