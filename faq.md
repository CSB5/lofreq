---
layout: page
title: FAQ
---

# FAQ for LoFreq-Star

##### How to Cite

Please cite: [Wilm et a. LoFreq: A sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets. _Nucleic Acids Res._ 2012; 40(22):11189-201.](http://www.ncbi.nlm.nih.gov/pubmed/23066108)

##### Where to get help

Join the LoFreq
[mailing list](https://sourceforge.net/p/lofreq/mailman/) or send us
an email directly to us (see [Contact](/#contact))

##### Do I need to filter LoFreq predictions?

You usually don't. Predicted variants are already filtered using
default parameters (which include coverage, strand-bias, snv-quality
etc). If you however, need extra conservative settings you can use
`lofreq filter` to add some more filter criteria.

##### What happens if i didn't recalibrate base-qualities?

We usually recommend to follow GATK's best practices on
post-processing of you BAM file which include base quality
recalibration. You can run LoFreq on non-recalibrated data, which
might however results in a higher rate of false positive calls. We
would recommend GATK base-recalibration even for non-human data or
targeted sequencing. GATK usually requires the input of known variant
sites (which is a circular problem) which are however not known for
most organisms non-human data (use dbsnp for human data). One option
is to use initial LoFreq predictions as vcf input to GATK and then run
LoFreq again.

##### Regions for Targeted Sequencing

If you want to run LoFreq only on certain regions use the appropriate
bed-file as input with `-l regions.bed`. Not doing so will negatively
affect LoFreq sensitivity and it might calls variants outside of the
desired regions.

##### PCR amplified data (amplicons)

If your data was heavily PCR-amplified LoFreq will likely call SNVs in
primer regions, due to ambiguous primer positions, primer impurities
etc. You should ignore primer positions after SNV calling. The best
way to achieve this is to create a bed-file that only lists non-primer
regions and use that as input to Lofreq with `-l regions.bed`.

Another problem with heavily PCRed input is that PCR artifacts will
show up as low-frequency SNVs. Their allele frequency will usually be
low, unless a mis-amplification happened in early cycles.
Computational tools will be unable to distinguish these from true
low-frequency variants, since the mis-incorporated bases look real to
the sequencing machine. To get rid of these you will either have to
run your samples in duplicates (before amplification) or remove SNVs
below the expected PCR error frequency.

##### I get fewer SNVs than expected when running LoFreq on a BAM file with Bowtie/BWA-SW

BWA-SW assigns very low mapping qualities to mapped reads. Consider
disabling the use of mapping quality on `lofreq call` with `-J`. The
same is true for Bowtie, depending on the parameters you used. Bowtie
it will sometimes produce alignments with reads that have either only
mapping quality 0 (non unique) or 255 (not available). It is unclear
to us why that happens, or what exactly that means (see also
[this thread](http://seqanswers.com/forums/showthread.php?t=3142)).
LoFreq uses mapping quality by default for the prediction of SNVs,
which doesn't help in this scenario. We would suggest to switch off
the use of mapping quality (see above) in such cases.

##### `call-parallel` fails with a cryptic error message 

This is easiest to debug by running the single-threaded version `call`
first with the same parameters, check the error message and fix the
command accordingly. Then execute `call-parallel` again with fixed
parameters.
