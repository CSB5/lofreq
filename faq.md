---
layout: page
title: FAQ
---

### What is LoFreq's software license?

LoFreq's source code is released under [the MIT license](http://opensource.org/licenses/MIT).


### Do I need to filter LoFreq predictions?

Predicted variants are already filtered using
default parameters (which include coverage, strand-bias, snv-quality
etc). If you however need extra conservative settings you can use
`lofreq filter` to add some more filter criteria.


### Troubleshooting missing variants (where did my favourite variant go?)

There are a number of reasons why a variant might not be predicted by LoFreq, 
even though it's "visible" in the pileup. These reasons include low(ish) variant
quality or high strand-bias, high alignment error possibilties etc.

Note that LoFreq by default applies filtering on variant quality and also runs
`lofreq filter` after each run. The latter  removes low coverage variants and
those with strand bias (see settings there).

A simple troubleshooting approach is to rerun LoFreq with very permissive options
on the region of interest. If the variants "reappears" then you can remove the
permissive options one by one and revert back to the default, to see which option
was responsible.

To run LoFreq in a very permissive way use:

    lofreq call --no-default-filter -A -B -a 1 -b 1 -r yourchr:yourstart-yourend -f yourref.fa your.bam

With this `lofreq filter` won't be called (`--no-default-filter`), alignment
qualities are ignored (`-A -B`) and multiple testing correction and p-value
threshold filtering (on variant quality) is switched off (`-a 1 -b 1`).

If the above calls the "missing" variant, check whether the strand bias value
(SB) is high or coverage (DP) is low.  If they are, then the default filter is
responsible for removing that variant. If those values are not problematic, then
remove `-A -B`. If the variant is afterwards not predicted or the variant
quality drops, then there is likely a base or indel-alignment issue that
reduces the variant quality to a value that makes LoFreq discard the variant.


### Should I add customer filters to `lofreq call`?

We would discourage you to filter on base or mapping quality etc. in `lofreq call`.
The reason is that these qualities are built into LoFreq's calling model, i.e. they 
are dealt with properly, which is reflected in the resulting variant quality. 
Using too many base filters can bias results.



### What if I didn't/can't recalibrate base-qualities?

You can run LoFreq on non-recalibrated data, even though we recommend
to follow GATK's best practices for pre-processing of your BAM file
for Illumina data. Non-recalibrated Illumina data might lead to a
slight increase in false positive calls, but having said this we've
rarely seen really bad examples. 

Newer GATK versions don't support recalibration of quality scores anymore.
An easy alternative is to use `lofreq indelqual` and `alnqual` which add
indel qualities (BI/BD) and alignment qualities (for indels and bases). 
The added advantage is that you don't need a dbSNP file to do this, i.e. 
you can use these subcommands on all species.

### Targeted Resequencing, Exome etc.

If you want to run LoFreq only on certain regions use the appropriate
bed-file as input with `-l regions.bed`. Not doing so will negatively
affect LoFreq's sensitivity (because the automatically applied
Bonferroni correction will be too harsh) and it also might call variants
outside of the desired regions (where some reads will map wrongly).

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

And yet another problem with highly PCR amplified is the high strand bias.
There is no one size fits all filtering setting for this.

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

### I can't see any indels

If you're using the `--call-indels` option and don't see any indels
in the output then your BAM file likely didn't contain  indel
qualities. This can be achieved by running GATK's BQSR or `lofreq
indelqual` (use option `--dindel` for Illumina).

### What about Rock 'n' Roll?

[LoFreq rocks](http://www.last.fm/music/Lofreq) and there are [variants](http://www.last.fm/music/Lo+Freq).
