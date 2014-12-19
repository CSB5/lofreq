---
layout: page
title: Commands
---

# [Basics](#basics)

LoFreq comes with a variety of subcommands. By just typing `lofreq`
you will get a [list of available commands]({{page.url}}index.html#cmdlist). The two most
important ones are

- `lofreq call`: simply call variants ([see documentation]({{page.url}}index.html#call))

and

- `lofreq somatic`: call somatic variants  ([see documentation]({{page.url}}index.html#somatic))

Calling lofreq with one of these subcommands and without further
arguments will display  corresponding usage information.

A few things to note:

- We suggest that you preprocess your BAM files by following
  [GATK's best practice protocol](http://www.broadinstitute.org/gatk/guide/best-practices),
  i.e. that you mark duplicates, realign indels and recalibrate base
  qualities with GATK.
- If you are working with exome data or targeted resequencing, don't
  forget to use the corresponding bed-file with `-l region.bed`
  (otherwise the automatically applied Bonferroni correction will
  reduce your sensitivity)
- If you also want to call indels add `--call-indels` to the parameter
  list (note, this requires special preprocessing of your data,
  explained elsewhere).
- Predicted variants are already filtered using default parameters
  (which include coverage, strand-bias, SNV-quality etc). There will
  usually be no need for you to filter the output again.
- Output format for variant calls is [vcf](http://samtools.github.io/hts-specs/VCFv4.1.pdf).

<!-- FIXME preprocessing needs separate article -->



---

# <a name="call">Calling variants: `lofreq call`</a>

Assuming your BAM file (`aln.bam`) contains reads aligned
against the sequence/s in `ref.fa` and you want to predict variants and
save them to `vars.vcf`, you would use the following
 
    lofreq call -f ref.fa -o vars.vcf aln.bam

If you also want to call indels add `--call-indels` to the parameter
list (note, this requires special preprocessing of your data).

<!-- FIXME preprocessing needs separate article-->

If you want to make use of multiple processing use `lofreq
call-parallel` instead of `lofreq call` and add `--pp-threads
THREADS`, where `THREADS` is the number of threads you want to use.
All other parameters stay the same, e.g.:

    lofreq call-parallel --pp-threads 8 -f ref.fa -o vars.vcf aln.bam


If you are dealing with human samples (or large genomes  in general) we
recommend the use if `-s` (source quality) in combination with `-S
dbsnp.vcf.gz` (make sure the dbsnp version matches your reference version) to get rid of
some mapping problems. Source quality is automatically enabled by default in the
somatic SNV calling subcommand.


# <a name="somatic">Calling somatic variants: `lofreq somatic`</a>

Assuming you have your normal reads mapped in a file called `normal.bam` and the tumor
reads in `tumor.bam`, both of which are mapped against `hg19.fa`
then you would use the following to call somatic SNVs, using 8 threads
and store the results to files with the prefix `out_`:

    lofreq somatic --threads 8 -n normal.bam -t tumor.bam -f hg19.fa -o out_ [-d dbsnp.vcf.gz]

The use of dbsnp is optional but highly recommended. It will remove
possibly undetected germline variants from the final output. Ideally
you have remove somatic variants from dbsnp (those matching SAO=2 or SAO=3).
`lofreq somatic` produces several output files. Depending on whether
you also enabled indels and whether you used dbsnp or not, they are
called:

SNVs before and after dbsnp removal:

- `out_somatic_final.snvs.vcf.gz`
- `out_somatic_final_minus-dbsnp.snvs.vcf.gz`

And if you enabled indel calls:

- `out_somatic_final.indels.vcf.gz`
- `out_somatic_final_minus-dbsnp.indels.vcf.gz`


---

# <a name="cmdlist">List of all commands</a>

## Main commands

##### `call`: Call variants

The main command for calling variants. 

##### `call-parallel`: Call variants in parallel

A wrapper around the `call` command that executes several instances of
`lofreq` (use `--pp-threads` to specify number of threads to use).

##### `somatic` : Call somatic variants

Calls somatic variants in matched tumor/normal pairs


---

## Preprocessing commands 

#### `viterbi`: Viterbi realignment

Probabilistic realignment of your already mapped reads, which corrects
mapping errors (run right after mapping).

#### `indelqual`: Insert indel qualities

Inserts indel qualities into your BAM file. Can be used instead of
GATK's BQSR or for non Illumina data.

#### `alnqual`: Insert base and indel alignment qualities

Rarely needed because computed on the fly with `lofreq call`. 


---

## Other commands


##### checkref: Check that reference fasta and BAM file match

Rarely needed if you created the BAM file yourself.

##### `filter`: Filter variants

This will rarely be needed as LoFreq calls this command automatically
with default parameters after predicting variants.


##### `uniq`: Test whether variants predicted really can't be called in other

Variants are sometimes only predicted from one sample, but not a
second, related because of coverage issues, borderline SNV-pvalues
etc. This command will tell you whether a SNV was really not possible
to be called in another. You will rarely need to use this command. It
is built into the somatic SNV calling pipeline. Also note, this is not
designed for very high coverage samples (>>1000X).


##### `plpsummary`: Print pileup summary per position

Mainly useful for debugging.


##### `vcfset`: VCF set operations

vcfset operations like intersection and complement (similar to
`bedtools intersect` and `bedtools subtract` but allele-aware by
default).


##### `version`: Print version info

Print LoFreq version.

---

##  Samtools clones

For convenience only

-    `faidx`         : Create index for fasta file
-    `index`         : Create index for BAM file
-    `idxstats`      : Print stats for indexed BAM file


---
##  Extra tools

These are optionally installed tools, that you might find useful. They
require [PyVCF](https://github.com/jamescasbon/PyVCF),
[Scipy](http://www.scipy.org/) & [Numpy](http://www.numpy.org/) as
well as [Matplotlib](http://matplotlib.org/) installed.
 
##### `vcfplot`: Plot VCF statistics

Plot properties of variants in vcf file

##### `cluster`: Cluster variants in VCF file

Clusters variants based on their frequency confidence into minimal
haplotype groups. Will give a lower haplotype estimate. Note, this is
designed for viral samples.

