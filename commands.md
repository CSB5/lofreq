---
layout: page
title: Commands
---

# [Basics](#basics)

LoFreq comes with a variety of subcommands. By just typing `lofreq`
you will get a [list of available commands](#cmdlist). The two most
important ones are

- `lofreq call`: simply call variants ([see documentation](#call))

and

- `lofreq somatic`: call somatic variants  ([see documentation](#somatic))

Calling lofreq with one of these subcommands and without further
arguments will display  corresponding usage information.

A few things to note:

- For Illumina data, we suggest that you preprocess your BAM files by
  following
  [GATK's best practice protocol](http://www.broadinstitute.org/gatk/guide/best-practices),
  i.e. that you mark duplicates (not for very high coverage data though),
  realign indels and recalibrate base qualities with GATK (BQSR). The
  latter also used to add indel qualities, which is needed for indel
  calling. Newer versions of GATK's BQSR (known for at least 4.1.4) do not add
  indel qualities (BI/BD tags) anymore. Please check with `samtools view your.bam | grep 'B[ID]:Z:'`. If these tags don't show, please run `lofreq indelqual` after
  BQSR (if you want to call indels later).
- If you are working with exome or targeted resequencing data, don't
  forget to provide LoFreq with the corresponding bed-file (`-l
  region.bed`). Otherwise the automatically applied Bonferroni
  correction will reduce your sensitivity
- If you also want to call indels add `--call-indels` to the parameter
  list. Note, this will only work if your BAM file contains indel
  qualities. This will be the case if you ran GATK's BQSR.
  Alternatively use `lofreq indelqual`
- Predicted variants are already filtered using default parameters
  (which include coverage, strand-bias, SNV-quality etc). There will
  usually be no need for you to filter the output again
- Output format for variant calls is [vcf](http://samtools.github.io/hts-specs/VCFv4.1.pdf)

<!-- FIXME preprocessing needs separate article -->



---

# <a name="call">Calling variants: `lofreq call`</a>

Assuming your BAM file (`aln.bam`) contains reads aligned
against the sequence/s in `ref.fa` and you want to predict variants and
save them to `vars.vcf`, you would use the following
 
    lofreq call -f ref.fa -o vars.vcf aln.bam

If you also want to call indels add `--call-indels` to the parameter
list. Note, this requires that your BAM file contains indel qualities
(automatically added by GATK's BQSR or with `lofreq indelqual`)

<!-- FIXME preprocessing needs separate article-->

If you want to make use of multiple processors, simply use `lofreq
call-parallel` instead of `lofreq call` and add `--pp-threads
THREADS`, where `THREADS` is the number of threads you want to use.
All other parameters stay the same, e.g.:

    lofreq call-parallel --pp-threads 8 -f ref.fa -o vars.vcf aln.bam


If you are dealing with human samples (or large genomes  in general) we
recommend the use of `-s` (source quality) in combination with `-S
dbsnp.vcf.gz` (make sure the dbsnp version matches your reference version) to get rid of
some mapping problems. Source quality is automatically enabled by default in the
somatic SNV calling subcommand.


# <a name="somatic">Calling somatic variants: `lofreq somatic`</a>


<!-- FIXME preprocessing needs separate article-->

Assuming you have preprocessed your BAM files nicely, e.g. following
GATK best practices, and your mapped reads of the normal sample are in
a file called `normal.bam` and the tumor reads in `tumor.bam`, both of
which are mapped against `hg19.fa` then you would use the following to
call somatic SNVs, using 8 threads and store the results to files with
the prefix `out_`:

    lofreq somatic -n normal.bam -t tumor.bam -f hg19.fa \
        --threads 8 -o out_ [-d dbsnp.vcf.gz]

The use of dbsnp is optional but **highly** recommended if you are
dealing with human samples. It will help remove possibly undetected
germline variants from the final output. Ideally you have removed
somatic variants from dbsnp (those matching SAO=2 or SAO=3). LoFreq
expects dbsnp to be tabix indexed for fast random access. Indexing can
be achieved by running `bgzip` and `tabix` on the dbsnp vcf file.

`lofreq somatic` produces several output files, most of which you can
ignore. Depending on whether you also enabled indels and whether you
used dbsnp or not, the final output files are called as follows:

For SNVs before and after dbsnp removal (if you used `-o out_`):

- `out_somatic_final.snvs.vcf.gz`
- `out_somatic_final_minus-dbsnp.snvs.vcf.gz`

For indels (if enabled and BAM file was properly preprocessed):

- `out_somatic_final.indels.vcf.gz`
- `out_somatic_final_minus-dbsnp.indels.vcf.gz`

If you need more sensitive calls increase the value for
`--tumor-mtc-alpha` (and `--indel-tumor-mtc-alpha` resp.)

---

# <a name="cmdlist">List of all commands</a>

## Main commands

##### `call`: Call variants

The main command for calling variants

##### `call-parallel`: Call variants in parallel

A wrapper for the `call` command that executes several instances of
`lofreq call`. Use `--pp-threads` to specify number of threads to use

##### `somatic` : Call somatic variants

Calls somatic variants in matched tumor/normal pairs


---

## Preprocessing commands 

#### `viterbi`: Viterbi realignment

Probabilistic realignment of your already mapped reads, which corrects
mapping errors (run right after mapping). Not recommended for
non-Illumina data.

#### `indelqual`: Insert indel qualities

Inserts indel qualities into your BAM file. Can be used instead of
GATK's BQSR or on non-Illumina data. If you have Illumina data and
don't want to use GATK's BQSR then the easist thing is to use 
the `--dindel` option. If you have non-Illumina data and have a good
guess of the indel error rate then use the `--uniform` option which
adds uniform indel qualities.

#### `alnqual`: Insert base and indel alignment qualities

Rarely needed because computed on the fly in `lofreq call`. 


---

## Other commands


##### `checkref`: Check that reference fasta file and BAM file match

Rarely needed if you created the BAM file yourself

##### `filter`: Filter variants

Rarely needed directly as LoFreq calls this command automatically with
default parameters after predicting variants.


##### `uniq`: Test whether variants are unique to one sample

Variants are sometimes only predicted in one sample, but not another.
This could be because of coverage issues, borderline p-values etc.
This command will tell you whether a SNV was really not possible to be
called in another. You will rarely need to use this command directly.
It is built into the somatic SNV calling pipeline. Also note, this is
not designed for very high coverage samples (>1000X).


##### `plpsummary`: Print pileup summary per position

Only useful for debugging.


##### `vcfset`: VCF set operations

Set operations (intersection and complement) on vcf files. Similar to
`bedtools intersect` and `bedtools subtract` but allele-aware (by
default) and with tabix support.


##### `version`: Print version info

Print LoFreq version.

---

##  Samtools clones

Convenience clones:

-    `faidx`: Create index for fasta file
-    `index`: Create index for BAM file
-    `idxstats`: Print stats for indexed BAM file


---
##  Extra Python Tools

These are optionally installed tools. They
require [PyVCF](https://github.com/jamescasbon/PyVCF),
[Scipy](http://www.scipy.org/) & [Numpy](http://www.numpy.org/) as
well as [Matplotlib](http://matplotlib.org/) installed.
 
##### `vcfplot`: Plot VCF statistics

Plot properties of variants in vcf file

##### `cluster`: Cluster variants in VCF file

Clusters variants based on their frequency confidence interval into
minimal haplotype groups. Will give a lower haplotype estimate. Note,
this is designed for viral samples.

