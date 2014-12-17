---
layout: page
title: Usage
---

# Usage for LoFreq-Star 


---
## Basics 

LoFreq comes with a variety of subcommands, all which are accessed by simply calling lofreq itself. By just typing `lofreq` you will get a list of available commands. 

The two most important commands are

- `lofreq call` or `lofreq call-parallel`: for simply calling variants

and

- `lofreq somatic`: for calling somatic variants in paired normal/tumor samples

Calling lofreq with one of these subcommands and without further arguments will display usage information for that command.

A few things to note:

- For target sequencing or exomes, you should restrict the analysis to the corresponding regions with `-l regions.bed`.
- Predicted variants are already filtered using default parameters (which include coverage, strand-bias, SNV-quality etc). There will usually be no need for you to filter the output again.
- We strongly suggest that you pre-process your BAM files by following [GATK's best practice protocol](http://www.broadinstitute.org/gatk/guide/best-practices), i.e. that you realign indels and recalibrate base qualities with GATK.
- Output format for variant calls is [vcf](http://samtools.github.io/hts-specs/VCFv4.1.pdf).

---
## Examples


##### Calling variants:

Assuming your BAM file (aligned_reads.bam) contains reads aligned against the sequence/s in ref.fa and you want to predict variants and save them to vars.vcf, you would use the following
 
    lofreq call -f refe.fa aligned_reads.bam -o vars.vcf

If you are dealing with human samples (or large genomes in general) we recommend the use if `-S` (source quality) in combination with `-V dbsnp.vcf` to get rid of some mapping problems (source quality is automatically used in the somatic SNV calling subcommand).


##### Calling somatic variants:

Assuming you have your normal reads mapped in normal.bam and the tumor reads in tumor.bam, both of which are mapped against hg19 (hg19.fa) then you would use the following to call somatic SNVs, using 8 threads and store the results to files with the prefix `somatic`:

    lofreq somatic --threads 8 -n normal.bam -t tumor.bam -f hg19.fa -o somatic [-d dbsnp.vcf.gz]

The use of dbsnp is optional but recommended. It will remove possibly undetected germline variants from the final output.
    
If you are working with Exome data, don't forget to provide the somatic command with the corresponding bed-file (`-l region.bed`)

---
## Main commands 

##### `call`: Call variants

The main command for calling variants. 

##### `call-parallel`: Parallel calling of  variants

A wrapper around the call command that executes several instances of lofreq by working on multiple regions. 

##### `somatic` : Call somatic variants in matched tumor/normal pairs

Runs the somatic SNV calling pipeline.


---
### Other commands 


You will rarely need to run any of these commands. Most of them are either automatically used by lofreq itself.


##### `filter`: Filter variants in (LoFreq) VCF file

This will rarely be needed as LoFreq calls this command automatically with default parameters after predicting variants.


##### `uniq`: Test whether variants predicted really can't be called in other

Variants are sometimes only predicted from one sample, but not a second, related because of coverage issues, borderline SNV-pvalues etc. This command will tell you whether a SNV was really not possible to be called in another. You will rarely need to use this command. It is built into the somatic SNV calling pipeline. Also note, this is not designed for very high coverage samples (>>1000X).


##### `plpsummary`: Print pileup summary per position

Mainly useful for debugging.


##### `vcfset`: VCF set operations
     
vcfset operations like intersection and complement (similar to `bedtools intersect` and `bedtools subtract` but base-aware by default).


##### `version`: Print version info


---
##  Python tools

These are optionally installed tools, that you might find useful. They require [PyVCF](https://github.com/jamescasbon/PyVCF), [Scipy](http://www.scipy.org/) & [Numpy](http://www.numpy.org/) as well as [Matplotlib](http://matplotlib.org/) installed.
 
##### `vcfplot`: Plot VCF statistics

Summarize properties of variants listed in vcf file. Part of the optional LoFreq Python tools. 

##### `cluster`: Cluster variants in VCF file

Clusters variants based on their frequency confidence into minimal haplotype groups. Will give a lower haplotype estimate. Note, this  is designed for viral samples.



---
##  Samtools clones:

Convenience clones of samtools subcommands.

##### `index`: Create index for BAM file

A clone of samtools index

##### `idxstats`: Print stats for indexed BAM file

A clone of samtools idxstats
