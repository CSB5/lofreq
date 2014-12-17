---
layout: page
title: Usage for version 0.3.1
---

## SNV calling with LoFreq ##


LoFreq takes a read mapping as input (see the end of this document for notes on short-read mapping). It's a good idea to be as stringent as possible with your mapping and to recalibrate base-call qualities as well. A simple LoFreq call would look like this:

    lofreq_snpcaller.py -f ref.fa -b mapping.bam -o raw-snv-output-file

In almost all cases you will want to post-process the predicted SNV calls by applying some filtering criteria, for which you should use `lofreq_filter.py` (see below). 

Please note:

* We highly recommend to set a Bonferroni correction factor, which allows you to control the false positive rate, will speed the computation up as well and makes subsequent quality-based SNV-filtering obsolete. Future versions will do this automatically. See *Multiple testing correction* under LoFreq Options.
* The default output format ('snp') is a simple csv-file: the 1st column gives the chromosome, 2nd column: SNV position, 3rd column: SNV-type, 4th column: frequency and the 5th column contains additional information about the SNV call. Alternatively, you can get output in vcf-format (`--format vcf`), however, downstream scripts in LoFreq do not work with this format at the moment. We are currently migrating to vcf as default.
* LoFreq will distinguish between consensus-variants (`type:consensus-var`) and low-frequency variants (`type:low-freq-var`). Consensus variants are majority/consensus changes with respect to the reference and do not have quality values assigned (LoFreq is not meant to be a genotyping program).


### LoFreq Options ###

`lofreq_snpcaller.py -h` will print the full help.
Some of the more important options are described in the following:

* *Multiple testing correction*: LoFreq calculates a p-value for each SNV call. To perform multiple testing correction on the fly (which will speed calculations for ultra-high coverage data) you can change the Bonferroni-factor (`--bonf`). A conservative setting would be genome-size multiplied by three. You can use the helper script `lofreq_bonf.py` to compute this value automatically.
* *Regions*: If you want to limit the analysis to certain regions then you can pass a [bed-file](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) describing those regions to LoFreq with the option `-l`.  
* *Base-call qualities*:
    * LoFreq uses samtools to read your data. Samtools can correct base-call qualities if indel-errors are likely. This makes sense if you allowed for indels during mapping (default in e.g. [BWA](http://bio-bwa.sourceforge.net/)). LoFreq uses sensitive BAQ (`-E`) by default. You can influence this with the `--baq` option.
    * To ignore any base below a certain quality use the option `-Q`. The default is 3, which is in accordance with Illumina guidelines.
    * If your BAM file does not contain base-call quality values or if they are meaningless, you can try LoFreq's quality agnostic SNV calling module (option: `--lofreq-nq-on --lofreq-q-off`).

---

## SNV filtering ##

Use `lofreq_filter.py` to filter SNV predictions produced by LoFreq. The three most common filter options would be

* minimum coverage (`--min-cov`) and
* SNV quality (`--snp-phred`; unnecessary if you used automatic Bonferroni above)
* strand-bias (recommended: `--strandbias-holmbonf`),

An example call looks like this:

    lofreq_filter.py --strandbias-holmbonf --min-cov 10 --snp-phred 60 \
        -i raw-snv-file -o filtered-snv-file
 
---

## Somatic SNV Calls / Unique SNV Calls in Sample Pairs ##

SNVs only called in one sample (e.g. cancer) but not in another paired sample (e.g. blood), can either be biologically interesting or simply be due to low coverage in one sample. You can use `lofreq_uniq.py` to find out whether a call made only in one sample cannot be simply explained by the low coverage in the other (e.g. blood). `lofreq_uniq.py` takes as minimal input a file listing SNVs predicted in only one sample (see also `lofreq_diff.py`) and the other sample's BAM file.  LoFreq comes with a script that automatically calls SNVs, filters them and finally derives unique SNVs (`lofreq_uniq_pipeline.py`). An example call looks like this:

    lofreq_uniq_pipeline.py --bam1 first.bam --bam2 second.bam \
        --ref ref-fasta --bed regions-bedfile -o output-dir 

This pipeline requires a bed-file describing the regions of interest to calculate a Bonferroni factor automatically. You can derive a template for such a file using `lofreq_regionbed.py`. Output files can be found in `output-dir`.
