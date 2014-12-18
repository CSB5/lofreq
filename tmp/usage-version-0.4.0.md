---
layout: page
title: Usage for version 0.4.0
---


# Installation #

You will need a C compiler and Python 2.6 or 2.7 (including the developer files, i.e. headers etc). [Download the source for the latest distribution](https://sourceforge.net/projects/lofreq/files/) (using the development version in GIT is not recommended), unpack it and change the working directory to the newly created directory. Assuming you have admin rights, use the following to install LoFreq:

    python setup.py install

If you dont't have admin rights or want to install to a non-standard directory use the `--prefix` flag, e.g.

    python setup.py install --prefix $HOME/local/

and make sure the corresponding installation sub-directory is in your `PYTHONPATH`.

The will list all available options:

    python setup.py --help

Note, that you will also need to install [samtools](http://samtools.sourceforge.net/) to be able to actually use LoFreq.
 
---

# Usage #

The following describes the usage of the most recent version of LoFreq, which is 0.3.2. Older versions of this document are available:

* [usage-version-0.3.1]


## SNV calling with LoFreq ##


LoFreq takes a read mapping as input (see the end of this document for notes on short-read mapping). It's a good idea to be as stringent as possible with your mapping and to recalibrate base-call qualities as well. A simple LoFreq call would look like this:

    lofreq_snpcaller.py -f ref.fa -b mapping.bam -o raw-snv-output-file

If you want to limit the analysis to certain regions and have those described in a bed-file, then add `-l bed-file`. In almost all cases you will want to post-process the predicted SNV calls by applying some filtering criteria, for which you should use `lofreq_filter.py` (see below). 

Please note:

* The default output format ('snp') is a simple csv-file: the 1st column gives the chromosome, 2nd column: SNV position, 3rd column: SNV-type, 4th column: frequency and the 5th column contains additional information about the SNV call. Alternatively, you produce output in [vcf-format](http://vcftools.sourceforge.net/specs.html) (`--format vcf`), however, some LoFreq scripts do not work with this format at the moment. We are currently migrating to vcf as default.
* LoFreq will distinguish between consensus-variants (`type:consensus-var`) and low-frequency variants (`type:low-freq-var`). Consensus variants are majority/consensus changes with respect to the reference and do not have quality values assigned (LoFreq is not meant to be a genotyping program). Low-frequency variants arise from subpopulations or non-dominant variants and by definition have an abundance/frequency of <50%.

### Options for SNV Calling ###

`lofreq_snpcaller.py -h` prints the full help.
Some of the more important options are described in the following:

* *Regions*: If you want to limit the analysis to certain regions then you can pass a [bed-file](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) describing those regions to LoFreq with the option `-l`.  
* *Multiple testing correction*: LoFreq calculates a p-value for each SNV call. Multiple testing correction (Bonferroni) is performed automatically by default (using genome-size by 3; see `--bonf` option). If you expect lots of zero-coverage columns and don't have a bed-file to filter those, then use `--bonf auto-ign-zero-cov` (available from version 0.4.0 onwards)
* *Base-call qualities*:
    * LoFreq uses samtools to read your data. Samtools can correct base-call qualities if indel-errors are likely. This makes sense if you allowed for indels during mapping (default in e.g. [BWA](http://bio-bwa.sourceforge.net/)). LoFreq uses sensitive BAQ (`-E`) by default. You can influence this with the `--baq` option.
    * To ignore any base below a certain quality use the option `-Q`. The default is 3, which is in accordance with Illumina guidelines.
    * If your BAM file does not contain base-call quality values or if they are meaningless, you can try LoFreq's quality agnostic SNV calling module (option: `--lofreq-nq-on --lofreq-q-off`).

---

## SNV filtering ##

Use `lofreq_filter.py` to filter SNV predictions produced by LoFreq. The three most common filter options would be:

* minimum coverage (`--min-cov`) and
* strand-bias (either `--strandbias-bonf` or `--strandbias-holmbonf`; the latter is recommended),
* SNV quality (`--snp-phred`)

An example call would look like this:

    lofreq_filter.py --strandbias-holmbonf --min-cov 10 \
        -i raw-snv-file -o filtered-snv-file

Note that SNV quality filtering is largely unnecessary if you used the automatic Bonferroni correction during SNV calling (which is the default).

 
---

## Somatic SNV Calls / Unique SNV Calls in Sample Pairs ##

SNVs only called in one sample (e.g. cancer) but not in another paired sample (e.g. blood), can either be biologically interesting or simply due to low coverage in one sample. You can use `lofreq_uniq.py` to find out whether a call made only in one sample cannot be simply explained by the low coverage in the other (e.g. blood). `lofreq_uniq.py` takes as minimal input a file listing SNVs predicted in only one sample (see also `lofreq_diff.py`) and the other sample's BAM file.  LoFreq comes with a script that automatically calls SNVs, filters them and finally derives unique SNVs (`lofreq_uniq_pipeline.py`). An example call looks like this:

    lofreq_uniq_pipeline.py --bam1 first.bam --bam2 second.bam \
        --ref ref.fa --bed regions.bed -o output-dir 

This pipeline requires a bed-file describing the regions of interest to calculate a Bonferroni factor automatically. You can derive a template for such a file using `lofreq_regionbed.py`. Output files can be found in `output-dir`.
