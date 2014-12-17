---
layout: page
title: Usage for version 0.5 to 0.6.1
---



# Introduction #

LoFreq is a fast and sensitive variant-caller for inferring single-nucleotide variants (SNVs) from high-throughput sequencing data. It is designed to robustly call low-frequency variants by exploiting base-call quality values. LoFreq has been used to call rare variants in viral and bacterial sequencing datasets and can be used to study mitochondrial heteroplasmy and rare somatic mutations in heterogeneous tumors.

LoFreq makes full use of base-call qualities (and versions >=0.5 also use read mapping qualities) which are usually ignored by other methods or only used for filtering. It is very sensitive; most notably, it is able to predict variants below the average base-call quality (i.e. sequencing error rate). Each SNV call is assigned a p-value which allows for rigorous false positive control. Even though it uses no approximations or heuristics, it is very efficient due to several runtime optimizations. LoFreq is generic and fast enough to be applied to high-coverage data and large genomes. It takes a minute to analyze Dengue genome sequencing data with nearly 4000X coverage, roughly one hour to call SNVs on a 600X coverage *E.coli* genome and 1.5 hours to run on a 100X coverage human exome dataset.

For more details see: Andreas Wilm, Pauline Poh Kim Aw, Denis Bertrand, Grace Hui Ting Yeo,
Swee Hoe Ong, Chang Hua Wong, Chiea Chuen Khor, Rosemary Petric, Martin Lloyd Hibberd and Niranjan Nagarajan. [LoFreq: A sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets.](http://www.ncbi.nlm.nih.gov/pubmed/23066108) _Nucleic Acids Res._ 2012; 40(22):11189-201.

---

# Installation #

You will need a C compiler, Python 2.7 (including the developer files, i.e. headers etc.) and the zlib developer files (all of which are very likely already installed on your system). [Download the source for the latest LoFreq distribution](https://sourceforge.net/projects/lofreq/files/) (using the development version in GIT is not recommended), unpack it and change the working directory to the newly created directory. Assuming you have admin rights, use the following to compile and install LoFreq:

    ./configure
    make install

If you dont't have admin rights or want to install LoFreq to a non-standard directory use the `--prefix` flag to configure, e.g.

    ./configure --prefix $HOME/local/
    make install

If you installed LoFreq to a non-system directory (e.g. your home-directory), you will have to make sure that the corresponding installation sub-directory is part of your `PATH` and `PYTHONPATH` environment variables. The corresponding changes needed for the prefix setting mentioned above would be:

    export PATH=$HOME/local/bin
    export PYTHONPATH=$HOME/local/lib/python2.7/site-packages/

The first line ensures that you can call the LoFreq scripts simply by typing their names (do not call LoFreq by using `python full-path-to-installation/lofreq_snpcaller.py`; simply use `lofreq_snpcaller.py`). The second makes sure that the Python libraries are found. 

If you need special compiler or linker options you can pass them either to configure or make by simply appending the following to either configure or make:

    CFLAGS="your-c-flags" LDFLAGS="your-ld-flags"

For example, on Mac systems using Xcode you might have to run the following to get LoFreq to comile properly

    ./configure CFLAGS="-arch x86_64" LDFLAGS="-arch x86_64"


## Troubleshooting

Many users don't have admin rights, which can make the installation of software sometimes a bit troublesome. This is mainly because binaries and libraries are installed to non-standard directories, which your system might not be set up to use. Please feel free to contact us if you have problems with the installation. In the following which compiled a list of more common problems and their solutions.

### Problem 1: PATH not updated
You get the following error:

    lofreq_snpcaller.py: command not found

`Reason and Solution`: LoFreq was installed into a directory which is not part of your PATH environment variable. Add the installation directory, or more exactly the bin directory to your PATH

### Problem 2: PYTHONPATH not updated
You get the following error:

    Traceback (most recent call last):
      File "./lofreq_snpcaller.py", line 47, in <module>
        from lofreq import conf
    ImportError: No module named lofreq

`Reason and Solution`: LoFreq can't find the installed libraries, because the installation directory, or more exactly the python-lib subdirectory, is not part of your PYTHONPATH environment variable. Add it to PYTHONPATH and the problem should be fixed.

### Problem 3: PYTHONPATH pointing to old installation
You get the following error:

    File "some-prefix/bin/lofreq_snpcaller.py", line 685, in main
        opts.baq, opts.max_depth, opts.region_bed, opts.join_mapq_and_baseq)
    TypeError: generate_pileup() takes at most 4 arguments (5 given)

`Reason and Solution`: You have an old installation of LoFreq and your PYTHONPATH still points to it. Remove the old installation directory from PYTHONPATH and add the new one.

---

# Usage #

The section below describes the most recent version of LoFreq, which is 0.5.0. The following older versions of this document are available:

* [usage-version-0.3.1]
* [usage-version-0.4.0]


## SNV calling with LoFreq ##


LoFreq takes a read mapping in BAM format as input. It's a good idea to be as stringent as possible with your mapping, to realign around indels and recalibrate base-call qualities as well. See the end of this document for recipes and scripts to achieve this (section 'Best practices...'). A simple LoFreq call would look like this:

    lofreq_snpcaller.py -f ref.fa -b mapping.bam -o raw-snv-output-file

If you want to limit the analysis to certain regions and have those described in a bed-file, then add `-l bed-file`. In almost all cases you will want to **post-process the predicted SNV calls by applying some filtering criteria**, for which you should use `lofreq_filter.py` (see section on filtering below). 

Newer versions of LoFreq will incorporate mapping qualities. This reduces the number of false positive calls on longer genomes. If you are analyzing very small genomes with no repetitive regions (e.g. RNA viruses) you can assume that there is no mapping uncertainty and increase LoFreq's sensitivity by using the option `--dont-join-mapq-and-baseq` for `lofreq_snpcaller.py`.

Please note:

* The default output format ('snp') is a simple csv-file: the 1st column gives the chromosome, 2nd column: SNV position, 3rd column: SNV-type, 4th column: frequency and the 5th column contains additional information about the SNV call. Alternatively, you produce output in [vcf-format](http://vcftools.sourceforge.net/specs.html) (`--format vcf`). However, some LoFreq scripts - most importantly `lofreq_filter.py` - do not work with this format at the moment. We are currently migrating to vcf as default.
* LoFreq will distinguish between consensus-variants (`type:consensus-var`) and low-frequency variants (`type:low-freq-var`). Consensus variants are majority/consensus changes with respect to the reference and do not have quality values assigned (LoFreq is not designed to be a genotype caller). Low-frequency variants arise from subpopulations or non-dominant variants and by definition have an abundance/frequency of <50%. These are the ones LoFreq was designed to predict.

### Options for SNV Calling ###

`lofreq_snpcaller.py -h` prints the full help if you need more advanced control.
Some of the more important options are described in the following:

* *Regions*: If you want to limit the analysis to certain regions then you can pass a [bed-file](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) describing those regions to LoFreq with the option `-l`. This is useful for example for Exome sequencing. 
* *Multiple testing correction*: LoFreq calculates a p-value for each SNV call. Multiple testing correction (Bonferroni) is performed automatically by default (see `--bonf` option), so there should be no need for phred-value/SNP-quality based filtering afterwards.
* If you get the error "Pileup was empty" and you are sure that the reads in your BAM are fine, you can try the option `-A` or `--anomalous-pairs-allowed` (analogous to `samtools mpileup -A`; supported from version 0.6.1 on), which takes anomalous reads into consideration, i.e. where a mate-pair is not mapped. Although generally not recommended there might be situations where this is needed (e.g. the mate-pair mapping region is known to be heavily mutated).
* *Base-call qualities*:
    * To ignore any base below a certain quality use the option `-Q`. The default is 3, which is in accordance with Illumina guidelines.
    * If your BAM file does not contain base-call quality values or if they are meaningless, you can try LoFreq's quality agnostic SNV calling module (option: `--lofreq-nq-on --lofreq-q-off`).
    * LoFreq makes use of samtools' base-call quality correction (BAQ), which is useful if indel-errors are likely. LoFreq uses sensitive BAQ (`-E`) by default. You can influence this with the `--baq` option (not recommended, unless you know exactly what you are doing).

 
LoFreq versions >= 0.6.0 have a new option `-c` or `--cons-as-ref`. This controls what base LoFreq uses as reference base. Normally this will come from the reference sequence provided by the user. If you however want to compute your own consensus per position and are only interested in low-frequency mutations (by definition <50%) use `--cons-as-ref`.

---

## SNV filtering ##

Use `lofreq_filter.py` to filter SNV predictions produced by LoFreq. The two most highly recommended filter options are:

* minimum coverage (`--min-cov`) and
* strand-bias (either `--strandbias-bonf` or `--strandbias-holmbonf`; the latter is recommended),

An example call with recommended settings would look like this:

    lofreq_filter.py --strandbias-holmbonf --min-cov 10 \
        -i raw-snv-file -o filtered-snv-file

Note that SNV quality filtering (`--snp-phred`) is largely unnecessary if the default and automatic Bonferroni correction was used during SNV calling.

 
---

## Somatic SNV Calls / Unique SNV Calls in Sample Pairs ##

SNVs only called in one sample (e.g. cancer) but not in another paired sample (e.g. blood), can either be biologically interesting or simply due to low coverage in one sample. You can use `lofreq_uniq.py` to find out whether a call made only in one sample cannot be simply explained by the low coverage in the other (e.g. blood). `lofreq_uniq.py` takes as minimal input a file listing SNVs predicted in only one sample (see also `lofreq_diff.py`) and the other sample's BAM file.  LoFreq comes with a script that automatically calls SNVs, filters them and finally derives unique SNVs (`lofreq_uniq_pipeline.py`). An example call looks like this:

    lofreq_uniq_pipeline.py --bam1 first.bam --bam2 second.bam \
        --ref ref.fa --bed regions.bed -o output-dir 

This pipeline requires a bed-file describing the regions of interest to calculate a Bonferroni factor automatically. You can derive a template for such a file using `lofreq_regionbed.py`. Output files can be found in `output-dir`.

Note, this is a add-on that we did not thoroughly benchmark.


---

# Best practices for creating a BAM file for LoFreq ##

Especially for low frequency SNV calling it's best to use very stringent mapping criteria. Only keep properly aligned reads and only allow uniquely mapped reads. LoFreq comes with a script called `bwa_unique.sh` which will help you to create a unique mapping of your reads to a genome with [BWA (Lee & Durbin, 2009)](http://bio-bwa.sourceforge.net/) 

We also highly recommend to perform a base-call quality calibration on your input. Spurious SNVs will otherwise be likely. For recalibration you can use [GATK (McKenna et al., 2010)](http://www.broadinstitute.org/gatk/). Since the GATK's usage is rather cumbersome especially for non-human data, LoFreq provides you with a wrapper script called `base_qual_calib_wrapper.sh` (based on GATK version 2). This will create a fake-vcf file of 'known' SNVs (needed for GATK) and execute the necessary GATK2 commands for recalibration. For human data you are better off providing the script with a 'real' vcf file from e.g. dbSNP. It also makes sense to use GATK's indel realigner (before base call quality recalibration). Note, that this step is not yet covered by the above mentioned script.

---

### Caveats ###

If your data was heavily PCR-amplified LoFreq will likely call SNVs in the primer regions as well, due to ambiguous primer positions, primer impurities etc. You should ignore primer positions after SNV calling.

Another problem with heavily PCRed input is that PCR artifacts might show up as low-frequency SNVs, especially if a mis-amplification happened in early cycles. Computational tools typically will be unable to distinguish these from true low-frequency variants, since the mis-incorporated bases look real to the sequencing machine.

---


### Related Tools ###

* [Breseq](http://barricklab.org/breseq): full-blown pipeline of which polymorphism-prediction is just a part
* [SNVer](http://snver.sourceforge.net/): low frequency SNV calling based on a frequentist approach
* [V-Phaser](http://www.broadinstitute.org/scientific-community/science/projects/viral-genomics/v-phaser): based on similar ideas. Uses phasing for enhanced sensitivity. Focus on viral 454 data.

---


### About ###

* Publication: Please cite [Wilm et al. LoFreq: A sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets. _Nucleic Acids Res._ 2012; 40(22):11189-201](http://www.ncbi.nlm.nih.gov/pubmed/23066108).


* LoFreq was developed in the [Genome Institute of Singapore](http://www.gis.a-star.edu.sg/)
* Sourceforge Admins: 
[[project_admins]]


Please feel free to contact us if you find bugs, have suggestions, need help etc. Use the discussion forum, the mailing-list or simply mail us directly.
