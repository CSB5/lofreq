---
layout: default
title: Home
---


# Introduction #

LoFreq-Star (also known as LoFreq&#42; or LoFreq version 2) is a fast and sensitive variant-caller for inferring single-nucleotide variants (SNVs) from high-throughput sequencing data. It is designed to robustly call low-frequency variants in next-gen sequencing data-sets. LoFreq has been used to call rare variants in viral and bacterial sequencing datasets and can be used to study mitochondrial heteroplasmy and rare somatic mutations in heterogeneous tumors.

LoFreq-Star makes full use of base-call qualities and other sources of errors in next-gen sequencing, e.g. mapping qualities, which are usually ignored by other methods or only used for filtering. It is very sensitive; most notably, it is able to predict variants below the average base-call quality (i.e. sequencing error rate). Each SNV call is assigned a p-value which allows for rigorous false positive control. Even though it uses no approximations or heuristics, it is very efficient due to several runtime optimizations. LoFreq is generic and fast enough to be applied to high-coverage data and large genomes. It takes a minute to analyze Dengue genome sequencing data with nearly 4000X coverage, roughly one hour to call SNVs on a 600X coverage *E.coli* genome and 1.5 hours to run on a 100X coverage human exome dataset.

For more details on the original version of LoFreq see: Andreas Wilm, Pauline Poh Kim Aw, Denis Bertrand, Grace Hui Ting Yeo, Swee Hoe Ong, Chang Hua Wong, Chiea Chuen Khor, Rosemary Petric, Martin Lloyd Hibberd and Niranjan Nagarajan. [LoFreq: A sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets.](http://www.ncbi.nlm.nih.gov/pubmed/23066108) _Nucleic Acids Res._ 2012; 40(22):11189-201.


---

# Citations etc.

- [Full list of publications citing LoFreq](http://scholar.google.com.sg/scholar?oi=bibs&hl=en&cites=12020456701536684432) via Google Scholar
- Coverage in [Featured Research of Science Daily](http://www.sciencedaily.com/releases/2013/03/130307145744.htm)
- [A-Star press release](http://www.research.a-star.edu.sg/research/6661)
- [Pacific Biosciences' minor variants calling protocol](http://files.pacb.com/Training/SMRTAnalysisv22Overview/story.html) is based on LoFreq

---

# [Contact](#contact)

Please contact us if you find bugs, have suggestions, need help etc.
You can either use our
[mailing list](https://sourceforge.net/p/lofreq/mailman/) or contact
us directly:

- [Andreas Wilm](mailto:wilma@gis.a-star.edu.sg)
- [Niranjan Nagarajan](mailto:nagarajann@gis.a-star.edu.sg)

LoFreq is developed in the [Genome Institute of Singapore](http://www.gis.a-star.edu.sg/)

