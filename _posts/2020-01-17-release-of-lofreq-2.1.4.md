---
layout: post
title: Release of 2.1.4
---

A week ago we published version 2.1.4 of LoFreq. It’s been two and a half years since version 2.1.3 was published. What are the changes and what took so long? And what’s the future of LoFreq?

## Version 2.1.4: keeping up with HTSlib

As [announced on Twitter](https://twitter.com/me_myself_andY/status/1213493490853634049), LoFreq was released early January 2020. The release contains one major new “feature”, which was contributed by the community, namely [John Marshall](https://github.com/jmarshall) from the University of Glasgow. We would like to repeat our thanks here again! [John’s patch](https://github.com/CSB5/lofreq/pull/86) introduced the use of the HTSlib API and got rid of the legacy samtools API. This was long on our list and it’s great to see that a samtools expert introduced these changes. What’s the big deal? The old version of LoFreq didn’t allow us to keep up with HTSlib development. Now, you can link LoFreq against newer HTSlib versions and furthermore CRAM support came for free. The build process also got easier because of this: instead of having to set environment variables SAMTOOLS and HTSLIB for `configure`, you can now use `--with-htslib` and point it, for example, to your conda installed HTSlib installation. This extends LoFreq’s shelf life substantially and makes the life of package maintainers much easier. A belated shout-out at this stage to all folks that made LoFreq available through common package managers over the years:

- The [Debian Med team](https://wiki.debian.org/DebianMed)
- [Bioconda](https://bioconda.github.io/) (John Greeley?)
- [Brew](https://brew.sh/) (Torsten Seemann)

Thanks a lot to all of you.

## Status of LoFreq version 2

We are more than delighted to see that LoFreq is still used by so many (we are hitting 400 citations according to [Microsoft Academic](https://academic.microsoft.com/paper/1964807436/) and Google Scholar). Because of its widespread use, we still receive a fair amount of questions via email and Github issues. Most of these fall into two categories:

1. Questions on how LoFreq works internally
1. Questions around filtering ("missing" variants etc.)

You might ask why questions about the internal workings are necessary. After all the tool has been around for a while and there is a publication. Well yes, but the publication only describes version 1. Version 2 brought substantial changes and additions, for example
the idea that you can join quality (like mapping, alignment and base qualities) into one unified error probability,
a quality-aware Viterbi-based (base) realigner, indel alignment qualities (think BAQ for indels), indel calling, etc. This work was mostly done by Hui Ting Grace Yeo (now at MIT).
We also developed LoFreq Somatic (see e.g. [Ewing et al., 2015](https://www.nature.com/articles/nmeth.3407)), which opened a can of worms (LoFreq’s main development goal was on viral, bacterial and high coverage data in general). LoFreq Somatic is based on a straightforward, subtractive approach that sort of works, but at least for human whole genome data there are much better callers out there to be fair.
We furthermore used LoFreq on long read data (PacBio), showing how its core idea and quality awareness hold for other types of data. The draft of the paper unfortunately remained a draft, mostly because I (Andreas) preferred starting new things over finishing old ones and eventually I changed jobs.

The second set of questions (the ones around filtering) are more difficult to deal with. We built default (variant) filtering in, because one of the earliest papers reviewing the earliest version of LoFreq didn’t use any filtering. So we forced this upon users. But even though the defaults work well under most circumstances, they don't always work (especially when it comes to strand bias) and additionally users like to tinker with base quality filtering settings for some reason, which makes things worse. Some of the filtering routines have side effects or are [buggy](https://github.com/CSB5/lofreq/issues/80) under certain circumstances and they are hard to change. Having said all this: the code works just fine as-is, but it has acquired a lot of technical depth.

Please keep in mind: you shouldn’t apply excessive base quality filtering for LoFreq in general, because it’s designed to be quality aware and build from the ground up to deal with base quality variations. Base quality filtering will bias results!

## LoFreq version 3

When [Filip Sodic](https://github.com/sodic), a highly capable student, joined as an intern, we started yet another experiment: rewriting the code  (a mix of Python and C) from scratch in a modern, efficient, expressive and elegant language ([Nim](https://nim-lang.org/)!). [LoFreq version 3](https://github.com/andreas-wilm/lofreq3) was born, but I (Andreas) left academia before the reimplementation was completed. LoFreq 3 splits pileup (quality joining) entirely from variant calling, which avoids all sorts of problems and makes the process more transparent. The current version can already be used as drop-in replacement for version 2's `lofreq call`. A release with binary will be published within weeks. In the meantime you can build from source.

##  Going forward

We will continue to support LoFreq version 2 as best as we can (answering questions, fixing critical bugs etc.). Since I (Andreas) left academia, I can only spend so much time on this, and most of it will go into the development of version 3. As mentioned, this already functions as a drop-in replacement for `lofreq call`. Implementations of other sub-commands, like `viterbi`, will follow.

And we have not yet given up on the idea to publish the paper! If you are or know a talented student, consider joining [Niranjan’s lab](http://csb5.github.io/) to work on this (and other cool things :))!
