---
layout: post
title: Release of LoFreq 2.1
---

Users rejoice: we've just released LoFreq 2.1!


The most important changes are the following:

- LoFreq can now call indels. Indel calling depends on a good alignment (use BWA-MEM and refine with `lofreq viterbi`), indel qualities (use `lofreq indelqual` or GATK's BQSR) and indel alignment qualities (computed internally). The feature is off by default, because it's still considered a bit suboptimal. To enable it, use `--call-indels` as extra argument for `lofreq call` and `lofreq somatic`.
- Parameters for somatic calls are now thoroughly tested 
- Release of MIT licensed source code
- LoFreq is now compiled against an external version of samtools and htslib (1.1) 
- bgzip and tabix support (all output from `somatic` is now bgzipped)
- Run-time issue with bed-files containing thousands of regions (e.g. human exome) has been resolved
- Lots of other tiny bug-fixes and improvements
 
This release comes as a binary package for Linux and MacOSX and now
also with MIT-licensed source-code. Please have a look at the
[files section](https://sourceforge.net/projects/lofreq/files/) on sourceforge.


Usage has changed slightly, but should be obvious from the commandline help of LoFreq's subcommands. We will be updating the online documentation soon after we've moved source code and website to github and github-pages (will be announced here) so check back soon!

Andreas


