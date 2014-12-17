---
layout: post
title: Performance issues when using bed-file with many regions
---
We noticed performance issues with LoFreq 2.0.0 when running the somatic and
call-parallel sub-command with bed-files containing many (i.e. thousands of)
regions, which is for example the case in human exome analysis. This will be
addressed in the upcoming release.
Andreas
PS: The next release will also come with source-code.
