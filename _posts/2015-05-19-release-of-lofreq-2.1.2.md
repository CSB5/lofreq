---
layout: post
title: Release of 2.1.2
---

A long overdue release with many, many smaller bug fixes and
improvements  (almost would have justified a version bump to 2.2). 
Main focus were the indel and somatic calling routines. By far the
 biggest visible change is that we got rid of the consensus variants
 (CONSVAR) concept and now assign qualities to all variants (in the
 past  CONSVARs could not be filtered based on quality)

See https://github.com/CSB5/lofreq/blob/master/Changelog for full details.

Andreas
