---
layout: post
title: GATK recalibration problems&#58; low coverage and few SNVs
---
We've had reports about very strange results mostly on viral genomes. In such
cases very few SNVs will be reported, with odd frequencies and very low
coverage. We could trace this back to GATK's base call recalibration step, in
which for yet unknown reasons, most base call qualities are set to the lowest
possible value by GATK. We are investigating this and ask users to run LoFreq
on the uncalibrated data in such cases (keeping in mind that spurious SNV calls
are possible).
