# LoFreq*: A sequence-quality aware, ultra-sensitive variant caller for NGS data



## Note

Most users will want to use either the binary or the
source-code package which is distributed via [LoFreq's
Sourceforge site](https://sourceforge.net/projects/lofreq/files/).
The source hosted here on github is mainly for developers.

## Status

[![Build Status](https://travis-ci.org/CSB5/lofreq.svg?branch=master)](https://travis-ci.org/CSB5/lofreq)


## Building the source

To build the LoFreq source you will need

- a C compiler (we used gcc and clang routinely)
- a Python 2.7 interpreter
- zlib developer files
- a compiled version of samtools (>=1.1)
- a compiled version of htslib (>= 1.1; usually part of the above; see
 e.g. [here](http://sourceforge.net/projects/samtools/files/samtools/1.1/samtools-1.1.tar.bz2/download)

Then, 

- either clone the repot or download the current master as zip packageand unpack
- run `./bootstrap` which will set up the required automake files
  - If you get an error like this `configure.ac:72: error: required file './ltmain.sh' not found`,  run libtoolize (or glibtoolize) first and then `bootstrap` again
- Run `./configure` with the **absolute** path to samtools and htslib: `./configure SAMTOOLS=/path-to-samtools HTSLIB=/path-to-htslib [--prefix=custom-path]`
- Run `make`
- At this point you can already start using lofreq: see `./bin/lofreq` for help
- To properly install the package type: `make install`.
- Depending on the used prefix you might need to adjust your PATH (and PYTHONPATH).


## Using LoFreq

See http://csb5.github.io/lofreq/ for full documentation.

On the command line simply type `lofreq` to get help.




