# LoFreq*: A sequence-quality aware, ultra-sensitive variant caller for NGS data

## Note

Most users will want to install LoFreq via [BioConda](https://bioconda.github.io/)!
The source code hosted here on Github is mainly for developers!

LoFreq was published 10 years ago, is considered very stable and has
almost 1000 citations at the time of writing this.  I (Andreas) have
long left academia and tried to maintain the code as best as I could
until now, but find that I have no more time to do this. For bugs,
suggestions, ideas, collaborations please contact [Niranjan Nagarajan](mailto:nagarajann@gis.a-star.edu.sg).


## Building the Source

### Current Build Status

[![Build Status](https://travis-ci.org/CSB5/lofreq.svg?branch=master)](https://travis-ci.org/CSB5/lofreq)

### Prerequisites

You will need:

- a C compiler (e.g. gcc or clang)
- a Python 3 interpreter
- zlib developer filesi (zlib1g-dev on Ubuntu)
- a compiled version of [HTSlib 1.4 or later](https://github.com/samtools/htslib)

### Compilation

- Clone the repo (or download the current master as zip package and unpack)
- Run `./bootstrap` to set up the required automake files
  - If you get an error like `required file './ltmain.sh'
    not found`, run `libtoolize` (or `glibtoolize`) first and then
    `bootstrap` again
  - Subsequent pulls won't require rerunning `./bootstrap`. This is
    only necesary when changing `configure.ac` or any of the `Makefile.am`
- Run `./configure` with the **absolute** path to HTSlib
  (e.g. `./configure --with-htslib=/home/user/miniconda [--prefix=inst-path]`)
- Run `make`
  - At this point you can already start using lofreq: `./bin/lofreq`
- Run `make install` to properly install the package
  - Default is `/usr/local/`. If `--prefix` was given to `configure`,
    the corresponding argument is used
  - Depending on the used prefix you might need to adjust your PATH (and PYTHONPATH).

## Documentation

- Simply calling `lofreq` on the command line will display a list of
subcommands
- `lofreq cmd` will then display help for `cmd`
- See [LoFreq's website](http://csb5.github.io/lofreq/) for full documentation

## License

LoFreq is licensed under the MIT License (see LICENSE).

Licenses for third party software that is part of the source:
- cdflib90 (see src/cdflib90.README)
- uthash (see src/uthash/LICENSE)
