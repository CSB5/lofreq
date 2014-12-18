# LoFreq*: A sequence-quality aware, ultra-sensitive variant caller for high-throughput sequencing data

See http://csb5.github.io/lofreq/ for more info.

## Building the source

To build the LoFreq source you will need

- a C compiler (e.g. gcc or clang)
- a Python 2.7 interpreter
- zlib developer files
- a compiled version of samtools (>=1.1) and htslib version >= 1.1
  (Note, [samtools 1.1](http://sourceforge.net/projects/samtools/files/samtools/1.1/samtools-1.1.tar.bz2/download)
  already contains htslib 1.1)

Then, 

- either clone the source or download the current master as zip and unpack
- run `./bootstrap`
- If you get an error like this `configure.ac:72: error: required file './ltmain.sh' not found`,  run libtoolize (or glibtoolize) first and then `bootstrap` again
- `./configure SAMTOOLS=/path-to-samtools HTSLIB=/path-to-htslib [--prefix=custom-path]`
- `make`
- You can already start using lofreq at this stage: see `./bin/lofreq` for help
- To install the package properly type: `make install`.  After that simply type `lofreq` to get more help
- Depending on the used prefix you might need to adjust your PATH (and PYTHONPATH).





