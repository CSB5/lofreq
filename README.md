# LoFreq*: A sequence-quality aware, ultra-sensitive variant caller for high-throughput sequencing data

See https://sourceforge.net/p/lofreq/wiki/LoFreq-Star/ for more info

## Installing/Compiling a repository clone

After cloning the repository, run `bootstrap`, which will set up things such
that you can use `./configure`

After that you can proceed as in 'Installing/Compiling the source distribution'

## Installing/Compiling the source distribution

Download and compile samtools and htslib version >= 1.1, e.g. [samtools version
1.1 (which also contains htslib)](http://sourceforge.net/projects/samtools/files/samtools/1.1/samtools-1.1.tar.bz2/download)

Then compile LoFreq as follows:
```
./configure SAMTOOLS=/full-path/to/samtools HTSLIB=/full-path/to/htslib [--prefix instdir]
make install
```

Depending on the used prefix you might need to adjust your PATH (and PYTHONPATH).

After that simply type `lofreq` to get more help


