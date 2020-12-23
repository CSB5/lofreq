---
layout: page
title: Installation
---


LoFreq comes bundled as two packages: [(i) a source package](#source)
or [(ii) a binary/precompiled package](#binary). See below for specific
instructions for each of the above listed options.

LoFreq is also available via [Bioconda](https://bioconda.github.io/) and [Homebrew](https://docs.brew.sh/Homebrew-on-Linux), which might be the most
convenient installation option for most.

# Prerequisites

From version 2.1.6 onwards LoFreq (`lofreq call-parallel` to be exact) requires `bcftools` to be installed.

# <a name="binary">Installation from Binary Packages</a>

Download the binary LoFreq distribution matching your system (e.g.
Linux or MacOSX) [from the dist folder]({{ site.github.dist }}) on Github and unpack it. LoFreq can
then simply be called with `./lofreq_star-2.1.0/bin/lofreq` (assuming
you downloaded version 2.1.0). You can move the folder anywhere you
like (just preserve its structure) or copy its contents to a
system-wide installation path, e.g. `/usr/local`

    cp -rv ./lofreq_star-2.1.0/* /usr/local/

After that you should be able to call LoFreq by simply typing `lofreq`
(assuming `/usr/local/bin/` is in your `PATH` which is the case on
most systems).




### Binary packages of versions before 2.1

Binary packages before 2.1 were installed differently: `cd` to the
unpacked directory and use the
following to install LoFreq (to `/usr/local/`):

    bash binary_installer.sh 

If you don't have admin rights or want to install LoFreq to a non-standard directory use the `--prefix` argument. e.g.

    bash binary_installer.sh --prefix $HOME/local/

Then follow the instructions given below under
[Installation to a non-standard directory](#prefix).

# <a name="source">Installation from Source</a>


To build the LoFreq source you will need

- a C compiler (e.g. gcc or clang)
- a Python 2.7 or Python 3 interpreter (2.6 won't work for sure)
- zlib developer files
- [HTSlib 1.4 or later](https://github.com/samtools/htslib)

If those requirements are met, download LoFreq's source package
[from the dist folder]({{ site.github.dist }}) on Github,
unpack it and `cd` to the newly created directory. Assuming you have
admin rights, use the following to compile and install LoFreq:

    ./configure --with-htslib=/path-to-htslib 
    make install

If you don't have admin rights or want to install LoFreq to a
non-standard directory use the `--prefix` argument. In that case
you will also have to follow the instructions given below under
[Installation to a non-standard directory](#prefix).

Note, you will need to give the absolute path to samtools and htslib,
not a relative one!

### Installation from Github

If you downloaded
[LoFreq's source from Github]({{ site.github.repo }})
(by cloning or downloading master as zip), you will need to do run
`./bootstrap`  (once) before proceeding with the instructions above.
- If you get an error like this `configure.ac:72: error: required file
'./ltmain.sh' not found`, run libtoolize (or glibtoolize) first and
then `bootstrap` again. After that you can proceed as described above
with the GNU triplejump:

- `configure`
- `make`
- `make install`


# <a name="prefix">Installation to a Non-Standard Directory</a>

If you installed LoFreq to a non-system directory (e.g. your
home-directory), you will have to make sure that the corresponding
bin sub-directory is part of your `PATH` environment
variable.

If you also want to use the non-essential Python tools that come with
LoFreq, please also make sure to add the corresponding
sub-directory to your `PYTHONPATH`.

For example, for the prefix setting mentioned above (i.e.
`$HOME/local/`), you would need the following:

    export PATH=$HOME/local/bin/
    export PYTHONPATH=$HOME/local/lib/python2.7/site-packages/

The first line ensures that you can call LoFreq by simply typing
`lofreq`. The second makes sure that the optionally installed Python tools work.

# Further help

If the documentation above is not sufficient and you have trouble
installing LoFreq, please ask for help on the mailing list or email us
directly (see <a href="{{ site.baseurl }}/contact">Contact</a>).

