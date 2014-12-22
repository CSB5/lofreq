---
layout: page
title: Installation
---


LoFreq comes bundled as two packages: [(i) a source package](#source)
or [(ii) a binary/precompiled package](#binary). Installing the binary
package will be the simplest for most users. See below for specific
instructions for each of the above listed options.

# <a name="binary">Installation from Binary Packages</a>

Download the binary LoFreq distribution matching your system (e.g.
Linux or MacOSX) [from the Files/Download section]({{ site.sourceforge.download }}) on sourceforge and unpack it. LoFreq can
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
- a Python 2.7 interpreter (no, 2.6 won't work and neither would >=3.0)
- zlib developer files
- a compiled version of samtools (>=1.1) and htslib version >= 1.1
  (note: [samtools 1.1](http://sourceforge.net/projects/samtools/files/samtools/1.1/samtools-1.1.tar.bz2/download)
  already contains htslib 1.1)

If those requirements are met, download LoFreq's source package
[from the Files/Download section]({{ site.sourceforge.download }}),
unpack it and `cd` to the newly created directory. Assuming you have
admin rights, use the following to compile and install LoFreq:

    ./configure SAMTOOLS=/path-to-samtools HTSLIB=/path-to-htslib
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

