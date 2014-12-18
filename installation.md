---
layout: page
title: Installation
---


LoFreq comes in two versions: (i) a source package and (ii) a
binary/precompiled package.
Installing the binary package will be the simplest for most users. See below for
specific instructions for each of the above listed options.

# Installation from Binary Packages

Download the binary LoFreq distribution matching your system
[from the Files/Download section](https://sourceforge.net/projects/lofreq/files/)
on sourceforge and unpack it. LoFreq can then simply be called with
`./lofreq_star-2.1.0/bin/lofreq` (assuming you downloaded
version 2.1.0).

You can move the folder anywhere you like (just
preserve its structure) or copy its contents to a system-wide
installation path, e.g. `/usr/local`

    cp -rv ./lofreq_star-2.1.0/* /usr/local/

After that you should be able to call LoFreq by simply typing `lofreq`
(assuming `/usr/local/bin/` is in your `PATH` as should be the case on
most systems).




### Binary packages of versions before 2.1

Binary packages before 2.1 used another way installation. There,
you would `cd` to the newly created directory and (assuming you have
admin rights) use the following to install LoFreq:

    bash binary_installer.sh 

If you don't have admin rights or want to install LoFreq to a non-standard directory use the `--prefix` argument. e.g.

    bash binary_installer.sh --prefix $HOME/local/

and follow the instruction given below under "Installation to a non-standard directory"

# Installation from Source


To build the LoFreq source you will need

- a C compiler (e.g. gcc or clang)
- a Python 2.7 interpreter
- zlib developer files
- a compiled version of samtools (>=1.1) and htslib version >= 1.1
  (Note, [samtools 1.1](http://sourceforge.net/projects/samtools/files/samtools/1.1/samtools-1.1.tar.bz2/download)
  already contains htslib 1.1)

If those requirements are met, download LoFreq's source package
[from the Files/Download section](https://sourceforge.net/projects/lofreq/files/),
unpack it and `cd` to the newly created directory. Assuming you have
admin rights, use the following to compile and install LoFreq:

    ./configure SAMTOOLS=/path-to-samtools HTSLIB=/path-to-htslib
    make install

If you don't have admin rights or want to install LoFreq to a non-standard directory use the `--prefix` argument, e.g.
In that case you will also have to follow the instruction given below under "Installation to a non-standard directory"

Note, you will need to give the absolute path to samtools and htslib,
not a relative one.

### Installation from Github


If you downloaded the source from Github instead of using a source
package, you will need to do the following (once) before proceeding with the
instructions above:

- Either clone the source of download the current master as zip and
  unpack
- run `./bootstrap`
- If you get an error like this `configure.ac:72: error: required file
'./ltmain.sh' not found`, 
run libtoolize (or glibtoolize) first and then `bootstrap` again

Now you can proceed as described above with `configure`, `make` and `make install`.


 Source
 
# Installation to a non-standard directory

If you installed LoFreq to a non-system directory (e.g. your
home-directory), you will have to make sure that the corresponding
installation sub-directory is part of your `PATH` environment
variable. If you also want to use the non-essential Python tools that
come with LoFreq, please make sure to add the corresponding
installation sub-directory to your `PYTHONPATH`. For example, for the
prefix setting mentioned above ($HOME/local/), you would need the
following:

    export PATH=$HOME/local/bin/
    export PYTHONPATH=$HOME/local/lib/python2.7/site-packages/

The first line ensures that you can call LoFreq by simply typing
`lofreq`. The second makes sure that the extra (optional) Python
scripts work.

# Further help

If the documentation above is not sufficient and you have trouble
installing LoFreq, please ask for help on the mailing list or email us
directly (See <a href="{{ site.baseurl }}/contact">Contact</a>).

