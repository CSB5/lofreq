---
layout: page
title: Installation
---

# Installation of LoFreq-Star

LoFreq comes in two versions: (i) a source package and (ii) a binary/precompiled package. Installing the binary package will be easier in most cases and has fewer dependencies. See below for specific instructions for each.

## Binary Installation

Download the binary LoFreq distribution matching your system [from the Files/Download section](https://sourceforge.net/projects/lofreq/files/), unpack it and cd to the newly created directory. Assuming you have admin rights, use the following to install LoFreq:

    bash binary_installer.sh 

If you don't have admin rights or want to install LoFreq to a non-standard directory use the `--prefix` argument. e.g.

    bash binary_installer.sh --prefix $HOME/local/

and follow the instruction given below under "Installation to a non-standard directory"

## Installation from Source

Please note: we are planning to release the code here within the next couple of months, i.e. it is not available yet. If you think you really need it now, please send us an email.

You will need a C compiler, a Python 2.7 interpreter as well as the zlib developer files installed (all of which are likely already present on your system). Download the source package of the latest LoFreq distribution [from the Files/Download section](https://sourceforge.net/projects/lofreq/files/), unpack it and cd to the newly created directory. Assuming you have admin rights, use the following to compile and install LoFreq:

    ./configure
    make install

If you don't have admin rights or want to install LoFreq to a non-standard directory use the `--prefix` argument, e.g.

    ./configure --prefix $HOME/local/
    make install

In that case you will also have to follow the instruction given below under "Installation to a non-standard directory"


# Installation to a non-standard directory

If you installed LoFreq to a non-system directory (e.g. your home-directory), you will have to make sure that the corresponding installation sub-directory is part of your `PATH` environment variable. If you also want to use the non-essential Python tools that come with LoFreq, please make sure to add the corresponding installation sub-directory to your `PYTHONPATH`. For example, for the prefix setting mentioned above ($HOME/local/), you would need the following:

    export PATH=$HOME/local/bin/
    export PYTHONPATH=$HOME/local/lib/python2.7/site-packages/

The first line ensures that you can call the LoFreq simply by typing their names. The second makes sure that the extra (optional) Python scripts work. 

# Further help

If the documentation above is not sufficient and you have trouble installing LoFreq, please ask for help on the mailing or email us directly.

