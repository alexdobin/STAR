Fixed a problem that can arise for very small genomes while using --alignIntronMax 1.

STAR 2.5
========
Spliced Transcripts Alignment to a Reference
© Alexander Dobin, 2009-2016
https://www.ncbi.nlm.nih.gov/pubmed/23104886

AUTHOR/SUPPORT
==============
Alex Dobin, dobin@cshl.edu
https://groups.google.com/d/forum/rna-star

MANUAL
======
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

[RELEASEnotes](RELEASEnotes.md) contains detailed information about the latest major release

DIRECTORY CONTENTS
==================
  * source: all source files required for compilation
  * bin: pre-compiled executables for Linux and Mac OS X
  * doc: documentation
  * extras: miscellaneous files and scripts
  * STAR-Fusion: fusion detection developed by Brian Haas, see https://github.com/STAR-Fusion/STAR-Fusion for details.
             To populate this submodule, clone STAR with `git clone --recursive https://github.com/alexdobin/STAR`
  * STAR-Fusion-x.x.x: latest release of the STAR-Fusion


COMPILING FROM SOURCE
=====================

To compile STAR from source, you must first download the latest [release](release) and uncompress it and then build it.

Linux
-----

```bash
# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz
tar -xzf 2.5.3a.tar.gz
cd STAR-2.5.3a

# Alternatively, get STAR source using git
git clone https://github.com/alexdobin/STAR.git
cd STAR/source

# Build STAR
make STAR

# To include STAR-Fusion
git submodule update --init --recursive

# If you have a TeX environment, you may like to build the documentation
make manual
```

Mac OS X
--------

```bash
# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz
tar -xzf 2.5.2b.tar.gz
cd STAR-2.5.2b

# Alternatively, get STAR source using git
git clone https://github.com/alexdobin/STAR.git
cd STAR

# To include STAR-Fusion
git submodule update --init --recursive

# Build STAR
cd source
make STARforMacStatic
```

All platforms
-------------

If g++ compiler (true g++, not Clang sym-link) is not on the path, you will need to tell `make` where to find it:

```bash
# Build STAR
cd source
make STARforMacStatic CXX=/path/to/gcc
```

If employing STAR only on a single machine or a homogeneously setup cluster, you may aim at helping the compiler to optimize in way that is tailored to your platform. The flags LDFLAGSextra and CXXFLAGSextra are appended to the default optimizations specified in source/Makefile.

```
# platform-specific optimization for gcc/g++
make CXXFLAGSextra=-march=native
# together with link-time optimization
make LDFLAGSextra=-flto CXXFLAGSextra="-flto -march=native"
```


Developers
==========

STAR developers with write access to https://github.com/alexdobin/STAR can update the `STAR-Fusion`
submodule to a specific tag by following these steps:

```bash
git clone --recursive https://github.com/alexdobin/STAR.git
cd STAR
# or:
#
# git clone //github.com/alexdobin/STAR.git
# cd STAR
# git git submodule update --init --recursive

# checkout a specific tag for the submodule
cd STAR-Fusion
git checkout v0.3.1

# Commit the change
cd ../
git add STAR-Fusion
git commit -m "Updated STAR-Fusion to v0.3.1"

# Push the change to GitHub
git push
```


HARDWARE/SOFTWARE REQUIREMENTS
==============================
  * x86-64 compatible processors
  * 64 bit Linux or Mac OS X 
  * 30GB of RAM for human genome 


LIMITATIONS
===========
This release was tested with the default parameters for human and mouse genomes.
Please contact the author for a list of recommended parameters for much larger or much smaller genomes.

