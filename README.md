STAR 2.6
========
Spliced Transcripts Alignment to a Reference
© Alexander Dobin, 2009-2018
https://www.ncbi.nlm.nih.gov/pubmed/23104886

AUTHOR/SUPPORT
==============
Alex Dobin, dobin@cshl.edu
https://groups.google.com/d/forum/rna-star

HARDWARE/SOFTWARE REQUIREMENTS
==============================
  * x86-64 compatible processors
  * 64 bit Linux or Mac OS X 

MANUAL
======
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

[RELEASEnotes](https://github.com/alexdobin/STAR/blob/master/RELEASEnotes.md) contains detailed information about the latest major release

DIRECTORY CONTENTS
==================
  * source: all source files required for compilation
  * bin: pre-compiled executables for Linux and Mac OS X
  * doc: documentation
  * extras: miscellaneous files and scripts

COMPILING FROM SOURCE
=====================

Download the latest [release from](https://github.com/alexdobin/STAR/releases) and uncompress it
--------------------------------------------------------

```bash
# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz
tar -xzf 2.6.0a.tar.gz
cd STAR-2.6.0a

# Alternatively, get STAR source using git
git clone https://github.com/alexdobin/STAR.git

Compile under Linux
-------------------

```bash
# Compile
cd STAR/source
make STAR
```

Compile under Mac OS X
----------------------

```bash
# Compile
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


LIMITATIONS
===========
This release was tested with the default parameters for human and mouse genomes.
Mammal genomes require at least 16GB of RAM, ideally 32GB.
Please contact the author for a list of recommended parameters for much larger or much smaller genomes.


FUNDING
=======
The developmenr of STAR is supported by the National Human Genome Research Institute of
the National Institutes of Health under Award Number R01HG009318. 
The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
