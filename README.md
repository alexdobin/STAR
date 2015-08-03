[![Build Status](https://travis-ci.org/nathanhaigh/STAR.svg?branch=master)](https://travis-ci.org/nathanhaigh/STAR)

STAR 2.4
########
Spliced Transcripts Alignment to a Reference
© Alexander Dobin, 2009-2015

AUTHOR/SUPPORT
##############
Alex Dobin, dobin@cshl.edu
https://groups.google.com/d/forum/rna-star

MANUAL
######
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

(RELEASEnotes)[RELEASEnotes] contains detailed information about the latest major release

DIRECTORY CONTENTS
##################
source: all source files required for compilation
bin: pre-compiled executables for Linux and Mac OS X
doc: documentation
extras: miscellaneous files and scripts

STAR-Fusion: fusion detection developed by Brian Haas, see https://github.com/STAR-Fusion/STAR-Fusion for details.
             To populate this submodule, clone STAR with `git clone --recursive https://github.com/alexdobin/STAR`
STAR-Fusion-x.x.x: latest release of the STAR-Fusion


COMPILING FROM SOURCE
#####################
Unzip and "cd" into source/ subdirectory.
Linux:    run `make STAR`
Mac OS X: run `make STARforMacStatic`
          If g++ compiler (true g++, not Clang sym-link) is not on the path, run `make STARforMacStatic CXX=/path/to/gcc`


HARDWARE/SOFTWARE REQUIREMENTS
##############################
x86-64 compatible processors
64 bit Linux or Mac OS X 
30GB of RAM for human genome 


LIMITATIONS
###########
This release was tested with the default parameters for human and mouse genomes.
Please contact the author for a list of recommended parameters for much larger or much smaller genomes.

