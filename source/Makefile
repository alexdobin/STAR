# user may define these whole flags
# LDFLAGS 
# CXXFLAGS

# or these user-set flags that will be added to standard flags
LDFLAGSextra ?=
CXXFLAGSextra ?=

# user may define the compiler
CXX ?=g++

# pre-defined flags
LDFLAGS_shared := -pthread -Lhtslib -Bstatic -lhts -Bdynamic -lz -lrt
LDFLAGS_static := -static -static-libgcc -pthread -Lhtslib -lhts -lz
LDFLAGS_Mac :=-pthread -lz htslib/libhts.a
LDFLAGS_Mac_static :=-pthread -lz -static-libgcc htslib/libhts.a
LDFLAGS_gdb := $(LDFLAGS_shared)

COMPTIMEPLACE := -D'COMPILATION_TIME_PLACE="$(shell echo `date` $(HOSTNAME):`pwd`)"'

CXXFLAGS_common := -pipe -std=c++0x -Wall -Wextra -fopenmp $(COMPTIMEPLACE)
CXXFLAGS_main := -O3 $(CXXFLAGS_common)
CXXFLAGS_gdb :=  -O0 -g $(CXXFLAGS_common)



##########################################################################################################

OBJECTS = SharedMemory.o PackedArray.o SuffixArrayFuns.o STAR.o Parameters.o InOutStreams.o SequenceFuns.o Genome.o Stats.o \
        Transcript.o Transcript_alignScore.o \
        ReadAlign.o ReadAlign_storeAligns.o ReadAlign_stitchPieces.o ReadAlign_multMapSelect.o ReadAlign_mapOneRead.o readLoad.o \
	ReadAlignChunk.o ReadAlignChunk_processChunks.o ReadAlignChunk_mapChunk.o \
	OutSJ.o outputSJ.o blocksOverlap.o ThreadControl.o sysRemoveDir.o \
        ReadAlign_maxMappableLength2strands.o binarySearch2.o\
	ReadAlign_outputAlignments.o  \
	ReadAlign_outputTranscriptSAM.o ReadAlign_outputTranscriptSJ.o ReadAlign_outputTranscriptCIGARp.o \
        ReadAlign_createExtendWindowsWithAlign.o ReadAlign_assignAlignToWindow.o ReadAlign_oneRead.o \
	ReadAlign_stitchWindowSeeds.o ReadAlign_chimericDetection.o \
        stitchWindowAligns.o extendAlign.o stitchAlignToTranscript.o alignSmithWaterman.o \
        genomeGenerate.o genomeParametersWrite.o \
	TimeFunctions.o ErrorWarning.o loadGTF.o streamFuns.o stringSubstituteAll.o \
        Transcriptome.o Transcriptome_quantAlign.o ReadAlign_quantTranscriptome.o Quantifications.o Transcriptome_geneCountsAddAlign.o \
        sjdbLoadFromFiles.o sjdbLoadFromStream.o sjdbPrepare.o sjdbBuildIndex.o sjdbInsertJunctions.o mapThreadsSpawn.o \
        Parameters_openReadsFiles.cpp Parameters_closeReadsFiles.cpp \
        BAMoutput.o BAMfunctions.o ReadAlign_alignBAM.o BAMbinSortByCoordinate.o signalFromBAM.o bamRemoveDuplicates.o BAMbinSortUnmapped.o \
        bam_cat.o \
        GlobalVariables.cpp

SOURCES := $(wildcard *.cpp) $(wildcard *.c)


%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $<

%.o : %.c
	$(CXX) -c $(CXXFLAGS) $<

all: STAR

.PHONY: clean
clean:
	rm -f *.o STAR STARstatic STARlong Depend.list

.PHONY: CLEAN
CLEAN:
	rm -f *.o STAR Depend.list
	$(MAKE) -C htslib clean

.PHONY: cleanRelease
cleanRelease:
	rm -f *.o Depend.list
	$(MAKE) -C htslib clean


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),cleanRelease)
ifneq ($(MAKECMDGOALS),CLEAN)
ifneq ($(MAKECMDGOALS),STARforMac)
ifneq ($(MAKECMDGOALS),STARforMacGDB)
Depend.list: $(SOURCES) parametersDefault.xxd htslib
	echo $(SOURCES)
	/bin/rm -f ./Depend.list
	$(CXX) $(CXXFLAGS_main) -MM $^ >> Depend.list
include Depend.list
endif
endif
endif
endif
endif

htslib : htslib/libhts.a

htslib/libhts.a :
	$(MAKE) -C htslib lib-static

parametersDefault.xxd: parametersDefault
	xxd -i parametersDefault > parametersDefault.xxd

STAR : CXXFLAGS ?= $(CXXFLAGSextra) $(CXXFLAGS_main)
STAR : LDFLAGS ?= $(LDFLAGSextra) $(LDFLAGS_shared)
STAR : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

POSIXSHARED : CXXFLAGS ?= $(CXXFLAGSextra) $(CXXFLAGS_main) -DPOSIX_SHARED_MEM
POSIXSHARED : LDFLAGS ?= $(LDFLAGSextra) $(LDFLAGS_shared)
POSIXSHARED : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARstatic : CXXFLAGS ?= $(CXXFLAGSextra) $(CXXFLAGS_main)
STARstatic : LDFLAGS ?= $(LDFLAGSextra) $(LDFLAGS_static)
STARstatic : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARlong : CXXFLAGS ?= $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_LONG_READS'
STARlong : LDFLAGS ?= $(LDFLAGSextra) $(LDFLAGS_shared)
STARlong : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STARlong $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARlongStatic : CXXFLAGS ?= $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_LONG_READS'
STARlongStatic : LDFLAGS ?= $(LDFLAGSextra) $(LDFLAGS_static)
STARlongStatic : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STARlong $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

gdb : CXXFLAGS ?= $(CXXFLAGSextra) $(CXXFLAGS_gdb)
gdb : LDFLAGS ?= $(LDFLAGSextra) $(LDFLAGS_gdb)
gdb : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

gdb-long : CXXFLAGS ?= $(CXXFLAGSextra) $(CXXFLAGS_gdb) -D'COMPILE_FOR_LONG_READS'
gdb-long : LDFLAGS ?= $(LDFLAGSextra) $(LDFLAGS_gdb)
gdb-long : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STARlong $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARforMacStatic : CXXFLAGS ?= $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_MAC'
STARforMacStatic : LDFLAGS ?= $(LDFLAGSextra) $(LDFLAGS_Mac_static)
STARforMacStatic : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

######################################################### all trargets below are not supported and not recommended!

STARforMac : CXXFLAGS=-D'COMPILE_FOR_MAC' -I ./Mac_Include/ $(CXXFLAGS_main)
STARforMac : parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CXXFLAGS) $(LDFLAGS_Mac) $(OBJECTS)


STARlongForMacStatic : CXXFLAGS=-D'COMPILE_FOR_LONG_READS' -D'COMPILE_FOR_MAC' -I ./Mac_Include/ $(CXXFLAGS_main)
STARlongForMacStatic : parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STARlong $(CXXFLAGS) $(LDFLAGS_Mac_static) $(OBJECTS)

#
STARforMacGDB : CXXFLAGS=-D'COMPILE_FOR_MAC' -I ./Mac_Include/ $(CXXFLAGS_gdb)
STARforMacGDB : parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CXXFLAGS_gdb) $(OBJECTS) $(LDFLAGS_gdb)

localChains : CXXFLAGS=-D'OUTPUT_localChains' $(CXXFLAGS_main)
localChains : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CXXFLAGS) $(LDFLAGS) $(OBJECTS)


