# user may define these whole flags
# LDFLAGS
# CPPFLAGS
# CXXFLAGS
# CFLAGS

# or these user-set flags that will be added to standard flags
LDFLAGSextra ?=
CXXFLAGSextra ?=

# user may define the compiler
CXX ?= g++

# pre-defined flags
LDFLAGS_shared := -pthread -Lhtslib -Bstatic -lhts -Bdynamic -lz
LDFLAGS_static := -static -static-libgcc -pthread -Lhtslib -lhts -lz
LDFLAGS_Mac :=-pthread -lz htslib/libhts.a
LDFLAGS_Mac_static :=-pthread -lz -static-libgcc htslib/libhts.a
LDFLAGS_gdb := $(LDFLAGS_shared)

COMPTIMEPLACE := -D'COMPILATION_TIME_PLACE="$(shell echo `date` $(HOSTNAME):`pwd`)"'

CXXFLAGS_common := -pipe -std=c++11 -Wall -Wextra -fopenmp $(COMPTIMEPLACE)
CXXFLAGS_main := -O3 $(CXXFLAGS_common)
CXXFLAGS_gdb :=  -O0 -g $(CXXFLAGS_common)

CFLAGS := -O3 -pipe -Wall -Wextra $(CFLAGS)


##########################################################################################################

OBJECTS = SharedMemory.o PackedArray.o SuffixArrayFuns.o STAR.o Parameters.o InOutStreams.o SequenceFuns.o Genome.o Stats.o \
        Transcript.o Transcript_alignScore.o Transcript_generateCigarP.o Chain.o \
        Transcript_variationAdjust.o Variation.o ReadAlign_waspMap.o \
        ReadAlign.o ReadAlign_storeAligns.o ReadAlign_stitchPieces.o ReadAlign_multMapSelect.o ReadAlign_mapOneRead.o readLoad.o \
	ReadAlignChunk.o ReadAlignChunk_processChunks.o ReadAlignChunk_mapChunk.o \
	OutSJ.o outputSJ.o blocksOverlap.o ThreadControl.o sysRemoveDir.o \
        ReadAlign_maxMappableLength2strands.o binarySearch2.o\
	ReadAlign_outputAlignments.o  \
	ReadAlign_outputTranscriptSAM.o ReadAlign_outputTranscriptSJ.o ReadAlign_outputTranscriptCIGARp.o ReadAlign_calcCIGAR.cpp \
        ReadAlign_createExtendWindowsWithAlign.o ReadAlign_assignAlignToWindow.o ReadAlign_oneRead.o \
	ReadAlign_stitchWindowSeeds.o \
        ReadAlign_peOverlapMergeMap.o ReadAlign_mappedFilter.o \
        ReadAlign_chimericDetection.o ReadAlign_chimericDetectionOld.o ReadAlign_chimericDetectionOldOutput.o\
        ChimericDetection.o ChimericDetection_chimericDetectionMult.o ReadAlign_chimericDetectionPEmerged.o \
        ChimericSegment.cpp ChimericAlign.cpp ChimericAlign_chimericJunctionOutput.o ChimericAlign_chimericStitching.o \
        stitchWindowAligns.o extendAlign.o stitchAlignToTranscript.o alignSmithWaterman.o \
        Genome_genomeGenerate.o genomeParametersWrite.o genomeScanFastaFiles.o genomeSAindex.o \
        Genome_insertSequences.o insertSeqSA.o funCompareUintAndSuffixes.o funCompareUintAndSuffixesMemcmp.o \
	TimeFunctions.o ErrorWarning.o loadGTF.o streamFuns.o stringSubstituteAll.o \
        Transcriptome.o Transcriptome_quantAlign.o ReadAlign_quantTranscriptome.o Quantifications.o Transcriptome_geneCountsAddAlign.o \
        sjdbLoadFromFiles.o sjdbLoadFromStream.o sjdbPrepare.o sjdbBuildIndex.o sjdbInsertJunctions.o mapThreadsSpawn.o \
        Parameters_openReadsFiles.cpp Parameters_closeReadsFiles.cpp Parameters_readSAMheader.o \
        BAMoutput.o BAMfunctions.o ReadAlign_alignBAM.o BAMbinSortByCoordinate.o signalFromBAM.o bamRemoveDuplicates.o BAMbinSortUnmapped.o \
        bam_cat.o \
        serviceFuns.o \
        GlobalVariables.cpp

SOURCES := $(wildcard *.cpp) $(wildcard *.c)


%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $<

%.o : %.c
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $<

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

.PHONY: install
install:
	'mv' STAR STARlong ../bin

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),cleanRelease)
ifneq ($(MAKECMDGOALS),CLEAN)
ifneq ($(MAKECMDGOALS),STARforMac)
ifneq ($(MAKECMDGOALS),STARforMacGDB)
Depend.list: $(SOURCES) parametersDefault.xxd htslib
	echo $(SOURCES)
	'rm' -f ./Depend.list
	$(CXX) $(CXXFLAGS_common) -MM $^ >> Depend.list
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

STAR : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
STAR : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
STAR : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

POSIXSHARED : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) -DPOSIX_SHARED_MEM $(CXXFLAGS)
POSIXSHARED : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
POSIXSHARED : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARstatic : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
STARstatic : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_static) $(LDFLAGS)
STARstatic : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARlong : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_LONG_READS' $(CXXFLAGS)
STARlong : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
STARlong : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STARlong $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARlongStatic : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_LONG_READS' $(CXXFLAGS)
STARlongStatic : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_static) $(LDFLAGS)
STARlongStatic : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STARlong $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

gdb : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_gdb) $(CXXFLAGS)
gdb : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_gdb) $(LDFLAGS)
gdb : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

gdb-long : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_gdb) -D'COMPILE_FOR_LONG_READS' $(CXXFLAGS)
gdb-long : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_gdb) $(LDFLAGS)
gdb-long : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STARlong $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARforMacStatic : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_MAC' $(CXXFLAGS)
STARforMacStatic : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_Mac_static) $(LDFLAGS)
STARforMacStatic : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARlongForMacStatic : CXXFLAGS := -D'COMPILE_FOR_LONG_READS' $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_MAC' $(CXXFLAGS)
STARlongForMacStatic : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_Mac_static) $(LDFLAGS)
STARlongForMacStatic : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STARlong $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)


