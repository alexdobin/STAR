OBJECTS = PackedArray.o SuffixArrayFuns.o STAR.o Parameters.o InOutStreams.o SequenceFuns.o Genome.o Stats.o \
        Transcript.o Transcript_alignScore.o \
        ReadAlign.o ReadAlign_storeAligns.o ReadAlign_stitchPieces.o ReadAlign_multMapSelect.o ReadAlign_mapOneRead.o readLoad.o \
	ReadAlignChunk.o ReadAlignChunk_processChunks.o ReadAlignChunk_mapChunk.o \
	OutSJ.o outputSJ.o blocksOverlap.o ThreadControl.o sysRemoveDir.o \
        ReadAlign_maxMappableLength2strands.o binarySearch2.o\
	ReadAlign_outputAlignments.o  \
	ReadAlign_outputTranscriptSAM.o ReadAlign_outputTranscriptSJ.o ReadAlign_outputTranscriptCIGARp.o \
        ReadAlign_createExtendWindowsWithAlign.o ReadAlign_assignAlignToWindow.o ReadAlign_oneRead.o \
	ReadAlign_stitchWindowSeeds.o ReadAlign_chimericDetection.o \
        stitchWindowAligns.o extendAlign.o stitchAlignToTranscript.o alignSmithWaterman.o genomeGenerate.o \
	TimeFunctions.o ErrorWarning.o loadGTF.o streamFuns.o stringSubstituteAll.o \
        Transcriptome.o Transcriptome_quantAlign.o ReadAlign_quantTranscriptome.o Quantifications.o Transcriptome_geneCountsAddAlign.o \
        sjdbLoadFromFiles.o sjdbLoadFromStream.o sjdbPrepare.o sjdbBuildIndex.o mapThreadsSpawn.o \
        Parameters_openReadsFiles.cpp Parameters_closeReadsFiles.cpp \
        BAMoutput.o BAMfunctions.o ReadAlign_alignBAM.o BAMbinSortByCoordinate.o signalFromBAM.o bamRemoveDuplicates.o BAMbinSortUnmapped.o \
        bam_cat.o
SOURCES := $(wildcard *.cpp) $(wildcard *.c)

LDFLAGS := -pthread -Lhtslib -Bstatic -lhts -Bdynamic -lz
LDFLAGS_static := -static -static-libgcc -pthread -Lhtslib -lhts -lz
LDFLAGS_Mac :=-pthread -lz htslib/libhts.a
LDFLAGS_Mac_static :=-pthread -lz -static-libgcc htslib/libhts.a

LDFLAGS_gdb := $(LDFLAGS)

COMPTIMEPLACE := -D'COMPILATION_TIME_PLACE="$(shell echo `date` $(HOSTNAME):`pwd`)"'
EXTRAFLAGS := 

CCFLAGS_common := -pipe -std=c++0x -Wall -Wextra -fopenmp $(COMPTIMEPLACE) $(OPTIMFLAGS) $(OPTIMFLAGS1) $(EXTRAFLAGS)
CCFLAGS_main := -O3 $(CCFLAGS_common)
CCFLAGS_gdb :=  -O0 -g $(CCFLAGS_common)

CXX :=g++


%.o : %.cpp
	$(CXX) -c $(CCFLAGS) $<

%.o : %.c
	$(CXX) -c $(CCFLAGS) $<

all: STAR

.PHONY: clean
clean:
	rm -f *.o STAR STARstatic Depend.list

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
	$(CXX) $(CCFLAGS_common) -MM $^ >> Depend.list
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

STAR : CCFLAGS=$(CCFLAGS_main)
STAR : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CCFLAGS) $(OBJECTS) $(LDFLAGS)

STARstatic : CCFLAGS=$(CCFLAGS_main)
STARstatic : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STARstatic $(OBJECTS) $(CCFLAGS) $(LDFLAGS_static)

STARlong : CCFLAGS=-D'COMPILE_FOR_LONG_READS' $(CCFLAGS_main)
STARlong : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CCFLAGS) $(OBJECTS) $(LDFLAGS)

STARlongStatic : CCFLAGS=-D'COMPILE_FOR_LONG_READS' $(CCFLAGS_main)
STARlongStatic : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STARstatic $(CCFLAGS) $(LDFLAGS_static) $(OBJECTS)


STARforMac : CCFLAGS=-D'COMPILE_FOR_MAC' -I ./Mac_Include/ $(CCFLAGS_main)
STARforMac : parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CCFLAGS) $(LDFLAGS_Mac) $(OBJECTS)

STARforMacStatic : CCFLAGS=-D'COMPILE_FOR_MAC' -I ./Mac_Include/ $(CCFLAGS_main)
STARforMacStatic : parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CCFLAGS) $(LDFLAGS_Mac_static) $(OBJECTS)


STARforMacGDB : CCFLAGS=-D'COMPILE_FOR_MAC' -I ./Mac_Include/ $(CCFLAGS_gdb)
STARforMacGDB : parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CCFLAGS_gdb) $(OBJECTS) $(LDFLAGS_gdb)

gdb : CCFLAGS= $(CCFLAGS_gdb)
gdb : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CCFLAGS_gdb) $(OBJECTS) $(LDFLAGS_gdb) 

gdb-long : CCFLAGS= -D'COMPILE_FOR_LONG_READS' $(CCFLAGS_gdb)
gdb-long : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CCFLAGS_gdb) $(OBJECTS) $(LDFLAGS_gdb) 

localChains : CCFLAGS=-D'OUTPUT_localChains' $(CCFLAGS_main)
localChains : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CCFLAGS) $(LDFLAGS) $(OBJECTS)


