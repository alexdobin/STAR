# user may define the whole sets of flags
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

DATE_FMT = --iso-8601=seconds
ifdef SOURCE_DATE_EPOCH
    BUILD_DATE ?= $(shell date -u -d "@$(SOURCE_DATE_EPOCH)" "$(DATE_FMT)" 2>/dev/null || date -u -r "$(SOURCE_DATE_EPOCH)" "$(DATE_FMT)" 2>/dev/null || date -u "$(DATE_FMT)")
else
    BUILD_DATE ?= $(shell date "$(DATE_FMT)")
endif

BUILD_PLACE ?= $(HOSTNAME):$(shell pwd)

COMPTIMEPLACE := -D'COMPILATION_TIME_PLACE="$(BUILD_DATE) $(BUILD_PLACE)"'


GIT_CHECK := $(shell git status 1> /dev/null 2> /dev/null && echo 0)
ifeq ($(GIT_CHECK),0)
	GIT_DIFF := $(shell git diff --name-only | tr '\n' ' ')
    GIT_BRANCH_COMMIT_DIFF := $(shell git status | head -n 1) ; \
                                                    $(shell git log | head -n 1) ; \
                                                    diff files: $(GIT_DIFF)
endif
GIT_BRANCH_COMMIT_DIFF := -D'GIT_BRANCH_COMMIT_DIFF="$(GIT_BRANCH_COMMIT_DIFF)"'

# Defaults, can be overridden by make arguments or environment
CXXFLAGS ?= -pipe -Wall -Wextra
CFLAGS ?= -pipe -Wall -Wextra -O3
CXXFLAGS_SIMD ?= -mavx2

# Unconditionally set essential flags and optimization options
CXXFLAGS_common := -std=c++11 -fopenmp $(COMPTIMEPLACE) $(GIT_BRANCH_COMMIT_DIFF)
CXXFLAGS_main := -O3 $(CXXFLAGS_common)
CXXFLAGS_gdb := -O0 -g3 $(CXXFLAGS_common)

##########################################################################################################
OBJECTS = systemFunctions.o funPrimaryAlignMark.o \
	SoloFeature_collapseUMI_Graph.o SoloFeature_collapseUMIall.o ParametersClip_initialize.o ClipMate_clip.o ClipCR4.o opal/opal.o ClipMate_clipChunk.o ClipMate_initialize.o \
	SoloFeature_loadRawMatrix.o SoloFeature_emptyDrops_CR.o soloInputFeatureUMI.o SoloFeature_countSmartSeq.o SoloFeature_redistributeReadsByCB.o \
	SoloFeature_quantTranscript.o SoloFeature_sumThreads.o SoloFeature_countVelocyto.o SoloFeature_countCBgeneUMI.o \
	Transcriptome_classifyAlign.o Transcriptome_geneFullAlignOverlap_ExonOverIntron.o Transcriptome_alignExonOverlap.cpp \
	SoloFeature_cellFiltering.o \
	SoloFeature_statsOutput.o bamSortByCoordinate.o SoloBarcode.o \
	ParametersSolo.o SoloRead.o SoloRead_record.o \
	SoloReadBarcode.o SoloReadBarcode_getCBandUMI.o SoloBarcode_extractBarcode.o \
	SoloReadFeature.o SoloReadFeature_record.o SoloReadFeature_inputRecords.o \
	Solo.o SoloFeature.o SoloFeature_outputResults.o SoloFeature_processRecords.o SoloFeature_addBAMtags.o \
	ReadAlign_transformGenome.o Genome_transformGenome.o Transcript_convertGenomeCigar.o \
	twoPassRunPass1.o samHeaders.o Genome_genomeLoad.o Genome_genomeOutLoad.o Transcript_transformGenome.o ReadAlign_outputSpliceGraphSAM.o \
	ReadAlign_mapOneReadSpliceGraph.o SpliceGraph.o SpliceGraph_swScoreSpliced.o SpliceGraph_swTraceBack.o \
	SpliceGraph_findSuperTr.o sjAlignSplit.o \
	GTF.o GTF_transcriptGeneSJ.o GTF_superTranscript.o SuperTranscriptome.o \
	ReadAlign_outputAlignments.o  \
	ReadAlign.o STAR.o \
	SharedMemory.o PackedArray.o SuffixArrayFuns.o Parameters.o Parameters_samAttributes.o InOutStreams.o SequenceFuns.o Genome.o ParametersGenome.o Stats.o \
	Transcript.o Transcript_alignScore.o Transcript_generateCigarP.o Chain.o \
	Transcript_variationAdjust.o Variation.o ReadAlign_waspMap.o \
	ReadAlign_storeAligns.o ReadAlign_stitchPieces.o ReadAlign_multMapSelect.o ReadAlign_mapOneRead.o readLoad.o \
	ReadAlignChunk.o ReadAlignChunk_processChunks.o ReadAlignChunk_mapChunk.o \
	OutSJ.o outputSJ.o blocksOverlap.o ThreadControl.o sysRemoveDir.o \
	ReadAlign_maxMappableLength2strands.o binarySearch2.o\
	ReadAlign_outputTranscriptSAM.o ReadAlign_outputTranscriptSJ.o ReadAlign_outputTranscriptCIGARp.o ReadAlign_calcCIGAR.cpp \
	ReadAlign_createExtendWindowsWithAlign.o ReadAlign_assignAlignToWindow.o ReadAlign_oneRead.o \
	ReadAlign_stitchWindowSeeds.o \
	ReadAlign_peOverlapMergeMap.o ReadAlign_mappedFilter.o \
	ParametersChimeric_initialize.o ReadAlign_chimericDetection.o ReadAlign_chimericDetectionOld.o ReadAlign_chimericDetectionOldOutput.o\
	ChimericDetection.o ChimericDetection_chimericDetectionMult.o ReadAlign_chimericDetectionPEmerged.o \
	stitchWindowAligns.o extendAlign.o stitchAlignToTranscript.o \
	ChimericSegment.cpp ChimericAlign.cpp ChimericAlign_chimericJunctionOutput.o ChimericAlign_chimericBAMoutput.o ChimericAlign_chimericStitching.o \
	Genome_genomeGenerate.o genomeParametersWrite.o genomeScanFastaFiles.o genomeSAindex.o \
	Genome_insertSequences.o insertSeqSA.o funCompareUintAndSuffixes.o funCompareUintAndSuffixesMemcmp.o \
	TimeFunctions.o ErrorWarning.o streamFuns.o stringSubstituteAll.o \
	Transcriptome.o Transcriptome_quantAlign.o Transcriptome_geneFullAlignOverlap.o \
	ReadAlign_quantTranscriptome.o Quantifications.o Transcriptome_geneCountsAddAlign.o \
	sjdbLoadFromFiles.o sjdbLoadFromStream.o sjdbPrepare.o sjdbBuildIndex.o sjdbInsertJunctions.o mapThreadsSpawn.o \
	Parameters_readFilesInit.o Parameters_openReadsFiles.cpp Parameters_closeReadsFiles.cpp Parameters_readSAMheader.o \
	bam_cat.o serviceFuns.o GlobalVariables.cpp \
	BAMoutput.o BAMfunctions.o ReadAlign_alignBAM.o BAMbinSortByCoordinate.o signalFromBAM.o bamRemoveDuplicates.o BAMbinSortUnmapped.o

SOURCES := $(wildcard *.cpp) $(wildcard *.c)


%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $<

%.o : %.c
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $<

all: cleanCompileInfo STAR$(SFX)

opal/opal.o : opal/opal.cpp opal/opal.h
	cd opal && \
	$(CXX) -c -I./ -std=c++11 $(CPPFLAGS) $(CXXFLAGS) $(CXXFLAGSextra) $(CXXFLAGS_SIMD) opal.cpp

.PHONY: clean
clean:
	'rm' -f *.o opal/opal.o STAR STARstatic STARlong Depend.list parametersDefault.xxd

.PHONY: CLEAN
CLEAN: clean
	$(MAKE) -C htslib clean


.PHONY: clean_solo
clean_solo:
	'rm' -f Solo*.o

.PHONY: cleanCompileInfo
cleanCompileInfo:
	'rm' -f STAR.o Parameters.o

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),CLEAN)
ifneq ($(MAKECMDGOALS),CLEAN)
ifneq ($(MAKECMDGOALS),clean_solo)
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
endif

htslib : htslib/libhts.a

htslib/libhts.a :
	$(MAKE) -C htslib lib-static

parametersDefault.xxd: parametersDefault
	xxd -i parametersDefault > parametersDefault.xxd

STAR$(SFX) : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
STAR$(SFX) : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
STAR$(SFX) : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR$(SFX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARstatic$(SFX) : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
STARstatic$(SFX) : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_static) $(LDFLAGS)
STARstatic$(SFX) : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR$(SFX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARlong$(SFX) : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_LONG_READS' $(CXXFLAGS)
STARlong$(SFX) : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
STARlong$(SFX) : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STARlong$(SFX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARlongStatic$(SFX) : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_LONG_READS' $(CXXFLAGS)
STARlongStatic$(SFX) : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_static) $(LDFLAGS)
STARlongStatic$(SFX) : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STARlong$(SFX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)



POSIXSHARED : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) -DPOSIX_SHARED_MEM $(CXXFLAGS)
POSIXSHARED : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
POSIXSHARED : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR$(SFX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

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
