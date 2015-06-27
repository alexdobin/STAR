COMPTIMEPLACE := -D'COMPILATION_TIME_PLACE="$(shell echo `date` $(HOSTNAME):`pwd`)"'

AM_CXXFLAGS = -std=c++0x -Wall -Wextra -pthread -fopenmp $(COMPTIMEPLACE)

bin_PROGRAMS = STAR

STAR_SOURCES = $(STAR_CXX_FILES)
STAR_CPPFLAGS = $(HTSLIB_CFLAGS)
STAR_LDFLAGS = -pthread $(HTSLIB_LIBS)

BUILT_SOURCES = parametersDefault.xxd

parametersDefault.xxd: parametersDefault
	ln -s $(top_srcdir)/source/parametersDefault
	xxd -i parametersDefault > parametersDefault.xxd

CLEANFILES = parametersDefault parametersDefault.xxd

STAR_CXX_FILES = \
	alignSmithWaterman.cpp \
	BAMbinSortByCoordinate.cpp \
	BAMbinSortUnmapped.cpp \
	bam_cat.cpp \
	BAMfunctions.cpp \
	BAMoutput.cpp \
	bamRemoveDuplicates.cpp \
	binarySearch2.cpp \
	blocksOverlap.cpp \
	ErrorWarning.cpp \
	extendAlign.cpp \
	Genome.cpp \
	genomeGenerate.cpp \
	genomeParametersWrite.cpp \
	InOutStreams.cpp \
	loadGTF.cpp \
	mapThreadsSpawn.cpp \
	outputSJ.cpp \
	OutSJ.cpp \
	PackedArray.cpp \
	Parameters_closeReadsFiles.cpp \
	Parameters.cpp \
	Parameters_openReadsFiles.cpp \
	Quantifications.cpp \
	ReadAlign_alignBAM.cpp \
	ReadAlign_assignAlignToWindow.cpp \
	ReadAlign_chimericDetection.cpp \
	ReadAlignChunk.cpp \
	ReadAlignChunk_mapChunk.cpp \
	ReadAlignChunk_processChunks.cpp \
	ReadAlign.cpp \
	ReadAlign_createExtendWindowsWithAlign.cpp \
	ReadAlign_mapOneRead.cpp \
	ReadAlign_maxMappableLength2strands.cpp \
	ReadAlign_multMapSelect.cpp \
	ReadAlign_oneRead.cpp \
	ReadAlign_outputAlignments.cpp \
	ReadAlign_outputTranscriptCIGARp.cpp \
	ReadAlign_outputTranscriptSAM.cpp \
	ReadAlign_outputTranscriptSJ.cpp \
	ReadAlign_quantTranscriptome.cpp \
	ReadAlign_stitchPieces.cpp \
	ReadAlign_stitchWindowSeeds.cpp \
	ReadAlign_storeAligns.cpp \
	readLoad.cpp \
	SequenceFuns.cpp \
	SharedMemory.cpp \
	signalFromBAM.cpp \
	sjdbBuildIndex.cpp \
	sjdbInsertJunctions.cpp \
	sjdbLoadFromFiles.cpp \
	sjdbLoadFromStream.cpp \
	sjdbPrepare.cpp \
	STAR.cpp \
	Stats.cpp \
	stitchAlignToTranscript.cpp \
	stitchWindowAligns.cpp \
	streamFuns.cpp \
	stringSubstituteAll.cpp \
	SuffixArrayFuns.cpp \
	sysRemoveDir.cpp \
	ThreadControl.cpp \
	TimeFunctions.cpp \
	Transcript_alignScore.cpp \
	Transcript.cpp \
	Transcriptome.cpp \
	Transcriptome_geneCountsAddAlign.cpp \
	Transcriptome_quantAlign.cpp
