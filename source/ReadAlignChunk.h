#ifndef READ_ALIGN_CHUNK_DEF
#define READ_ALIGN_CHUNK_DEF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "ReadAlign.h"
#include "OutSJ.h"
#include "Transcriptome.h"
#include "BAMoutput.h"
#include "Quantifications.h"

#if !defined(_WIN32) && defined(USE_PTHREAD)
#include <pthred>
#else
#include <mutex>
#endif

class ReadAlignChunk {//chunk of reads and alignments
public:
    Parameters* P;
    ReadAlign* RA;

    Transcriptome *chunkTr;
    
    char **chunkIn; //space for the chunk of input reads
    char *chunkOutBAM, *chunkOutBAM1;//space for the chunk of output SAM
    OutSJ *chunkOutSJ, *chunkOutSJ1;

    BAMoutput *chunkOutBAMcoord, *chunkOutBAMunsorted, *chunkOutBAMquant;
    Quantifications *chunkQuants;
    
    istringstream** readInStream;
    ostringstream*  chunkOutBAMstream;
    ofstream chunkOutBAMfile;
    string chunkOutBAMfileName;
    
    bool noReadsLeft;
	bool foundEnd; 
    uint iChunkIn; //current chunk # as read from .fastq
    uint iChunkOutSAM; //current chunk # writtedn to Aligned.out.sam
    int iThread; //current thread
    uint chunkOutBAMtotal, chunkOutBAMtotal1; //total number of bytes in the write buffer
            
    ReadAlignChunk(Parameters* Pin, Genome &genomeIn, Transcriptome *TrIn, int iChunk);
    void processChunks();
    void mapChunk();
    void chunkFstreamOpen(string filePrefix, int iChunk, fstream &fstreamOut);

#if !defined(_WIN32) && defined(USE_PTHREAD)
	void chunkFstreamCat(fstream &chunkOut, ofstream &allOut, bool mutexFlag, pthread_mutex_t &mutexVal);
#else
	void chunkFstreamCat(fstream &chunkOut, ofstream &allOut, bool mutexFlag, std::mutex &mutexVal);
#endif

    void chunkFilesCat(ostream *allOut, string filePrefix, uint &iC);
};
#endif
