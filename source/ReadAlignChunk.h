#ifndef CODE_ReadAlignChunk
#define CODE_ReadAlignChunk

#include "IncludeDefine.h"
#include "Parameters.h"
#include "ReadAlign.h"
#include "OutSJ.h"
#include "Transcriptome.h"
#include "BAMoutput.h"
#include "Quantifications.h"

class ReadAlignChunk {//chunk of reads and alignments
public:
    Parameters& P;
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
    uint iChunkIn; //current chunk # as read from .fastq
    uint iChunkOutSAM; //current chunk # writtedn to Aligned.out.sam
    int iThread; //current thread
    uint chunkOutBAMtotal; //total number of bytes in the write buffer

    ReadAlignChunk(Parameters& Pin, Genome &genomeIn, Transcriptome *TrIn, int iChunk);
    void processChunks();
    void mapChunk();
    void chunkFstreamOpen(string filePrefix, int iChunk, fstream &fstreamOut);
    void chunkFstreamCat (fstream &chunkOut, ofstream &allOut, bool mutexFlag, pthread_mutex_t &mutexVal);
    void chunkFilesCat(ostream *allOut, string filePrefix, uint &iC);
    
    Genome &mapGen;    
private:
};
#endif
