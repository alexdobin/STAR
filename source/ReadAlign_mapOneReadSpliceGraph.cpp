#include "ReadAlign.h"
#include "SequenceFuns.h"
#include "Stats.h"
#include "serviceFuns.cpp"

void ReadAlign::mapOneReadSpliceGraph() 
{
    #ifdef OFF_BEFORE_SEEDING
        #warning OFF_BEFORE_SEEDING
        nW=0;
        return;
    #endif

    if (Lread<P.outFilterMatchNmin) {//read is too short (trimmed too much?)
        mapMarker=MARKER_READ_TOO_SHORT;
        trBest->rLength=0; //min good piece length
        nW=0;
        return;
    };       
        
    resetN(); //reset aligns counters to 0

    //reset/initialize a transcript
    trInit->reset();
    trInit->Chr=0;    trInit->Str=0; trInit->roStr=0;    trInit->cStart=0;     trInit->gLength=0; //to generate nice output of 0 for non-mapped reads
    trInit->iRead=iRead;
    trInit->Lread=Lread;
    trInit->nExons=0;
    trInit->readLengthOriginal=readLengthOriginal;
    trInit->readLengthPairOriginal=readLengthPairOriginal;
    trInit->readLength=readLength;
    trInit->readNmates=readNmates; //not readNends: this is alignment
    trInit->readName=readName;

    trBest=trInit;

    splGraph->findSuperTr(Read1[0], Read1[2], Lread, readName, mapGen);
    
    return;
};
