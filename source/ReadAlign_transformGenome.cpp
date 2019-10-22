#include "ReadAlign.h"
#include "ErrorWarning.h"

void ReadAlign::transformGenome() 
{//convert to new genome
    if (!mapGen.genomeOut.convYes || mapGen.pGe.transform.type==0)
        return;
    
    uint32 nTr1=0;
    for (uint32 iTr=0; iTr<nTr; iTr++) {//convert output transcripts into new genome
        *trMultOut[nTr1]=*trMult[iTr];//copy information before conversion
        if (trMult[iTr]->transformGenome(*mapGen.genomeOut.g, *trMultOut[nTr1])) {
            ++nTr1;
            trMult[nTr1-1] = trMultOut[nTr1-1]; //point to new transcsript
        };
    };
    nTr=nTr1;
};
