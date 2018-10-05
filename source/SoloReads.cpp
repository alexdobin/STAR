#include "SoloReads.h"
SoloReads::SoloReads(Parameters Pin) : P(Pin), pSolo(P.pSolo) {
    if (pSolo==0)
        return;
    //allocate arrays
//     reads.nMax=P.limitNreadsSoft/P.runNthreads;
//     read.cb=new uint32[nMax];
//     read.umi=new uint32[nMax];
//     read.gene=new uint32[nMax];
//     reads.N=0;
    
};