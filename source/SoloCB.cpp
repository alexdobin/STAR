#include "SoloCB.h"
#include "streamFuns.h"

SoloCB::SoloCB(Parameters &Pin, int iChunk) : P(Pin), pSolo(P.pSolo) {
    if (pSolo.type==0)
        return;
    
    cbReadCount = new uint32[pSolo.cbWL.size()];
    for (uint32 ii=0; ii<pSolo.cbWL.size(); ii++)
        cbReadCount[ii]=0;
    
    strU_0 = &ofstrOpen(P.outFileTmp+"/soloCB_U_0_"+std::to_string(iChunk),ERROR_OUT, P);
    strU_1 = &ofstrOpen(P.outFileTmp+"/soloCB_U_1_"+std::to_string(iChunk),ERROR_OUT, P);
    
};