#include "SoloCB.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"

void SoloCB::readCB(const string &readNameExtra, const uint nTr, const vector<int32> &readGene) {
    if (pSolo.type==0)
        return;
    
//     if (nTr > 1) //unique mappers only for now
// maps        return;
    
    int32 rG=readGene.at(P.pReads.strand);
    
//     if (rG < 0) //only non-ambiguous genes for now
//         return;
    
    uint32 cbB;
    if (!convertNuclStrToInt32(readNameExtra.substr(0,pSolo.cbL),cbB))
        return; //non-ACGT symbols
    
    int32 cbI=binarySearchExact(cbB,pSolo.cbWL.data(),pSolo.cbWL.size());
    
    if (cbI>=0) {
        *strU_0 << cbI <<' '<<rG<< '\n';
    } else {
    };   
};