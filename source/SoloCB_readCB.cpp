#include "SoloCB.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"

void SoloCB::readCB(const string &readNameExtra, const uint nTr, const vector<int32> &readGene) {
    if (pSolo.type==0)
        return;
    
    int32 rG=readGene.at(P.pReads.strand);

//     if (nTr > 1) //unique mappers only for now
//        return;

    if (nTr > 1)
        rG=-3;

    
//     if (rG < 0) //only non-ambiguous genes for now
//         return;
    
    uint32 cbB, umiB;
    int32 cbI=-2;
    if (convertNuclStrToInt32(readNameExtra.substr(0,pSolo.cbL),cbB) \
        && convertNuclStrToInt32(readNameExtra.substr(pSolo.cbL,pSolo.umiL),umiB)) {
        cbI=binarySearchExact(cbB,pSolo.cbWL.data(),pSolo.cbWL.size());
    } else {
        return; //TODO stats for non-recorded
    };

    if (cbI!=-1) {
        *strU_0 << cbI <<' '<< rG <<' '<< umiB <<'\n';
    } else {
        for (uint32 ii=0; ii<2*pSolo.cbL; ii+=2) {
            for (uint32 jj=1; jj<4; jj++) {
                cbI=binarySearchExact(cbB^(jj<<ii),pSolo.cbWL.data(),pSolo.cbWL.size());
                if (cbI>=0) {
                    *strU_1 << cbI <<' '<< readNameExtra.at(pSolo.cbL+pSolo.umiL+1+ii/2) <<' ';
                };
            };
        };
        *strU_1 << rG <<' '<< umiB <<'\n';
    };   
};
