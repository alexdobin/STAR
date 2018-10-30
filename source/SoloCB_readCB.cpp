#include "SoloCB.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"

void SoloCB::readCB(const string &readNameExtra, const uint nTr, const vector<int32> &readGenes, vector<uint32> &readTranscripts, set<uint32> &readTrGenes) {
    if (pSolo.type==0)
        return;
    
//     int32 rG=readGene.at(P.pReads.strand);
//     if (nTr > 1) //unique mappers only for now
//        return;

    if (readTrGenes.size()==0) {
        stats.V[stats.nNoGene]++;
        return;
    } else if (readTrGenes.size()>1) {
        stats.V[stats.nAmbigGene]++;
        if (nTr>1)
                stats.V[stats.nAmbigGeneMultimap]++;
        return;
    };
    
    //debug
    if (readGenes.size()==1 && readGenes.at(0)!=*readTrGenes.begin())
        cout <<readGenes.at(0) <<' '<< *readTrGenes.begin() <<'\n';
    //
    
    uint32 cbB, umiB;
    int32 cbI=-2;
    if (convertNuclStrToInt32(readNameExtra.substr(0,pSolo.cbL),cbB) \
        && convertNuclStrToInt32(readNameExtra.substr(pSolo.cbL,pSolo.umiL),umiB)) {
        cbI=binarySearchExact(cbB,pSolo.cbWL.data(),pSolo.cbWL.size());
    } else {
        stats.V[stats.nNinBarcode]++;
        return;
    };

    if (cbI!=-1) {
        *strU_0 << cbI <<' '<< *readTrGenes.begin() <<' '<< umiB <<'\n';
    } else {
        for (uint32 ii=0; ii<2*pSolo.cbL; ii+=2) {
            for (uint32 jj=1; jj<4; jj++) {
                cbI=binarySearchExact(cbB^(jj<<ii),pSolo.cbWL.data(),pSolo.cbWL.size());
                if (cbI>=0) {
                    *strU_1 << cbI <<' '<< readNameExtra.at(pSolo.cbL+pSolo.umiL+1+ii/2) <<' ';
                };
            };
        };
        *strU_1 << *readTrGenes.begin() <<' '<< umiB <<'\n';
    };   
};
