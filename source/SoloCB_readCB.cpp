#include "SoloCB.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"

void SoloCB::readCB(const uint64 &iReadAll, const string &readNameExtra, const uint nTr, const vector<int32> &readGenes, vector<uint32> &readTranscripts, set<uint32> &readTrGenes) {
    if (pSolo.type==0)
        return;
    
//     int32 rG=readGene.at(P.pReads.strand);

    if (readTrGenes.size()==0) {
        stats.V[stats.nNoGene]++;
        return;
    } else if (readTrGenes.size()>1) {
        stats.V[stats.nAmbigGene]++;
        if (nTr>1)
            stats.V[stats.nAmbigGeneMultimap]++;
        return;
    };
    
    /*debug
    if (readGenes.size()==1 && readGenes.at(0)!=*readTrGenes.begin())
        cout <<readGenes.at(0) <<' '<< *readTrGenes.begin() <<'\n';
    */
    
    uint32 cbB, umiB;
    int32 cbI=-2;
    if (!(convertNuclStrToInt32(readNameExtra.substr(0,pSolo.cbL),cbB) \
        && convertNuclStrToInt32(readNameExtra.substr(pSolo.cbL,pSolo.umiL),umiB))) {//convert and check for Ns
        stats.V[stats.nNinBarcode]++;
        return;
    };
    
    if (umiB==homoPolymer[0] || umiB==homoPolymer[1] || umiB==homoPolymer[2] || umiB==homoPolymer[3]) {
        stats.V[stats.nUMIhomopolymer]++;
        return;
    };
    
    cbI=binarySearchExact(cbB,pSolo.cbWL.data(),pSolo.cbWL.size());


//     if (cbI!=-1) {
//         *strU_0 << cbI <<' '<< *readTrGenes.begin() <<' '<< umiB <<' '<< iReadAll  <<'\n';
//         cbReadCount[cbI]++;
//     } else {
//         for (uint32 ii=0; ii<2*pSolo.cbL; ii+=2) {
//             for (uint32 jj=1; jj<4; jj++) {
//                 cbI=binarySearchExact(cbB^(jj<<ii),pSolo.cbWL.data(),pSolo.cbWL.size());
//                 if (cbI>=0) {                        
//                     *strU_1 << cbI <<' '<< readNameExtra.at(pSolo.cbL+pSolo.umiL+1+ii/2) <<' ';
//                     cbReadCount[cbI]++;
//                 };
//             };
//         };
//         *strU_1 << *readTrGenes.begin() <<' '<< umiB <<' '<< iReadAll <<'\n';
//     };   
    
    //simple procedure: accept only if one WL CB exists with 1MM
    if (cbI>=0) {
        stats.V[stats.nExactMatch]++;
    } else {
        for (uint32 ii=0; ii<2*pSolo.cbL; ii+=2) {
            for (uint32 jj=1; jj<4; jj++) {
                int32 cbI1=binarySearchExact(cbB^(jj<<ii),pSolo.cbWL.data(),pSolo.cbWL.size());
                if (cbI1>=0) {                        
                    if (cbI>=0) {//had another match already
                        stats.V[stats.nTooMany]++;
                        cout << iReadAll << endl;
                        return;
                    };
                    cbI=cbI1;
                };
            };
        };
    };
    if (cbI<0) {
        stats.V[stats.nNoMatch]++;
        return;
    };
    
    stats.V[stats.nMatch]++;
    //output to file
    cbReadCount[cbI]++;
    *strU_0 << cbI <<' '<< *readTrGenes.begin() <<' '<< umiB <<'\n';
    //*strU_0 << cbI <<' '<< *readTrGenes.begin() <<' '<< umiB <<' '<< iReadAll  <<'\n';
};
