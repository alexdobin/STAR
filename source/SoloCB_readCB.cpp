#include "SoloCB.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"

uint32 outputReadCB(fstream *streamOut, int32 featureType, uint32 umiB, uint32 gene, vector<array<uint64,2>> &readSJs, string stringCB)
{
    if (featureType==0) {//genes
        *streamOut << umiB <<' '<< gene <<' '<< stringCB <<'\n';
        return 1;
    } else if (featureType==1) {//sjs
        for (auto &sj : readSJs) {
            *streamOut << umiB <<' '<< sj[0] <<' '<< sj[1] <<' '<< stringCB <<'\n';
        };
        return readSJs.size();
    };
    
    return 0; //this should not happen
};

void SoloCB::readCB(uint64 &iReadAll, string &readNameExtra, uint nTr, set<uint32> &readTrGenes, Transcript *alignOut) 
{
    if (pSolo.type==0  || !pSolo.featureYes[featureType])
        return;
    
    int32 cbI=-999;
    uint32 cbB, umiB;
    
    cbSeq=readNameExtra.substr(0,pSolo.cbL);
    umiSeq=readNameExtra.substr(pSolo.cbL,pSolo.umiL);
    cbQual=readNameExtra.substr(pSolo.cbL+pSolo.umiL+1,pSolo.cbL);
    umiQual=readNameExtra.substr(pSolo.cbL+pSolo.umiL+1+pSolo.cbL,pSolo.umiL);
    
    int32 posN=convertNuclStrToInt32(cbSeq,cbB);
    if (posN==-2) {//>2 Ns, might already be filtered by Illumina
        stats.V[stats.nNinBarcode]++;
        return;
    } else if (posN==-1) {//no Ns
        cbI=binarySearchExact(cbB,pSolo.cbWL.data(),pSolo.cbWL.size());
        if (cbI>=0) {//exact match
            cbReadCountExact[cbI]++;//note that this simply counts reads per exact CB, no checks of genes or UMIs
        };
    };
    
    if (nTr==0) {
        stats.V[stats.nUnmapped]++;
        return;
    };
    
    vector<array<uint64,2>> readSJs;
    if (featureType==0) {//genes
        //check genes, return if no gene
        if (readTrGenes.size()==0) {
            stats.V[stats.nNoFeature]++;
            return;
        } else if (readTrGenes.size()>1) {
            stats.V[stats.nAmbigFeature]++;
            if (nTr>1)
                stats.V[stats.nAmbigFeatureMultimap]++;
            return;
        };
    } else if (featureType==1) {//SJs
        alignOut->extractSpliceJunctions(readSJs);
        if (readSJs.empty()) {
            stats.V[stats.nNoFeature]++;        
        } else if (nTr>1) {//record SJs from the read
            stats.V[stats.nAmbigFeature]++;
            stats.V[stats.nAmbigFeatureMultimap]++;
        };
    };

    //check UMIs, return if bad UMIs
    if (convertNuclStrToInt32(umiSeq,umiB)!=-1) {//convert and check for Ns
        stats.V[stats.nNinBarcode]++;//UMIs are not allowed to have MMs
        return;
    };
    if (umiB==homoPolymer[0] || umiB==homoPolymer[1] || umiB==homoPolymer[2] || umiB==homoPolymer[3]) {
        stats.V[stats.nUMIhomopolymer]++;
        return;
    };
    
    if (cbI>=0) {//exact match
        cbReadCount[cbI] += outputReadCB(strU_0, featureType, umiB, *readTrGenes.begin(), readSJs, to_string(cbI));
        return;
    };
        
    if (posN>=0) {//one N
        uint32 posNshift=2*(pSolo.cbL-1-posN);//shift bits for posN
        for (uint32 jj=0; jj<4; jj++) {
            uint32 cbB1=cbB^(jj<<posNshift);
            int32 cbI1=binarySearchExact(cbB1,pSolo.cbWL.data(),pSolo.cbWL.size());
            if (cbI1>=0) {                        
                if (cbI>=0) {//had another match already
                    stats.V[stats.nTooMany]++;
                    return;//with N in CB, do not allow matching >1 in WL
                };
                cbI=cbI1;
            };
        };
        if (cbI>=0) {
            cbReadCount[cbI]+= outputReadCB(strU_1, featureType, umiB, *readTrGenes.begin(), readSJs, to_string(cbI));
        } else {//no match
            stats.V[stats.nNoMatch]++;
        };
        return;
    } else {//look for 1MM, posN==-1, no Ns
        string cbOutString("");
        uint32 ncbOut1=0;
        for (uint32 ii=0; ii<pSolo.cbL; ii++) {
            for (uint32 jj=1; jj<4; jj++) {
                int32 cbI1=binarySearchExact(cbB^(jj<<(ii*2)),pSolo.cbWL.data(),pSolo.cbWL.size());
                if (cbI1>=0) {//found match
                    //output all
                    cbI=cbI1;
                    ++ncbOut1;
                    cbOutString += ' ' +to_string(cbI) + ' ' + cbQual.at(pSolo.cbL-1-ii);
                    cbReadCount[cbI]++;//this read may be counted in many CBs. This is safe for allocating arrays later. The final number of reads per CB will be calculated later.                   
                };
            };
        };

        if (ncbOut1==0) {//no matches    
            stats.V[stats.nNoMatch]++;
            return;
        } else if (ncbOut1==1){//only one match
            cbReadCount[cbI] += outputReadCB(strU_1, featureType, umiB, *readTrGenes.begin(), readSJs, to_string(cbI));
        } else {//more than one match
            cbReadCount[cbI] += outputReadCB(strU_2, featureType, umiB, *readTrGenes.begin(), readSJs, to_string(ncbOut1) + to_string(cbI));
        };
    };

};
