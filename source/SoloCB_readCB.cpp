#include "SoloCB.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"
#include "SoloBarcodes.h"

uint32 outputReadCB(fstream *streamOut, int32 featureType, uint32 umiB, uint32 gene, vector<array<uint64,2>> &readSJs, string &stringCB)
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

bool SoloCB::readCB(SoloBarcodes &soloBar, uint nTr, set<uint32> &readTrGenes, Transcript *alignOut) 
{
    if (pSolo.type==0  || !pSolo.featureYes[featureType] || soloB.cbMatch<0)
        return false;

    //unmapped
    if (nTr==0) {
        stats.V[stats.nUnmapped]++;
        return false;
    };
    
    vector<array<uint64,2>> readSJs;
    if (featureType==0) {//genes
        //check genes, return if no gene
        if (readTrGenes.size()==0) {
            stats.V[stats.nNoFeature]++;
            return false;
        } else if (readTrGenes.size()>1) {
            stats.V[stats.nAmbigFeature]++;
            if (nTr>1)
                stats.V[stats.nAmbigFeatureMultimap]++;
            return false;
        };
    } else if (featureType==1) {//SJs
        alignOut->extractSpliceJunctions(readSJs);
        if (readSJs.empty()) {
            stats.V[stats.nNoFeature]++; 
            return false;
        } else if (nTr>1) {//record SJs from the read
            stats.V[stats.nAmbigFeature]++;
            stats.V[stats.nAmbigFeatureMultimap]++;
            return false;
        };
    };
    
    if (soloB.cbMatch==0) {//exact match
        cbReadCount[cbI] += outputReadCB(strU_0, featureType, soloB.umiB, *readTrGenes.begin(), readSJs, to_string(soloB.cbI));
        return true;
    } else if (soloB.cbMatch==1) {//1 match with 1MM
        cbReadCount[cbI]+= outputReadCB(strU_1, featureType, umiB, *readTrGenes.begin(), readSJs, to_string(cbI));
        return true;
    } else {
        cbReadCount[cbI] += outputReadCB(strU_2, featureType, umiB, *readTrGenes.begin(), readSJs, to_string(soloB.cbMatch) + cbOutString);
        return true;        
    };
};
