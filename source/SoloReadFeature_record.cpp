#include "SoloReadFeature.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"

uint32 outputReadCB(fstream *streamOut, int32 featureType, uint32 umiB, uint32 gene, vector<array<uint64,2>> &readSJs, const string &stringCB)
{
    if (featureType==0 || featureType==2) {//genes
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

void SoloReadFeature::record(SoloReadBarcode &soloBar, uint nTr, set<uint32> &readGene, set<uint32> &readGeneFull, Transcript *alignOut)
{
    if (pSolo.type==0 || soloBar.cbMatch<0)
        return;

    //unmapped
    if (nTr==0) {
        stats.V[stats.nUnmapped]++;
        return;
    };

    vector<array<uint64,2>> readSJs;

    set<uint32> *readGe;
    if (featureType==0) {
        readGe = &readGene;
    } else if (featureType==2) {
        readGe = &readGeneFull;
    };

    if (featureType==0 || featureType==2) {//genes
        //check genes, return if no gene of multimapping
        if (readGe->size()==0) {
            stats.V[stats.nNoFeature]++;
            return;
        };
        if (readGe->size()>1) {
            stats.V[stats.nAmbigFeature]++;
            if (nTr>1)
                stats.V[stats.nAmbigFeatureMultimap]++;
            return;
        };
    } else if (featureType==1) {//SJs
        if (nTr>1) {//reject all multimapping junctions
            stats.V[stats.nAmbigFeatureMultimap]++;
            return;
        };

        //for SJs, still check genes, return if multi-gene
        if (readGene.size()>1) {
            stats.V[stats.nAmbigFeature]++;
            return;
        };
        bool sjAnnot;
        alignOut->extractSpliceJunctions(readSJs, sjAnnot);
        if ( readSJs.empty() || (sjAnnot && readGene.size()==0) ) {//no junctions, or annotated junction buy no gene (i.e. read does not fully match transcript)
            stats.V[stats.nNoFeature]++;
            return;
        };
    };

    if (soloBar.cbMatch==0) {//exact match
        uint32 n1 = outputReadCB(strU_0, featureType, soloBar.umiB, *readGe->begin(), readSJs, to_string(soloBar.cbI));
        if (pSolo.cbWL.size()>0) {//WL
            cbReadCount[soloBar.cbI] += n1;
        } else {//no WL
            cbReadCountMap[soloBar.cbI] += n1;
        };
        return;
    } else if (soloBar.cbMatch==1) {//1 match with 1MM
        cbReadCount[soloBar.cbI]+= outputReadCB(strU_1, featureType, soloBar.umiB, *readGe->begin(), readSJs, to_string(soloBar.cbI));
        return;
    } else {//>1 matches
        uint32 nfeat=outputReadCB(strU_2, featureType, soloBar.umiB, *readGe->begin(), readSJs, to_string(soloBar.cbMatch) + soloBar.cbMatchString);
        for (auto &cbi : soloBar.cbMatchInd)
            cbReadCount[cbi] += nfeat;
        return;
    };
};
