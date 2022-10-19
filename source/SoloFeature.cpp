#include "SoloFeature.h"
#include "streamFuns.h"

SoloFeature::SoloFeature(Parameters &Pin, ReadAlignChunk **RAchunk, Transcriptome &inTrans, int32 feTy, SoloReadBarcode *readBarSumIn, SoloFeature **soloFeatAll)
            : P(Pin), RAchunk(RAchunk), Trans(inTrans), featureType(feTy), soloFeatAll(soloFeatAll), pSolo(P.pSolo), readBarSum(readBarSumIn)
{
    if (featureType>=0) {//otherwise we do not need these arrays - e.g. with --runMode soloCellFiltering 
        readFeatSum = new SoloReadFeature(featureType,P,-1);
        readFeatAll = new SoloReadFeature*[P.runThreadN];
    };
    
    //number of features
    switch (featureType) {
        case SoloFeatureTypes::Gene :
        case SoloFeatureTypes::GeneFull :
        case SoloFeatureTypes::GeneFull_Ex50pAS :
        case SoloFeatureTypes::GeneFull_ExonOverIntron :
        case SoloFeatureTypes::Velocyto :
            featuresNumber=Trans.nGe;
            break;
        case SoloFeatureTypes::SJ :
            featuresNumber=P.sjAll[0].size();
            break;
        default:
            featuresNumber = -1; //undefined
    };    
};

void SoloFeature::clearLarge()
{
    cbFeatureUMImap.clear();
    cbFeatureUMImap.shrink_to_fit();
    countCellGeneUMI.clear();
    countCellGeneUMI.shrink_to_fit();
    countCellGeneUMIindex.clear();
    countCellGeneUMIindex.shrink_to_fit();
    countMatMult.i.clear();
    countMatMult.i.shrink_to_fit();
    countMatMult.m.clear();
    countMatMult.m.shrink_to_fit();
    //indCB.clear(); //needed for Velocyto
    //indCB.shrink_to_fit();
    indCBwl.clear();
    indCBwl.shrink_to_fit();
    nGenePerCB.clear();
    nGenePerCB.shrink_to_fit();
    nGenePerCBmulti.clear();
    nGenePerCBmulti.shrink_to_fit();
    nReadPerCB.clear();
    nReadPerCB.shrink_to_fit();
    nReadPerCBtotal.clear();
    nReadPerCBtotal.shrink_to_fit();
    nReadPerCBunique.clear();
    nReadPerCBunique.shrink_to_fit();
    nUMIperCB.clear();
    nUMIperCB.shrink_to_fit();
    nUMIperCBmulti.clear();
    nUMIperCBmulti.shrink_to_fit();
    nUMIperCBsorted.clear();
    nUMIperCBsorted.shrink_to_fit();
    sjAll[0].clear();
    sjAll[0].shrink_to_fit();
    sjAll[1].clear();
    sjAll[1].shrink_to_fit();
};
