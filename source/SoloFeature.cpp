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
