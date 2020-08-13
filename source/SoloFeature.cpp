#include "SoloFeature.h"
#include "streamFuns.h"

SoloFeature::SoloFeature(int32 feTy, Parameters &Pin, Transcriptome &inTrans, SoloReadBarcode *readBarSumIn, SoloFeature **soloFeatAll)
            :featureType(feTy), P(Pin), Trans(inTrans), soloFeatAll(soloFeatAll), pSolo(P.pSolo), readBarSum(readBarSumIn)
{

    readFeatSum = new SoloReadFeature(featureType,P,-1);
    readFeatAll = new SoloReadFeature*[P.runThreadN];

    if (pSolo.type==0)
        return;
};
