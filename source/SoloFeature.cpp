#include "SoloFeature.h"
#include "streamFuns.h"

SoloFeature::SoloFeature(int feTy, Parameters &Pin, Transcriptome &inTrans) 
          :  featureType(feTy), P(Pin), pSolo(P.pSolo), Trans(inTrans)
{
    if (pSolo.type==0)
        return;

    readFeatSum = new SoloReadFeature(featureType,P,-1);
    readBarSum = new SoloReadBarcode(P);    
    readFeatAll = new SoloReadFeature*[P.runThreadN];
    readBarAll = new SoloReadBarcode*[P.runThreadN];

    statsStream = &ofstrOpen(P.outFileNamePrefix+pSolo.outFileNames[0]+pSolo.featureNames[featureType]+".stats",ERROR_OUT, P);
};
