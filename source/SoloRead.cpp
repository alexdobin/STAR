#include "SoloRead.h"

SoloRead::SoloRead(Parameters &Pin, int32 iChunkIn) :  iChunk(iChunkIn), P(Pin), pSolo(P.pSolo)
{
    readBar = new SoloReadBarcode(P);
    readFeat = new SoloReadFeature*[pSolo.nFeatures];
    
    if (pSolo.type==0)
        return;

    for (uint32 ii=0; ii<pSolo.nFeatures; ii++)
        readFeat[ii] = new SoloReadFeature(pSolo.features[ii], P, iChunk);
    
};
