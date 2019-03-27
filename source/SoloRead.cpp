#include "SoloRead.h"

SoloRead::SoloRead(Parameters &Pin, int32 iChunkIn) :  iChunk(iChunkIn), P(Pin), pSolo(P.pSolo)
{
    if (pSolo.type==0)
        return;    
    
    readBar = new SoloReadBarcode(P);
    readFeat = new SoloReadFeature*[pSolo.nFeatures];

    for (uint32 ii=0; ii<pSolo.nFeatures; ii++)
        readFeat[ii] = new SoloReadFeature(pSolo.features[ii], P, iChunk);

};
