#include "Solo.h"

Solo::Solo(ReadAlignChunk **RAchunkIn, Parameters &Pin, Transcriptome &inTrans) 
          :  RAchunk(RAchunkIn), P(Pin), pSolo(P.pSolo), Trans(inTrans)
{
    if (pSolo.type==0 )
        return;

    soloFeat = new SoloFeature*[pSolo.nFeatures];
    for (uint32 ii=0; ii<pSolo.nFeatures; ii++)
        soloFeat[ii] = new SoloFeature(pSolo.features[ii], P, Trans);
};

////////////////////////////////////////////////////////////////////////////////////
void Solo::processAndOutput()
{
    if (pSolo.type==0 )
    return;
    
    for (uint32 ii=0; ii<pSolo.nFeatures; ii++) {
        soloFeat[ii]->processRecords(RAchunk);
    };
};
