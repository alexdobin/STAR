#include "Solo.h"
#include "TimeFunctions.h"

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

    *P.inOut->logStdOut << timeMonthDayTime() << " ..... started Solo counting\n" <<flush;
    P.inOut->logMain    << timeMonthDayTime() << " ..... started Solo counting\n" <<flush;
        
    for (uint32 ii=0; ii<pSolo.nFeatures; ii++) {
        soloFeat[ii]->processRecords(RAchunk);
    };
    
    *P.inOut->logStdOut << timeMonthDayTime() << " ..... finished Solo counting\n" <<flush;
    P.inOut->logMain    << timeMonthDayTime() << " ..... finished Solo counting\n" <<flush;
    
};
