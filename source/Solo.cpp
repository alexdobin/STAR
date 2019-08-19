#include "Solo.h"
#include "TimeFunctions.h"

Solo::Solo(ReadAlignChunk **RAchunkIn, Parameters &Pin, Transcriptome &inTrans)
          :  RAchunk(RAchunkIn), P(Pin), Trans(inTrans), pSolo(P.pSolo)
{
    if ( pSolo.type==0 )
        return;
    if (pSolo.type==pSolo.SoloTypes::CB_samTagOut)
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
    if (pSolo.type==pSolo.SoloTypes::CB_samTagOut)
        return;

    *P.inOut->logStdOut << timeMonthDayTime() << " ..... started Solo counting\n" <<flush;
    P.inOut->logMain    << timeMonthDayTime() << " ..... started Solo counting\n" <<flush;
        
    for (uint32 ii=0; ii<pSolo.nFeatures; ii++) {
        soloFeat[ii]->processRecords(RAchunk);
    };
    
    *P.inOut->logStdOut << timeMonthDayTime() << " ..... finished Solo counting\n" <<flush;
    P.inOut->logMain    << timeMonthDayTime() << " ..... finished Solo counting\n" <<flush;
};
