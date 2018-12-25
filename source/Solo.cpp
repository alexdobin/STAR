#include "Solo.h"
#include "streamFuns.h"

Solo::Solo(int feTy, Parameters &Pin, Transcriptome &inTrans) 
          :  featureType(feTy), P(Pin), pSolo(P.pSolo), Trans(inTrans)
{
    if (pSolo.type==0  || !pSolo.featureYes[featureType])
        return;

    soloCBsum = new SoloCB(featureType,P,-1);
    soloCBall = new SoloCB*[P.runThreadN];

    soloStatsStream = &ofstrOpen(P.outFileNamePrefix+pSolo.outFileNames[0]+pSolo.featureNames[featureType]+".stats",ERROR_OUT, P);
};
