#include "Solo.h"
#include "streamFuns.h"

Solo::Solo(Parameters &Pin, Transcriptome &inTrans) : P(Pin), pSolo(P.pSolo), Trans(inTrans) {
    if (pSolo.type==0)
        return;

    soloCBsum = new SoloCB(P,-1);
    soloCBall = new SoloCB*[P.runThreadN];

    soloStatsStream = &ofstrOpen(P.outFileNamePrefix+pSolo.outFileNames[0]+"Solo.stats",ERROR_OUT, P);    
};
