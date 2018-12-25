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
    
    if (featureType==1) {//read the SJ array
        ifstream sjStream((P.outFileNamePrefix+"SJ.out.tab").c_str());
        uint64 sj[2];
        while (sjStream >> sj[0] >> sj[1]) {
            for (int ii=0; ii<2; ii++)
                sjAll[ii].push_back(sj[ii]);
            sjStream.ignore ((uint32) (-1), '\n');
        };
        sjStream.close();
    };
};
