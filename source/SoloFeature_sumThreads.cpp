#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "Stats.h"
#include "GlobalVariables.h"

void SoloFeature::sumThreads(ReadAlignChunk **RAchunk)
{      
    ///////////////////////////// collect RAchunk->RA->soloRead->readFeat
    nReadsInput=g_statsAll.readN+1; //reserve 1 extra
            
    for (int ii=0; ii<P.runThreadN; ii++) {//point to
        readFeatAll[ii]= RAchunk[ii]->RA->soloRead->readFeat[pSolo.featureInd[featureType]];
    };
    
    for (int ii=0; ii<P.runThreadN; ii++) {
        readFeatSum->addCounts(*readFeatAll[ii]);
    };

    if (!pSolo.cbWLyes) {//now we can define WL and counts ??? we do not need to do it for every feature???
        pSolo.cbWLsize=readFeatSum->cbReadCountMap.size();
        pSolo.cbWL.resize(pSolo.cbWLsize);
        pSolo.cbWLstr.resize(pSolo.cbWLsize);
        uint64 ii=0;
        for (auto &cb : readFeatSum->cbReadCountMap) {
            pSolo.cbWL[ii] = cb.first;
            pSolo.cbWLstr[ii] = convertNuclInt64toString(pSolo.cbWL[ii],pSolo.cbL); 
            ii++;
        };
        readFeatSum->cbReadCount = new uint32[pSolo.cbWLsize];
        readBarSum->cbReadCountExact = new uint32[pSolo.cbWLsize];

        uint64 icb=0;
        for (auto ii=readFeatSum->cbReadCountMap.cbegin(); ii!=readFeatSum->cbReadCountMap.cend(); ++ii) {
            pSolo.cbWL[icb]=ii->first;
            readFeatSum->cbReadCount[icb]=ii->second;
            readBarSum->cbReadCountExact[icb]=ii->second;
            ++icb;
        };
    };

    //detected CBs
    nCB=0;nReadsMapped=0;
    for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
        if (readFeatSum->cbReadCount[ii]>0) {
            nCB++;
            nReadsMapped += readFeatSum->cbReadCount[ii];
        };
    };
    
    indCB = new uint32[nCB];
    nCB=0;//will count it again below
    for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
        if (readFeatSum->cbReadCount[ii]>0) {
            indCB[nCB]=ii;
            ++nCB;
        };
    };
    
    //pseudocounts
    if (pSolo.CBmatchWL.mm1_multi_pc) {
        for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
            readBarSum->cbReadCountExact[ii]++;//add one to exact counts
        };
    };
};
    