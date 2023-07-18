#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "Stats.h"
#include "GlobalVariables.h"

void SoloFeature::sumThreads()
{   
    //stats
    nReadsInput=g_statsAll.readN+1; //reserve 1 extra

    ///////////////////////////// collect RAchunk->RA->soloRead->readFeat            
    for (int ii=0; ii<P.runThreadN; ii++) {//point to
        readFeatAll[ii]= RAchunk[ii]->RA->soloRead->readFeat[pSolo.featureInd[featureType]];
        readFeatAll[ii]->streamReads->flush();
        readFeatSum->addCounts(*readFeatAll[ii]);        
    };       
    
    // if WL was not defined
    if (!pSolo.cbWLyes) {//now we can define WL and counts ??? we do not need to do it for every feature???
        pSolo.cbWLsize=readFeatSum->cbReadCountMap.size();
        pSolo.cbWL.resize(pSolo.cbWLsize);
        pSolo.cbWLstr.resize(pSolo.cbWLsize);
        readFeatSum->cbReadCount.resize(pSolo.cbWLsize);
        readBarSum->cbReadCountExact.resize(pSolo.cbWLsize);

        if (pSolo.CBtype.type==1) {//sequence cb
            uint64 icb=0;
            for (auto &cb : readFeatSum->cbReadCountMap) {
                pSolo.cbWL[icb] = cb.first;
                pSolo.cbWLstr[icb] = convertNuclInt64toString(pSolo.cbWL[icb],pSolo.cbL);
                readFeatSum->cbReadCount[icb]=cb.second;
                readBarSum->cbReadCountExact[icb]=cb.second;
                icb++;
            };
        } else if (pSolo.CBtype.type==2) {//string cb
            vector< std::unordered_map<string,uint32>::iterator > cbiter(pSolo.CBtype.strMap.size());
            for (auto cbi=pSolo.CBtype.strMap.begin(); cbi!=pSolo.CBtype.strMap.end(); cbi++)
                cbiter[cbi->second] = cbi;

            uint64 icb=0;
            for (auto &cb : readFeatSum->cbReadCountMap) {
                pSolo.cbWL[icb] = cb.first;
                pSolo.cbWLstr[icb] = cbiter[cb.first]->first;
                readFeatSum->cbReadCount[icb]=cb.second;
                readBarSum->cbReadCountExact[icb]=cb.second;
                icb++;
            };
        };

        //pseudocounts
        if (pSolo.CBmatchWL.mm1_multi_pc) {
            for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
                readBarSum->cbReadCountExact[ii]++;//add one to exact counts
            };
        };
    };

    // if restarting from _STARtmp/solo* file
    if (P.runRestart.type==1) {//this could happen if the run is restarted. Would be better to save/load cbReadCount, or recalculate it from
        for (int ii=0; ii<P.runThreadN; ii++) {
            readFeatAll[ii]->streamReads->clear(); //just in case EOF was reached in previous reading
            readFeatAll[ii]->streamReads->seekg(0,ios::beg);
            string line1;
            while (std::getline(*readFeatAll[ii]->streamReads, line1)) {
                istringstream line1stream(line1);
                uint64 cb1;            
                line1stream >> cb1 >> cb1 >> cb1;
                if (featureType==SoloFeatureTypes::SJ)
                    line1stream >> cb1;
                line1stream >> cb1;
                //if (cb1>readFeatSum->cbReadCount.size())
                //    continue;//this should not happen!
                readFeatSum->cbReadCount[cb1]++;
            };
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
    
    indCBwl.resize(pSolo.cbWLsize, (uint32) -1);
    indCB.resize(nCB);
    nCB=0;//will count it again below
    for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
        if (readFeatSum->cbReadCount[ii]>0) {
            indCB[nCB]=ii;
            indCBwl[ii]=nCB;
            ++nCB;
        };
    };
    
};
    
