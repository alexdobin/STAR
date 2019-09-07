#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"

void SoloFeature::countCBgeneUMI()
{    
    time_t rawTime;
    
    rguStride=2;
    if (pSolo.readInfoYes[featureType]) {
        rguStride=3; //to keep readI column
        readInfo.resize(nReadsInput,{(uint64)-1,(uint32)-1});
        time(&rawTime);
        P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Allocated and initialized readInfo array, nReadsInput = " << nReadsInput <<endl;        
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

    //allocate arrays to store CB/gene/UMIs for all reads
    nCB=0;nReadsMapped=0;
    for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
        if (readFeatSum->cbReadCount[ii]>0) {
            nCB++;
            nReadsMapped += readFeatSum->cbReadCount[ii];
        };
    };
    
    //pseudocounts
    if (pSolo.CBmatchWL.mm1_multi_pc) {
        for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
            readBarSum->cbReadCountExact[ii]++;
        };
    };

    rGeneUMI = new uint32[rguStride*nReadsMapped]; //big array for all CBs - each element is gene and UMI
    rCBp = new uint32*[nCB+1];
    uint32 **rCBpa = new uint32*[pSolo.cbWLsize+1];
    indCB = new uint32[nCB];

    uint32 nReadPerCBmax=0;
    rCBp[0]=rGeneUMI;
    rCBpa[0]=rGeneUMI;
    nCB=0;//will count it again below
    for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
        if (readFeatSum->cbReadCount[ii]>0) {
            indCB[nCB]=ii;
            rCBp[nCB+1] = rCBp[nCB] + rguStride*readFeatSum->cbReadCount[ii];
            ++nCB;
        };
        rCBpa[ii+1]=rCBp[nCB];
    };

    //read and store the CB/gene/UMI from files
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished allocating arrays for Solo " << nReadsMapped*rguStride*4.0/1024/1024/1024 <<" GB" <<endl;

    ///////////////////////////////////////////////////////////////////////////
    ////////////// Input records
    for (int ii=0; ii<P.runThreadN; ii++) {//TODO: this can be parallelized
        readFeatAll[ii]->inputRecords(rCBpa, rguStride, readBarSum->cbReadCountExact, streamTranscriptsOut, readInfo);
    };

    for (uint32 iCB=0; iCB<nCB; iCB++) {
        uint64 nr=(rCBpa[indCB[iCB]]-rCBp[iCB])/rguStride;  //number of reads that were matched to WL, rCBpa accumulated reference to the last element+1
        if (nr>nReadPerCBmax)
            nReadPerCBmax=nr;
        readFeatSum->stats.V[readFeatSum->stats.nMatch] += nr;
    };

    for (int ii=0; ii<P.runThreadN; ii++) {
        readFeatSum->addStats(*readFeatAll[ii]);
    };

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished reading reads from Solo files nCB="<<nCB <<", nReadPerCBmax="<<nReadPerCBmax;
    P.inOut->logMain <<", nMatch="<<readFeatSum->stats.V[readFeatSum->stats.nMatch]<<endl;
    
    if (featureType==SoloFeatureTypes::Transcript3p) {
        streamTranscriptsOut->flush();
        ofstream &outStr=ofstrOpen(P.outFileNamePrefix+pSolo.outFileNames[0]+"transcripts.tsv",ERROR_OUT, P);        
        for (uint32 ii=0; ii<Trans.nTr; ii++)
            outStr << Trans.trID[ii] <<"\t"<< Trans.trLen[ii] <<"\t"<< Trans.geName[Trans.trGene[ii]] << '\n';
        outStr.close();
        return; //the rest not implemented yet
    };
    
    //////////////////////////////////////////////////////////////////////////////
    /////////////////////////// collapse each CB
    nUMIperCB.resize(nCB);
    nGenePerCB.resize(nCB);
    nReadPerCB.resize(nCB);
    uint32 *umiArray = new uint32[nReadPerCBmax*umiArrayStride];
    
    for (uint32 icb=0; icb<nCB; icb++) {//main collapse cycle
        nReadPerCB[icb] =( rCBpa[indCB[icb]]-rCBp[icb] ) / rguStride; //number of reads that were matched to WL, rCBpa accumulated reference to the last element+1
        
        collapseUMI(rCBp[icb], nReadPerCB[icb], nGenePerCB[icb], nUMIperCB[icb], umiArray,indCB[icb]);
        
        readFeatSum->stats.V[readFeatSum->stats.nUMIs] += nUMIperCB[icb];
        if (nGenePerCB[icb]>0)
            ++readFeatSum->stats.V[readFeatSum->stats.nCellBarcodes];
    };

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished collapsing UMIs" <<endl;
};