#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"

void SoloFeature::countCBgeneUMI()
{    
    time_t rawTime;
    
    rguStride=2;
    if (pSolo.readIndexYes[featureType])
        rguStride=3; //to keep readI column

    if (pSolo.readInfoYes[featureType]) {
        readInfo.resize(nReadsInput,{(uint64)-1,(uint32)-1});
        time(&rawTime);
        P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Allocated and initialized readInfo array, nReadsInput = " << nReadsInput <<endl;        
    };

    rGeneUMI = new uint32[rguStride*nReadsMapped]; //big array for all CBs - each element is gene and UMI
    rCBp = new uint32*[nCB+1];
    uint32 **rCBpa = new uint32*[pSolo.cbWLsize+1];
    
    rCBp[0]=rGeneUMI;
    rCBpa[0]=rGeneUMI;
    nCB=0;//will count it again below
    for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
        if (readFeatSum->cbReadCount[ii]>0) {
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
        readFeatAll[ii]->inputRecords(rCBpa, rguStride, readBarSum->cbReadCountExact, readInfo);
    };

    nReadPerCB.resize(nCB);
    uint32 nReadPerCBmax=0;
    for (uint32 iCB=0; iCB<nCB; iCB++) {
        nReadPerCB[iCB] = (rCBpa[indCB[iCB]]-rCBp[iCB])/rguStride;  //number of reads that were matched to WL, rCBpa accumulated reference to the last element+1
        nReadPerCBmax=max(nReadPerCBmax,nReadPerCB[iCB]);
        //readFeatSum->stats.V[readFeatSum->stats.nMatch] += nReadPerCB[iCB];
    };    
    
    for (int ii=0; ii<P.runThreadN; ii++) {
        readFeatSum->addStats(*readFeatAll[ii]);
    };
    
    //readFeatSum->stats.calcUnique(pSolo.multiMap.yes.multi && (featureType==SoloFeatureTypes::Gene || featureType==SoloFeatureTypes::GeneFull));    

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished reading reads from Solo files nCB="<<nCB <<", nReadPerCBmax="<<nReadPerCBmax;
    P.inOut->logMain <<", nMatch="<<readFeatSum->stats.V[readFeatSum->stats.nMatch]<<endl;
    
    //////////////////////////////////////////////////////////////////////////////
    /////////////////////////// collapse each CB
    nUMIperCB.resize(nCB);
    nGenePerCB.resize(nCB);
    
    uint32 *umiArray = new uint32[nReadPerCBmax*umiArrayStride];//temp array for collapsing UMI
                     //dedup options        //gene ID
    countMatStride = pSolo.umiDedup.yes.N + 1;
    countCellGeneUMI.resize(nReadsMapped*countMatStride/5+16); //5 is heuristic, will be resized if needed
    countCellGeneUMIindex.resize(nCB+1, 0);
    
    if (pSolo.multiMap.yes.multi) {
                    //gene   //uniform  //rescue
        countMatMult.s = 1 + pSolo.multiMap.yes.N * pSolo.umiDedup.yes.N;
        countMatMult.m.resize(nReadsMapped*countMatMult.s/5+16);
        countMatMult.i.resize(nCB+1, 0);
    };
    
    nReadPerCBtotal.resize(nCB);
    nReadPerCBunique.resize(nCB);
    for (uint32 icb=0; icb<nCB; icb++) {//main collapse cycle
        
        collapseUMIall(icb, umiArray);
        
        readFeatSum->stats.V[readFeatSum->stats.nUMIs] += nUMIperCB[icb];
        if (nGenePerCB[icb]>0) //nGenePerCB contains only unique
            ++readFeatSum->stats.V[readFeatSum->stats.nCellBarcodes];
        
        readFeatSum->stats.V[readFeatSum->stats.nMatch] += nReadPerCBtotal[icb];        
        readFeatSum->stats.V[readFeatSum->stats.nMatchUnique ] += nReadPerCBunique[icb];        
    };
        
    delete[] rGeneUMI;
    //delete[] rCBp;
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished collapsing UMIs" <<endl;
};
