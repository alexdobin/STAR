#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "systemFunctions.h"

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
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished allocating arrays for Solo " << nReadsMapped*rguStride*4.0/1024/1024/1024 <<" GiB" <<endl;

    ///////////////////////////////////////////////////////////////////////////
    ////////////// Input records
    readFlagCounts.flagCounts.reserve(nCB*3/2);
    readFlagCounts.flagCountsNoCB = {};
    vector<uint32> nReadPerCBunique1(pSolo.cbWLsize), nReadPerCBmulti1(pSolo.cbWLsize); //temp arrays to record # of reads for all cells in the WL
    if (P.runThreadN==1) {
        for (int ii=0; ii<P.runThreadN; ii++) {
            readFeatAll[ii]->inputRecords(rCBpa, rguStride, readBarSum->cbReadCountExact, readInfo, readFlagCounts, nReadPerCBunique1, nReadPerCBmulti1);
            readFeatSum->addStats(*readFeatAll[ii]);//sum stats: has to be done after inputRecords, since the stats values are updated there
        };
    } else {//parse the files concurrently, then make the shared writes in file order, as the serial loop does
        vector<vector<SoloRecordBuf>> recs(P.runThreadN);
        vector<SoloReadFlagClass> flagCountsFile(P.runThreadN);
        vector<vector<uintCB>> flagCountsOrder(P.runThreadN);

        #pragma omp parallel for num_threads(P.runThreadN) schedule(dynamic,1)
        for (int ii=0; ii<P.runThreadN; ii++)
            readFeatAll[ii]->inputRecordsBuffered(readBarSum->cbReadCountExact, readInfo, flagCountsFile[ii], recs[ii], flagCountsOrder[ii]);

        for (int ii=0; ii<P.runThreadN; ii++) {
            readFeatAll[ii]->replayRecords(rCBpa, rguStride, recs[ii], nReadPerCBunique1, nReadPerCBmulti1);
            vector<SoloRecordBuf>().swap(recs[ii]);

            //CBs enter flagCounts in the serial loop's order, so the unordered_map is built the same way
            for (auto cb : flagCountsOrder[ii]) {
                auto &cIn = flagCountsFile[ii].flagCounts.find(cb)->second;
                auto &cOut = readFlagCounts.flagCounts.insert({cb, {} }).first->second;
                for (uint32 ibit=0; ibit<readFlagCounts.nBits; ibit++)
                    cOut[ibit] += cIn[ibit];
            };
            readFlagCounts.countsAddNoCBarray(flagCountsFile[ii].flagCountsNoCB);

            readFeatSum->addStats(*readFeatAll[ii]);
        };
    };
    readFlagCounts.countsAddNoCBarray(readFeatSum->readFlag.flagCountsNoCB);//add no-CB counts calculated in SoloReadFeature_record.cpp and not recorded to temp Solo files

    nReadPerCBtotal.resize(nCB);
    nReadPerCBunique.resize(nCB);
    for (uint32 icb=0; icb<nCB; icb++) {
        nReadPerCBunique[icb] = nReadPerCBunique1[indCB[icb]];
        nReadPerCBtotal[icb] = nReadPerCBunique[icb] + nReadPerCBmulti1[indCB[icb]];
    };

    //debug
    /*{
        uint64 n1=0,n2=0;
        for (uint32 icb=0; icb<nCB; icb++) {
            n1 += nReadPerCBtotal[icb];
        };

        for (auto &c: readFlagCounts.flagCounts) {
            n2+= c.second[readFlagCounts.countedU] + c.second[readFlagCounts.countedM];
        };

        cout << "n1,2=" << n1<<" "<<n2<<endl;
    };*/

    nReadPerCB.resize(nCB);
    nReadPerCBmax=0;
    for (uint32 iCB=0; iCB<nCB; iCB++) {
        nReadPerCB[iCB] = (rCBpa[indCB[iCB]]-rCBp[iCB])/rguStride;  //number of reads that were matched to WL, rCBpa accumulated reference to the last element+1
                                                                    //for multimappers this is the number of all alignments > number of reads
        nReadPerCBmax=max(nReadPerCBmax,nReadPerCB[iCB]);
        //readFeatSum->stats.V[readFeatSum->stats.yesWLmatch] += nReadPerCB[iCB];
    };    
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished reading reads from Solo files nCB="<<nCB <<", nReadPerCBmax="<<nReadPerCBmax;
    P.inOut->logMain <<", yesWLmatch="<<readFeatSum->stats.V[readFeatSum->stats.yesWLmatch]<<endl;
    
    //////////////////////////////////////////////////////////////////////////////
    /////////////////////////// collapse each CB
    nUMIperCB.resize(nCB);
    nGenePerCB.resize(nCB);
    
                     //dedup options        //gene ID
    countMatStride = pSolo.umiDedup.yes.N + 1;
    //countCellGeneUMI and countMatMult.m are sized in collapseUMIall
    countCellGeneUMIindex.resize(nCB+1, 0);
    
    if (pSolo.multiMap.yes.multi) {
                    //gene   //uniform  //rescue
        countMatMult.s = 1 + pSolo.multiMap.yes.N * pSolo.umiDedup.yes.N;
        countMatMult.i.resize(nCB+1, 0);
    };

    collapseUMIall();
        
    P.inOut->logMain << "RAM for solo feature "<< SoloFeatureTypes::Names[featureType] <<"\n"
                     <<  linuxProcMemory() << flush;        
    delete[] rGeneUMI;
    delete[] rCBp;
    delete[] rCBpa;
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished collapsing UMIs" <<endl;
};
