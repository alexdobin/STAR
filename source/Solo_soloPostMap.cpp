#include "Solo.h"
#include "streamFuns.h"
#include "TimeFunctions.h"

void Solo::soloPostMap(ReadAlignChunk **RAchunk) {
    
    for (int ii=0; ii<P.runThreadN; ii++) {
        soloCBall[ii]=RAchunk[ii]->RA->soloCB;
    };

    //summary statistics
    for (int ii=0; ii<P.runThreadN; ii++) {
        soloCBsum->addSoloCBcounts(*soloCBall[ii]);
    };
       
    //allocate arrays to store CB/gene/UMIs for all reads
    nCB=0;nReadsMapped=0;
    for (uint32 ii=0; ii<pSolo.cbWL.size(); ii++) {
        if (soloCBsum->cbReadCountExact[ii]>0) {
            nCB++;
            nReadsMapped += soloCBsum->cbReadCount[ii];
        };
    };

    rGeneUMI = new uint32[2*nReadsMapped]; //big array for all CBs - each element is gene and UMI
    rCBp = new uint32*[nCB+1];
    uint32 **rCBpa = new uint32*[pSolo.cbWL.size()+1];
    indCB = new uint32[nCB];
    
    uint32 nReadPerCBmax=0;
    rCBp[0]=rGeneUMI;
    rCBpa[0]=rGeneUMI;
    nCB=0;//will count it again below
    for (uint32 ii=0; ii<pSolo.cbWL.size(); ii++) {
        if (soloCBsum->cbReadCountExact[ii]>0) {//if no exact matches, this CB is not present
            indCB[nCB]=ii;
            rCBp[nCB+1] = rCBp[nCB] + 2*soloCBsum->cbReadCount[ii];
            ++nCB;
        };
        rCBpa[ii+1]=rCBp[nCB];
    };
    
    //read and store the CB/gene/UMI from files
    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished allocating arrays for Solo " << nReadsMapped*2.0*8/1024/1024/1024 <<" GB" <<endl;

    for (int ii=0; ii<P.runThreadN; ii++) {//TODO: this can be parallelized
        soloCBall[ii]->readCBgeneUMIfromFiles(rCBpa,soloCBsum->cbReadCountExact);
    };
    
    for (uint32 iCB=0; iCB<nCB; iCB++) {
        uint64 nr=(rCBpa[indCB[iCB]]-rCBp[iCB])/2;  //number of reads that were matched to WL, rCBpa accumulated reference to the last element+1
        if (nr>nReadPerCBmax)
            nReadPerCBmax=nr;
        soloCBsum->stats.V[soloCBsum->stats.nMatch] += nr;        
    };
    
    for (int ii=0; ii<P.runThreadN; ii++) {
        soloCBsum->addSoloCBstats(*soloCBall[ii]);
    };    
        
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished reading reads from Solo files nCB="<<nCB <<", nReadPerCBmax="<<nReadPerCBmax;
    P.inOut->logMain <<", nMatch="<<soloCBsum->stats.V[soloCBsum->stats.nMatch]<<endl;
    
    //collapse each CB
    nUperCB = new uint32[nCB];
    nGperCB = new uint32[nCB];
    uint32 umiArray[nReadPerCBmax*umiArrayStride];
    
    for (uint32 iCB=0; iCB<nCB; iCB++) {
        nUperCB[2*iCB+1]=iCB;
        uint64 nr=(rCBpa[indCB[iCB]]-rCBp[iCB])/2; //number of reads that were matched to WL, rCBpa accumulated reference to the last element+1
        collapseUMI(rCBp[iCB],nr,nGperCB[iCB],nUperCB[iCB],umiArray);
        soloCBsum->stats.V[soloCBsum->stats.nCellBarcodes] += nUperCB[iCB];
        nCellGeneEntries += nGperCB[iCB];        
    };
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished collapsing UMIs" <<endl;
    soloCBsum->stats.V[soloCBsum->stats.nCellBarcodes]=nCB;
    soloCBsum->statsOut(*soloStatsStream);

    //output nU per gene per CB
    outputNumUMIperGeneCB();

};
