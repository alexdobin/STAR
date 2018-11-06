#include "Solo.h"
#include "streamFuns.h"
#include "TimeFunctions.h"

void Solo::soloPostMap(ReadAlignChunk **RAchunk) {
    
    for (int ii=0; ii<P.runThreadN; ii++) {
        soloCBall[ii]=RAchunk[ii]->RA->soloCB;
    };

    //summary statistics
    for (int ii=0; ii<P.runThreadN; ii++) {
        soloCBsum->addSoloCB(*soloCBall[ii]);
    };
    soloCBsum->statsOut(*soloStatsStream);
       
    //allocate arrays to store CB/gene/UMIs for all reads
    nReadsMapped=soloCBsum->stats.V[soloCBsum->stats.nMatch];
    rGeneUMI = new uint32[2*nReadsMapped]; 

    nCB=0;
    for (uint32 ii=0; ii<pSolo.cbWL.size()-1; ii++)
        nCB += (soloCBsum->cbReadCount[ii]>0 ? 1 : 0 );

    rCBp = new uint32*[nCB];
    uint32 **rCBpa = new uint32*[pSolo.cbWL.size()];
    rCBn = new uint32[nCB];
    indCB = new uint32[nCB];
    rCBp[0]=rGeneUMI;
    rCBpa[0]=rGeneUMI;
    nCB=0;//will count it again below
    for (uint32 ii=0; ii<pSolo.cbWL.size()-1; ii++) {
        if (soloCBsum->cbReadCount[ii]>0) {
            indCB[nCB]=ii;
            rCBn[nCB]=soloCBsum->cbReadCount[ii];
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
        soloCBall[ii]->readCBgeneUMIfromFiles(rCBpa);
    };
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished reading reads from Solo files nCB="<<nCB <<endl;
    
    //collapse each CB
    nUperCB = new uint32[2*nCB];
    nGperCB = new uint32[nCB];
    for (uint32 iCB=0; iCB<nCB; iCB++) {
        nUperCB[2*iCB+1]=iCB;
        collapseUMI(iCB,nGperCB[iCB],nUperCB[2*iCB]);
    };
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished collapsing UMIs" <<endl;

    //output nU per gene per CB
    outputNumUMIperGeneCB();
    

    
//     statsStream.close();
};
