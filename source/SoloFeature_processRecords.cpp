#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"

void SoloFeature::processRecords(ReadAlignChunk **RAchunk) 
{
    if (pSolo.type==0)
        return;
 
    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Starting Solo post-map for " <<pSolo.featureNames[featureType] <<endl;    
    
    for (int ii=0; ii<P.runThreadN; ii++) {//point to 
        readFeatAll[ii]=RAchunk[ii]->RA->soloRead->readFeat[featureType];
        readBarAll[ii]=RAchunk[ii]->RA->soloRead->readBar;        
    };

    for (int ii=0; ii<P.runThreadN; ii++) {
        readFeatSum->addCounts(*readFeatAll[ii]);
        readBarSum->addCounts(*readBarAll[ii]);        
    };
       
    //allocate arrays to store CB/gene/UMIs for all reads
    nCB=0;nReadsMapped=0;
    for (uint32 ii=0; ii<pSolo.cbWL.size(); ii++) {
        if (readBarSum->cbReadCountExact[ii]>0) {
            nCB++;
            nReadsMapped += readFeatSum->cbReadCount[ii];
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
        if (readBarSum->cbReadCountExact[ii]>0) {//if no exact matches, this CB is not present
            indCB[nCB]=ii;
            rCBp[nCB+1] = rCBp[nCB] + 2*readFeatSum->cbReadCount[ii];
            ++nCB;
        };
        rCBpa[ii+1]=rCBp[nCB];
    };
    
    //read and store the CB/gene/UMI from files
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished allocating arrays for Solo " << nReadsMapped*2.0*8/1024/1024/1024 <<" GB" <<endl;

    for (int ii=0; ii<P.runThreadN; ii++) {//TODO: this can be parallelized
        readFeatAll[ii]->inputRecords(rCBpa,readBarSum->cbReadCountExact);
    };
    
    for (uint32 iCB=0; iCB<nCB; iCB++) {
        uint64 nr=(rCBpa[indCB[iCB]]-rCBp[iCB])/2;  //number of reads that were matched to WL, rCBpa accumulated reference to the last element+1
        if (nr>nReadPerCBmax)
            nReadPerCBmax=nr;
        readFeatSum->stats.V[readFeatSum->stats.nMatch] += nr;        
    };
    
    for (int ii=0; ii<P.runThreadN; ii++) {
        readFeatSum->addStats(*readFeatAll[ii]);
        readBarSum->addStats(*readBarAll[ii]);
    };    
        
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished reading reads from Solo files nCB="<<nCB <<", nReadPerCBmax="<<nReadPerCBmax;
    P.inOut->logMain <<", nMatch="<<readFeatSum->stats.V[readFeatSum->stats.nMatch]<<endl;
    
    //collapse each CB
    nUperCB = new uint32[nCB];//record pair: nUMIs per CB and iCB, for sorting if needed
    nGperCB = new uint32[nCB];
    uint32 *umiArray = new uint32[nReadPerCBmax*umiArrayStride];
    nCellGeneEntries=0;
    
    for (uint32 iCB=0; iCB<nCB; iCB++) {
        uint64 nr=(rCBpa[indCB[iCB]]-rCBp[iCB])/2; //number of reads that were matched to WL, rCBpa accumulated reference to the last element+1
        collapseUMI(rCBp[iCB],nr,nGperCB[iCB],nUperCB[iCB],umiArray);
        readFeatSum->stats.V[readFeatSum->stats.nUMIs] += nUperCB[iCB];
        if (nGperCB[iCB]>0)
            ++readFeatSum->stats.V[readFeatSum->stats.nCellBarcodes];
        nCellGeneEntries += nGperCB[iCB];        
    };
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished collapsing UMIs" <<endl;
    
    *statsStream << setw(50)<< "Barcodes:\n";
    readBarSum->statsOut(*statsStream);
    *statsStream << setw(50)<< pSolo.featureNames[featureType] <<":\n";
    readFeatSum->statsOut(*statsStream);

    //output nU per gene per CB
    outputResults();

};
