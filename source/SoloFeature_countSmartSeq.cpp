#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "serviceFuns.cpp"

void SoloFeature::countSmartSeq()
{    
    time_t rawTime;
    
    nCB=pSolo.cbWLsize; //all cells are recorded
    nReadPerCB.resize(nCB);
    cbFeatureUMImap.resize(nCB);
    ///////////////////////////////////////////////////////////////////////////
    ////////////// Input records
    for (int ii=0; ii<P.runThreadN; ii++) {//TODO: this can be parallelized
        readFeatAll[ii]->inputRecords(NULL, rguStride, readBarSum->cbReadCountExact, readInfo, this);
    };

    uint32 nReadPerCBmax=0;
    indCB = new uint32[nCB];
    for (uint32 iCB=0; iCB<nCB; iCB++) {
        nReadPerCBmax=max(nReadPerCBmax,nReadPerCB[iCB]);
        readFeatSum->stats.V[readFeatSum->stats.nMatch] += nReadPerCB[iCB];
        indCB[iCB]=iCB;//one to one correspondence, keep all cells
    };

    for (int ii=0; ii<P.runThreadN; ii++) {
        readFeatSum->addStats(*readFeatAll[ii]);
    };

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished reading reads from Solo files nCB="<<nCB <<", nReadPerCBmax="<<nReadPerCBmax;
    P.inOut->logMain <<", nMatch="<<readFeatSum->stats.V[readFeatSum->stats.nMatch]<<endl;
    
    //////////////////////////////////////////////////////////////////////////////
    /////////////////////////// collapse each CB
    nUMIperCB.resize(nCB);
    nGenePerCB.resize(nCB);
      
    countMatStride=2; //hard-coded for now, works for both Gene*/SJ and Velocyto
    
    uint64 ccgN=0;//total number of entries in the countCellGeneUMI
    for (auto &cbf :  cbFeatureUMImap) {
        ccgN+=cbf.size();
    };
    countCellGeneUMI.resize(ccgN*countMatStride);
    countCellGeneUMIindex.resize(nCB+1);
    countCellGeneUMIindex[0]=0;
    
    auto ind1=countCellGeneUMIindex[0];
    for (uint32 icb=0; icb<nCB; icb++) {//count UMIs
                
//         for ( uint32 ig=0; ig<Trans.nGe; ig++) {
//             if (cbFeatureUMImap[icb].count(ig)==0)
//                 continue;
//             nGenePerCB[icb]++;
//             uint64 numi=cbFeatureUMImap[icb][ig].size();
//             nUMIperCB[icb]+=numi;
//             
//             countCellGeneUMI[ind1] = ig;
//             countCellGeneUMI[ind1+1] = numi;
//             ind1+=countMatStride;
//         };
        
        for ( auto cbf=cbFeatureUMImap[icb].begin(); cbf!=cbFeatureUMImap[icb].end(); ++cbf ) {
            nGenePerCB[icb]++;
            uint64 numi=cbf->second.size();
            nUMIperCB[icb]+=numi;
            
            countCellGeneUMI[ind1] = cbf->first;
            countCellGeneUMI[ind1+1] = numi;
            ind1+=countMatStride;
        };
        
        uint32 *s1 = countCellGeneUMI.data()+countCellGeneUMIindex[icb]; //pointer to the start of this icb
        qsort((void*) s1, cbFeatureUMImap[icb].size(), countMatStride*sizeof(uint32), funCompareNumbers<uint32>);
        
        countCellGeneUMIindex[icb+1]=ind1;
                
        readFeatSum->stats.V[readFeatSum->stats.nUMIs] += nUMIperCB[icb];
        if (nGenePerCB[icb]>0)
            ++readFeatSum->stats.V[readFeatSum->stats.nCellBarcodes];
    };
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished collapsing UMIs" <<endl;
};
