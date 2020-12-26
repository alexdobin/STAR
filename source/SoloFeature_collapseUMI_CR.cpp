#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "serviceFuns.cpp"
#include <unordered_map>
#include "SoloCommon.h"
#include <bitset>

#define def_MarkNoColor  (uint32) -1

inline int funCompareSolo1 (const void *a, const void *b) {
    uint32 *va= (uint32*) a;
    uint32 *vb= (uint32*) b;

    if (va[1]>vb[1]) {
        return 1;
    } else if (va[1]<vb[1]) {
        return -1;
    } else if (va[0]>vb[0]){
        return 1;
    } else if (va[0]<vb[0]){
        return -1;
    } else {
        return 0;
    };
};

void SoloFeature::collapseUMI_CR(uint32 iCB, uint32 *umiArray) 
{
                                 
    uint32 *rGU=rCBp[iCB];
    uint32 rN=nReadPerCB[iCB]; 
    
    qsort(rGU,rN,rguStride*sizeof(uint32),funCompareNumbers<uint32>); //sort by gene index

    //compact reads per gene
    uint32 gid1=-1;//current gID
    uint32 nGenes=0; //number of genes
    uint32 *gID = new uint32[min(featuresNumber,rN)+1]; //gene IDs
    uint32 *gReadS = new uint32[min(featuresNumber,rN)+1]; //start of gene reads TODO: allocate this array in the 2nd half of rGU
    for (uint32 iR=0; iR<rN*rguStride; iR+=rguStride) {
        if (rGU[iR+rguG]!=gid1) {//record gene boundary
            gReadS[nGenes]=iR;
            gid1=rGU[iR+rguG];
            gID[nGenes]=gid1;
            ++nGenes;
        };
    };
    gReadS[nGenes]=rguStride*rN;//so that gReadS[nGenes]-gReadS[nGenes-1] is the number of reads for nGenes, see below in qsort

    //unordered_map<uint32, uint32> umiMaxGeneCount;//for each umi, max counts of reads per gene
    
    unordered_map <uintUMI, unordered_map<uint32,uint32>> umiGeneHash, umiGeneHash0;
                   //UMI                 //Gene //Count
    
    vector<unordered_map <uintUMI,uintUMI>> umiCorrected(nGenes);

    if (countCellGeneUMI.size() < countCellGeneUMIindex[iCB] + nGenes*countMatStride)
        countCellGeneUMI.resize((countCellGeneUMI.size() + nGenes*countMatStride )*2);//allocated vector too small
    
    nGenePerCB[iCB]=0;
    nUMIperCB[iCB]=0;
    countCellGeneUMIindex[iCB+1]=countCellGeneUMIindex[iCB];
    /////////////////////////////////////////////
    /////////// main cycle over genes
    for (uint32 iG=0; iG<nGenes; iG++) {//collapse UMIs for each gene
        uint32 *rGU1=rGU+gReadS[iG];

        if (gReadS[iG+1]==gReadS[iG])
            continue; //no reads

        qsort(rGU1, (gReadS[iG+1]-gReadS[iG])/rguStride, rguStride*sizeof(uint32), funCompareTypeShift<uint32,rguU>);
        
        //exact collapse
        uint32 iR1=-umiArrayStride; //number of distinct UMIs for this gene
        uint32 u1=-1;
        for (uint32 iR=rguU; iR<gReadS[iG+1]-gReadS[iG]; iR+=rguStride) {//count and collapse identical UMIs
            
            if (rGU1[iR]!=u1) {
                iR1 += umiArrayStride;
                u1=rGU1[iR];
                umiArray[iR1]=u1;
                umiArray[iR1+1]=0;
            };
            umiArray[iR1+1]++;
            //if ( umiArray[iR1+1]>nRumiMax) nRumiMax=umiArray[iR1+1];
        };

        uint32 nU0=(iR1+umiArrayStride)/umiArrayStride;
        uint32 nU1=nU0;//2 types of 1MM collapsing
       
        if (pSolo.umiFiltering.MultiGeneUMI) {
            for (uint64 iu=0; iu<nU0*umiArrayStride; iu+=umiArrayStride) {
                umiGeneHash0[umiArray[iu+0]][iG]+=umiArray[iu+1];//this sums read counts over UMIs that were collapsed
            };
        };       
       
        qsort(umiArray, nU0, umiArrayStride*sizeof(uint32), funCompareSolo1);
        for (uint64 iu=0; iu<(nU0-1)*umiArrayStride; iu+=umiArrayStride) {
            uint64 iuu;
            for (iuu=(nU0-1)*umiArrayStride; iuu>iu; iuu-=umiArrayStride) {

                uint32 uuXor=umiArray[iu+0] ^ umiArray[iuu+0];

                if ( (uuXor >> (__builtin_ctz(uuXor)/2)*2) <= 3 ) {//1MM
                    if (readInfo.size()>0) {//record corrections
                        umiCorrected[iG][umiArray[iu+0]]=umiArray[iuu+0];//replace iu with iuu
                    };                    
                    umiArray[iu+0]=umiArray[iuu+0];//replace iu with iuu
                    break;
                };
            };
        };

        if (pSolo.umiFiltering.MultiGeneUMI) {
            for (uint64 iu=0; iu<nU0*umiArrayStride; iu+=umiArrayStride) {
                umiGeneHash[umiArray[iu+0]][iG]+=umiArray[iu+1];//this sums read counts over UMIs that were collapsed
            };

        } else {//no MultiGeneUMI
            qsort(umiArray, nU0, umiArrayStride*sizeof(uint32), funCompareNumbers<uint32>);
            nU1=1;
            for (uint64 iu=umiArrayStride; iu<nU0*umiArrayStride; iu+=umiArrayStride) {
                if (umiArray[iu+0]!=umiArray[iu+0-umiArrayStride]) {
                    nU1++;
                };
            };
            nGenePerCB[iCB]++;
            nUMIperCB[iCB]+=nU1;
            countCellGeneUMI[countCellGeneUMIindex[iCB+1] + 0] = gID[iG];
            countCellGeneUMI[countCellGeneUMIindex[iCB+1] + pSolo.umiDedup.countInd.CR] = nU1;
            countCellGeneUMIindex[iCB+1] = countCellGeneUMIindex[iCB+1] + countMatStride;//iCB+1 accumulates the index

            if (readInfo.size()>0) {//record cb/umi for each read
                for (uint32 iR=0; iR<gReadS[iG+1]-gReadS[iG]; iR+=rguStride) {//cycle over reads
                    uint64 iread1 = rGU1[iR+rguR];
                    readInfo[iread1].cb = indCB[iCB] ;
                    uint32 umi=rGU1[iR+rguU];
                    
                    if (umiCorrected[iG].count(umi)>0)
                        umi=umiCorrected[iG][umi]; //correct UMI
                    readInfo[iread1].umi=umi;
                };      
            };            
        };

    };

    if (pSolo.umiFiltering.MultiGeneUMI) {
        
        vector<uint32> geneCounts(nGenes,0);
        
        vector<unordered_set<uintUMI>> geneUmiHash;
        if (readInfo.size()>0)
            geneUmiHash.resize(nGenes);
        
        for (auto &iu: umiGeneHash) {//loop over UMIs for all genes
                       
            uint32 maxu=0, maxg=-1;
            for (auto &ig : iu.second) {
                if (ig.second>maxu) {
                    maxu=ig.second;
                    maxg=ig.first;
                } else if (ig.second==maxu) {
                    maxg=-1;
                };
            };

            if ( maxg+1==0 )
                continue; //this umi is not counted for any gene, because two genes have the same read count for this UMI
            
            for (auto &ig : umiGeneHash0[iu.first]) {//check that this umi/gene had also top count for uncorrected umis
                if (ig.second>umiGeneHash0[iu.first][maxg]) {
                    maxg=-1;
                    break;
                };
            };

            if ( maxg+1!=0 ) {//this UMI is counted
                geneCounts[maxg]++;
                if (readInfo.size()>0)
                    geneUmiHash[maxg].insert(iu.first);
            };
        };

        for (uint32 ig=0; ig<nGenes; ig++) {
            if (geneCounts[ig] == 0)
                continue; //no counts for this gene
            nGenePerCB[iCB]++;
            nUMIperCB[iCB] += geneCounts[ig];
            countCellGeneUMI[countCellGeneUMIindex[iCB+1] + 0] = gID[ig];
            countCellGeneUMI[countCellGeneUMIindex[iCB+1] + pSolo.umiDedup.countInd.CR] = geneCounts[ig];
            countCellGeneUMIindex[iCB+1] = countCellGeneUMIindex[iCB+1] + countMatStride;//iCB+1 accumulates the index
        };
        
        if (readInfo.size()>0) {//record cb/umi for each read
            for (uint32 iG=0; iG<nGenes; iG++) {//cycle over genes
                uint32 *rGU1=rGU+gReadS[iG];            
                for (uint32 iR=0; iR<gReadS[iG+1]-gReadS[iG]; iR+=rguStride) {//cycle over reads
                    uint64 iread1 = rGU1[iR+rguR];
                    readInfo[iread1].cb = indCB[iCB] ;
                    uint32 umi=rGU1[iR+rguU];
                    
                    if (umiCorrected[iG].count(umi)>0)
                        umi=umiCorrected[iG][umi]; //correct UMI

                    //cout << iG << "-" << iR << " " <<flush ;
                    if (geneUmiHash[iG].count(umi)>0) {
                        readInfo[iread1].umi=umi;
                    } else {
                        readInfo[iread1].umi=(uintUMI) -1;
                    };
                };
            };
        };
    };
};
