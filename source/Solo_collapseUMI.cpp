#include "Solo.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "serviceFuns.cpp"

// uint32 searchSortedForChange(uint32 *A, uint32 nA, uint32 strideA, uint32 stepMax) {//find the next element in the sorted list
//     uint32 iP2=0;
//     while (iP2<nA-1 && A[iP2]==A[0] ) //rough search
//         iP2 += strideA*stepMax;
//     iP2=min(iP2,nA-1);
//     
//     uint32 iP1=iP-stepMax;
//     uint32 iP;
//     while (iP1<iP2+1) {//binary search
//         iP=(iP1+iP2)/2;
//         if (A[iP*strideA]==A[0]) {
//             iP1=iP;
//         } else {
//             iP2=iP;
//         };
//     };
//     
//     return iP2; //returns the start of the next element, or the last element in the array
// };

void Solo::collapseUMI(uint32 iCB, uint32 &nGenes, uint32 &nUtot) {//iCB = CB to collapse, nReads=number of reads for this CB
    
    uint32 *rGU=rCBp[iCB];
    uint32 rN=rCBn[iCB];
    //sort 
    qsort(rGU,rN,2*sizeof(uint32),funCompareNumbers<uint32>); //sort by gene number
    
    //compact reads per gene
    uint32 gid1=-1;//current gID
    nGenes=0; //number of genes
    uint32 *gID = new uint32[min(Trans.nGe,rN)+1]; //gene IDS
    uint32 *gReadS = new uint32[min(Trans.nGe,rN)+1]; //start of gene reads TODO: allocate this array in the 2nd half of rGU
    for (uint32 iR=0; iR<rN; iR++) {
        if (rGU[iR*2]!=gid1) {//record gene boundary
            gReadS[nGenes]=iR;
            gid1=rGU[iR*2];
            gID[nGenes]=gid1;            
            ++nGenes;
        };
        rGU[iR]=rGU[iR*2+1]; //shift UMIs        
    };
    gReadS[nGenes]=rN;//so that gReadS[nGenes]-gReadS[nGenes-1] is the number of reads for nGenes
 
    uint32 *nUg = new uint32[nGenes*3];//3 types of counts
    nUtot=0;
    for (uint32 iG=0; iG<nGenes; iG++) {//collapse UMIs for each gene
        qsort(rGU+gReadS[iG],gReadS[iG+1]-gReadS[iG],sizeof(uint32),funCompareNumbers<uint32>);

        //exact collapse
        uint32 *nrU = new uint32[gReadS[iG+1]-gReadS[iG]]; //counts of reads per UMI
        uint32 nU0=-1; //number of distinct UMIs for this gene
        uint32 u1=-1;
        for (uint32 ir=gReadS[iG]; ir<gReadS[iG+1]; ir++) {//count and collapse identical UMIs
            if (rGU[ir]!=u1) {
                ++nU0;
                u1=rGU[ir];                
                rGU[nU0]=u1;
                nrU[nU0]=0;                
            };
            nrU[nU0]++;             
        };
        ++nU0;//nU0 was the index of the last added U
        
        //collapse with 1MM
        uint8 *dupU = new uint8 [nU0]; //TODO use bits in existing array
        for (uint32 iu=0; iu<nU0; iu++)
            dupU[iu]=0;
        uint32 nU1=nU0, nU2=nU0;//2 types of 1MM collapsing

        for (uint32 iu=0; iu<nU0; iu++) {//each UMI
            uint32 iuu=iu+1;
            for (; iuu<nU0; iuu++) {//compare to all UMIs down

                if (dupU[iuu]==3)
                    continue;//this one was already found duplicated for both collapse types

                uint32 uuXor=rGU[iu] ^ rGU[iuu];
                
                if (uuXor==0)
                    exit(1); //debug
                
                if ( uuXor > pSolo.umiMaskLow)
                    break; //upper half is different

                if (uuXor >> (__builtin_ctz(uuXor)/2)*2 > 3) //shift by even number of trailing zeros
                    continue;//>1MM
                
                //1MM UMI
                if ( (dupU[iuu] & 1) == 0) {
                    dupU[iuu] |=1;
                    --nU1;//subtract the duplicated UMIs
                };
                if ( (dupU[iuu] & 2) == 0 && nrU[iu]>(2*nrU[iuu]+1) ) {
                    dupU[iuu] |=2;
                    --nU2;//subtract the duplicated UMIs
                };
                if ( (dupU[iu] & 2) == 0 && nrU[iuu]>(2*nrU[iu]+1) ) {
                    dupU[iu] |=2;
                    --nU2;//subtract the duplicated UMIs
                };            
            };
        };
                
        nUg[3*iG]=nU0;
        nUg[3*iG+1]=nU1;
        nUg[3*iG+2]=nU2;
        nUtot+=nU0;
    };

    uint32 *rGUp=rGU;
    for (uint32 iG=0; iG<nGenes; iG++) {//output for all genes
        rGUp[0]=gID[iG];
        rGUp[1]=nUg[3*iG];
        if (nUg[3*iG]>1) {//record 2 more counts
            rGUp[2]=nUg[3*iG+1];
            rGUp[3]=nUg[3*iG+2];
            rGUp += 4;
        } else {//only one count recorded, save space
            rGUp += 2;
        };
    };
};
