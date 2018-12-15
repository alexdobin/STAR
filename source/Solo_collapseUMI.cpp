#include "Solo.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "serviceFuns.cpp"

void collapseUMIwith1MMlowHalf(uint32 *rGU, uint32 umiMaskLow, uint32 nU0, uint32 &nU1, uint32 &nU2) {
    
    const uint32 bitTop=1<<31;
    const uint32 bitTop1=1<<30;
    const uint32 bitMaskTop2bits=(~bitTop) & (~bitTop1);
    
    for (uint32 iu=0; iu<2*nU0; iu+=2) {//each UMI
        uint32 iuu=iu+2;
        for (; iuu<2*nU0; iuu+=2) {//compare to all UMIs down

            if ( (rGU[iuu+1] & bitTop & bitTop1) > 0)
                continue;//this one was already found duplicated for both collapse types

            uint32 uuXor=rGU[iu] ^ rGU[iuu];

//             if (uuXor==0)
//                 exit(1); //debug

            if ( uuXor > umiMaskLow)
                break; //upper half is different

            if (uuXor >> (__builtin_ctz(uuXor)/2)*2 > 3) //shift by even number of trailing zeros
                continue;//>1MM

            //1MM UMI
            if ( (rGU[iuu+1] & bitTop) == 0) {
                rGU[iuu+1] |= bitTop;
                --nU1;//subtract the duplicated UMIs
            } else if ( (rGU[iu+1] & bitTop) == 0) {
                rGU[iu+1] |= bitTop;
                --nU1;//subtract the duplicated UMIs                
            };
            if ( (rGU[iuu+1] & bitTop1) == 0 && (rGU[iu+1] & bitMaskTop2bits)>(2*(rGU[iuu+1] & bitMaskTop2bits)-1) ) {//iuu is duplicate of iu
                rGU[iuu+1] |= bitTop1;
                --nU2;//subtract the duplicated UMIs
            } else if ( (rGU[iu+1] & bitTop1) == 0 && (rGU[iuu+1] & bitMaskTop2bits)>(2*(rGU[iu+1] & bitMaskTop2bits)-1) ) {//iu is duplicate of iuu
                rGU[iu+1] |= bitTop1;
                --nU2;//subtract the duplicated UMIs
            };            
        };
    };
};
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
    for (uint32 iR=0; iR<2*rN; iR+=2) {
        if (rGU[iR]!=gid1) {//record gene boundary
            gReadS[nGenes]=iR;
            gid1=rGU[iR];
            gID[nGenes]=gid1;            
            ++nGenes;
        };
        rGU[iR]=rGU[iR+1]; //shift UMIs
        //rGU[iR+1] storage this will be used later for counting
    };
    gReadS[nGenes]=2*rN;//so that gReadS[nGenes]-gReadS[nGenes-1] is the number of reads for nGenes
 
    uint32 *nUg = new uint32[nGenes*3];//3 types of counts
    nUtot=0;
    for (uint32 iG=0; iG<nGenes; iG++) {//collapse UMIs for each gene
        uint32 *rGU1=rGU+gReadS[iG];
        
        qsort(rGU1, (gReadS[iG+1]-gReadS[iG])/2, 2*sizeof(uint32), funCompareNumbers<uint32>);

        //exact collapse
        uint32 iR1=-2; //number of distinct UMIs for this gene
        uint32 u1=-1;
        for (uint32 iR=0; iR<gReadS[iG+1]-gReadS[iG]; iR+=2) {//count and collapse identical UMIs
            if (rGU1[iR]!=u1) {
                iR1 += 2;
                u1=rGU1[iR];                
                rGU1[iR1]=u1;
                rGU1[iR1+1]=0;                
            };
            rGU1[iR1+1]++;             
        };
        uint32 nU0=(iR1+2)/2;
        
        //collapse with 1MM
        uint32 nU1=nU0, nU2=nU0;//2 types of 1MM collapsing

        collapseUMIwith1MMlowHalf(rGU1,pSolo.umiMaskLow, nU0, nU1, nU2);
        
        //exchange low and high half of UMIs, re-sort, and look for 1MM again
        for (uint32 iu=0; iu<2*nU0; iu+=2) {
            uint32 high=rGU1[iu]>>(pSolo.umiL);
            rGU1[iu] &= pSolo.umiMaskLow; //remove high
            rGU1[iu] <<= (pSolo.umiL); //move low to high
            rGU1[iu] |= high; //add high
            rGU1[iu+1] &= ~(3<<30); //remove the duplication flags (top 2 bits)
        };
        qsort(rGU1, nU0, 2*sizeof(uint32), funCompareNumbers<uint32>);
        collapseUMIwith1MMlowHalf(rGU1,pSolo.umiMaskLow, nU0, nU1, nU2);
                
        nUg[3*iG]=nU0;
        nUg[3*iG+1]=nU1;
        nUg[3*iG+2]=nU2;
        nUtot+=nU1;//TODO user makes the choice 
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
