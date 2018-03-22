#include "ReadAlign.h"
#include "BAMfunctions.h"

void ReadAlign::chimericDetectionPEmerged(ReadAlign &seRA) {

    chimRecord=false;
    if (P.pCh.segmentMin==0) {//no chimeric detection requested
        return;
    };
    
    seRA.multMapSelect(); //this needs to be done for ChimericDetectionOld, may not need it for the new algorithm
    seRA.mappedFilter();
    
    chimRecord=seRA.chimericDetectionOld();
    if (!chimRecord) {
        return;
    };

    statsRA.chimericAll++;

//     seRA.chimericDetectionPEmergedTrim();
    
    //convert merged into PE   
    for (uint ii=0; ii<2; ii++) {
        trChim[ii]=*trInit;
        trChim[ii].peOverlapSEtoPE(peOv.mateStart,seRA.trChim[ii]);
    };

    uint segLen[2][2], segEx[2];
    uint segLmin=-1LLU, i1,i2;    
    for (uint ii=0; ii<2; ii++) {
        segLen[ii][0]=0;
        segLen[ii][1]=0;
        for (uint iex=0; iex<trChim[ii].nExons; iex++) {
            if (trChim[ii].exons[iex][EX_iFrag]==trChim[ii].exons[0][EX_iFrag]) {
                segLen[ii][0]+=trChim[ii].exons[iex][EX_L];
                segEx[ii]=iex;
            } else {
                segLen[ii][1]+=trChim[ii].exons[iex][EX_L];
            };
        };
        for (uint jj=0; jj<2; jj++) {
            if (segLen[ii][jj]<segLmin) {
                segLmin=segLen[ii][jj];
                i1=ii;
                i2=jj;
            };
        };
    };
    
    if (i2==1) {
        trChim[i1].nExons=segEx[i1]+1;
    } else {
        for (uint iex=segEx[i1]+1; iex<trChim[i1].nExons; iex++) {
            for (uint ii=0; ii<EX_SIZE; ii++) {
                trChim[i1].exons[iex-segEx[i1]-1][ii]=trChim[i1].exons[iex][ii];
            };
        };
        trChim[i1].nExons=trChim[i1].nExons-segEx[i1]-1;
    };
    
    chimericDetectionOldOutput();
    
    return;
};


// void ReadAlign::chimericDetectionPEmergedTrim() {
//     
//    
// //     uint roSE[2][2];
// //     
// //     for (uint ii=0; ii<2; ii++) {
// //         if (trChim[ii].Str==0) {
// //             roSE[ii][0]=trChim[ii].exons[0][EX_R];
// //             roSE[ii][1]=trChim[ii].exons[trChim[ii].nExons-1][EX_R]+trChim[ii].exons[trChim[ii].nExons-1][EX_L]-1;  
// //         } else {
// //             roSE[ii][1]=Lread-trChim[ii].exons[0][EX_R]-1;
// //             roSE[ii][0]=Lread-(trChim[ii].exons[trChim[ii].nExons-1][EX_R]+trChim[ii].exons[trChim[ii].nExons-1][EX_L]);  
// //         };
// //     };
// //    
// //     uint chR = roSE[0][1]+1;//length of the left segment on the original read   
// //     if (chR<(Lread-peOv.nOv)) {//no trimming necessay, chim-junctions is outside of overlap
// //         return;
// //     };
// //     
// //     uint trim0=roSE[0][0];
// //     uint trim1=Lread-roSE[1][1]-1;
// //     
// //     uint seg0=min(chR-trim0, peOv.ovS+peOv.nOv)
//     
//     
//                 
// };
