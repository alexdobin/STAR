#include "ReadAlign.h"
#include "BAMfunctions.h"

void ReadAlign::chimericDetectionPEmerged(ReadAlign &seRA) {

    chimRecord=false;
    if (P.pCh.segmentMin==0) {//no chimeric detection requested
        return;
    };

    if (P.pCh.multimapNmax==0) {

        seRA.multMapSelect(); //this needs to be done for ChimericDetectionOld, may not need it for the new algorithm
        seRA.mappedFilter();

        chimRecord=seRA.chimericDetectionOld();
        if (!chimRecord) {
            return;
        };

        //convert merged into PE   
        for (uint ii=0; ii<2; ii++) {
            trChim[ii]=*trInit;
            trChim[ii].peOverlapSEtoPE(peOv.mateStart,seRA.trChim[ii]);
        };

        uint segLen[2][2]; //segment length [trChim][mate]
        uint segEx[2];//last exon of the mate0 [trChim]
        uint segLmin=-1LLU, i1=0,i2=0;    
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
                    i1=ii;//trChim of the shortest segment length
                    i2=jj;//mate of the shortest segment length
                };
            };
        };

        if (i2==1) {//eliminate mate1: simply cut the exons that belong to mate1
            trChim[i1].nExons=segEx[i1]+1;
        } else {//eliminate mate 0: shift mate1 exon to the beginning
            for (uint iex=0; iex<trChim[i1].nExons; iex++) {
                uint iex1=iex+segEx[i1]+1;
                for (uint ii=0; ii<EX_SIZE; ii++) {
                    trChim[i1].exons[iex][ii]=trChim[i1].exons[iex1][ii];
                };
                trChim[i1].canonSJ[iex]=trChim[i1].canonSJ[iex1];
                trChim[i1].sjAnnot[iex]=trChim[i1].sjAnnot[iex1];
                trChim[i1].sjStr[iex]=trChim[i1].sjStr[iex1];
                trChim[i1].shiftSJ[iex][0]=trChim[i1].shiftSJ[iex1][0];
                trChim[i1].shiftSJ[iex][1]=trChim[i1].shiftSJ[iex1][1];
            };
            trChim[i1].nExons=trChim[i1].nExons-segEx[i1]-1;
        };

        chimericDetectionOldOutput();
    
    } else if (trBest->maxScore <= (int) (readLength[0]+readLength[1]) - (int) P.pCh.nonchimScoreDropMin) {//require big enough drop in the best score
        chimRecord=seRA.chimDet->chimericDetectionMult(seRA.nW, seRA.readLength);
    };    
    
    if ( chimRecord ) {
        statsRA.chimericAll++;    
    };    
    
    return;
};

