#include "IncludeDefine.h"
#include "Parameters.h"

int stitchGapIndel (uint rAend, uint gAend, uint rBstart, uint gBstart, uint L, uint gapStart, uint gapEnd, char* R, char* G, Parameters& P,\
                    uint &iRbest, uint &nMM){//returns stitch score

    uint gapLength = gapEnd-gapStart+1;
    sint inDel= (sint) (gBstart-gAend-1) - (sint) gapLength - (sint) (rBstart-rAend-1); //>0: deletion; <0: insertion

    if (inDel==0) {//this should not happen, it should have been caught in the first stitching
        return -1;
    };
    int score2best;
    int score2;

    if (inDel>0) {//
        score2=0;
        score2best=-1;
        iRbest=0;
        for (uint iR=1; iR<rBstart-rAend; iR++) {//scan to find the best position iR of the indel
            uint iG1=gAend+iR;
            uint iG2=iG1+(uint) inDel;
            if (iG1>=gapStart) iG1 += gapLength;//exclude gap
            if (iG2>=gapStart) iG2 += gapLength;

            if  ( R[rAend+iR]==G[iG1] && R[rAend+iR]!=G[iG2] ) {
                score2++;
            } else if  ( R[rAend+iR]!=G[iG1] && R[rAend+iR]==G[iG2] ) {
                score2--;
            };

            if (score2>score2best) {
                score2best=score2;
                iRbest=iR;
            };
        };

        //score the alignment with inDel at iRbest
        nMM=0;
        score2= L - inDel*P.scoreDelBase - P.scoreDelOpen; //score B and deletion
        for (uint iR=1; iR<rBstart-rAend; iR++) {//scan to find the best position iR of the indel
            uint iG=gAend+iR;
            if (iR>iRbest) iG += (uint) inDel;
            if (iG>=gapStart) iG += gapLength;//exclude gap

            if  ( R[rAend+iR]==G[iG] ) {
                score2++;
            } else if (R[rAend+iR]!=G[iG] && R[rAend+iR]<4 && G[iG]<4) {//only penalize mismatches for non-N bases
                score2--;
                nMM++;
            };
        };

    } else {
        return -1;
    };

    return score2;
};
