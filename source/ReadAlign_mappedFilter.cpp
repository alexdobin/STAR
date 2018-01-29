#include "ReadAlign.h"

void ReadAlign::mappedFilter() {//filter mapped read, add to stats
    unmapType=-1;//mark as mapped
    if ( nW==0 ) {//no good windows
        statsRA.unmappedOther++;
        unmapType=0;
    } else if ( (trBest->maxScore < P.outFilterScoreMin) || (trBest->maxScore < (intScore) (P.outFilterScoreMinOverLread*(Lread-1))) \
              || (trBest->nMatch < P.outFilterMatchNmin)  || (trBest->nMatch < (uint) (P.outFilterMatchNminOverLread*(Lread-1))) ) {//too short
        statsRA.unmappedShort++;
        unmapType=1;
    } else if ( (trBest->nMM > outFilterMismatchNmaxTotal) || (double(trBest->nMM)/double(trBest->rLength)>P.outFilterMismatchNoverLmax) ) {//too many mismatches
        statsRA.unmappedMismatch++;
        unmapType=2;
    } else if (nTr > P.outFilterMultimapNmax){//too multi
        statsRA.unmappedMulti++;
        unmapType=3;
    };
    
    return;        
};