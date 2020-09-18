#include "ReadAlign.h"
#include "BAMfunctions.h"


void ReadAlign::chimericDetectionPEmerged(ReadAlign &seRA) {

    chimRecord=false;
    if (P.pCh.segmentMin==0) {//no chimeric detection requested
        return;
    };

    if (P.pCh.multimapNmax==0) {

        // runs old chimeric detection routines.

        seRA.multMapSelect(); //this needs to be done for ChimericDetectionOld, may not need it for the new algorithm
        seRA.mappedFilter();

        chimRecord=seRA.chimericDetectionOld();
        if (!chimRecord) {
            return;
        };

        peOverlapChimericSEtoPE(&seRA.trChim[0], &seRA.trChim[1], &trChim[0], &trChim[1]);
        chimericDetectionOldOutput();

    } else if (trBest->maxScore <= (int) (readLength[0]+readLength[1]) - (int) P.pCh.nonchimScoreDropMin) {//require big enough drop in the best score

        // new chimeric detection routine

        chimRecord=seRA.chimDet->chimericDetectionMult(seRA.nW, seRA.readLength, seRA.trBest->maxScore, this);
    };

    if ( chimRecord ) {
        statsRA.chimericAll++;
    };

    return;
};
