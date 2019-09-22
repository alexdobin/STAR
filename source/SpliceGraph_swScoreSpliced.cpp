/*
 * Created by Fahimeh Mirhaj on 6/18/19.
*/

#include "SpliceGraph.h"
// #include <climits> //INT_MIN
// #include <iostream>
// #include <string>
// #include <algorithm> //find()
// #include <stdlib.h> //abs()
// #include <utility>

SpliceGraph::typeAlignScore SpliceGraph::swScoreSpliced(char *readSeq, uint32 readLength, uint32 suTrInd, array<SpliceGraph::typeSeqLen, 2> &alignEnds) 
{//Smith-Waterman alignment  
    uint32 superTrLen = superTr.length[suTrInd];
    
    //fill 0th row and column
    for(uint i = 0; i <= readLength; i++) {
        scoringMatrix[0][i] = 0;
    };
    for(uint j = 0; j <= superTrLen; j++) {
        scoringMatrix[j][0] = 0;
    };

    typeAlignScore maxScore = 0;
    for(uint64 col=1; col<=superTrLen; col++) {//main cycle over columns
        vector<uint32> sjColumn;
        for(uint64 sj1 = 0; sj1 < superTr.sjC[suTrInd].size(); sj1++) {
            if( col == superTr.sjC[suTrInd][sj1][sjE] ) {
                sjColumn.push_back(superTr.sjC[suTrInd][sj1][sjS]);
            };
        };
        
        for(uint row = 1; row <= readLength; row++) {
            vector<typeAlignScore> scoreVec(3 + (2 * sjColumn.size()));
            scoreVec[0] = scoringMatrix[col][row-1] + gapPenalty;
            scoreVec[1] = scoringMatrix[col-1][row] + gapPenalty;
            scoreVec[2] = scoringMatrix[col-1][row-1] + (readSeq[row-1] == superTr.seqp[suTrInd][col-1] ? matchScore : misMatchPenalty);
            
            for(uint32 i = 0; i < sjColumn.size(); i++) {
                auto indx = sjColumn[i];
                scoreVec[3 + i * 2] = scoringMatrix[indx + 1][row] + gapPenalty;
                scoreVec[3 + i * 2 + 1] = scoringMatrix[indx + 1][row - 1] + (readSeq[row - 1] == superTr.seqp[suTrInd][indx+1] ? matchScore : misMatchPenalty);
            };
            
            int kMax = 0;
            typeAlignScore localMaxScore = 0;
            for(uint k = 0; k < scoreVec.size(); k++) {
                if(localMaxScore < scoreVec[k]) {
                    kMax = k;
                    localMaxScore = scoreVec[k];
                };
            };
            scoringMatrix[col][row] = localMaxScore;
            if(maxScore < scoringMatrix[col][row]) {
                maxScore = scoringMatrix[col][row];
                alignEnds[0] = row;                
                alignEnds[1] = col;
            };
            
            if(kMax == 0) {
                directionMatrix[col][row].first = col;
                directionMatrix[col][row].second = row-1;
            } else if (kMax == 1) {
                directionMatrix[col][row].first = col-1;
                directionMatrix[col][row].second = row;
            } else if (kMax == 2) {
                directionMatrix[col][row].first = col-1;
                directionMatrix[col][row].second = row-1;
            } else {
                directionMatrix[col][row].first = sjColumn[(kMax - 3) /2];
                directionMatrix[col][row].second = ((kMax - 3) %2) == 0 ? row: row - 1;
            };
        }; // row for loop
    }; // col for loop
    return maxScore;
};

