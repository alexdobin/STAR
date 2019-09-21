//
//  SWComputeScore.cpp
//  
//
//  Created by Fahimeh Mirhaj on 6/18/19.
//

#include "SWComputeScore.hpp"
#include "SmithWatermanAlignment.h"
#include "GTF.h"
#include <climits> //INT_MIN
#include <iostream>
#include <string>
#include <algorithm> //find()
#include <stdlib.h> //abs()
#include <utility>
SmithWatermanAlignment::scoreType SmithWatermanAlignment::computeScoreWithSpliceJunctions(uint32 transId, vector<uint8> read, array<SmithWatermanAlignment::seqLenType, 2> &indexToAbsMaxScore) {
    vector<uint8> &superTrancript = sequenceOfSuperTranscripts[transId];
    uint32 readLength = read.size();
    uint32 STLength = superTrancript.size();
    vector<array<uint32, 2>> &eachSTSP = superTranscriptsSpliceJunctionsCollapsed[transId];
    
    for(uint i = 0; i <= readLength; i++) {
        scoringMatrix[0][i] = 0;
    }
    for(uint j = 0; j <= STLength; j++) {
        scoringMatrix[j][0] = 0;
    }
    // Compute scores including splice junctions and finding absolute maximum

    scoreType maxScore = 0;
    for(uint col = 1; col <= STLength; col++) {
        vector<uint32> splicesForEachColumn;
        for(uint splices = 0; splices < eachSTSP.size(); splices++) {
            if(col == eachSTSP[splices][sjE]) {
                splicesForEachColumn.push_back(eachSTSP[splices][sjS]);
            }
        }
        for(uint row = 1; row <= readLength; row++) {
            scoreType spLocalMax = 0;
            scoreType prevSPLocalMax = -1;
            pair<seqLenType, seqLenType> spCellIndexes{-1, -1};
            pair<seqLenType, seqLenType> cellIndexes{-1, -1};
            
            vector<scoreType> allPossibleScores(3 + (2 * splicesForEachColumn.size()));
            allPossibleScores[0] = scoringMatrix[col][row-1] + gapPenalty;
            allPossibleScores[1] = scoringMatrix[col-1][row] + gapPenalty;
            allPossibleScores[2] = scoringMatrix[col-1][row-1] + (read[row-1] == superTrancript[col-1] ? matchScore : misMatchPenalty);
            
            scoreType multiMaxOnSPs = 0;
            for(uint32 i = 0; i < splicesForEachColumn.size(); i++) {
                auto indx = splicesForEachColumn[i];
                allPossibleScores[3 + i * 2] = scoringMatrix[indx + 1][row] + gapPenalty;
                allPossibleScores[3 + i * 2 + 1] = scoringMatrix[indx + 1][row - 1] + (read[row - 1] == superTrancript[indx+1] ? matchScore : misMatchPenalty);
            }
        
            //================== new approach =========
            // Finding maximum & absolute maximun indexes
            int kMax = 0;
            scoreType localMaxScore = 0;
            for(uint k = 0; k < allPossibleScores.size(); k++) {
                if(localMaxScore < allPossibleScores[k]) {
                    kMax = k;
                    localMaxScore = allPossibleScores[k];
                }
            }
            scoringMatrix[col][row] = localMaxScore;
            if(maxScore < scoringMatrix[col][row]) {
                maxScore = scoringMatrix[col][row];
                indexToAbsMaxScore[0] = col;
                indexToAbsMaxScore[1] = row;
            }
                //======= Backward Move ==================
            if(kMax == 0) {
                directionMatrix[col][row].first = col;
                directionMatrix[col][row].second = row-1;
            }
            else if(kMax == 1) {
                directionMatrix[col][row].first = col-1;
                directionMatrix[col][row].second = row;
            }
            else if(kMax == 2) {
                directionMatrix[col][row].first = col-1;
                directionMatrix[col][row].second = row-1;
            }
            else {
                directionMatrix[col][row].first = splicesForEachColumn[(kMax - 3) /2];
                directionMatrix[col][row].second = ((kMax - 3) %2) == 0 ? row: row - 1;
            }
        } // row for loop
    } // col for loop
    return maxScore;
}
SmithWatermanAlignment::scoreType SmithWatermanAlignment::computeScore(uint32 transId, vector<uint8> read, array<SmithWatermanAlignment::seqLenType, 2> &indexToAbsMaxScore) {
    vector<uint8> &superTrancript = sequenceOfSuperTranscripts[transId];
    uint32 readLength = read.size();
    uint32 STLength = superTrancript.size();
    
    for(uint i = 0; i < readLength; i++) {
        scoringMatrix[0][i] = 0;
    }
    for(uint j = 0; j < STLength; j++) {
        scoringMatrix[j][0] = 0;
    }
    // Compute scores and finding absolute maximum
    scoreType maxScore = 0;
    
    for(uint col = 1; col <= STLength; col++) {
        for(uint row = 1; row <= readLength; row++) {
            scoringMatrix[col][row]= max(0, max(max(scoringMatrix[col][row-1] + gapPenalty, scoringMatrix[col-1][row] + gapPenalty), scoringMatrix[col-1][row-1] + (read[row-1] == superTrancript[col-1] ? matchScore : misMatchPenalty)));
            
            if(maxScore < scoringMatrix[col][row]) {
                maxScore = scoringMatrix[col][row];
                indexToAbsMaxScore[0] = col;
                indexToAbsMaxScore[1] = row;
            }
        }
    }
    return maxScore;
}
void SmithWatermanAlignment::traceBack(uint32 transId, array<seqLenType, 2> &maxScoreIndexes, seqLenType &SuperTransStart) {
    eachReadStats[0] = 0;
    eachReadStats[1] = 0;
    eachReadStats[2] = 0;
    eachReadStats[3] = 0;
    
    int col = maxScoreIndexes[0];
    int row = maxScoreIndexes[1];
    int rowT = 0;
    int colT = 0;
    while(col > 0 && row > 0 && scoringMatrix[col][row] != 0) {
        rowT = directionMatrix[col][row].second;
        colT = directionMatrix[col][row].first;

        row = rowT;
        col = colT;
    }
    SuperTransStart = col;
}
void SmithWatermanAlignment::splicJunctionsTransformation() {
    superTranscriptsSpliceJunctions.resize(sequenceOfSuperTranscripts.size());
    uint32 prevTransId;
    map <uint32, set<array<uint32, 2>>> transIdToStartEndMap;
    bool firtInsert = true;
    for(uint i = 0; i < spliceJunctions.size(); i++) {
        uint32 transId = spliceJunctions[i][sjSu];
        array<uint32, 3> sp;
        sp[0] = spliceJunctions[i][sjS];
        sp[1] = spliceJunctions[i][sjE];
        sp[2] = spliceJunctions[i][sjT];
        superTranscriptsSpliceJunctions[transId].push_back(sp); // the transcript id calculation might be wrong!
        array<uint32, 2> startEndArray = {sp[0], sp[1]};
        transIdToStartEndMap[transId].insert(startEndArray);
    }
    
    superTranscriptsSpliceJunctionsCollapsed.resize(sequenceOfSuperTranscripts.size());
    for(auto k = transIdToStartEndMap.begin(); k != transIdToStartEndMap.end(); k++) {
        for(auto& l: k->second) {
            superTranscriptsSpliceJunctionsCollapsed[k->first].push_back(l);
        }
    }
    for(uint j = 0; j < superTranscriptsSpliceJunctionsCollapsed.size(); j++) {
        sort(superTranscriptsSpliceJunctionsCollapsed[j].begin(), superTranscriptsSpliceJunctionsCollapsed[j].end(),
             [](const array<uint32, 2>& e1, const array<uint32, 2>& e2) {
                 return e1[0] < e2[0];
             });
    }
}
