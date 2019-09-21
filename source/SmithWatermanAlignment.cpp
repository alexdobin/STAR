//
//  SmithWatermanAlignment.cpp
//  
//
//  Created by Fahimeh Mirhaj on 6/10/19.
//

#include "SmithWatermanAlignment.h"
#include "SWReadSequenceGeneration.hpp"
#include "SWComputeScore.hpp"
#include "GTF.h"
#include <iostream>
using namespace std;
SmithWatermanAlignment::SmithWatermanAlignment (vector<vector<uint8>> &sequenceOfSTs, vector<array<uint32,sjL>> &spJunctions, GTF &gtfIn, vector<pair<uint64, uint64>> &normalTransIntervalInST, vector<vector<uint8>> &seqOfNormalT, vector<uint64> &normalTSuperTInd): sequenceOfSuperTranscripts(sequenceOfSTs), spliceJunctions(spJunctions), gtf(gtfIn), normalTranscriptIntervalsInST(normalTransIntervalInST), sequenceOfNormalTranscripts(seqOfNormalT), normalTranscriptSuperTindex(normalTSuperTInd) {
    scoringMatrix = new scoreType*[maxSeqLength];
    directionMatrix = new pair<seqLenType,seqLenType>*[maxSeqLength];
    superTIntervals = new pair<uint64, uint64>[gtfIn.superTranscriptIntervals.size()];
    for(uint i = 0; i < maxSeqLength; ++i) {
        directionMatrix[i] = new pair<seqLenType,seqLenType>[10000];
    }
    for(uint i = 0; i < maxSeqLength; i++) {
        scoringMatrix[i] = new scoreType[10000];
    }
    splicJunctionsTransformation();
};
//Destructor
SmithWatermanAlignment::~SmithWatermanAlignment() {
    for(uint i = 0; i < maxSeqLength; i++) {
        delete[] scoringMatrix[i];
    }
    delete[] scoringMatrix;
    for(uint i = 0; i < maxSeqLength; i++) {
        delete[] directionMatrix[i];
    }
    delete[] directionMatrix;
}
