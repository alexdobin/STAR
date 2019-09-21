//
//  Test.cpp
//  
//
//  Created by Fahimeh Mirhaj on 6/18/19.
//

#include "Test.hpp"
#include "SWComputeScore.hpp"
#include "SmithWatermanAlignment.h"
#include "GTF.h"
#include <time.h>
#include <cstdlib> //rand()
#include <chrono>
#include <iostream>
#include <fstream>
void SmithWatermanAlignment::testSmithWatermanScoreComp() {
    uint countSizeOverflowingCases = 0;
    for(uint j = 73203; j < sequenceOfNormalTranscripts.size(); j++) {
        if(sequenceOfNormalTranscripts[j].size() > 10000 || sequenceOfSuperTranscripts[normalTranscriptSuperTindex[j]].size() > 100000) { // Filtering normal transcripts with length longer than 10000 and super transcripts longer than 100000
            countSizeOverflowingCases++;
            continue;
        }
        uint readStart = 0;
        uint readEnd = 0;
        //uint seed = time(NULL);
        uint seed = 0;
        srand(seed);
        int trueScore;
        array<seqLenType, 2> maxScoreIndexes = {0, 0};
        seqLenType superTStart = 0;
        vector<uint8> read = generateReadSeq(sequenceOfNormalTranscripts[j], trueScore, readStart, readEnd); // Generating simulated reads
        uint normalTSeqReadEnd = normalTranscriptIntervalsInST[j].second - (sequenceOfNormalTranscripts[j].size() - readEnd); //Computing read end
        
        auto computeScoreWithSPStart = chrono::steady_clock::now();
        scoreType maxScoreWithSPs = computeScoreWithSpliceJunctions(normalTranscriptSuperTindex[j], read, maxScoreIndexes);
        auto computeScoreWithSPEnd = chrono::steady_clock::now();
        
        auto tracebackStart = chrono::steady_clock::now();
        traceBack(normalTranscriptSuperTindex[j], maxScoreIndexes, superTStart);
        auto tracebackEnd = chrono::steady_clock::now();
        // Computing time taken by score computation and traceback functions
        auto computeScoreWithSPsTimeTaken = chrono::duration_cast<chrono::nanoseconds>(computeScoreWithSPEnd - computeScoreWithSPStart).count();
        auto traceBackComputationTimeTaken = chrono::duration_cast<chrono::nanoseconds> (tracebackEnd - tracebackStart).count();
    
    //if(abs(normalTranscriptIntervalsInST[j].second - maxScoreIndexes[0] + 1) > 4) {
        cout << j << " " << normalTranscriptSuperTindex[j] << " " << sequenceOfNormalTranscripts[j].size() << " " << trueScore << " " << /*maxScore << " " << */maxScoreWithSPs << " " << normalTranscriptIntervalsInST[j].second << " " << maxScoreIndexes[0] - 1 << " " << normalTranscriptIntervalsInST[j].first << " " << superTStart << " " << computeScoreWithSPsTimeTaken << " " << sequenceOfNormalTranscripts[j].size() * sequenceOfSuperTranscripts[normalTranscriptSuperTindex[j]].size() << " " << traceBackComputationTimeTaken << " " << /*eachReadStats[0] << " " << eachReadStats[1] << " " << eachReadStats[2] << " " << eachReadStats[3] << " " <<*/ superTranscriptsSpliceJunctionsCollapsed[normalTranscriptSuperTindex[j]].size() << endl;
    //}
    }
    cout << "number of cases with overflowing size:" << countSizeOverflowingCases << endl;
}
