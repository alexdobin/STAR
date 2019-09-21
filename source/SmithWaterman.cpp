/*
 * Created by Fahimeh Mirhaj on 6/10/19.
*/

#include "SmithWaterman.h"
#include "GTF.h"
using namespace std;
SmithWaterman::SmithWaterman (vector<vector<uint8>> &sequenceOfSTs, vector<array<uint32,sjL>> &spJunctions, GTF &gtfIn, vector<pair<uint64, uint64>> &normalTransIntervalInST, vector<vector<uint8>> &seqOfNormalT, vector<uint64> &normalTSuperTInd): sequenceOfSuperTranscripts(sequenceOfSTs), spliceJunctions(spJunctions), gtf(gtfIn), normalTranscriptIntervalsInST(normalTransIntervalInST), sequenceOfNormalTranscripts(seqOfNormalT), normalTranscriptSuperTindex(normalTSuperTInd) {
    scoringMatrix = new scoreType*[maxSeqLength];
    directionMatrix = new pair<seqLenType,seqLenType>*[maxSeqLength];
    superTIntervals = new pair<uint64, uint64>[gtfIn.superTranscriptIntervals.size()];
    for(uint i = 0; i < maxSeqLength; ++i) {
        directionMatrix[i] = new pair<seqLenType,seqLenType>[10000];
    };
    for(uint i = 0; i < maxSeqLength; i++) {
        scoringMatrix[i] = new scoreType[10000];
    };
    splicJunctionsTransformation();
};

SmithWaterman::~SmithWaterman() {
    for(uint i = 0; i < maxSeqLength; i++) {
        delete[] scoringMatrix[i];
    };
    delete[] scoringMatrix;
    for(uint i = 0; i < maxSeqLength; i++) {
        delete[] directionMatrix[i];
    };
    delete[] directionMatrix;
}

//Unslpiced score, for testing - not used in main code
// SmithWaterman::scoreType SmithWaterman::computeScore(uint32 transId, vector<uint8> read, array<SmithWaterman::seqLenType, 2> &indexToAbsMaxScore) {
//     vector<uint8> &superTrancript = sequenceOfSuperTranscripts[transId];
//     uint32 readLength = read.size();
//     uint32 STLength = superTrancript.size();
//     
//     for(uint i = 0; i < readLength; i++) {
//         scoringMatrix[0][i] = 0;
//     }
//     for(uint j = 0; j < STLength; j++) {
//         scoringMatrix[j][0] = 0;
//     }
//     // Compute scores and finding absolute maximum
//     scoreType maxScore = 0;
//     
//     for(uint col = 1; col <= STLength; col++) {
//         for(uint row = 1; row <= readLength; row++) {
//             scoringMatrix[col][row]= max(0, max(max(scoringMatrix[col][row-1] + gapPenalty, scoringMatrix[col-1][row] + gapPenalty), scoringMatrix[col-1][row-1] + (read[row-1] == superTrancript[col-1] ? matchScore : misMatchPenalty)));
//             
//             if(maxScore < scoringMatrix[col][row]) {
//                 maxScore = scoringMatrix[col][row];
//                 indexToAbsMaxScore[0] = col;
//                 indexToAbsMaxScore[1] = row;
//             }
//         }
//     }
//     return maxScore;
// }
