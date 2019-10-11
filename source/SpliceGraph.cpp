/*
 * Created by Fahimeh Mirhaj on 6/10/19.
*/
using namespace std;

#include "SpliceGraph.h"
#include "GTF.h"
SpliceGraph::SpliceGraph (SuperTranscript &superTr, Parameters &P) : superTr(superTr), P(P)
{
    //find candidate superTr
    superTrSeedCount = new typeSuperTrSeedCount[2*superTr.N];//TODO: for stranded data, do not need 2nd strand
    
    //Smith-Waterman
    scoringMatrix = new typeAlignScore*[superTr.sjNmax+2];
    directionMatrix = new pair<typeSeqLen,typeSeqLen>*[maxSeqLength];
    //superTIntervals = new pair<uint64, uint64>[gtfIn.superTr.startEndInFullGenome.size()];
    for(uint i = 0; i < maxSeqLength; ++i) {
        directionMatrix[i] = new pair<typeSeqLen,typeSeqLen>[10000];
    };
    for(uint i = 0; i < maxSeqLength; i++) {
        scoringMatrix[i] = new typeAlignScore[10000];
    };
};

SpliceGraph::~SpliceGraph() {
    for(uint i = 0; i < maxSeqLength; i++) {
        delete[] scoringMatrix[i];
    };
    delete[] scoringMatrix;
    for(uint i = 0; i < maxSeqLength; i++) {
        delete[] directionMatrix[i];
    };
    delete[] directionMatrix;
};

//Unslpiced score, for testing - not used in main code
// SpliceGraph::typeAlignScore SpliceGraph::computeScore(uint32 transId, vector<uint8> read, array<SpliceGraph::typeSeqLen, 2> &indexToAbsMaxScore) {
//     vector<uint8> &superTrancript = superTr.seq[transId];
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
//     typeAlignScore maxScore = 0;
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
