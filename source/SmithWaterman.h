/*
 * Created by Fahimeh Mirhaj on 6/14/19.
*/

#ifndef H_SmithWaterman
#define H_SmithWaterman

#include "GTF.h"
class SmithWaterman {
public:
    typedef int32 scoreType;
    typedef uint32 seqLenType;
    vector<array<seqLenType, 2>> readAndSuperTranscript;
    const static seqLenType maxSeqLength = 100000;
    scoreType** scoringMatrix;
    pair<seqLenType,seqLenType>** directionMatrix;
    uint8** cellStats;
    array<uint, 4> eachReadStats;
    array<bool, 3> stats;
    pair<uint64, uint64>* superTIntervals;
    int8_t gapPenalty = -1;
    int8_t matchScore = 1;
    int8_t misMatchPenalty = -1;
    enum {sjS,sjE,sjT,sjSu,sjL};
    
    SmithWaterman(vector<vector<uint8>> &sequenceOfSTs, vector<array<uint32,sjL>> &spJunctions, GTF &gtfIn, vector<pair<uint64, uint64>> &normalTransIntervalInST, vector<vector<uint8>> &seqOfNormalT, vector<uint64> &normalTSuperTInd);
    ~SmithWaterman();

    void splicJunctionsTransformation();
    //scoreType computeScore(uint32 transId, vector<uint8> read, array<seqLenType, 2> &indexToAbsMaxScore);
    scoreType scoreSpliced(uint32 transId, vector<uint8> read, array<seqLenType, 2> &indexToAbsMaxScore);
    void traceBack(uint32 transId, array<seqLenType, 2> &indexToAbsMaxScore, seqLenType &SuperTransStart);
    
private:
    
    void clearDPTables();
    
    vector<vector<uint8>> &superTrSeq;
    vector<array<uint32,sjL>> &spliceJunctions;
    vector<vector<array<uint32, 3>>> superTranscriptsSpliceJunctions;
    vector<vector<array<uint32, 2>>> superTranscriptsSpliceJunctionsCollapsed;
    vector<pair<uint64, uint64>> &transcriptSuperTrStartEnd;
    vector<vector<uint8>> &transcriptSeq;
    vector<uint64> &transcriptSuperTrIndex;
    GTF &gtf;
};

#endif
