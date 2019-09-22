/*
 * Created by Fahimeh Mirhaj on 6/14/19.
*/

#ifndef H_SpliceGraph
#define H_SpliceGraph

#include "GTF.h"
class SpliceGraph {
public:
    typedef int32 typeAlignScore;
    typedef uint32 typeSeqLen;
    
    SuperTranscript &superTr;
    
    //vector<array<typeSeqLen, 2>> readAndSuperTranscript;
    const static typeSeqLen maxSeqLength = 100000;
    typeAlignScore** scoringMatrix;
    pair<typeSeqLen,typeSeqLen>** directionMatrix;
    
    //uint8** cellStats;
    //array<uint, 4> eachReadStats;
    //array<bool, 3> stats;
    
    //pair<uint64, uint64>* superTIntervals;
    
    int8_t gapPenalty = -1;
    int8_t matchScore = 1;
    int8_t misMatchPenalty = -1;
    enum {sjS,sjE,sjT,sjSu,sjL};//same as in SuperTranscript class, need to make more consitent - maybe convert to structure?
    
    SpliceGraph(SuperTranscript &superTr);
    ~SpliceGraph();

    void splicJunctionsTransformation();
    //typeAlignScore computeScore(uint32 transId, vector<uint8> read, array<typeSeqLen, 2> &indexToAbsMaxScore);
    typeAlignScore swScoreSpliced(char *readSeq, uint32 readLength, uint32 suTrInd, array<SpliceGraph::typeSeqLen, 2> &alignEnds);
    void swTraceBack(array<typeSeqLen, 2> &alignEnds, array<typeSeqLen, 2> &alignStarts);

private:  
//    void clearDPTables();
// 
//     vector<vector<uint8>> &superTr.seq;
//     vector<array<uint32,sjL>> &superTr.sj;
//     vector<vector<array<uint32, 3>>> superTranscriptsSpliceJunctions;
//     vector<vector<array<uint32, 2>>> superTranscriptsSpliceJunctionsCollapsed;
//     vector<pair<uint64, uint64>> &superTr.trStartEnd;
//     vector<vector<uint8>> &transcriptSeq;
//     vector<uint64> &superTr.trIndex;
//     GTF &gtf;
};

#endif
