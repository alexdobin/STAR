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
    
    int8_t gapPenalty = -1;
    int8_t matchScore = 1;
    int8_t misMatchPenalty = -1;
    
    //seed-and-rank
    typedef uint16 typeSuperTrSeedCount;
    typeSuperTrSeedCount *superTrSeedCount;
    
    SpliceGraph(SuperTranscript &superTr, Parameters &P);
    ~SpliceGraph();

    typeAlignScore swScoreSpliced(const char *readSeq, const uint32 readLength, const uint32 suTrInd, array<SpliceGraph::typeSeqLen, 2> &alignEnds);
    void swTraceBack(array<typeSeqLen, 2> &alignEnds, array<typeSeqLen, 2> &alignStarts);
    void findSuperTr(const char *readSeq, const char *readSeqRevCompl, const uint32 readLen, const string &readName, Genome &mapGen);
    
private:
    Parameters &P;
};

#endif
