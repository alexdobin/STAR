/*
 * Created by Fahimeh Mirhaj on 6/14/19.
*/

#ifndef H_SpliceGraph
#define H_SpliceGraph
#include "GTF.h"
#include "Transcript.h"

class ReadAlign;

class SpliceGraph {
public:
    typedef int32 typeAlignScore;
    typedef uint32 typeSeqLen;
    
    SuperTranscriptome &superTrome;
    
    //vector<array<typeSeqLen, 2>> readAndSuperTranscript;
    const static typeSeqLen maxSeqLength = 100000;//make user parameter?
    typeAlignScore **scoringMatrix, *scoreTwoColumns[2];
    const static typeSeqLen maxMatrixSize = 1000000000; //1B elements in the matrix for now
    uint8 directionMatrix[maxMatrixSize];
    uint32 *sjDindex;
    
    int8_t gapPenalty = -1;
    int8_t matchScore = 1;
    int8_t misMatchPenalty = -1;
    
    //seed-and-rank
    typedef uint16 typeSuperTrSeedCount;
    typeSuperTrSeedCount *superTrSeedCount;
    
    //output
    struct {
        uint32 nMap, nMM, nI, nD, nSJ;
        array<SpliceGraph::typeSeqLen, 2> aStart, aEnd;
    } alignInfo;
    //vector<array<int32,4>> blockCoord;
    //vector<int32> blockSJ;
    //vector<int32> rowCol;
    //vector<array<int32,2>> rowSJ;
    
    
    
    SpliceGraph(SuperTranscriptome &superTrome, Parameters &P, ReadAlign *RA);
    ~SpliceGraph();

    typeAlignScore swScoreSpliced(const char *readSeq, const uint32 readLen, const SuperTranscript &superTr, vector<array<uint32,2>> &cigar);
    //void swTraceBack(array<typeSeqLen, 2> &alignEnds, array<typeSeqLen, 2> &alignStarts);
    void findSuperTr(const char *readSeq, const char *readSeqRevCompl, const uint32 readLen, const string &readName, Genome &mapGen);
    
private:
    Parameters &P;
    ReadAlign *RA;
};

#endif
