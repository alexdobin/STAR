/*
 * Created by Fahimeh Mirhaj on 6/10/19.
*/
using namespace std;

#include "SpliceGraph.h"
#include "GTF.h"
SpliceGraph::SpliceGraph (SuperTranscriptome &superTrome, Parameters &P, ReadAlign *RA) : superTrome(superTrome), P(P), RA(RA)
{
    //find candidate superTr
    superTrSeedCount = new typeSuperTrSeedCount[2*superTrome.N];//TODO: for stranded data, do not need 2nd strand
    
    //Smith-Waterman
    scoringMatrix = new typeAlignScore*[superTrome.sjDonorNmax+2];
    scoreTwoColumns[0] = new typeAlignScore[maxSeqLength];
    scoreTwoColumns[1] = new typeAlignScore[maxSeqLength];
    for(uint32 ii = 0; ii < superTrome.sjDonorNmax+2; ii++) {
        scoringMatrix[ii] = new typeAlignScore[maxSeqLength];//TODO make it a user parameter
    };
    sjDindex = new uint32[superTrome.sjDonorNmax];
    
    //rowCol.reserve(100000);
    //rowSJ.reserve(100000);
    //blockCoord.reserve(100000);
    //blockSJ.reserve(10000);
};

SpliceGraph::~SpliceGraph() {
    for(uint i = 0; i < maxSeqLength; i++) {
        delete[] scoringMatrix[i];
    };
    delete[] scoringMatrix;
};
