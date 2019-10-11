/*
 * Created by Fahimeh Mirhaj on 6/18/19.
*/

#include "SpliceGraph.h"

#define macro_CompareScore(score1,scoreMax,dirInd,dirIndMax)	if(score1>scoreMax){dirIndMax=dirInd;scoreMax=score1;}


SpliceGraph::typeAlignScore SpliceGraph::swScoreSpliced(const char *readSeq, const uint32 readLength, const uint32 suTrInd, array<SpliceGraph::typeSeqLen, 2> &alignEnds)
{//Smith-Waterman alignment with splices
    
    uint32 superTrLen = superTr.length[suTrInd];

    typeAlignScore scoreMaxGlobal = 0;
    
    uint32 iDonor=0; //current donor in the donor list
    auto pColPrev=scoringMatrix[0];//pointer to prev column
    auto pCol=scoringMatrix[0]+1;//pointer to current column
    
    for(uint i = 0; i <= readLength; i++)//fill 0th column
        pColPrev[i] = 0;
    
	uint32 matIndex=readLength;//current matIndex: starting from 1st column (not 0th)
	
    for(uint64 col=1; col<=superTrLen; col++) {//main cycle over columns
		pCol[0] = 0;
		
        vector<uint32> sjColumn;
        for(uint64 sj1 = 0; sj1 < superTr.sjC[suTrInd].size(); sj1++) {
            if( col == superTr.sjC[suTrInd][sj1][1] ) {
                sjColumn.push_back(superTr.sjC[suTrInd][sj1][0]);
            };
        };
        
        for(uint row = 1; row <= readLength; row++) {
			//matIndex++;
			
			uint8 dirIndMax = 0 ; //direction index
			typeAlignScore scoreMax = 0; //max score for this cell
			typeAlignScore score1; //calculated score
			
			score1 = pCol[row-1] + gapPenalty;
			macro_CompareScore(score1,scoreMax,1,dirIndMax);
			
			score1 = pColPrev[row] + gapPenalty;
			macro_CompareScore(score1,scoreMax,2,dirIndMax);
			
			score1 = pColPrev[row-1] + (readSeq[row-1] == superTr.seqp[suTrInd][col-1] ? matchScore : misMatchPenalty);
			macro_CompareScore(score1,scoreMax,3,dirIndMax);
			
            for(uint32 ii = 0; ii < sjColumn.size(); ii++) {
                auto sjCol = sjColumn[ii];
				
                score1 = scoringMatrix[sjCol + 1][row] + gapPenalty;
				macro_CompareScore(score1,scoreMax,4+ii*2,dirIndMax);
				
                score1 = scoringMatrix[sjCol + 1][row - 1] + (readSeq[row - 1] == superTr.seqp[suTrInd][sjCol+1] ? matchScore : misMatchPenalty);
				macro_CompareScore(score1,scoreMax,5+ii*2,dirIndMax);				
            };
			
// 			directionMatrix[matIndex]=dirIndMax;			
            pCol[row] = scoreMax;
            if(scoreMaxGlobal < scoreMax) {
                scoreMaxGlobal = scoreMax;
                alignEnds[0] = row;                
                alignEnds[1] = col;
            };
        }; // row for loop
        
        if (col-1==superTr.sjDonor[suTrInd][iDonor]) {//need to save previous column
            pColPrev=pCol;//current -> previous        
            ++pCol;//advance by one
            ++iDonor;
        } else {//no need to save previous column, swap it with current and refill
            swap(pCol,pColPrev);
        };
        
    }; // col for loop
    return scoreMaxGlobal;
};

