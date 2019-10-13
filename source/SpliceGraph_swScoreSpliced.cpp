/*
 * Created by Fahimeh Mirhaj on 6/18/19.
*/

#include "SpliceGraph.h"

#define macro_CompareScore(score1,scoreMax,dirInd,dirIndMax)	if(score1>scoreMax){dirIndMax=dirInd;scoreMax=score1;}

SpliceGraph::typeAlignScore SpliceGraph::swScoreSpliced(const char *readSeq, const uint32 readLen, const uint32 suTrInd, array<SpliceGraph::typeSeqLen, 2> &alignStarts, array<SpliceGraph::typeSeqLen, 2> &alignEnds)
{//Smith-Waterman alignment with splices
    
    uint32 superTrLen = superTr.length[suTrInd];

    typeAlignScore scoreMaxGlobal = 0;
    
    int32 iDonor=0; //current donor in the donor list
    int32 iAcceptor=0;//current acceptor in the sjC list (sorted by acceptors
    bool sjYes=superTr.sjC[suTrInd].size() > 0; //spliced superTr
    
    typeAlignScore *scoreColumn, *scoreColumnPrev;//pointer to column and prev column    
    uint32 iTwoColumn = 0;//selects the column from scoreTwoColumn 2-col matrix
    scoreColumn = scoreTwoColumns[iTwoColumn];
    for(uint32 ii = 0; ii <= readLen; ii++){//fill 0th column
        scoreColumn[ii] = 0;
    };
    
	uint32 matIndex=0;//current matIndex: starting from 1st column (not 0th)
    
    for(uint32 col=0; col<superTrLen; col++) {//main cycle over columns: note that columns are counted from 0!
        scoreColumnPrev = scoreColumn;//prev column pointer
        iTwoColumn = iTwoColumn==0 ? 1 : 0; //switch columns      
        scoreColumn = scoreTwoColumns[iTwoColumn];
        
        //find the splice junction connected donor columns, if any
        uint32 sjN=0;
        if (sjYes) {//col matches acceptor column
            while (col==superTr.sjC[suTrInd][iAcceptor][1]) {//col matches acceptor column: find all donors
                sjDindex[sjN]=superTr.sjC[suTrInd][iAcceptor][2];//index of donor
                ++sjN;
                ++iAcceptor;
            };
            
            if (col==superTr.sjDonor[suTrInd][iDonor]) {//donor column, has to be recorded
                scoreColumn=scoringMatrix[iDonor];//point to the stored columns
                ++iDonor; //advance for the next donor column
            };
        };
        
        scoreColumn[0] = 0; //initialize top row
        for(uint row = 1; row <= readLen; row++) {//rows are counte from 1, 0th row is filled with 0
			
			uint8 dirIndMax = 0 ; //direction index
			typeAlignScore scoreMax = 0; //max score for this cell
			typeAlignScore score1; //calculated score
			
            //down
			score1 = scoreColumn[row-1] + gapPenalty;
			macro_CompareScore(score1,scoreMax,1,dirIndMax);
			
            //right
			score1 = scoreColumnPrev[row] + gapPenalty;
			macro_CompareScore(score1,scoreMax,2,dirIndMax);
			
            //diagonal
			score1 = scoreColumnPrev[row-1] + (readSeq[row-1] == superTr.seqp[suTrInd][col] ? matchScore : misMatchPenalty);
			macro_CompareScore(score1,scoreMax,3,dirIndMax);
			
            for(uint32 ii = 0; ii < sjN; ii++) {
                score1 = scoringMatrix[sjDindex[ii]][row] + gapPenalty;
				macro_CompareScore(score1,scoreMax,4+ii*2,dirIndMax);
				
                score1 = scoringMatrix[sjDindex[ii]][row - 1] + (readSeq[row - 1] == superTr.seqp[suTrInd][col] ? matchScore : misMatchPenalty);
				macro_CompareScore(score1,scoreMax,5+ii*2,dirIndMax);				
            };
			
			directionMatrix[matIndex]=dirIndMax;	
            matIndex++;//matIndex stride is readLen (not +1)

            scoreColumn[row] = scoreMax;
            if(scoreMaxGlobal < scoreMax) {
                scoreMaxGlobal = scoreMax;
                alignEnds[0] = row;                
                alignEnds[1] = col;
            };
        }; // row for loop
    }; // col for loop
    
    ///////////traceback
    int32 row = alignEnds[0];    
    int32 col = alignEnds[1];
    matIndex=row+col*readLen;
    --iAcceptor; //= last junction
    while(col >= 0 && row > 0) {
        uint32 dir1= (uint32) directionMatrix[row+col*readLen];
        if (dir1==0) //reached scoringMatrix==0
            break;
        switch (dir1) 
        {
            case 1:
                --row;
                break;
            case 2:
                --col;
                break;
            case 3:
                --row;
                --col;
                break;
            default: //junctions
                while (iAcceptor+1!=0 && col <= (int32)superTr.sjC[suTrInd][iAcceptor][1]) {//find (acceptor-1) that matches this column
                    --iAcceptor;
                };                
                col=superTr.sjDonor[suTrInd][ superTr.sjC[suTrInd][iAcceptor+1+(dir1-4)/2][2] ];
                if ((dir1-4)%2 == 1) {
                    --row;
                };
        };
    };
    alignStarts[0]=(uint32)max(0,row-1);
    alignStarts[1]=(uint32)max(0,col);
    alignEnds[0]--; //substract 1 from row
    return scoreMaxGlobal;
};

