/*
 * Created by Fahimeh Mirhaj on 6/18/19.
*/

#include "SpliceGraph.h"
#define macro_CompareScore(score1,scoreMax,dirInd,dirIndMax)	if(score1>scoreMax){dirIndMax=dirInd;scoreMax=score1;}

SpliceGraph::typeAlignScore SpliceGraph::swScoreSpliced(const char *readSeq, const uint32 readLen, const SuperTranscript &superTr, 
                                                        array<SpliceGraph::typeSeqLen, 2> &alignStarts, array<SpliceGraph::typeSeqLen, 2> &alignEnds)
{//Smith-Waterman alignment with splices
    
    uint32 superTrLen = superTr.length;

    typeAlignScore scoreMaxGlobal = 0;
    
    int32 iDonor=0; //current donor in the donor list
    int32 iAcceptor=0;//current acceptor in the sjC list (sorted by acceptors
    bool sjYes=superTr.sjC.size() > 0; //spliced superTr
    
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
            while (iAcceptor < (int32) superTr.sjC.size() && col==superTr.sjC[iAcceptor][1]) {//col matches acceptor column: find all donors
                sjDindex[sjN]=superTr.sjC[iAcceptor][2];//index of donor
                ++sjN;
                ++iAcceptor;
            };
            
            if (iDonor < (int32) superTr.sjDonor.size() && col==superTr.sjDonor[iDonor]) {//donor column, has to be recorded
                scoreColumn=scoringMatrix[iDonor];//point to the stored columns
                ++iDonor; //advance for the next donor column
            };
        };
        
        scoreColumn[0] = 0; //initialize top row
        ++matIndex; //since we are skippin the 1st row in the loop
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
			score1 = scoreColumnPrev[row-1] + (readSeq[row-1] == superTr.seqP[col] ? matchScore : misMatchPenalty);
			macro_CompareScore(score1,scoreMax,3,dirIndMax);
			
            for(uint32 ii = 0; ii < sjN; ii++) {
                score1 = scoringMatrix[sjDindex[ii]][row] + gapPenalty;
				macro_CompareScore(score1,scoreMax,4+ii*2,dirIndMax);
				
                score1 = scoringMatrix[sjDindex[ii]][row - 1] + (readSeq[row - 1] == superTr.seqP[col] ? matchScore : misMatchPenalty);
				macro_CompareScore(score1,scoreMax,5+ii*2,dirIndMax);				
            };
			
			directionMatrix[matIndex]=dirIndMax;
            ++matIndex;//matIndex stride is (readLen+1)

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
    
    //blockSJ.clear();//index of junction blocks, recorded acceptors, then converted to donors
    blockCoord.clear();//bR,bG,bL,type recorded ends, then converted to starts
    array<int32,4> block1;
    int32 bType=16; //0: continuing block; 1: D, 2: I, 3: DI, 4: SJ
    
    --iAcceptor; //= last junction
    while(col >= 0 && row > 0) {
        uint32 dir1= (uint32) directionMatrix[row+col*(readLen+1)];
        if (dir1==0) //reached scoringMatrix==0
            break;
        switch (dir1) 
        {
            case 1:
                --row;
                bType=bType | 2;//insertion
                break;
            case 2:
                --col;
                bType=bType | 1;//deletion
                break;
            case 3:
                if (bType!=0) {//start new block
                    block1={row-1, col, 0, bType};
                    bType=0;//marks continuing block
                };
                --row;
                --col;                  
                block1[2]++;//increase block length
                break;
            default: //junction jump
                bType=4;//marks junction for now TODO: motif, annotation
                //blockSJ.push_back(blockCoord.size());
                while (iAcceptor+1!=0 && col <= (int32)superTr.sjC[iAcceptor][1]) {//find (acceptor-1) that matches this column
                    --iAcceptor;
                };                
                col=superTr.sjDonor[ superTr.sjC[iAcceptor+1+(dir1-4)/2][2] ];
                if ((dir1-4)%2 == 1) {
                    --row;
                    block1[2]++;//increase block length
                };
        };
        if (bType!=0 && block1[2]>0) {//block has ended: record the old one
            blockCoord.push_back(block1);
            block1[2]=0;
        };
    };   
    
    alignStarts[0]=(uint32)max(0,row-1);
    alignStarts[1]=(uint32)max(0,col);
    alignEnds[0]--; //substract 1 from row
    
    if (block1[2]>0) {
        blockCoord.push_back(block1);//last block
    };
    std::reverse(blockCoord.begin(), blockCoord.end());
    for (auto &b : blockCoord) {
        b[0] -= b[2]-1;
        b[1] -= b[2]-1;
    };
//     std::reverse(blockSJ.begin(), blockSJ.end());
//     for (auto &b : blockSJ) {
//         b = blockCoord.size()-1-b-1;
//     };
    return scoreMaxGlobal;
};

