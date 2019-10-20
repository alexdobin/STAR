/*
 * Created by Fahimeh Mirhaj on 6/18/19.
 */

#include "SpliceGraph.h"
#define macro_CompareScore(score1,scoreMax,dirInd,dirIndMax)	if(score1>scoreMax){dirIndMax=dirInd;scoreMax=score1;}

SpliceGraph::typeAlignScore SpliceGraph::swScoreSpliced
                            (const char *readSeq, const uint32 readLen, const SuperTranscript &superTr, vector<array<uint32,2>> &cigar)
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
            ++matIndex;//matIndex stride is (readLen)

            scoreColumn[row] = scoreMax;
            if(scoreMaxGlobal < scoreMax) {
                scoreMaxGlobal = scoreMax;
                alignInfo.aEnd[0] = row;                
                alignInfo.aEnd[1] = col;
            };
        }; // row for loop
    }; // col for loop
    alignInfo.aEnd[0]--;//true row
    
    ///////////traceback
//    int32 row = alignInfo.aEnd[0];//true row
//    int32 col = alignInfo.aEnd[1];
    
//    uint32 nMapped=0, nMM=0, nI=0, nD=0, nSJ=0;
//     blockSJ.clear();//index of junction blocks, recorded acceptors, then converted to donors
//     blockCoord.clear();//bR,bG,bL,type recorded ends, then converted to starts
//     vector<uint32> rowCol(readLen), rowSJ(readLen,0); //records col vs row TODO define outside for speed
//     --iAcceptor; //= last junction
//     
//     rowCol.clear();
//     rowSJ.clear();
//     rowCol.resize(readLen,-1);
//     rowSJ.resize(readLen+1,{-1,-1});//one extra element since we are going to check row+1
//     
//     while(col >= 0 && row >= 0) {
//         uint32 dir1= (uint32) directionMatrix[row+col*readLen];
//         if (dir1==0) //reached scoringMatrix==0
//             break;
//         switch (dir1) 
//         {
//             case 1:
//                 --row;
//                 ++nI;
//                 break;
//             case 2:
//                 --col;
//                 ++nD;
//                 break;
//             case 3:
//                 ++nMapped;
//                 nMM+=(uint32)(readSeq[row]!=superTr.seqP[col]);
//                 rowCol[row]=col;//for each row, point to the column where blocj starts
//                 --row;
//                 --col;                  
//                 break;
//             default: //junction jump
//                 ++nSJ;
//                 while (iAcceptor+1!=0 && col <= (int32)superTr.sjC[iAcceptor][1]) {//find (acceptor-1) that matches this column
//                     --iAcceptor;
//                 };                
//                 //row: same or -1
//                 if ((dir1-4)%2 == 1) {//diagonal-like
//                     ++nMapped;
//                     nMM+=(uint32)(readSeq[row]!=superTr.seqP[col]);
//                     rowCol[row]=col;//for each row, point to the column where block starts                     
//                     --row;
//                 } else {
//                     //++nD; //?
//                 };
//                 rowSJ[row][1]=col;//acceptor
//                 //column: jump to donor
//                 col=superTr.sjDonor[ superTr.sjC[iAcceptor+1+(dir1-4)/2][2] ];
//                 rowSJ[row][0]=col;//donor
//         };
//     };   
//     
//     row=max(0,row);
//     while (rowCol[row]<0)
//         row++; //find first col>=0 - this is the alignment start. 
//    alignInfo.aStart[0]=row;
//    alignInfo.aStart[1]=rowCol[row];
    
    cigar.clear();
    cigar.reserve(readLen);
    int32 row = alignInfo.aEnd[0];//true row
    int32 col = alignInfo.aEnd[1];
    
    alignInfo.nMap=0;
    alignInfo.nMM=0;
    alignInfo.nI=0;
    alignInfo.nD=0;
    alignInfo.nSJ=0;
    iAcceptor=superTr.sjC.size()-1; //= last junction
    
    if (row!=(int32)readLen-1) //soft-clip
        cigar.push_back({BAM_CIGAR_S, readLen-1-row}); 
    
    uint32 cigarOp=0, cigarLen=0, cigarOpPrev=(uint32)-1;
    uint32 sjGap=0;
    while(col >= 0 && row >= 0) {
        uint32 dir1= (uint32) directionMatrix[row+col*readLen];
        
        if (dir1==0) //reached scoringMatrix==0
            break;
        
        switch (dir1) 
        {
            case 1:
                --row;
                alignInfo.nI++;
                cigarOp=BAM_CIGAR_I;
                break;
            case 2:
                --col;
                alignInfo.nD++;
                cigarOp=BAM_CIGAR_D;
                break;
            case 3:
                alignInfo.nMap++;
                alignInfo.nMM += (uint32)(readSeq[row]!=superTr.seqP[col]);
                cigarOp=BAM_CIGAR_M;
                --row;
                --col;                  
                break;
            default: //junction jump
                alignInfo.nSJ++;
                while (iAcceptor+1!=0 && col <= (int32)superTr.sjC[iAcceptor][1]) {//find (acceptor-1) that matches this column
                    --iAcceptor;
                };
                
                //row: same or -1
                if ((dir1-4)%2 == 1) {//diagonal-like
                    alignInfo.nMap++;
                    alignInfo.nMM+=(uint32)(readSeq[row]!=superTr.seqP[col]);
                    --row;
                    cigarOp=BAM_CIGAR_M;
                } else {
                    alignInfo.nD++;
                    cigarOp=BAM_CIGAR_D;
                };
                
                //column: jump to donor
                sjGap=col;
                col=superTr.sjDonor[ superTr.sjC[iAcceptor+1+(dir1-4)/2][2] ];
                sjGap=sjGap-col-1;
        };
        
        if (cigarOp!=cigarOpPrev) {//changed direction - record new cigar op
            if (cigarLen>0)
                cigar.push_back({cigarOpPrev, cigarLen});
            cigarLen=0;
            cigarOpPrev=cigarOp;
        };
        ++cigarLen;
        if (sjGap>0) {
            cigar.push_back({cigarOp, cigarLen});//record previous opeartion
            cigar.push_back({BAM_CIGAR_N, sjGap});//record N
            cigarLen=0; //keep cigarLen=0: after SJ, the blocks start anew
            cigarOpPrev=(uint32)-1;
            sjGap=0;
        };
    };
    
    if (cigarLen>0)
        cigar.push_back({cigarOp, cigarLen});
        
//     row=max(0,row);
//     col=max(0,col);
        
    ++row;
    ++col;
    alignInfo.aStart[0]=row;
    alignInfo.aStart[1]=col;        
    if (row>0)
        cigar.push_back({BAM_CIGAR_S, (uint32)row}); 
        
    std::reverse(cigar.begin(), cigar.end());
    return scoreMaxGlobal;
};

