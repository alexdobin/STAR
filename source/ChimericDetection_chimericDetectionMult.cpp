//#include "blocksOverlap.h"
#include "ChimericDetection.h"
#include "ChimericSegment.h"

int chimericAlignScore (ChimericSegment & seg1, ChimericSegment & seg2)
{
    int chimScore=0;
    uint chimOverlap = seg2.roS>seg1.roS ?  (seg2.roS>seg1.roE ? 0 : seg1.roE-seg2.roS+1) : (seg2.roE<seg1.roS ? 0 : seg2.roE-seg1.roS+1);
    bool diffMates=(seg1.roE < seg1.align.readLength[0] && seg2.roS >= seg1.align.readLength[0]) || (seg2.roE < seg1.align.readLength[0] && seg1.roS >= seg1.align.readLength[0]);

    //segment lengths && (different mates || small gap between segments)
    if (seg1.roE > seg1.P.pCh.segmentMin + seg1.roS + chimOverlap && seg2.roE > seg1.P.pCh.segmentMin + seg2.roS + chimOverlap  \
        && ( diffMates || ( (seg1.roE + seg1.P.pCh.segmentReadGapMax + 1) >= seg2.roS && (seg2.roE + seg1.P.pCh.segmentReadGapMax + 1) >= seg1.roS ) ) ) 
    {
       chimScore = seg1.align.maxScore + seg2.align.maxScore - (int)chimOverlap; //subtract overlap to avoid double counting
    };
    
    return chimScore;
};

/////////////////////////////////////////////////////////////
bool ChimericDetection::chimericDetectionMult(uint nW, uint *readLength) {

    chimRecord=false;
    vecAligns.clear();
    chimAligns.clear();
    chimScoreBest=0;
                
    for (uint iW1=0; iW1<nW; iW1++) {//cycle windows 
        for (uint iA1=0; iA1<nWinTr[iW1]; iA1++) {//cycle aligns in the window
            
            ChimericSegment seg1(P,*trAll[iW1][iA1]);
    
            if (!seg1.segmentCheck())
                continue; //seg1 is bad - skip
           
            for (uint iW2=iW1; iW2<nW; iW2++) {//check all windows for chimeras
                for (uint iA2=(iW1!=iW2 ? 0 : iA1+1); iA2<nWinTr[iW2]; iA2++) {//cycle over aligns in the window
                    //for the same window, start iA2 at iA1+1 to avoid duplicating
                    
                    ChimericSegment seg2(P,*trAll[iW2][iA2]);

                    if (!seg2.segmentCheck())
                        continue; //seg2 is bad - skip

                    if (seg1.str!=0 && seg2.str!=0 && seg2.str!=seg1.str) 
                        continue; //chimeric segments have to have consistent strands. TODO: make this optional                

                    int chimScore=chimericAlignScore(seg1,seg2);

                    if  (chimScore>0)
                    {//candidate chimera
                        ChimericAlign chAl(seg1, seg2, chimScore, vecAligns);
                        
                        if (!chAl.chimericCheck())
                            continue; //check chimeric alignment
       
                        if (chimScore>=chimScoreBest-(int)P.pCh.multimapScoreRange) 
                            chimAligns.push_back(chAl);//add this chimeric alignment
                            
                        if ( chimScore > chimScoreBest && chimScore >= P.pCh.scoreMin && chimScore >= (int)(readLength[0]+readLength[1]) - P.pCh.scoreDropMax ) {
                            chimAligns.back().chimericStitching(outGen.G, Read1[0]);
                            if (chimAligns.back().chimScore > chimScoreBest)
                                chimScoreBest=chimAligns.back().chimScore;
                        };

                    };
                };//cycle over window2 aligns
            };//cycle over window2
        };//cycle over window1 aligns
    };//cycle over window1                

    if (chimScoreBest==0)
        return chimRecord;
    
    chimN=0;   
    for (uint ic=0; ic<chimAligns.size(); ic++) {//scan all chimeras, find the number within score range
        if (chimAligns.at(ic).chimScore >= chimScoreBest - (int)P.pCh.multimapScoreRange)
            ++chimN;
    };
    if (chimN > 2*P.pCh.multimapNmax) //too many loci (considering 2* more candidates for stitching below)
        return chimRecord;
    
    chimN=0;   
    for (uint ic=0; ic<chimAligns.size(); ic++) {//re-scan all chimeras: stitch and re-check the score
        if (chimAligns.at(ic).chimScore >= chimScoreBest-(int)P.pCh.multimapScoreRange) {
            chimAligns.at(ic).chimericStitching(outGen.G, Read1[0]);
            if (chimAligns.at(ic).chimScore >= chimScoreBest - (int)P.pCh.multimapScoreRange)
                ++chimN;
        };
    };
    if (chimN > P.pCh.multimapNmax) //too many loci
        return chimRecord;
    
    for (uint ic=0; ic<chimAligns.size(); ic++) {//output chimeras within score range
        if (chimAligns.at(ic).chimScore >= chimScoreBest-(int)P.pCh.multimapScoreRange)
            chimAligns.at(ic).chimericJunctionOutput(ostreamChimJunction, chimN);
    };

    if (chimN>0)
        chimRecord=true;
    
    return chimRecord;
};//END
