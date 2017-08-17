#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ReadAlign.h"
#include "BAMfunctions.h"
#include "blocksOverlap.h"
#include "ChimericSegment.h"
#include "ChimericAlign.h"

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
bool ReadAlign::chimericDetectionMult() {

    bool chimRecord=false;
    
    //////////////////// chimeras
    //stich windows => chimeras
    //stich only the best window with one of the lower score ones for now - do not stich 2 lower score windows
    //stitch only one window on each end of the read
    
    if (nTr>P.pCh.mainSegmentMultNmax && nTr!=2)
    {//multimapping main segment, nTr==2 is a special case to be checked later
        return chimRecord;
    };
    
    vector <ChimericAlign> chimAligns;
    int chimScoreBest=0,chimScoreNext=0;
                
    for (uint iW1=0; iW1<nW; iW1++)
    {//cycle windows 
        for (uint iA1=0; iA1<nWinTr[iW1]; iA1++)
        {//cycle aligns in the window
            
            ChimericSegment seg1(P,*trAll[iW1][iA1]);
    
            bool seg1yes = true;
            seg1yes = seg1yes && seg1.align.rLength >= P.pCh.segmentMin; //mapped length >= chim segmentMin
            seg1yes = seg1yes && (seg1.align.exons[seg1.align.nExons-1][EX_R] + seg1.align.exons[seg1.align.nExons-1][EX_L] + P.pCh.segmentMin <= Lread \
                      || seg1.align.exons[0][EX_R] >= P.pCh.segmentMin); //uncovered by seg1 read length is <= segmentMin
            seg1yes = seg1yes && seg1.align.intronMotifs[0]==0; //no non-canonical juncions. TODO: allow non-canonical anotated
            seg1yes = seg1yes && (seg1.align.intronMotifs[1]==0 || seg1.align.intronMotifs[2]==0); //consistent intron motifs. TODO: switch to strands
            
            if (!seg1yes)
                continue; //seg1 is bad
            
            for (uint iW2=0; iW2<nW; iW2++) 
            {//check all windows for chimeras
                for (uint iA2=0; iA2<nWinTr[iW2]; iA2++)
                {//cycle over aligns in the window
                    if (trBest!=trAll[iW2][0] && iA2>0) break; //check only best transcripts in each window 2(i.e. iA2=0), for all windows except that of the trBest
                    if (trBest==trAll[iW2][0] && iA2==0) continue;//for trBest window, check all iA2>0

                    ChimericSegment seg2(P,*trAll[iW2][iA2]);

                    if (seg2.align.intronMotifs[0]>0) 
                        continue; //do not stitch to a window with non-canonical junctions
                    if (seg1.str!=0 && seg2.str!=0 && seg2.str!=seg1.str) 
                        continue; //chimeric segments have to have consitent strands                

                    int chimScore=chimericAlignScore(seg1,seg2);

                    if ( &seg1.align!=trBest && &seg2.align!=trBest )
                        continue; //debug

                    if  (chimScore>0)
                    {//candidate chimera                       

                        ChimericSegment *s1=&seg1,*s2=&seg2;
                        if (seg1.align.roStart > seg2.align.roStart) swap (s1,s2);
                        uint e1 = s1->align.Str==1 ? 0 : s1->align.nExons-1;
                        uint e2 = s2->align.Str==0 ? 0 : s2->align.nExons-1;   
                        if ( s1->align.exons[e1][EX_iFrag] > s2->align.exons[e2][EX_iFrag] )  
                            continue; //strange configuration
                        
                        //         if ( trChim[0].exons[e0][EX_iFrag] > trChim[1].exons[e1][EX_iFrag] ) {//strange configuration, rare, similar to the next one
                        //             chimN=0;//reject such chimeras
                        //             //good test example:
                        //             //CTTAGCTAGCAGCGTCTTCCCAGTGCCTGGAGGGCCAGTGAGAATGGCACCCTCTGGGATTTTTGCTCCTAGGTCT
                        //             //TTGAGGTGAAGTTCAAAGATGTGGCTGGCTGTGAGGAGGCCGAGCTAGAGATCATGGAATTTGTGAATTTCTTGAA
                        //         } else 
                        
                        if (s1->align.exons[e1][EX_L] < P.pCh.junctionOverhangMin &&  s2->align.exons[e2][EX_L] < P.pCh.junctionOverhangMin)                       
                            continue; //junction overhangs too short

                        uint overlap1=0;
                        if (iA2>0 && chimScoreBest>0)
                        {//overlap between chimeric candidate segment and the best chimeric segment so far. Maybe non-zero only if both are in the same window.
                            overlap1=blocksOverlap(chimAligns.back().seg2.align,seg2.align);
                        };
                        
                        if (chimScore > chimScoreBest && chimScore >= P.pCh.scoreMin && chimScore+P.pCh.scoreDropMax >= (int) (readLength[0]+readLength[1]) ) 
                        {
                            chimAligns.clear();
                            chimAligns.push_back(ChimericAlign(seg1, seg2, chimScore));
                            if (overlap1==0)
                            {
                                chimScoreNext=chimScoreBest;
                            };
                            chimScoreBest=chimScore;

                        } else if (chimScore>chimScoreNext && overlap1==0) {//replace the nextscore if it's not the best one and is higher than the previous one
                            chimScoreNext=chimScore;
                        };
                    };
                };//cycle over window2 aligns
            };//cycle over window2
        };//cycle over window1 aligns
    };//cycle over window1                

    if (nTr>P.pCh.mainSegmentMultNmax)
    {//check main segment for multi-aligns
     //this is nTr==2 - a special case: chimeras are allowed only if the 2nd chimeric segment is the next best alignment
        if ( &chimAligns[0].seg2.align!=trMult[0] && &chimAligns[0].seg2.align!=trMult[1] )
        {
            return chimRecord;
        };
    };

    chimN=0;   
    for (uint ic=0; ic<chimAligns.size(); ic++) {//scan all chimeras, find the number within score range
        if (chimAligns[ic].chimScore >= chimScoreBest - (int)P.pCh.multimapScoreRange)
            ++chimN;
    };

    if (chimN > 2*P.pCh.multimapNmax) //too many loci (consider 2* more candidates for stitching below)
        return chimRecord;
    
    //TODO deal with the case chimScoreBest chimera is eliminated
    chimN=0;
    for (uint ic=0; ic<chimAligns.size(); ic++) {//re-scan all chimeras: stitch and re-check the score
        if (chimAligns[ic].chimScore >= chimScoreBest-(int)P.pCh.multimapScoreRange) {
            chimAligns[ic].chimericStitching(mapGen.G, Read1[0]);
            if (chimAligns[ic].chimScore >= chimScoreBest-(int)P.pCh.multimapScoreRange)
                ++chimN;
        };
    };
    
    if (chimN > P.pCh.multimapNmax) //too many loci
        return chimRecord;
    
    if (chimScoreNext >= chimScoreBest - P.pCh.scoreSeparation)
        return chimRecord;
    
    for (uint ic=0; ic<chimAligns.size(); ic++) {//output chimeras within score range
        if (chimAligns[ic].chimScore >= chimScoreBest-(int)P.pCh.multimapScoreRange)
            chimAligns[ic].chimericJunctionOutput(chunkOutChimJunction);
    };
//     if (chimScoreNext + P.pCh.scoreSeparation < chimScoreBest) {//report only if chimera is unique
//             chimAligns[0].chimericStitching(mapGen.G, Read1[0]);
//             if (chimAligns[0].chimScore>0)
//                 chimAligns[0].chimericJunctionOutput(chunkOutChimJunction);
//     };
    if (chimN>0)
        chimRecord=true;
    return chimRecord;
};//END