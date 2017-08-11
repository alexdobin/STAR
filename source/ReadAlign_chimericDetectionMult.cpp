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
    bool diffMates=(seg1.roE < seg1.readLength[0] && seg2.roS >= seg1.readLength[0]) || (seg2.roE < seg1.readLength[0] && seg1.roS >= seg1.readLength[0]);

    //segment lengths && (different mates || small gap between segments)
    if (seg1.roE > seg1.P.pCh.segmentMin + seg1.roS + chimOverlap && seg2.roE > seg1.P.pCh.segmentMin + seg2.roS + chimOverlap  \
        && ( diffMates || ( (seg1.roE + seg1.P.pCh.segmentReadGapMax + 1) >= seg2.roS && (seg2.roE + seg1.P.pCh.segmentReadGapMax + 1) >= seg1.roS ) ) ) 
    {
       chimScore=chimScore= seg1.align.maxScore + seg2.align.maxScore - (int)chimOverlap; //subtract overlap to avoid double counting
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
            
            ChimericSegment seg1(P,*trAll[iW1][iA1],Lread,readLength);
    
            bool seg1yes = true;
            seg1yes = seg1yes && seg1.align.rLength >= P.pCh.segmentMin; //mapped length >= chim segmentMin
            seg1yes = seg1yes && (seg1.align.exons[seg1.align.nExons-1][EX_R] + seg1.align.exons[seg1.align.nExons-1][EX_L] + P.pCh.segmentMin <= Lread) \
                      || (seg1.align.exons[0][EX_R] >= P.pCh.segmentMin); //uncovered by seg1 read length is <= segmentMin
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

                    ChimericSegment seg2(P,*trAll[iW2][iA2],Lread,readLength);

                    if (seg2.align.intronMotifs[0]>0) 
                        continue; //do not stitch to a window with non-canonical junctions
                    if (seg1.str!=0 && seg2.str!=0 && seg2.str!=seg1.str) 
                        continue; //chimeric segments have to have consitent strands                

                    int chimScore=chimericAlignScore(seg1,seg2);

                    if ( &seg1.align!=trBest && &seg2.align!=trBest )
                        continue; //debug

                    if  (chimScore>0)
                    {
                        uint overlap1=0;
                        if (iA2>0 && chimScoreBest>0)
                        {//overlap between chimeric candidate segment and the best chimeric segment so far. Maybe non-zero only if both are in the same window.
                            overlap1=blocksOverlap(chimAligns.back().seg2.align,seg2.align);
                        };

                        if (chimScore > chimScoreBest && chimScore >= P.pCh.scoreMin && chimScore+P.pCh.scoreDropMax >= (int) (readLength[0]+readLength[1]) ) 
                        {
                            chimAligns.clear();
                            chimAligns.push_back(ChimericAlign(seg1,seg2));
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

    if (chimAligns.size()==0)
        return chimRecord;

    trChim[0]=chimAligns[0].seg1.align;
    trChim[1]=chimAligns[0].seg2.align;
//     trChim[1].roStart = trChim[1].roStr ==0 ? trChim[1].rStart : Lread - trChim[1].rStart - trChim[1].rLength;
//     trChim[1].cStart  = trChim[1].gStart - P.chrStart[trChim[1].Chr];        

    chimStr = max(chimAligns[0].seg1.str,chimAligns[0].seg2.str); //segment strands are either equal, or one is zero - select the non-zero strand
    
    chimN=0;
    if (chimScoreNext + P.pCh.scoreSeparation < chimScoreBest) {//report only if chimera is unique

        if (trChim[0].roStart > trChim[1].roStart) swap (trChim[0],trChim[1]);

        uint e0 = trChim[0].Str==1 ? 0 : trChim[0].nExons-1;
        uint e1 = trChim[1].Str==0 ? 0 : trChim[1].nExons-1;

        uint chimRepeat0=0,chimRepeat1=0,chimJ0=0,chimJ1=0;
        int chimMotif=0;
        chimN=2;
        if ( trChim[0].exons[e0][EX_iFrag] > trChim[1].exons[e1][EX_iFrag] ) {//strange configuration, rare, similar to the next one
            chimN=0;//reject such chimeras
            //good test example:
            //CTTAGCTAGCAGCGTCTTCCCAGTGCCTGGAGGGCCAGTGAGAATGGCACCCTCTGGGATTTTTGCTCCTAGGTCT
            //TTGAGGTGAAGTTCAAAGATGTGGCTGGCTGTGAGGAGGCCGAGCTAGAGATCATGGAATTTGTGAATTTCTTGAA
        } else if ( trChim[0].exons[e0][EX_iFrag] < trChim[1].exons[e1][EX_iFrag] ) {//mates bracket the chimeric junction
            chimN=2;
            chimRepeat=0;
            chimMotif=-1;
            if (trChim[0].Str==1) {//negative strand
                chimJ0=trChim[0].exons[e0][EX_G]-1;
            } else {
                chimJ0=trChim[0].exons[e0][EX_G]+trChim[0].exons[e0][EX_L];
            };
            if (trChim[1].Str==0) {//positive strand
                chimJ1=trChim[1].exons[e1][EX_G]-1;
            } else {
                chimJ1=trChim[1].exons[e1][EX_G]+trChim[1].exons[e1][EX_L];
            };
        } else {//chimeric junctions is within one of the mates, check and shift chimeric junction if necessary
            if (trChim[0].exons[e0][EX_L]>=P.pCh.junctionOverhangMin && trChim[1].exons[e1][EX_L]>=P.pCh.junctionOverhangMin ) {//large enough overhang required
                uint roStart0 = trChim[0].Str==0 ? trChim[0].exons[e0][EX_R] : Lread - trChim[0].exons[e0][EX_R] - trChim[0].exons[e0][EX_L];
                uint roStart1 = trChim[1].Str==0 ? trChim[1].exons[e1][EX_R] : Lread - trChim[1].exons[e1][EX_R] - trChim[1].exons[e1][EX_L];

                uint jR, jRbest=0;
                int jScore=0,jMotif=0,jScoreBest=-999999,jScoreJ=0;
                for (jR=0; jR<roStart1+trChim[1].exons[e1][EX_L]-roStart0-1; jR++) {//scan through the exons to find a canonical junction, and check for mismatches

                    if (jR==readLength[0]) jR++; //skip the inter-mate base

                    char bR=Read1[0][roStart0+jR];

                    char b0,b1;
                    if (trChim[0].Str==0) {
                        b0=mapGen.G[trChim[0].exons[e0][EX_G]+jR];
                    } else {
                        b0=mapGen.G[trChim[0].exons[e0][EX_G]+trChim[0].exons[e0][EX_L]-1-jR];
                        if (b0<4) b0=3-b0;
                    };

                    if (trChim[1].Str==0) {
                        b1=mapGen.G[trChim[1].exons[e1][EX_G]-roStart1+roStart0+jR];
                    } else {
                        b1=mapGen.G[trChim[1].exons[e1][EX_G]+trChim[1].exons[e1][EX_L]-1+roStart1-roStart0-jR];
                        if (b1<4) b1=3-b1;
                    };

                    if ( ( P.pCh.filter.genomicN && (b0>3 || b1>3) ) || bR>3) {//chimera is not called if there are Ns in the genome or in the read
                        chimN=0;
                        break;
                    };

                    char b01,b02,b11,b12;
                    if (trChim[0].Str==0) {
                        b01=mapGen.G[trChim[0].exons[e0][EX_G]+jR+1];
                        b02=mapGen.G[trChim[0].exons[e0][EX_G]+jR+2];
                    } else {
                        b01=mapGen.G[trChim[0].exons[e0][EX_G]+trChim[0].exons[e0][EX_L]-1-jR-1];
                        if (b01<4) b01=3-b01;
                        b02=mapGen.G[trChim[0].exons[e0][EX_G]+trChim[0].exons[e0][EX_L]-1-jR-2];
                        if (b02<4) b02=3-b02;
                    };
                    if (trChim[1].Str==0) {
                        b11=mapGen.G[trChim[1].exons[e1][EX_G]-roStart1+roStart0+jR-1];
                        b12=mapGen.G[trChim[1].exons[e1][EX_G]-roStart1+roStart0+jR];
                    } else {
                        b11=mapGen.G[trChim[1].exons[e1][EX_G]+trChim[1].exons[e1][EX_L]-1+roStart1-roStart0-jR+1];
                        if (b11<4) b11=3-b11;
                        b12=mapGen.G[trChim[1].exons[e1][EX_G]+trChim[1].exons[e1][EX_L]-1+roStart1-roStart0-jR];
                        if (b12<4) b12=3-b12;
                    };

                    jMotif=0;
                    if (b01==2 && b02==3 && b11==0 && b12==2) {//GTAG
                        if (chimStr!=2) {
                            jMotif=1;
                        };
                    } else if(b01==1 && b02==3 && b11==0 && b12==1) {//CTAC
                        if (chimStr!=1) {
                            jMotif=2;
                        };
                    };

                    if (bR==b0 && bR!=b1) {
                        jScore++;
                    } else if (bR!=b0 && bR==b1) {
                        jScore--;
                    };

                    jScoreJ =jMotif==0 ? jScore +  P.pCh.scoreJunctionNonGTAG : jScore ;

                    if ( jScoreJ > jScoreBest || (jScoreJ == jScoreBest && jMotif>0) ) {
                        chimMotif=jMotif;
                        jRbest=jR;
                        jScoreBest=jScoreJ;
                    };
                };//jR cycle
                if (chimN>0) {//else the chimera was rejected because of mismatches

                    //shift junction in trChim
                    if (trChim[0].Str==1) {
                        trChim[0].exons[e0][EX_R] +=trChim[0].exons[e0][EX_L]-jRbest-1;
                        trChim[0].exons[e0][EX_G] +=trChim[0].exons[e0][EX_L]-jRbest-1;
                        trChim[0].exons[e0][EX_L]=jRbest+1;
                        chimJ0=trChim[0].exons[e0][EX_G]-1;
                    } else {
                        trChim[0].exons[e0][EX_L]=jRbest+1;
                        chimJ0=trChim[0].exons[e0][EX_G]+trChim[0].exons[e0][EX_L];
                    };

                    if (trChim[1].Str==0) {
                        trChim[1].exons[e1][EX_R] +=roStart0+jRbest+1-roStart1;
                        trChim[1].exons[e1][EX_G] +=roStart0+jRbest+1-roStart1;
                        trChim[1].exons[e1][EX_L]=roStart1+trChim[1].exons[e1][EX_L]-roStart0-jRbest-1;
                        chimJ1=trChim[1].exons[e1][EX_G]-1;
                    } else {
                        trChim[1].exons[e1][EX_L]=roStart1+trChim[1].exons[e1][EX_L]-roStart0-jRbest-1;
                        chimJ1=trChim[1].exons[e1][EX_G]+trChim[1].exons[e1][EX_L];
                    };
                    //find repeats
                    char b0,b1;
                    uint jR;
                    for (jR=0;jR<100;jR++) {//forward check
                        if (trChim[0].Str==0) {
                            b0=mapGen.G[chimJ0+jR];
                        } else {
                            b0=mapGen.G[chimJ0-jR];
                            if (b0<4) b0=3-b0;
                        };

                        if (trChim[1].Str==0) {
                            b1=mapGen.G[chimJ1+1+jR];
                        } else {
                            b1=mapGen.G[chimJ1-1-jR];
                            if (b1<4) b1=3-b1;
                        };
                        if (b0!=b1) break;
                    };
                    chimRepeat1=jR;
                    for (jR=0;jR<100;jR++) {//reverse check
                        if (trChim[0].Str==0) {
                            b0=mapGen.G[chimJ0-1-jR];
                        } else {
                            b0=mapGen.G[chimJ0+1+jR];
                            if (b0<4) b0=3-b0;
                        };

                        if (trChim[1].Str==0) {
                            b1=mapGen.G[chimJ1-jR];
                        } else {
                            b1=mapGen.G[chimJ1+jR];
                            if (b1<4) b1=3-b1;
                        };
                        if (b0!=b1) break;
                    };
                    chimRepeat0=jR;
                };//chimN>0
            };//large enough overhang
        };//chimeric junction is within a mate

        //chimeric alignments output
        if ( chimN==2 \
                && ( trChim[0].Str!=trChim[1].Str ||  trChim[0].Chr!=trChim[1].Chr \
                || (trChim[0].Str==0 ? chimJ1-chimJ0+1LLU : chimJ0-chimJ1+1LLU) > (chimMotif>=0 ? P.alignIntronMax :  P.alignMatesGapMax) ) )
        {//
            if (chimMotif>=0 && \
                 (trChim[0].exons[e0][EX_L]<P.pCh.junctionOverhangMin+chimRepeat0 || trChim[1].exons[e1][EX_L]<P.pCh.junctionOverhangMin+chimRepeat1) )
            {//filter out linear junctions that are very close to chimeric junction
                return false;
            };

            //junction + SAMp
            chunkOutChimJunction << P.chrName[trChim[0].Chr] <<"\t"<< chimJ0 - P.chrStart[trChim[0].Chr]+1 <<"\t"<< (trChim[0].Str==0 ? "+":"-") \
                    <<"\t"<< P.chrName[trChim[1].Chr] <<"\t"<< chimJ1 - P.chrStart[trChim[1].Chr]+1 <<"\t"<< (trChim[1].Str==0 ? "+":"-") \
                    <<"\t"<< chimMotif <<"\t"<< chimRepeat0  <<"\t"<< chimRepeat1 <<"\t"<< readName+1 \
                    <<"\t"<< trChim[0].exons[0][EX_G] - P.chrStart[trChim[0].Chr]+1 <<"\t"<< outputTranscriptCIGARp(trChim[0]) \
                    <<"\t"<< trChim[1].exons[0][EX_G] - P.chrStart[trChim[1].Chr]+1 <<"\t"<<  outputTranscriptCIGARp(trChim[1]) <<"\n"; //<<"\t"<< trChim[0].exons[0][EX_iFrag]+1 --- no need for that, since trChim[0] is always on the first mate
        };
    };//chimeric score
    return chimRecord;
};//END
