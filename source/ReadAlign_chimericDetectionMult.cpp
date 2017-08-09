#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ReadAlign.h"
#include "BAMfunctions.h"
#include "blocksOverlap.h"
#include "ChimericSegment.h"

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
    
    if (P.pCh.segmentMin>0 && trBest->rLength >= P.pCh.segmentMin \
            && ( trBest->exons[trBest->nExons-1][EX_R] + trBest->exons[trBest->nExons-1][EX_L] + P.pCh.segmentMin <= Lread \
              || trBest->exons[0][EX_R] >= P.pCh.segmentMin ) \
             && trBest->intronMotifs[0]==0 && (trBest->intronMotifs[1]==0 || trBest->intronMotifs[2]==0) ) {
            //there is unmapped space at the start/end, and the main window is not a multimapping window, and non non-canonical junctions, and consistend junction motif
        int chimScoreBest=0,chimScoreNext=0;
        uint chimStrBest=0;
        
        trChim[0]=*trBest;
        Transcript* trChim1=NULL;

        ChimericSegment seg1(P,trChim[0],Lread,readLength);
        chimStr=seg1.str;

        for (uint iW=0; iW<nW; iW++) {//check all other windows for chimeras
            for (uint iWt=0; iWt<nWinTr[iW]; iWt++){//cycl over transcripts in the window
                if (trBest!=trAll[iW][0] && iWt>0) break; //for all windows except that of the best transcript - hceck only iWt=0 (best trnascripts)
                if (trBest==trAll[iW][0] && iWt==0) continue;
                
                ChimericSegment seg2(P,*trAll[iW][iWt],Lread,readLength);
                if (seg2.align.intronMotifs[0]>0) continue; //do not stitch a window to itself, or to a window with non-canonical junctions
                
                if (seg1.str!=0 && seg2.str!=0 && seg2.str!=seg1.str) continue; //chimeric segments have to have consitent strands                

                int chimScore=chimericAlignScore(seg1,seg2);
                
                if  (chimScore>0)
                {
                    uint overlap1=0;
                    if (iWt>0 && chimScoreBest>0)
                    {//overlap between chimeric candidate segment and the best chimeric segment so far. Maybe non-zero only if both are in the same window.
                        overlap1=blocksOverlap(trChim[1],*trAll[iW][iWt]);
                    };

                    if (chimScore > chimScoreBest && chimScore >= P.pCh.scoreMin && chimScore+P.pCh.scoreDropMax >= (int) (readLength[0]+readLength[1]) ) {
                        trChim[1]=*trAll[iW][iWt];
                        trChim1=trAll[iW][iWt];
                        if (overlap1==0)
                        {
                            chimScoreNext=chimScoreBest;
                        };
                        chimScoreBest=chimScore;
                        trChim[1].roStart = trChim[1].roStr ==0 ? trChim[1].rStart : Lread - trChim[1].rStart - trChim[1].rLength;
                        trChim[1].cStart  = trChim[1].gStart - P.chrStart[trChim[1].Chr];
                        chimStrBest=seg1.str;
                    } else if (chimScore>chimScoreNext && overlap1==0) {//replace the nextscore if it's not the best one and is higher than the previous one
                        chimScoreNext=chimScore;
                    };
                };
            };//cycle over window transcripts
        };//cycle over windows

        if (nTr>P.pCh.mainSegmentMultNmax)
        {//check main segment for multi-aligns
         //this is nTr==2 - a special case: chimeras are allowed only if the 2nd chimeric segment is the next best alignment
            if ( trChim1!=trMult[0] && trChim1!=trMult[1] )
            {
                return chimRecord;
            };
        };

        
        if (chimStr==0) chimStr=chimStrBest;

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

                chimRecord=true; //chimeric alignment was recorded

                //re-calculate the score for chimeric transcripts
                trChim[0].alignScore(Read1, mapGen.G, P);
                trChim[1].alignScore(Read1, mapGen.G, P);

                int chimRepresent=-999, chimType=0;
                if (trChim[0].exons[0][EX_iFrag]!=trChim[0].exons[trChim[0].nExons-1][EX_iFrag]) {//tr0 has both mates
                    chimRepresent = 0;
                    chimType = 1;
                    trChim[0].primaryFlag=true;//paired portion is primary
                    trChim[1].primaryFlag=false;
                } else if (trChim[1].exons[0][EX_iFrag]!=trChim[1].exons[trChim[1].nExons-1][EX_iFrag]) {//tr1 has both mates
                    chimRepresent = 1;
                    chimType = 1;
                    trChim[1].primaryFlag=true;//paired portion is primary
                    trChim[0].primaryFlag=false;
                } else if (trChim[0].exons[0][EX_iFrag]!=trChim[1].exons[0][EX_iFrag]) {//tr0 and tr1 are single different mates
                    chimRepresent = -1;
                    chimType = 2;
                    trChim[0].primaryFlag=true;
                    trChim[1].primaryFlag=true;
                } else  {//two chimeric segments are on the same mate - this can only happen for single-end reads
                    chimRepresent = (trChim[0].maxScore > trChim[1].maxScore) ? 0 : 1;
                    chimType = 3;
                    trChim[chimRepresent].primaryFlag=true;
                    trChim[1-chimRepresent].primaryFlag=false;
                };

                if (P.pCh.out.bam) {//BAM output
                    int alignType, bamN=0, bamIsuppl=-1, bamIrepr=-1;
                    uint bamBytesTotal=0;//estimate of the total size of all bam records, for output buffering
                    uint mateChr,mateStart;
                    uint8_t mateStrand;
                    for (uint itr=0;itr<chimN;itr++) {//generate bam for all chimeric pieces
                        if (chimType==2) {//PE, encompassing
                            mateChr=trChim[1-itr].Chr;
                            mateStart=trChim[1-itr].exons[0][EX_G];
                            mateStrand=(uint8_t) (trChim[1-itr].Str!=trChim[1-itr].exons[0][EX_iFrag]);
                            alignType=-10;
                        } else {//spanning chimeric alignment, could be PE or SE
                            mateChr=-1;mateStart=-1;mateStrand=0;//no need fot mate info unless this is the supplementary alignment
                            if (chimRepresent==(int)itr) {
                                alignType=-10; //this is representative part of chimeric alignment, record is as normal; if encompassing chimeric junction, both are recorded as normal
                                bamIrepr=( (itr%2)==(trChim[itr].Str) && chimType!=3) ? bamN+1 : bamN;//this is the mate that is chimerically split
                            } else {//"supplementary" chimeric segment
                                alignType=P.pCh.out.bamHardClip ? ( ( itr%2==trChim[itr].Str ) ? -12 : -11) : -13 ; //right:left chimeric junction
                                bamIsuppl=bamN;
                                if (chimType==1) {//PE alignment, need mate info for the suppl
                                    uint iex=0;
                                    for (;iex<trChim[chimRepresent].nExons-1;iex++) {
                                        if (trChim[chimRepresent].exons[iex][EX_iFrag]!=trChim[itr].exons[0][EX_iFrag]) {
                                            break;
                                        };
                                    };
                                    mateChr=trChim[chimRepresent].Chr;
                                    mateStart=trChim[chimRepresent].exons[iex][EX_G];
                                    mateStrand=(uint8_t) (trChim[chimRepresent].Str!=trChim[chimRepresent].exons[iex][EX_iFrag]);
                                };
                            };

                        };

                        bamN+=alignBAM(trChim[itr], 1, 0, P.chrStart[trChim[itr].Chr],  mateChr, mateStart-P.chrStart[mateChr], mateStrand, \
                                        alignType, NULL, P.outSAMattrOrder, outBAMoneAlign+bamN, outBAMoneAlignNbytes+bamN);
                        bamBytesTotal+=outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1];//outBAMoneAlignNbytes[1] = 0 if SE is recorded
                    };

                    //write all bam lines
                    for (int ii=0; ii<bamN; ii++) {//output all pieces
                        int tagI=-1;
                        if (ii==bamIrepr) {
                          tagI=bamIsuppl;
                        } else if (ii==bamIsuppl) {
                          tagI=bamIrepr;
                        };
                        if (tagI>=0) {
                            bam1_t *b;
                            b=bam_init1();
                            bam_read1_fromArray(outBAMoneAlign[tagI], b);
                            uint8_t* auxp=bam_aux_get(b,"NM");
                            uint32_t auxv=bam_aux2i(auxp);
                            string tagSA1="SAZ"+P.chrName[b->core.tid]+','+to_string((uint)b->core.pos+1) +',' + ( (b->core.flag&0x10)==0 ? '+':'-') + \
                                    ',' + bam_cigarString(b) + ',' + to_string((uint)b->core.qual) + ',' + to_string((uint)auxv) + ';' ;

                             memcpy( (void*) (outBAMoneAlign[ii]+outBAMoneAlignNbytes[ii]), tagSA1.c_str(), tagSA1.size()+1);//copy string including \0 at the end
                             outBAMoneAlignNbytes[ii]+=tagSA1.size()+1;
                             * ( (uint32*) outBAMoneAlign[ii] ) = outBAMoneAlignNbytes[ii]-sizeof(uint32);
                        };

                        if (P.outBAMunsorted) outBAMunsorted->unsortedOneAlign(outBAMoneAlign[ii], outBAMoneAlignNbytes[ii], ii>0 ? 0 : bamBytesTotal);
                        if (P.outBAMcoord)    outBAMcoord->coordOneAlign(outBAMoneAlign[ii], outBAMoneAlignNbytes[ii], (iReadAll<<32) );
                    };
                };


                for (uint iTr=0;iTr<chimN;iTr++) 
                {//write all chimeric pieces to Chimeric.out.sam/junction
                    if (P.readNmates==2) {//PE: need mate info
                        uint iex=0;
                        if ( trChim[1-iTr].exons[0][EX_iFrag] != trChim[1-iTr].exons[trChim[1-iTr].nExons-1][EX_iFrag] )
                        {//the other segment has 2 mates, need to find the opposite mate
                            for (;iex<trChim[1-iTr].nExons;iex++) {
                                if (trChim[1-iTr].exons[iex][EX_iFrag]!=trChim[iTr].exons[0][EX_iFrag]) {
                                    break;
                                };
                            };
                        };
                        
                        uint mateChr=trChim[1-iTr].Chr;
                        uint mateStart=trChim[1-iTr].exons[iex][EX_G];
                        char mateStrand=(char) (trChim[1-iTr].Str!=trChim[1-iTr].exons[iex][EX_iFrag]);
                        
                        outputTranscriptSAM(trChim[iTr], chimN, iTr, mateChr, mateStart, mateStrand, -1, NULL, &chunkOutChimSAM);
                    } else 
                    {
                        outputTranscriptSAM(trChim[iTr], chimN, iTr, -1, -1, -1, -1, NULL, &chunkOutChimSAM);
                    };
                };
                //junction + SAMp
                chunkOutChimJunction << P.chrName[trChim[0].Chr] <<"\t"<< chimJ0 - P.chrStart[trChim[0].Chr]+1 <<"\t"<< (trChim[0].Str==0 ? "+":"-") \
                        <<"\t"<< P.chrName[trChim[1].Chr] <<"\t"<< chimJ1 - P.chrStart[trChim[1].Chr]+1 <<"\t"<< (trChim[1].Str==0 ? "+":"-") \
                        <<"\t"<< chimMotif <<"\t"<< chimRepeat0  <<"\t"<< chimRepeat1 <<"\t"<< readName+1 \
                        <<"\t"<< trChim[0].exons[0][EX_G] - P.chrStart[trChim[0].Chr]+1 <<"\t"<< outputTranscriptCIGARp(trChim[0]) \
                        <<"\t"<< trChim[1].exons[0][EX_G] - P.chrStart[trChim[1].Chr]+1 <<"\t"<<  outputTranscriptCIGARp(trChim[1]) <<"\n"; //<<"\t"<< trChim[0].exons[0][EX_iFrag]+1 --- no need for that, since trChim[0] is always on the first mate
            };
        };//chimeric score
    };//chimeric search
    return chimRecord;
};//END
