#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "extendAlign.h"
#include "binarySearch2.h"
// #include "stitchGapIndel.cpp"


intScore stitchAlignToTranscript(uint rAend, uint gAend, uint rBstart, uint gBstart, uint L, uint iFragB, uint sjAB, Parameters& P, char* R, Genome &mapGen, Transcript *trA, const uint outFilterMismatchNmaxTotal) {
    //stitch together A and B, extend in the gap, returns max score

    char *G=mapGen.G;
    int Score=0;
//     int score2;

    if (sjAB!=((uint) -1) && trA->exons[trA->nExons-1][EX_sjA]==sjAB \
            && trA->exons[trA->nExons-1][EX_iFrag]==iFragB && rBstart==rAend+1 && gAend+1<gBstart ) {//simple stitching if junction belongs to a database
        if (mapGen.sjdbMotif[sjAB]==0 && (L<=mapGen.sjdbShiftRight[sjAB] || trA->exons[trA->nExons-1][EX_L]<=mapGen.sjdbShiftLeft[sjAB]) ) {
            return -1000006; //too large repeats around non-canonical junction
        };
        trA->exons[trA->nExons][EX_L] = L; //new exon length
        trA->exons[trA->nExons][EX_R] = rBstart; //new exon r-start
        trA->exons[trA->nExons][EX_G] = gBstart; //new exon g-start
        trA->canonSJ[trA->nExons-1]=mapGen.sjdbMotif[sjAB]; //mark sj-db
        trA->shiftSJ[trA->nExons-1][0]=mapGen.sjdbShiftLeft[sjAB];
        trA->shiftSJ[trA->nExons-1][1]=mapGen.sjdbShiftRight[sjAB];
        trA->sjAnnot[trA->nExons-1]=1;
        trA->sjStr[trA->nExons-1]=mapGen.sjdbStrand[sjAB];;
        trA->nExons++;
        trA->nMatch+=L;
        for (uint ii=rBstart;ii<rBstart+L;ii++) Score+=scoreMatch; //add QS for mapped portions
        Score+=P.pGe.sjdbScore;
    } else {//general stitching
        trA->sjAnnot[trA->nExons-1]=0;
        trA->sjStr[trA->nExons-1]=0;

        if (trA->exons[trA->nExons-1][EX_iFrag]==iFragB) {//stitch aligns on the same fragment
            uint gBend=gBstart+L-1;
            uint rBend=rBstart+L-1;

//             {//debug
//                 if (sjAB!=((uint) -1) && trA->exons[trA->nExons-1][EX_sjA]!=((uint) -1) && rBend<=rAend) {//
//                     Score -= rAend-rBstart+1;
//                     gAend -= rAend-rBstart+1;
//                     rAend = rBstart-1;
//                     trA->exons[trA->nExons-1][EX_L] =rAend-trA->exons[trA->nExons-1][EX_R]+1;
//                 };
//             };

            //check if r-overlapping fully and exit
            if (rBend<=rAend) return -1000001;
            if (gBend<=gAend && trA->exons[trA->nExons-1][EX_iFrag]==iFragB) return -1000002;

            //shift the B 5' if overlaps A 3'
            if (rBstart<=rAend) {
                gBstart+=rAend-rBstart+1;
                rBstart=rAend+1;
                L=rBend-rBstart+1;
            };

            for (uint ii=rBstart;ii<=rBend;ii++) Score+=scoreMatch; //add QS for mapped portions

            int gGap=gBstart-gAend-1; //could be < 0 for insertions
            int rGap=rBstart-rAend-1;//>0 always since we removed overlap

            uint nMatch=L;
            uint nMM=0;
            uint Del=0, Ins=0;
            uint nIns=0, nDel=0;
            int jR=0; //junction location in R-space
            int jCan=999; //canonical junction type
            uint gBstart1=gBstart-rGap-1;//the last base of the intron if all read gap belongs to acceptor, i.e. jR=0


            // check all the different combinations of gGap and rGap
            if ( gGap==0 && rGap==0 ) {//just joined the pieces, w/o stiching or gaps
                //do nothing for now
            } else if ( gGap>0 && rGap>0 && rGap==gGap ) {//no gaps, just try to fill space
                //simple stitching, assuming no insertion in the read

                for (int ii=1;ii<=rGap;ii++) {
                    if (G[gAend+ii]<4 && R[rAend+ii]<4) {//only score genome bases that are not Ns
                        if ( R[rAend+ii]==G[gAend+ii] ) {
                            Score+=scoreMatch;
                            nMatch++;
                        } else {
                            Score-=scoreMatch;
                            nMM++;
                        };
                    };
                };

            } else if ( gGap>rGap ) {//genomic gap (Deletion)

                nDel=1;
                Del=gGap-rGap; //gGap>0 here

                if (Del>P.alignIntronMax && P.alignIntronMax>0) {
                    return -1000003; //large gaps not allowed
                };

                int Score1=0;
                int jR1=1; //junction location in R-space
                do { // 1. move left, until the score for MM is less than canonical advantage
                    jR1--;
                    if ( R[rAend+jR1]!=G[gBstart1+jR1] && G[gBstart1+jR1]<4 && R[rAend+jR1]==G[gAend+jR1]) Score1 -= scoreMatch;
                }  while ( Score1+P.scoreStitchSJshift >= 0 && int(trA->exons[trA->nExons-1][EX_L]) + jR1 > 1);//>=P.alignSJoverhangMin); //also check that we are still within the exon

                int maxScore2=-999999;
                Score1=0;
                int jPen=0;
                do { // 2. scan to the right to find the best junction locus
                    // ?TODO? if genome base is N, how to score?
                    if  ( R[rAend+jR1]==G[gAend+jR1] && R[rAend+jR1]!=G[gBstart1+jR1] )  Score1+=scoreMatch;
                    if  ( R[rAend+jR1]!=G[gAend+jR1] && R[rAend+jR1]==G[gBstart1+jR1] )  Score1-=scoreMatch;

                    int jCan1=-1; //this marks Deletion
                    int jPen1=0;
                    int Score2=Score1;

                    if (Del>=P.alignIntronMin) {//only check intron motif for large gaps= non-Dels
                        //check if the intron is canonical, or semi-canonical
                        if ( G[gAend+jR1+1]==2 && G[gAend+jR1+2]==3 && G[gBstart1+jR1-1]==0 && G[gBstart1+jR1]==2 ) {//GTAG
                            jCan1=1;
                        } else if ( G[gAend+jR1+1]==1 && G[gAend+jR1+2]==3 && G[gBstart1+jR1-1]==0 && G[gBstart1+jR1]==1 ) {//CTAC
                            jCan1=2;
                        } else if ( G[gAend+jR1+1]==2 && G[gAend+jR1+2]==1 && G[gBstart1+jR1-1]==0 && G[gBstart1+jR1]==2 ) {//GCAG
                            jCan1=3;
                            jPen1=P.scoreGapGCAG;
                        } else if ( G[gAend+jR1+1]==1 && G[gAend+jR1+2]==3 && G[gBstart1+jR1-1]==2 && G[gBstart1+jR1]==1 ) {//CTGC
                            jCan1=4;
                            jPen1=P.scoreGapGCAG;
                        } else if ( G[gAend+jR1+1]==0 && G[gAend+jR1+2]==3 && G[gBstart1+jR1-1]==0 && G[gBstart1+jR1]==1 ) {//ATAC
                            jCan1=5;
                            jPen1=P.scoreGapATAC;
                        } else if ( G[gAend+jR1+1]==2 && G[gAend+jR1+2]==3 && G[gBstart1+jR1-1]==0 && G[gBstart1+jR1]==3 ) {//GTAT
                            jCan1=6;
                            jPen1=P.scoreGapATAC;
                        } else {
                            jCan1=0;
                            jPen1=P.scoreGapNoncan;
                        };

                        Score2 += jPen1;
                    };

                    if (maxScore2 < Score2 ) {//check if the score is the highest. TODO: record the next highest score
                        maxScore2=Score2;
                        jR=jR1; //this is the last base of donor
                        jCan=jCan1;
                        jPen=jPen1;
                    };
                        jR1++;
                } while ( jR1 < int(rBend) - int(rAend) );// - int(P.alignSJoverhangMin) );//TODO: do not need to search the full B-transcript, can stop as soon as Score goes down by more than

                //repeat length: go back and forth around jR to find repeat length
                uint jjL=0,jjR=0;
                while ( gAend+jR>=jjL && G[gAend-jjL+jR]==G[gBstart1-jjL+jR] && G[gAend-jjL+jR]<4 && jjL<=MAX_SJ_REPEAT_SEARCH) {//go back
                    jjL++;
                };

                while ( gAend+jjR+jR+1<mapGen.nGenome && G[gAend+jjR+jR+1]==G[gBstart1+jjR+jR+1] && G[gAend+jjR+jR+1]<4 && jjR<=MAX_SJ_REPEAT_SEARCH) {//go forward
                    jjR++;
                };

                if (jCan<=0) {//flush deletions and non-canonical junction to the left
                    jR-=jjL;
                    if (int(trA->exons[trA->nExons-1][EX_L])+jR<1) return -1000005;
                    jjR+=jjL;
                    jjL=0;
                };

                //TODO check here if the internal exon length < minDa, if so exit w/o stitiching

                for (int ii=min(1,jR+1);ii<=max(rGap,jR);ii++) {//score donor and acceptor
                    uint g1=(ii<=jR) ? (gAend+ii):(gBstart1+ii);
                    if (G[g1]<4 && R[rAend+ii]<4) {//only penalize non-N bases
                        if ( R[rAend+ii]==G[g1] ) {
                            if (ii>=1 && ii <=rGap) {//only add +score and matches within the gap
                                Score+=scoreMatch;
                                nMatch++;
                            };
                        } else {//add -score and MM for all bases
                            Score-=scoreMatch;
                            nMM++;
                            if (ii<1 || ii>rGap) {//subtract previuosly presumed matches
                                Score-=scoreMatch;
                                nMatch--;
//                                 if (ii<=jR) nMM--;
                            };
                        };
                    };
                };

                //score the gap
                if (mapGen.sjdbN>0) {//check if the junction is annotated
                        uint jS=gAend+jR+1, jE=gBstart1+jR;//intron start/end
                        int sjdbInd=binarySearch2(jS,jE,mapGen.sjdbStart,mapGen.sjdbEnd,mapGen.sjdbN);
                        if (sjdbInd<0) {
                            if (Del>=P.alignIntronMin) {
                                Score += P.scoreGap + jPen; //genome gap penalty + non-canonical penalty
                            } else {//deletion
                                Score += Del*P.scoreDelBase + P.scoreDelOpen;
                                jCan=-1;
                                trA->sjAnnot[trA->nExons-1]=0;
//                                 jjR-=jjL;
//                                 jR-=jjL;
//                                 jjL=0;
//                                 trA->shiftSJ[trA->nExons-1][0]=0;
//                                 trA->shiftSJ[trA->nExons-1][1]=jjR;
                            };
                        } else {//annotated
                            jCan=mapGen.sjdbMotif[sjdbInd];
                            if (mapGen.sjdbMotif[sjdbInd]==0) {//shift to match annotations
                                if (L<=mapGen.sjdbShiftLeft[sjdbInd] || trA->exons[trA->nExons-1][EX_L]<=mapGen.sjdbShiftLeft[sjdbInd]) {
                                    return -1000006;
                                };
                                jR += (int) mapGen.sjdbShiftLeft[sjdbInd];
                                jjL=mapGen.sjdbShiftLeft[sjdbInd];
                                jjR=mapGen.sjdbShiftRight[sjdbInd];
                            };
                            trA->sjAnnot[trA->nExons-1]=1;
                            trA->sjStr[trA->nExons-1]=mapGen.sjdbStrand[sjdbInd];
                            Score += P.pGe.sjdbScore;
                        };
                } else {//no annotation
                    if (Del>=P.alignIntronMin) {//junction, not short deletion
                        Score += P.scoreGap + jPen;
                    } else {
                        Score += Del*P.scoreDelBase + P.scoreDelOpen;
                        jCan=-1;
                        trA->sjAnnot[trA->nExons-1]=0;
                    };
                };

                trA->shiftSJ[trA->nExons-1][0]=jjL;
                trA->shiftSJ[trA->nExons-1][1]=jjR;
                trA->canonSJ[trA->nExons-1]=jCan;

                if (trA->sjAnnot[trA->nExons-1]==0) {//strand for unannotated junctions
                    if (jCan>0) {
                         trA->sjStr[trA->nExons-1]=2-jCan%2; //1=+,2=-
                    } else {
                         trA->sjStr[trA->nExons-1]=0;
                    };
                };

            } else if ( rGap>gGap ) {//insertion: if also gGap>0, need to stitch
                Ins=rGap-gGap;
                nIns=1;
                if (gGap==0) {//simple insertion, no need to stitch
                    jR=0;
                } else if (gGap<0) {//overlapping seeds: reduce the score
                    jR=0;
                    for (int ii=0; ii<-gGap; ii++) {
                        Score -= scoreMatch;
                    };
                } else {//stitch: define the exon boundary jR
                    int Score1=0; int maxScore1=0;
                    for (int jR1=1;jR1<=gGap;jR1++) {//scan to the right to find the best score

                        if (G[gAend+jR1]<4) {//only penalize goog genome bases
                            Score1+=( R[rAend+jR1]==G[gAend+jR1] ) ? scoreMatch:-scoreMatch;
                            Score1+=( R[rAend+Ins+jR1]==G[gAend+jR1] ) ? -scoreMatch:+scoreMatch;
                        };

                        if (Score1>maxScore1 || (Score1==maxScore1 && P.alignInsertionFlush.flushRight)) {//equal sign (>=) flushes insertions to the right
                            maxScore1=Score1;
                            jR=jR1;
                        };
                    };
                    for (int ii=1;ii<=gGap;ii++) {//score donor and acceptor
                        uint r1=rAend+ii+(ii<=jR ? 0:Ins);
                        if (G[gAend+ii]<4 && R[r1]<4) {
                            if ( R[r1]==G[gAend+ii] ) {
                                Score+=scoreMatch;
                                nMatch++;
                            } else {//add -score and MM for all bases
                                Score-=scoreMatch;
                                nMM++;
                            };
                        };
                    };
                };

                if (P.alignInsertionFlush.flushRight) {
                    for (; jR<(int)rBend-(int)rAend-(int)Ins; jR++ ){//flush the indel to the right as much as possible
                        if (R[rAend+jR+1]!=G[gAend+jR+1] || G[gAend+jR+1]==4) {
                            break;
                        };
                    };
                    if (jR==(int)rBend-(int)rAend-(int)Ins) {//nothing left of the B-piece
                        return -1000009;
                    };
                };
                Score += Ins*P.scoreInsBase + P.scoreInsOpen;
                jCan=-3;
            }; //different types of gaps selection



        #ifdef COMPILE_FOR_LONG_READS
            if ( (trA->nMM + nMM)<=outFilterMismatchNmaxTotal )
//             if ( Score>0 && nMM<=200 )

        #else
            if ( (trA->nMM + nMM)<=outFilterMismatchNmaxTotal  \
                         && ( jCan<0 || (jCan<7 && nMM<= (uint) P.alignSJstitchMismatchNmax[(jCan+1)/2]) ) )
        #endif
            {//stitching worked only if there no mis-matches for non-GT/AG junctions
                trA->nMM += nMM;
                trA->nMatch += nMatch;

                if (Del>=P.alignIntronMin) {
                    trA->nGap += nDel;
                    trA->lGap += Del;
                } else {
                    trA->nDel += nDel;
                    trA->lDel += Del;
                };

                //modify exons
                if (Del==0 && Ins==0) {//no gap => no new exon, extend the boundary of the previous exon
                    trA->exons[trA->nExons-1][EX_L] += rBend-rAend;
                } else if (Del>0) { //deletion:ca only have Del> or Ins>0
                    trA->exons[trA->nExons-1][EX_L] += jR; //correct the previous exon boundary
                    trA->exons[trA->nExons][EX_L] = rBend-rAend-jR; //new exon length
                    trA->exons[trA->nExons][EX_R] = rAend+jR+1; //new exon r-start
                    trA->exons[trA->nExons][EX_G] = gBstart1+jR+1; //new exon g-start
                    trA->nExons++;
                } else if (Ins>0) { //Ins>0;
                    trA->nIns += nIns;
                    trA->lIns += Ins;
                    trA->exons[trA->nExons-1][EX_L] += jR; //correct the previous exon boundary
                    trA->exons[trA->nExons][EX_L] = rBend-rAend-jR-Ins; //new exon length
                    trA->exons[trA->nExons][EX_R] = rAend+jR+Ins+1; //new exon r-start
                    trA->exons[trA->nExons][EX_G] = gAend+1+jR; //new exon g-start
                    trA->canonSJ[trA->nExons-1]=-2; //mark insertion
                    trA->sjAnnot[trA->nExons-1]=0;
                    trA->nExons++;
                };
            } else {//stitching was not accepted
                return -1000007;
            };
        } else if (gBstart+trA->exons[0][EX_R]+P.alignEndsProtrude.nBasesMax >= trA->exons[0][EX_G] || trA->exons[0][EX_G] < trA->exons[0][EX_R]){//if (iFragA==iFragB) stitch aligns from different fragments
                                                                                                            //CHECK: this second confdition does not make sense
            if (P.alignMatesGapMax>0 && gBstart > trA->exons[trA->nExons-1][EX_G] + trA->exons[trA->nExons-1][EX_L] + P.alignMatesGapMax) {
                return -1000004; //gap between mates too large
            };
            //extend the fragments inside
            //note, that this always works, i.e. Score>0

            for (uint ii=rBstart;ii<rBstart+L;ii++) Score+=scoreMatch; //add QS for mapped portions

            Transcript trExtend;

            //TODO: compare extensions to the left and right, pick the best one to be performed first
            //otherwise if a large nMM is reached in the 2st extension, it will prevent the 2nd extension
            //use the following example:
            //>1
            //TTCTGTGTCTCCCCCTCCCCCACTGGCTACATGGAGACAGGGGGGGGGGGCCGGGCGGTTCCCGGGCAGAAAAAAA
            //>1
            //AATATTTGGAACACTTATGTGAAAAATGATTTGTTTTTCTGAAATTTACGTTTCTCTCTGAGTCCTGTAACTGTCC


            trExtend.reset();
            if ( extendAlign(R, G, rAend+1, gAend+1, 1, 1, DEF_readSeqLengthMax, trA->nMatch, trA->nMM, outFilterMismatchNmaxTotal, P.outFilterMismatchNoverLmax, \
                             P.alignEndsType.ext[trA->exons[trA->nExons-1][EX_iFrag]][1], &trExtend) ) {

                trA->add(&trExtend);
                Score += trExtend.maxScore;

                trA->exons[trA->nExons-1][EX_L] += trExtend.extendL;
            };// if extendAlign for read A

            trA->exons[trA->nExons][EX_R] = rBstart;
            trA->exons[trA->nExons][EX_G] = gBstart;
            trA->exons[trA->nExons][EX_L] = L;
            trA->nMatch += L;

            trExtend.reset();
            //if end extension needs to be forced, use large length. Otherwise, only extend until the beginning of the transcript
            uint extlen=P.alignEndsType.ext[iFragB][1] ? DEF_readSeqLengthMax : gBstart-trA->exons[0][EX_G]+trA->exons[0][EX_R];
            if ( extendAlign(R, G, rBstart-1, gBstart-1, -1, -1, extlen, trA->nMatch, trA->nMM, outFilterMismatchNmaxTotal, P.outFilterMismatchNoverLmax, \
                             P.alignEndsType.ext[iFragB][1], &trExtend) ) {

                trA->add(&trExtend);
                Score += trExtend.maxScore;

                trA->exons[trA->nExons][EX_R] -= trExtend.extendL;
                trA->exons[trA->nExons][EX_G] -= trExtend.extendL;
                trA->exons[trA->nExons][EX_L] += trExtend.extendL;
            }; //if extendAlign B

            trA->canonSJ[trA->nExons-1]=-3; //mark different fragments junction
            trA->sjAnnot[trA->nExons-1]=0;

            trA->nExons++;
        } else {//no stitching possible
            return -1000008;
        };
    };

    trA->exons[trA->nExons-1][EX_iFrag]=iFragB; //the new exon belongs to fragment iFragB
    trA->exons[trA->nExons-1][EX_sjA]=sjAB;

    return Score;
};
