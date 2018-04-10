#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ReadAlign.h"
#include "stitchAlignToTranscript.h"
#include "extendAlign.h"
#include <math.h>
#include "binarySearch2.h"
#include "ErrorWarning.h"

void ReadAlign::stitchWindowSeeds (uint iW, uint iWrec, bool *WAexcl, char *R) {//stitches all seeds in one window: iW

    for (uint iS1=0;iS1<nWA[iW];iS1++) {
        scoreSeedBest[iS1]=0;
        scoreSeedBestMM[iS1]=0;
        scoreSeedBestInd[iS1]=-1;
        if (WAexcl!=NULL && WAexcl[iS1]) continue; //do not include this seed
        intScore score2=0;
        for (uint iS2=0;iS2<=iS1;iS2++) {
            if (WAexcl!=NULL && WAexcl[iS1]) continue; //do not include this seed
            trA1=*trInit;//initialize trA1
            if (iS2<iS1) {
                trA1.nExons=1;
                trA1.nMM=scoreSeedBestMM[iS2];
                trA1.exons[0][EX_R] = WA[iW][iS2][WA_rStart];
                trA1.exons[0][EX_G] = WA[iW][iS2][WA_gStart];
                trA1.exons[0][EX_L] = WA[iW][iS2][WA_Length];
                trA1.exons[0][EX_iFrag]=WA[iW][iS2][WA_iFrag];
                trA1.exons[0][EX_sjA]=WA[iW][iS2][WA_sjA];
                score2=\
                    stitchAlignToTranscript(WA[iW][iS2][WA_rStart]+WA[iW][iS2][WA_Length]-1, WA[iW][iS2][WA_gStart]+WA[iW][iS2][WA_Length]-1,\
                                        WA[iW][iS1][WA_rStart], WA[iW][iS1][WA_gStart], WA[iW][iS1][WA_Length], WA[iW][iS1][WA_iFrag],  WA[iW][iS1][WA_sjA], \
                                        P, R, mapGen, &trA1, outFilterMismatchNmaxTotal);

                if (P.outFilterBySJoutStage==2 && trA1.nExons>1)
                {//junctions have to be present in the filtered set P.sjnovel
                    uint iex=0;
                    if (trA1.canonSJ[iex]>=0 && trA1.sjAnnot[iex]==0)
                    {
                        uint jS=trA1.exons[iex][EX_G]+trA1.exons[iex][EX_L];
                        uint jE=trA1.exons[iex+1][EX_G]-1;
                        if ( binarySearch2(jS,jE,P.sjNovelStart,P.sjNovelEnd,P.sjNovelN) < 0 ) return;
                    };
                };

                //check the length of the iS2 exon. TODO: build the transcripts vs iS1, check the actual exon length
                bool exonLongEnough = trA1.exons[0][EX_L] >= ( trA1.sjAnnot[0]==0 ? P.alignSJoverhangMin : P.alignSJDBoverhangMin );
                
                if (exonLongEnough && score2>0 && score2+scoreSeedBest[iS2] > scoreSeedBest[iS1] ) {
                    scoreSeedBest[iS1]=score2+scoreSeedBest[iS2];
                    scoreSeedBestMM[iS1]=trA1.nMM;
                    scoreSeedBestInd[iS1]=iS2;
                };
            } else {//extend to the left
                score2=WA[iW][iS1][WA_Length];
                if ( WA[iW][iS1][WA_rStart]>0 \
                     && extendAlign(R, mapGen.G, WA[iW][iS1][WA_rStart]-1, WA[iW][iS1][WA_gStart]-1, -1, -1, WA[iW][iS1][WA_rStart], 100000, 0, outFilterMismatchNmaxTotal, P.outFilterMismatchNoverLmax, \
                                    P.alignEndsType.ext[WA[iW][iS1][WA_iFrag]][trA.Str], &trA1) ) {//if could extend
                    score2 += trA1.maxScore;
                };
                
                bool exonLongEnough = (WA[iW][iS1][WA_Length]+trA1.extendL) >= P.alignSJoverhangMin; //TODO new parameter to control end exons length
                
                if (exonLongEnough && score2 > scoreSeedBest[iS1] ) {
                    scoreSeedBest[iS1]=score2;
                    scoreSeedBestInd[iS1]=iS1;
//                     scoreSeedBestMM[iS1]=trA1.nMM;
                };
            };
        };
    };

    intScore scoreBest=0;
    uint scoreBestInd=0;
    

    for (uint iS1=0;iS1<nWA[iW];iS1++) {//find the best alignment
        trA1=*trInit;//initialize trA1
        uint tR2=WA[iW][iS1][WA_rStart]+WA[iW][iS1][WA_Length];
        uint tG2=WA[iW][iS1][WA_gStart]+WA[iW][iS1][WA_Length];
        if ( tR2 < Lread-1 \
            && extendAlign(R, mapGen.G, tR2, tG2, +1, +1, Lread-tR2, 100000, scoreSeedBestMM[iS1], outFilterMismatchNmaxTotal, P.outFilterMismatchNoverLmax, \
                           P.alignEndsType.ext[WA[iW][iS1][WA_iFrag]][1-trA.Str], &trA1) )
        {//extend to the right
            scoreSeedBest[iS1]+=trA1.maxScore;
        };
        
        bool exonLongEnough = (WA[iW][iS1][WA_Length]+trA1.extendL) >= P.alignSJoverhangMin; //TODO new parameter to control end exons length
        
        if (exonLongEnough && scoreSeedBest[iS1]>scoreBest) {//record new best transcript
            scoreBest=scoreSeedBest[iS1];
            scoreBestInd=iS1;
        };
    };

    uint seedN=0;
    while (true) {//construct the sequence of seeds
        seedChain[seedN++]=scoreBestInd;
        WAincl[scoreBestInd]=true;
        if (scoreBestInd>scoreSeedBestInd[scoreBestInd]){//keep going
            scoreBestInd=scoreSeedBestInd[scoreBestInd];
        } else {//this seed is the first one
            break;
        };
    };

    int Score=0;
    {//build final transcript form seedChain
        {//initiate transcript

            uint iS1=seedChain[seedN-1];
            Score= WA[iW][iS1][WA_Length];
            trA.maxScore = Score;
            trA.nMatch = WA[iW][iS1][WA_Length]; //# of matches
            trA.nMM = 0;

            trA.exons[0][EX_R] = trA.rStart = WA[iW][iS1][WA_rStart];
            trA.exons[0][EX_G] = trA.gStart = WA[iW][iS1][WA_gStart];
            trA.exons[0][EX_L] = WA[iW][iS1][WA_Length];
            trA.exons[0][EX_iFrag]=WA[iW][iS1][WA_iFrag];
            trA.exons[0][EX_sjA]=WA[iW][iS1][WA_sjA];

            trA.nExons=1;

        };

        for (uint iSc=seedN-1; iSc>0; iSc--) {//stitch seeds from the chain
            uint iS1=seedChain[iSc], iS2=seedChain[iSc-1];
            int scoreStitch= stitchAlignToTranscript(WA[iW][iS1][WA_rStart]+WA[iW][iS1][WA_Length]-1, WA[iW][iS1][WA_gStart]+WA[iW][iS1][WA_Length]-1,\
                                        WA[iW][iS2][WA_rStart], WA[iW][iS2][WA_gStart], WA[iW][iS2][WA_Length], WA[iW][iS2][WA_iFrag],  WA[iW][iS2][WA_sjA], \
                                        P, R, mapGen, &trA, outFilterMismatchNmaxTotal);
//            if (scoreStitch>0) {
                Score+=scoreStitch;
//           } else {
//                 cout <<"BUG"<<endl;
//                return;//this should not happen
//            };
        };

        trA.maxScore=Score;

        {//extend to the left
            trA1=*trInit;
            if ( trA.exons[0][EX_R]>0 \
                 && extendAlign(R, mapGen.G, trA.exons[0][EX_R]-1, trA.exons[0][EX_G]-1, -1, -1, trA.exons[0][EX_R], 100000, 0, outFilterMismatchNmaxTotal, P.outFilterMismatchNoverLmax,
                    P.alignEndsType.ext[trA.exons[0][EX_iFrag]][trA.Str], &trA1) ) {//if could extend

                trA.add(&trA1);

                trA.exons[0][EX_R] -=  trA1.extendL;
                trA.exons[0][EX_G] -=  trA1.extendL;
                trA.exons[0][EX_L] +=  trA1.extendL;
                trA.rStart = trA.exons[0][EX_R];
                trA.gStart = trA.exons[0][EX_G];
            };
        };

        {//extend to the right
            uint iS1=seedChain[0];
            trA1=*trInit;//initialize trA1
            uint tR2=WA[iW][iS1][WA_rStart]+WA[iW][iS1][WA_Length];
            uint tG2=WA[iW][iS1][WA_gStart]+WA[iW][iS1][WA_Length];
            if ( tR2 < Lread \
                && extendAlign(R, mapGen.G, tR2, tG2, +1, +1, Lread-tR2, 100000, scoreSeedBestMM[iS1], outFilterMismatchNmaxTotal, P.outFilterMismatchNoverLmax, \
                    P.alignEndsType.ext[trA.exons[trA.nExons-1][EX_iFrag]][1-trA.Str], &trA1) ) {//if could extend
                    trA.add(&trA1);
                    trA.exons[trA.nExons-1][EX_L] += trA1.extendL;//extend the length of the last exon
            };
        };

    };

    //debug: recalculate the number of MM
//     {
//         uint nMM1=0;
//         for (uint iex=0;iex<trA.nExons;iex++) {
//             for (uint ii=0;ii<trA.exons[iex][EX_L];ii++) {
//                 if ( R[ii+trA.exons[iex][EX_R]]!=G[ii+trA.exons[iex][EX_G]] && mapGen.G[ii+trA.exons[iex][EX_G]]<4 && R[ii+trA.exons[iex][EX_R]]<4) {
//                     nMM1++;
//                 };
//             };
//         };
//         if (nMM1!=trA.nMM) {
//             cout <<nMM1<<"   "<<trA.nMM<<"    "<<readName<<"   "<<iRead<<endl;
//         };
//     };

    {//calculate some final values for the transcript
        trA.rLength = 0;
        for (uint isj=0;isj<trA.nExons;isj++) {
            trA.rLength += trA.exons[isj][EX_L];
        };
        trA.gLength = trA.exons[trA.nExons-1][EX_G]+1-trA.gStart;

        //calculate some final values for the transcript
        trA.roStart = (trA.roStr == 0) ? trA.rStart : Lread - trA.rStart - trA.rLength;

        if (trA.exons[0][EX_iFrag]==trA.exons[trA.nExons-1][EX_iFrag]) {//mark single fragment transcripts
            trA.iFrag=trA.exons[0][EX_iFrag];
            maxScoreMate[trA.iFrag] = max (maxScoreMate[trA.iFrag] , trA.maxScore);
        } else {
            trA.iFrag=-1;
        };

        if (P.scoreGenomicLengthLog2scale!=0) {//add gap length score
            trA.maxScore += int(ceil( log2( (double) ( trA.exons[trA.nExons-1][EX_G]+trA.exons[trA.nExons-1][EX_L] - trA.exons[0][EX_G]) ) \
                     * P.scoreGenomicLengthLog2scale - 0.5));
            trA.maxScore = max(0,trA.maxScore);
        };

        //filter strand consistency
        uint sjN=0;
        trA.intronMotifs[0]=0;trA.intronMotifs[1]=0;trA.intronMotifs[2]=0;        
        for (uint iex=0;iex<trA.nExons-1;iex++) {
            if (trA.canonSJ[iex]>=0) 
            {//junctions - others are indels
                sjN++;
                trA.intronMotifs[trA.sjStr[iex]]++;
            };
        };

        if (trA.intronMotifs[1]>0 && trA.intronMotifs[2]==0)
            trA.sjMotifStrand=1;
        else if (trA.intronMotifs[1]==0 && trA.intronMotifs[2]>0)
            trA.sjMotifStrand=2;
        else
            trA.sjMotifStrand=0;

        if (trA.intronMotifs[1]>0 && trA.intronMotifs[2]>0 && P.outFilterIntronStrands=="RemoveInconsistentStrands")
                return;        
        
        if (sjN>0 && trA.sjMotifStrand==0 && P.outSAMstrandField=="intronMotif") {//strand not defined for a junction
            return;
        };

        if (P.outFilterIntronMotifs=="None") {//no filtering

        } else if (P.outFilterIntronMotifs=="RemoveNoncanonical") {
            for (uint iex=0;iex<trA.nExons-1;iex++) {
                if (trA.canonSJ[iex]==0) return;
            };
        } else if (P.outFilterIntronMotifs=="RemoveNoncanonicalUnannotated") {
            for (uint iex=0;iex<trA.nExons-1;iex++) {
                if (trA.canonSJ[iex]==0 && trA.sjAnnot[iex]==0) return;
            };
        } else {
            ostringstream errOut;
            errOut << "EXITING because of FATAL INPUT error: unrecognized value of --outFilterIntronMotifs=" <<P.outFilterIntronMotifs <<"\n";
            errOut << "SOLUTION: re-run STAR with --outFilterIntronMotifs = None -OR- RemoveNoncanonical -OR- RemoveNoncanonicalUnannotated\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };

//         if (P.outFilterIntronMotifs=="KeepCanonical" && (trA.intronMotifs[0]>0 || (trA.intronMotifs[1]>0 && trA.intronMotifs[2]>0) ) ) {//keep only conistent canonical introns
//             return;
//         };


        //check exons lenghts including repeats, do not report a transcript with short exons
//        for (uint isj=0;isj<trA.nExons-1;isj++) {//check exons for min length, if they precede a junction
//            if ( trA.canonSJ[isj]>=0 &&
//               ( trA.exons[isj][EX_L] < P.alignSJoverhangMin + trA.shiftSJ[isj][0]
//              || trA.exons[isj+1][EX_L] < P.alignSJoverhangMin + trA.shiftSJ[isj][1]) ) {
//                return;//do not record this transcript in wTr
//            };
//        };
    };

    if (WAexcl==NULL)
    {//record the transcript TODO: allow for multiple transcripts in one window
        *(trAll[iWrec][0])=trA;
        nWinTr[iWrec]=1;
    } else
    {//record 2nd best alignment in this window
        *(trAll[iWrec][1])=trA;
        nWinTr[iWrec]=2;
    };
};
