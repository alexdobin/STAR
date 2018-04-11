#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ReadAlign.h"
#include "SequenceFuns.h"
#include "stitchWindowAligns.h"
#include "sjAlignSplit.cpp"
#include "PackedArray.h"
#include "alignSmithWaterman.h"
#include "GlobalVariables.h"
#include <time.h>

void ReadAlign::stitchPieces(char **R, uint Lread) {

    //zero-out winBin
    memset(winBin[0],255,sizeof(winBin[0][0])*P.winBinN);
    memset(winBin[1],255,sizeof(winBin[0][0])*P.winBinN);


//     for (uint iWin=0;iWin<nWall;iWin++) {//zero out winBin
//         if (WC[iWin][WC_gStart]<=WC[iWin][WC_gEnd]) {//otherwise the window is dead
//             memset(&(winBin[WC[iWin][WC_Str]][WC[iWin][WC_gStart]]),255,sizeof(winBin[0][0])*(WC[iWin][WC_gEnd]-WC[iWin][WC_gStart]+1));
//         };
// //         for (uint ii=C[iWin][WC_gStart]; ii<WC[iWin][WC_gEnd]; ii++) {
// //             winBin[WC[WC_Str]
// //         };
//     };

//     //debug
//     for (uint ii=0;ii<P.winBinN;ii++){
//         if (winBin[0][ii]!=uintWinBinMax || winBin[1][ii]!=uintWinBinMax) {
//             cerr<< "BUG in stitchPieces: ii="<<ii<<"   "<< winBin[0][ii] <<"   "<<winBin[1][ii] <<"   iRead="<<iRead<<"   nW="<<nW<<endl;
//             for (uint iWin=0;iWin<nW;iWin++) {
//                 cerr <<WC[iWin][WC_gStart]<<"   " <<WC[iWin][WC_gEnd] <<"   "<<WC[iWin][WC_Str] <<endl;
//             };
//             exit(1);
//         };
//     };


    nW=0; //number of windows
    for (uint iP=0; iP<nP; iP++) {//scan through all anchor pieces, create alignment windows

//          if (PC[iP][PC_Nrep]<=P.winAnchorMultimapNmax || PC[iP][PC_Length]>=readLength[PC[iP][PC_iFrag]] ) {//proceed if piece is an anchor, i.e. maps few times or is long enough
       if (PC[iP][PC_Nrep]<=P.winAnchorMultimapNmax ) {//proceed if piece is an anchor, i.e. maps few times

            uint aDir   = PC[iP][PC_Dir];
            uint aLength= PC[iP][PC_Length];

            for (uint iSA=PC[iP][PC_SAstart]; iSA<=PC[iP][PC_SAend]; iSA++) {//scan through all alignments of this piece

                uint a1 = mapGen.SA[iSA];
                uint aStr = a1 >> mapGen.GstrandBit;
                a1 &= mapGen.GstrandMask; //remove strand bit

                //convert to positive strand
                if (aDir==1 && aStr==0) {
                    aStr=1;
                } else if (aDir==0 && aStr==1) {
                    a1 = mapGen.nGenome - (aLength+a1);
                } else if (aDir==1 && aStr==1) {
                    aStr=0;
                    a1 = mapGen.nGenome - (aLength+a1);
                };

                //final strand
                if (revertStrand) { //modified strand according to user input CHECK!!!!
                    aStr=1-aStr;
                };

                if (a1>=mapGen.sjGstart) {//this is sj align
                    uint a1D, aLengthD, a1A, aLengthA, sj1;
                    if (sjAlignSplit(a1, aLength, mapGen, a1D, aLengthD, a1A, aLengthA, sj1)) {//align crosses the junction

                        int addStatus=createExtendWindowsWithAlign(a1D, aStr);//add donor piece
                        if (addStatus==EXIT_createExtendWindowsWithAlign_TOO_MANY_WINDOWS) {//too many windows
                            break;
                        };
                        addStatus=createExtendWindowsWithAlign(a1A, aStr);//add acceptor piece
                        if (addStatus==EXIT_createExtendWindowsWithAlign_TOO_MANY_WINDOWS) {//too many windows
                            break;
                        };
                    };
                } else {//this is a normal genomic read
                    int addStatus=createExtendWindowsWithAlign(a1, aStr);
                    if (addStatus==EXIT_createExtendWindowsWithAlign_TOO_MANY_WINDOWS) {//too many windows
                        break;
                    };
                };
            }; //for (uint iSA=PC[iP][PC_SAstart]; iSA<=PC[iP][PC_SAend]; iSA++) //scan through all alignments of this piece
        };//if (PC[iP][PC_Nrep]<=P.winAnchorMultimapNmax) //proceed if anchor
    };//for (uint iP=0; iP<nP; iP++) //scan through all anchor pieces, create alignment windows


    for (uint iWin=0;iWin<nW;iWin++) {//extend windows with flanks
        if (WC[iWin][WC_gStart]<=WC[iWin][WC_gEnd]) {//otherwise the window is dead

            uint wb=WC[iWin][WC_gStart];
            for (uint ii=0; ii<P.winFlankNbins && wb>0 && mapGen.chrBin[(wb-1) >> P.winBinChrNbits]==WC[iWin][WC_Chr];ii++) {
                wb--;
                winBin[ WC[iWin][WC_Str] ][ wb ]=(uintWinBin) iWin;
            };
            WC[iWin][WC_gStart] = wb;

            wb=WC[iWin][WC_gEnd];
            for (uint ii=0; ii<P.winFlankNbins && wb+1<P.winBinN && mapGen.chrBin[(wb+1) >> P.winBinChrNbits]==WC[iWin][WC_Chr];ii++) {
                wb++;
                winBin[ WC[iWin][WC_Str] ][ wb ]=(uintWinBin) iWin;
            };
            WC[iWin][WC_gEnd] = wb;


        };
        nWA[iWin]=0; //initialize nWA
        WALrec[iWin]=0; //initialize rec-length
        WlastAnchor[iWin]=-1;
    };

    nWall=nW;

    #ifdef OFF_BEFORE_SEEDdistribution
        #warning OFF_BEFORE_SEEDdistribution
        nW=0;
        nTr=0;
        return;
    #endif

    for (uint iP=0; iP<nP; iP++) {//scan through all pieces/aligns, add them to alignment windows, create alignment coordinates
        uint aNrep=PC[iP][PC_Nrep];
        uint aFrag=PC[iP][PC_iFrag];
        uint aLength=PC[iP][PC_Length];
        uint aDir=PC[iP][PC_Dir];

        bool aAnchor=(aNrep<=P.winAnchorMultimapNmax); //this align is an anchor or not

        for (uint ii=0;ii<nW;ii++) {//initialize nWAP
            nWAP[ii]=0;
        };

        for (uint iSA=PC[iP][PC_SAstart]; iSA<=PC[iP][PC_SAend]; iSA++) {//scan through all alignments

            uint a1 = mapGen.SA[iSA];
            uint aStr = a1 >> mapGen.GstrandBit;
            a1 &= mapGen.GstrandMask; //remove strand bit
            uint aRstart=PC[iP][PC_rStart];

            //convert to positive strand
            if (aDir==1 && aStr==0) {
                aStr=1;
                aRstart = Lread - (aLength+aRstart);
            } else if (aDir==0 && aStr==1) {
                aRstart = Lread - (aLength+aRstart);
                a1 = mapGen.nGenome - (aLength+a1);
            } else if (aDir==1 && aStr==1) {
                aStr=0;
                a1 = mapGen.nGenome - (aLength+a1);
            };

            //final strand
            if (revertStrand) { //modified strand according to user input CHECK!!!!
                aStr=1-aStr;
            };


            if (a1>=mapGen.sjGstart) {//this is sj read
                uint a1D, aLengthD, a1A, aLengthA, isj1;
                if (sjAlignSplit(a1, aLength, mapGen, a1D, aLengthD, a1A, aLengthA, isj1)) {//align crosses the junction

                        assignAlignToWindow(a1D, aLengthD, aStr, aNrep, aFrag, aRstart, aAnchor, isj1);
                        assignAlignToWindow(a1A, aLengthA, aStr, aNrep, aFrag, aRstart+aLengthD, aAnchor, isj1);

                    } else {//align does not cross the junction
                        continue; //do not check this align, continue to the next one
                    };

                } else {//this is a normal genomic read
                    assignAlignToWindow(a1, aLength, aStr, aNrep, aFrag, aRstart, aAnchor, -1);
                };
        };

//         for (uint ii=0;ii<nW;ii++) {//check of some pieces created too many aligns in some windows, and remove those from WA (ie shift nWA indices
//             if (nWAP[ii]>P.seedNoneLociPerWindow) nWA[ii] -= nWAP[ii];
//         };
    };

    //TODO remove windows that have too many alignments
    //aligns are still sorted by original read coordinates, change direction for negative strand
    // DOES NOT HELP!!!
//     for ( uint iW=0;iW<nW;iW++ ) {
//         if (WA[iW][0][WA_rStart]>WA[iW][nWA[iW]-1][WA_rStart]) {//swap
//             for (uint iA=0;iA<nWA[iW]/2;iA++) {
//                 for (uint ii=0;ii<WA_SIZE;ii++) {
//                     uint dummy=WA[iW][iA][ii];
//                     WA[iW][iA][ii]=WA[iW][nWA[iW]-1-iA][ii];
//                     WA[iW][nWA[iW]-1-iA][ii]=dummy;
//                 };
//             };
//         };
//     };

#ifdef COMPILE_FOR_LONG_READS
uint swWinCovMax=0;
for (uint iW=0;iW<nW;iW++) {//check each window
    swWinCov[iW]=0;
    if (nWA[iW]>0) {
        //select good windows by coverage
        uint rLast=0;

        for (uint ia=0; ia<nWA[iW]; ia++) {//calculate coverage from all aligns
            uint L1=WA[iW][ia][WA_Length];
            uint r1=WA[iW][ia][WA_rStart];

            if (r1+L1>rLast+1) {
                if (r1>rLast) {
                    swWinCov[iW] += L1;
                } else {
                    swWinCov[iW] += r1+L1-(rLast+1);
                };
                rLast=r1+L1-1;
            };
        };//for (uint ia=0; ia<nWA[iW]; ia++)

        if (swWinCov[iW]>swWinCovMax) swWinCovMax=swWinCov[iW];
    };//if (nWA[iW]>0)
};//for (uint iW=0;iW<nW;iW++)
for (uint iW=0;iW<nW;iW++) {
    if (swWinCov[iW]<swWinCovMax*P.winReadCoverageRelativeMin || swWinCov[iW]<P.winReadCoverageBasesMin) {//remove windows that are not good enough
        nWA[iW]=0;
    } else {//merge pieces that are adjacent in R- and G-spaces
        uint ia1=0;
        for (uint ia=1; ia<nWA[iW]; ia++) {
            if ( WA[iW][ia][WA_rStart] == (WA[iW][ia1][WA_rStart]+WA[iW][ia1][WA_Length]) \
              && WA[iW][ia][WA_gStart] == (WA[iW][ia1][WA_gStart]+WA[iW][ia1][WA_Length]) \
              && WA[iW][ia][WA_iFrag]  ==  WA[iW][ia1][WA_iFrag]     ) {//merge

                WA[iW][ia1][WA_Length] += WA[iW][ia][WA_Length];
                WA[iW][ia1][WA_Anchor]=max(WA[iW][ia1][WA_Anchor],WA[iW][ia][WA_Anchor]);
                //NOTE: I am not updating sjA and Nrep fields - this could cause trouble in some cases

            } else {//do not merge
                ia1++;
                if (ia1!=ia) {//move from ia to ia1
                    for (uint ii=0; ii<WA_SIZE; ii++) {
                        WA[iW][ia1][ii]=WA[iW][ia][ii];
                    };
                };
            };
        };
        nWA[iW]=ia1+1;
    };
};

//mapping time initialize
std::time(&timeStart);
#endif

    #ifdef OFF_BEFORE_STITCH
        #warning OFF_BEFORE_STITCH
        nW=0;
        return;
    #endif
    //generate transcript for each window, choose the best
    trBest =trInit; //initialize next/best
    uint iW1=0;//index of non-empty windows
    uint trNtotal=0; //total number of recorded transcripts

    for (uint iW=0; iW<nW; iW++) {//transcripts for all windows

        if (nWA[iW]==0) continue; //the window does not contain any aligns because it was merged with other windows

//         {//debug
//             if ( WA[iW][0][WA_iFrag]==WA[iW][nWA[iW]-1][WA_iFrag] ) continue;
//         };
//
        if (WlastAnchor[iW]<nWA[iW]) {
            WA[ iW ][ WlastAnchor[iW] ][ WA_Anchor]=2; //mark the last anchor
        };

        for (uint ii=0;ii<nWA[iW];ii++) WAincl[ii]=false; //initialize mask

        trA=*trInit; //that one is initialized
        trA.Chr = WC[iW][WC_Chr];
        trA.Str = WC[iW][WC_Str];
        trA.roStr = revertStrand ? 1-trA.Str : trA.Str; //original strand of the read
        trA.maxScore=0;

        trAll[iW1]=trArrayPointer+trNtotal;
        if (trNtotal+P.alignTranscriptsPerWindowNmax >= P.alignTranscriptsPerReadNmax) {
            P.inOut->logMain << "WARNING: not enough space allocated for transcript. Did not process all windows for read "<< readName+1 <<endl;
            P.inOut->logMain <<"   SOLUTION: increase alignTranscriptsPerReadNmax and re-run\n" << flush;
            break;
        };
        *(trAll[iW1][0])=trA;
        nWinTr[iW1]=0; //initialize number of transcripts per window


    #ifdef COMPILE_FOR_LONG_READS
        stitchWindowSeeds(iW, iW1, NULL, R[trA.roStr==0 ? 0:2]);
        if (P.pCh.segmentMin>0) {
            for (uint ia=0;ia<nWA[iW];ia++)
            {//mark all seeds that overlap the best (and only for now) transcript trA
                if (WAincl[ia]) continue;
                for (uint iex=0;iex<trA.nExons;iex++)
                {
                    if ( WA[iW][ia][WA_rStart] < (trA.exons[iex][EX_R]+trA.exons[iex][EX_L]) && \
                        (WA[iW][ia][WA_rStart]+WA[iW][ia][WA_Length]) > trA.exons[iex][EX_R] && \
                         WA[iW][ia][WA_gStart] < (trA.exons[iex][EX_G]+trA.exons[iex][EX_L]) && \
                        (WA[iW][ia][WA_gStart]+WA[iW][ia][WA_Length]) > trA.exons[iex][EX_G] )
                    {
                        WAincl[ia]=true;
                        break;
                    };
                    
                };
            };
            stitchWindowSeeds(iW, iW1, WAincl, R[trA.roStr==0 ? 0:2]);
        };
    #else
        stitchWindowAligns(0, nWA[iW], 0, WAincl, 0, 0, trA, Lread, WA[iW], R[trA.roStr==0 ? 0:2], mapGen, P, trAll[iW1], nWinTr+iW1, this);
    #endif

        if (nWinTr[iW1]==0) {
            continue;
        };

        if (trAll[iW1][0]->maxScore > trBest->maxScore || (trAll[iW1][0]->maxScore == trBest->maxScore && trAll[iW1][0]->gLength < trBest->gLength ) ) {
            trBest=trAll[iW1][0];
        };

        trNtotal += nWinTr[iW1];
        iW1++;
    };

    nW=iW1;//only count windows that had alignments

//     {//debug
//         std::time(&timeFinish);
//         double timeDiff=difftime(timeFinish,timeStart);
//         cout << "     "<< timeDiff << "     "<<trBest->maxScore*100/Lread<<"   "<<iRead<<endl;;
//     };

    if (trBest->maxScore==0) {//no window was aligned (could happen if for all windows too many reads are multiples)
        mapMarker = MARKER_NO_GOOD_WINDOW;
        nW=0;
        return;
    };

};//end of function

