#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ReadAlign.h"
#include "ErrorWarning.h"
#include <random>

void ReadAlign::multMapSelect() {//select multiple mappers from all transcripts of all windows

    nTr=0;
    if (nW==0) {//no good windows
        return;
    };

    maxScore=-10*Lread;
    for (uint iW=0; iW<nW; iW++) {//scan windows
        if (maxScore < trAll[iW][0]->maxScore) maxScore = trAll[iW][0]->maxScore;
    };

    if (maxScore!=trBest->maxScore) {
        ostringstream errOut;
        errOut  << "BUG: maxScore!=trBest->maxScore in multMapSelect";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);
    };

    for (uint iW=0; iW<nW; iW++) {//scan windows
        for (uint iTr=0; iTr<nWinTr[iW]; iTr++) {//scan transcripts
            if ( (trAll[iW][iTr]->maxScore + P.outFilterMultimapScoreRange) >= maxScore  ) {//record this alignment
                // if paired-end, record alignments from ALL windows
                if (nTr==MAX_N_MULTMAP) {//too many alignments for this read, do not record it
                    ostringstream errOut;
                    errOut  << "EXITING: Fatal ERROR: number of alignments exceeds MAX_N_MULTMAP, increase it and re-compile STAR";
                    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
                };

                trMult[nTr]=trAll[iW][iTr];
                trMult[nTr]->Chr = trAll[iW][0]->Chr;
                trMult[nTr]->Str = trAll[iW][0]->Str;
                trMult[nTr]->roStr = trAll[iW][0]->roStr;

                if ( (trAll[iW][iTr]->maxScore + P.outFilterMultimapScoreRange) >= maxScore) nTrMate++;

                nTr++;
            };
        };
    };

    if (nTr > P.outFilterMultimapNmax || nTr==0)
    {//too multi OR no alignments, no need for further processing, since it will be considered unmapped
        return;
    };

    for (uint iTr=0; iTr<nTr; iTr++)
    {
        trMult[iTr]->roStart = trMult[iTr]->roStr==0 ? trMult[iTr]->rStart : Lread - trMult[iTr]->rStart - trMult[iTr]->rLength;
        trMult[iTr]->cStart=trMult[iTr]->gStart - mapGen.chrStart[trMult[iTr]->Chr];
    };

//     if (P.outMultimapperOrder.sortCoord)
//     {//sort multimappers by coordinate
//         uint *s=new uint[nTr*2];
//         Transcript **t=new Transcript*[nTr];
//         for (uint itr=0; itr<nTr; itr++)
//         {//fill the array of starts and pointers
//             s[itr*2]=trMult[itr]->exons[0][EX_G];
//             s[itr*2+1]=itr;
//             t[itr]=trMult[itr];
//         };
//         qsort((void*) s, nTr, sizeof(uint)*2, funCompareUint1);
//         for (uint itr=0; itr<nTr; itr++)
//         {
//             trMult[itr]=t[s[itr*2+1]];
//         };
//         delete [] s;
//     };

    if (nTr==1)
    {//unique mappers
        trMult[0]->primaryFlag=true;
    } else
    {//multimappers
        int nbest=0;
        if (P.outMultimapperOrder.random || P.outSAMmultNmax != (uint) -1 )
        {//bring the best alignment to the top of the list. TODO sort alignments by the score?
            for (uint itr=0; itr<nTr; itr++)
            {//move the best aligns to the top of the list
                if ( trMult[itr]->maxScore == maxScore )
                {
                    swap(trMult[itr],trMult[nbest]);
                    ++nbest;
                };
            };
        };
        if (P.outMultimapperOrder.random)
        {//shuffle separately the best aligns, and the rest
            for (int itr=nbest-1; itr>=1; itr--)
            {//Fisher-Yates-Durstenfeld-Knuth shuffle
                int rand1=int (rngUniformReal0to1(rngMultOrder)*itr+0.5);
                swap(trMult[itr],trMult[rand1]);
            };
            for (int itr=nTr-nbest-1; itr>=1; itr--)
            {//Fisher-Yates-Durstenfeld-Knuth shuffle
                int rand1=int (rngUniformReal0to1(rngMultOrder)*itr+0.5);
                swap(trMult[nbest+itr],trMult[nbest+rand1]);
            };
        };

        if ( P.outSAMprimaryFlag=="AllBestScore" )
        {
            for (uint itr=0; itr<nTr; itr++)
            {//mark all best score aligns as primary
                if ( trMult[itr]->maxScore == maxScore ) trMult[itr]->primaryFlag=true;
            };
        } else if (P.outMultimapperOrder.random || P.outSAMmultNmax != (uint) -1)
        {
            trMult[0]->primaryFlag=true;//mark as primary the first one in the random ordered list: best scoring aligns are already in front of the list
    //         for (uint itr=0; itr<nTr; itr++)
    //         {
    //             if ( trMult[itr]->maxScore == maxScore ) trMult[itr]->primaryFlag=true;
    //             break;
    //         };
        } else
        {//old way
            trBest->primaryFlag=true;
        };
    };
};

