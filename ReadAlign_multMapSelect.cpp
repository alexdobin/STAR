#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ReadAlign.h"
#include "ErrorWarning.h"

void ReadAlign::multMapSelect() {//select multiple mappers from all transcripts of all windows
    
    maxScore=0;
    for (uint iW=0; iW<nW; iW++) {//scan windows
        if (maxScore < trAll[iW][0]->maxScore) maxScore = trAll[iW][0]->maxScore;
    };
    
    if (maxScore!=trBest->maxScore) {
        ostringstream errOut;
        errOut  << "BUG: maxScore!=trBest->maxScore in multMapSelect";
        exitWithError(errOut.str(), std::cerr, P->inOut->logMain, EXIT_CODE_BUG, *P);                    
    };
    
    trBest->primaryFlag=true;
    
    bool chimRecord = false;

    nTr=0;
    for (uint iW=0; iW<nW; iW++) {//scan windows
        for (uint iTr=0; iTr<nWinTr[iW]; iTr++) {//scan transcripts
            if ( ( (trAll[iW][iTr]->maxScore + P->outFilterMultimapScoreRange) >= maxScore ) || \
                 ( chimRecord && trAll[iW][iTr]->iFrag>=0 && (trAll[iW][iTr]->maxScore + P->outFilterMultimapScoreRange) >=maxScoreMate[trAll[iW][iTr]->iFrag] )  ) {//record this alignment
                // if paired-end, record alignments from ALL windows
                if (nTr==MAX_N_MULTMAP) {//too many alignments for this read, do not record it
                    ostringstream errOut;
                    errOut  << "EXITING: Fatal ERROR: number of alignments exceeds MAX_N_MULTMAP, increase it and re-compile STAR";
                    exitWithError(errOut.str(), std::cerr, P->inOut->logMain, EXIT_CODE_PARAMETER, *P);                    
                };
                    
                trMultScores[nTr]=trAll[iW][iTr]->maxScore;
                trMult[nTr]=trAll[iW][iTr];
                trMult[nTr]->Chr = trAll[iW][0]->Chr;
                trMult[nTr]->Str = trAll[iW][0]->Str;                  
                trMult[nTr]->roStr = trAll[iW][0]->roStr;                                 
                if ( P->outSAMprimaryFlag=="AllBestScore" && trMultScores[nTr] == maxScore ) trMult[nTr]->primaryFlag=true; 
                
                if ( (trAll[iW][iTr]->maxScore + P->outFilterMultimapScoreRange) >= maxScore) nTrMate++;
                
                nTr++;           
            };
        };
    };
    
    for (uint iTr=0; iTr<nTr; iTr++) {//scan transcripts                 
        trMult[iTr]->roStart = trMult[iTr]->roStr==0 ? trMult[iTr]->rStart : Lread - trMult[iTr]->rStart - trMult[iTr]->rLength;
        trMult[iTr]->cStart=trMult[iTr]->gStart - P->chrStart[trMult[iTr]->Chr];                        
    };

};

