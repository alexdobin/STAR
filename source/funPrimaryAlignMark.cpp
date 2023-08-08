#include "funPrimaryAlignMark.h"

void funPrimaryAlignMark(Transcript **trMult, uint64 nTr, 
                         Parameters &P, int maxScore, std::uniform_real_distribution<double> rngUniformReal0to1, std::mt19937 rngMultOrder)
{
    if (nTr==1){//unique mappers
        trMult[0]->primaryFlag=true;
    } else {//multimappers
        int nbest=0;
        if (P.outMultimapperOrder.random || P.outSAMmultNmax != (uint) -1 ) {//bring the best alignment to the top of the list. TODO sort alignments by the score?
            for (uint itr=0; itr<nTr; itr++) {//move the best aligns to the top of the list
                if ( trMult[itr]->maxScore == maxScore ) {
                    swap(trMult[itr],trMult[nbest]);
                    ++nbest;
                };
            };
        };

        if (P.outMultimapperOrder.random) {//shuffle separately the best aligns, and the rest
            for (int itr=nbest-1; itr>=1; itr--) {//Fisher-Yates-Durstenfeld-Knuth shuffle
                int rand1=int (rngUniformReal0to1(rngMultOrder)*itr+0.5);
                swap(trMult[itr],trMult[rand1]);
            };
            for (int itr=nTr-nbest-1; itr>=1; itr--) {//Fisher-Yates-Durstenfeld-Knuth shuffle
                int rand1=int (rngUniformReal0to1(rngMultOrder)*itr+0.5);
                swap(trMult[nbest+itr],trMult[nbest+rand1]);
            };
        };

        if ( P.outSAMprimaryFlag=="AllBestScore" ) {
            for (uint itr=0; itr<nTr; itr++)
            {//mark all best score aligns as primary
                if ( trMult[itr]->maxScore == maxScore ) trMult[itr]->primaryFlag=true;
            };
        } else if (P.outMultimapperOrder.random || P.outSAMmultNmax != (uint) -1) {
            trMult[0]->primaryFlag=true;//mark as primary the first one in the random ordered list: best scoring aligns are already in front of the list
        } else {//old way
            //trBest->primaryFlag=true; //cannot do it, trBest may not be defined
            trMult[0]->primaryFlag=true;
        };
    };
};

