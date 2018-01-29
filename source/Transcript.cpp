#include "Transcript.h"

Transcript::Transcript()
{
    reset();
};

void Transcript::reset() {
    extendL=0;

//     for (uint ii=0;ii<4;ii++) {
//         polyXlength[ii]=0;
//         polyXnMM[ii]=0;
//     };
    primaryFlag=false;

    rStart=0; roStart=0; rLength=0; gStart=0; gLength=0; //read and genomic coordinates

    maxScore=0;
    nMatch=0;
    nMM=0;

    nGap=0; lGap=0; lDel=0; lIns=0; nDel=0; nIns=0;

    nUnique=nAnchor=0;
};

void Transcript::add(Transcript *trIn) {
    maxScore+=trIn->maxScore;
    nMatch+=trIn->nMatch;
    nMM+=trIn->nMM;
    nGap+=trIn->nGap; lGap+=trIn->lGap;
    lDel+=trIn->lDel; nDel+=trIn->nDel;
    lIns+=trIn->lIns; nIns+=trIn->nIns;
    nUnique+=trIn->nUnique;
};


