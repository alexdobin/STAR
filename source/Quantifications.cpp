#include "Quantifications.h"

Quantifications::Quantifications (uint32 nGeIn) {

    geneCounts.nType=3;
    geneCounts.cAmbig = new uintQ[geneCounts.nType];
    geneCounts.cNone = new uintQ[geneCounts.nType];

    geneCounts.nGe=nGeIn;
    geneCounts.gCount = new uintQ* [geneCounts.nType];

    geneCounts.cMulti=0;
    for (int itype=0; itype<geneCounts.nType; itype++)
    {
        geneCounts.cAmbig[itype]=0;
        geneCounts.cNone[itype]=0;
        geneCounts.gCount[itype] = new uintQ [geneCounts.nGe];
        for (uint32 ii=0; ii<geneCounts.nGe; ii++)
        {
            geneCounts.gCount[itype][ii]=0;
        };
    };
};

void Quantifications::addQuants(const Quantifications & quantsIn)
{
    geneCounts.cMulti += quantsIn.geneCounts.cMulti;
    for (int itype=0; itype<geneCounts.nType; itype++)
    {
        geneCounts.cAmbig[itype] += quantsIn.geneCounts.cAmbig[itype];
        geneCounts.cNone[itype] += quantsIn.geneCounts.cNone[itype];
        for (uint32 ii=0; ii<geneCounts.nGe; ii++)
        {
            geneCounts.gCount[itype][ii] += quantsIn.geneCounts.gCount[itype][ii];
        };
    };
};