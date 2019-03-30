#include "SoloRead.h"

void SoloRead::record(uint64 nTr, set<uint32> &readGene, set<uint32> &readGeneFull, Transcript *alignOut)
{
    if (pSolo.type==0)
        return;

    for (uint32 ii=0; ii<pSolo.nFeatures; ii++)
        readFeat[ii]->record(*readBar, nTr, readGene, readGeneFull, alignOut);
};
