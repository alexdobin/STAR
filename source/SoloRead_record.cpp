#include "SoloRead.h"

void SoloRead::record(uint64 nTr, set<uint32> &readGene, set<uint32> &readGeneFull, Transcript *alignOut, uint64 iRead, const vector<array<uint32,2>> &readTranscripts)
{
    if (pSolo.type==0)
        return;

    for (uint32 ii=0; ii<pSolo.nFeatures; ii++)
        readFeat[ii]->record(*readBar, nTr, readGene, readGeneFull, alignOut, iRead, readTranscripts);
};
