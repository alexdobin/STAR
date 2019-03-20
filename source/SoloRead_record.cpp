#include "SoloRead.h"

void SoloRead::record(string &barcodeSeq, uint64 nTr, set<uint32> &readTrGenes, Transcript *alignOut)
{
    if (pSolo.type==0)
        return;

    readBar->getCBandUMI(barcodeSeq);

    for (uint32 ii=0; ii<pSolo.nFeatures; ii++)
        readFeat[ii]->record(*readBar, nTr, readTrGenes, alignOut);
};
