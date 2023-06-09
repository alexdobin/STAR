#include "SoloRead.h"

void SoloRead::record(uint64 nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot)
{
    if (pSolo.type==pSolo.SoloTypes::None)
        return;
    if (pSolo.type==pSolo.SoloTypes::CB_samTagOut)
        return;

    if (pSolo.readStats.yes)
        readFlagReset();

    for (uint32 ii=0; ii<pSolo.nFeatures; ii++)
        readFeat[ii]->record(*readBar, nTr, alignOut, iRead, readAnnot);
};
