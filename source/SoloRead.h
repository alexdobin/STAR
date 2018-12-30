#ifndef H_SoloRead
#define H_SoloRead

#include "SoloReadBarcode.h"
#include "SoloReadFeature.h"

class SoloRead {
public:
    SoloReadBarcode *readBar;
    SoloReadFeature **readFeat;
    
    SoloRead(Parameters &Pin, int32 iChunkIn);
    void record(string &barcodeSeq, uint64 nTr, set<uint32> &readTrGenes, Transcript *alignOut);
    
private:
    const int32 iChunk;
    Parameters &P;
    ParametersSolo &pSolo;
};

#endif
