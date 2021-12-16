#ifndef H_SoloRead
#define H_SoloRead

#include "SoloReadBarcode.h"
#include "SoloReadFeature.h"
#include "ReadAnnotations.h"

class SoloRead {
public:
    SoloReadBarcode *readBar;
    SoloReadFeature **readFeat;
    
    SoloRead(Parameters &Pin, int32 iChunkIn);
    void readFlagReset();
    void record(uint64 nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot);
    
private:
    const int32 iChunk;
    Parameters &P;
    ParametersSolo &pSolo;
};

#endif
