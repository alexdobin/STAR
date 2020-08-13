#ifndef H_Solo
#define H_Solo
#include "IncludeDefine.h"
#include "ReadAlignChunk.h"
#include "Transcriptome.h"
#include <fstream>

#include "SoloFeature.h"


class Solo {
private:
    ReadAlignChunk **RAchunk;
    Parameters &P;
    Transcriptome &Trans;

public:
    ParametersSolo &pSolo;
    SoloFeature **soloFeat;
    
    SoloReadBarcode *readBarSum;

    Solo(ReadAlignChunk **RAchunk, Parameters &Pin, Transcriptome &inTrans);
    
    Solo(Parameters &Pin, Transcriptome &inTrans);//for soloCellFiltering

    void processAndOutput();
};

#endif
