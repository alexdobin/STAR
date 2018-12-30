#ifndef H_Solo
#define H_Solo
#include "IncludeDefine.h"
#include "ReadAlignChunk.h"
#include "Transcriptome.h"
#include <fstream>

#include "SoloFeature.h"


class Solo {   
public:

    SoloFeature **soloFeat;
    
    Solo(ReadAlignChunk **RAchunk, Parameters &Pin, Transcriptome &inTrans);
    void processAndOutput();
    
private:   
    ReadAlignChunk **RAchunk;
    Parameters &P;
    ParametersSolo &pSolo;
    Transcriptome &Trans;
};

#endif
