#ifndef H_Solo
#define H_Solo
#include "IncludeDefine.h"
#include "ReadAlignChunk.h"
#include "SoloRead.h"
#include "Transcriptome.h"
#include <fstream>

class SoloFeature {   
public:

    SoloFeature *soloFeat;
  

private:   
    
    Parameters &P;
    ParametersSolo &pSolo;
    Transcriptome &Trans;   
};

#endif
