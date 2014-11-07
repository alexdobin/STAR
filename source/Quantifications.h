#ifndef CODE_Quantifications
#define CODE_Quantifications
#include "IncludeDefine.h"

#define uintQ unsigned long

class Quantifications {
    public:
        struct {//counting reads per gene, similar to HTseq
            uint32  nGe;      //number of genes    
            uintQ cMulti,cAmbig,cNone;     //count multimappers, ambigouous, no-feature
            uintQ **uStr;     // array of read counts per gene for two strands
        } geneCounts;

    Quantifications (uint32 nGeIn);
};

#endif
