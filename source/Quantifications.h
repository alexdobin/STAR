#ifndef CODE_Quantifications
#define CODE_Quantifications
#include "IncludeDefine.h"

#define uintQ unsigned long

class Quantifications {
    public:
        struct {//counting reads per gene, similar to HTseq
            uint32  nGe;      //number of genes
            int nType; //number of count types (columns)
            uintQ cMulti;     //count multimappers
            uintQ *cAmbig, *cNone;//ambigouous, no-feature
            uintQ **gCount;     // array of read counts per gene for two strands
        } geneCounts;

    Quantifications (uint32 nGeIn);

    void addQuants(const Quantifications & quantsIn); //adds quantsIn to the quants
};

#endif
