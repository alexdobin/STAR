#ifndef H_ReadAnnotations
#define H_ReadAnnotations

#include "IncludeDefine.h"

class ReadAnnotations {
public:
            set<uint32> geneFull, geneConcordant;
            vector<array<uint32,2>> transcriptConcordant;
            vector<int32> geneExonOverlap;
            array<uint32, 3> geneVelocyto;//first element is gene, then counts of transcript types
};

#endif


