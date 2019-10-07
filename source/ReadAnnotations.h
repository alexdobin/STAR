#ifndef H_ReadAnnotations
#define H_ReadAnnotations

#include "IncludeDefine.h"
#include "SoloCommon.h"

class ReadAnnotations {
public:
            set<uint32> geneFull, geneConcordant;
            vector<array<uint32,2>> transcriptConcordant;
            vector<int32> geneExonOverlap;
            array<uint32,2> geneVelocytoSimple;//first element is gene, then counts of transcript types
            vector<trTypeStruct> trVelocytoType;//first element is gene, then counts of transcript types
            
};

#endif


