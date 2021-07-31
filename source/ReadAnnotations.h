#ifndef H_ReadAnnotations
#define H_ReadAnnotations

#include "IncludeDefine.h"
#include "SoloCommon.h"

class ReadAnnotations {
public:
            set<uint32> geneFull, geneFull_CR, geneFull_ExonOverIntron, geneConcordant;
            uint32 geneFullTr, geneFull_CR_Tr, geneFull_ExonOverIntron_Tr, geneConcordantTr; //index of the annotated align - for multimappers that aligned to one gene only
            
            vector<array<uint32,2>> transcriptConcordant;
            vector<int32> geneExonOverlap;
            array<uint32,2> geneVelocytoSimple;//first element is gene, then counts of transcript types
            vector<trTypeStruct> trVelocytoType;//first element is gene, then counts of transcript types
            
            //vector<array<uint64,2>> sj;
            //bool sjAnnot;
            
};

#endif


