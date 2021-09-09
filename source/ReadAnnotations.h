#ifndef H_ReadAnnotations
#define H_ReadAnnotations

#include "IncludeDefine.h"
#include "SoloCommon.h"

class ReadAnnotations {
public:
            set<uint32> geneFull, geneFull_Ex50pAS, geneFull_ExonOverIntron, geneConcordant;
            uint32 geneFullTr, geneFull_Ex50pAS_Tr, geneFull_ExonOverIntron_Tr, geneConcordantTr; //index of the annotated align - for multimappers that aligned to one gene only
            
            vector<int32> geneFull_Al, geneFull_Ex50pAS_Al, geneFull_ExonOverIntron_Al, geneConcordant_Al; //gene for each align

            vector<array<uint32,2>> transcriptConcordant;
            vector<int32> geneExonOverlap;
            array<uint32,2> geneVelocytoSimple;//first element is gene, then counts of transcript types
            vector<trTypeStruct> trVelocytoType;//first element is gene, then counts of transcript types
            
            //vector<array<uint64,2>> sj;
            //bool sjAnnot;
            
};

#endif


