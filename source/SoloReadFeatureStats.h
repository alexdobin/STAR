#ifndef H_SoloReadFeatureStats
#define H_SoloReadFeatureStats
#include "IncludeDefine.h"

class SoloReadFeatureStats {
public:
    vector<string> names;
    enum {      noUnmapped,  noNoFeature,  noTooManyWLmatches,  noNoExactWLmatch,  yesWLmatch,  yessubWLmatchExact, yesWLmatch_UniqueFeature,  yesWLmatch_MultiFeature,  yessubWLmatch_MultiFeatureMultiGenomic,  yesCellBarcodes,  yesUMIs, nStats};
    uint64 V[nStats];    
    SoloReadFeatureStats() 
    {
        names={"noUnmapped","noNoFeature","noTooManyWLmatches","noNoExactWLmatch","yesWLmatch","yessubWLmatchExact","yesWLmatch_UniqueFeature","yesWLmatch_MultiFeature","yessubWLmatch_MultiFeatureMultiGenomic","yesCellBarcodes","yesUMIs"};
        for (uint32 ii=0; ii<nStats; ii++)
            V[ii]=0;
    };
    
    uint64 numInvalidBarcodes()
    {
        return V[noTooManyWLmatches]+V[noNoExactWLmatch];
    };
    
    uint64 numMappedToTranscriptome()
    {
        return V[yesWLmatch];
    };
    
    /*
    void calcUnique(bool multYes)
    {
        if ( multYes  ) {//multi-genic reads were counted in yesWLmatch
            V[yesWLmatch_UniqueFeature] = V[yesWLmatch] - V[yesWLmatch_MultiFeature];
        } else {//no multi-genic reads
            V[yesWLmatch_UniqueFeature] =  V[yesWLmatch];
        };
    };
    */
    
    uint64 numMappedToTranscriptomeUnique()
    {
        return V[yesWLmatch_UniqueFeature];
    };
    
    double numSequencingSaturation()
    {
        return 1.0 - double(V[yesUMIs])/V[yesWLmatch_UniqueFeature]; //yesUMIs is calculated for unqiue-gene reads
    };
};

#endif