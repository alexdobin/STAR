#ifndef H_SoloReadFeatureStats
#define H_SoloReadFeatureStats
#include "IncludeDefine.h"

class SoloReadFeatureStats {
public:
    vector<string> names;
    enum {      noUnmapped,  noNoFeature,  MultiFeature,  subMultiFeatureMultiGenomic, noTooManyWLmatches,  noMMtoWLwithoutExact,    yesWLmatch,  yessubWLmatchExact, yessubWLmatch_UniqueFeature,  yesCellBarcodes,  yesUMIs, nStats};
    uint64 V[nStats];    
    SoloReadFeatureStats() 
    {
        names={"noUnmapped","noNoFeature","MultiFeature","subMultiFeatureMultiGenomic","noTooManyWLmatches","noMMtoWLwithoutExact","yesWLmatch","yessubWLmatchExact","yessubWLmatch_UniqueFeature","yesCellBarcodes","yesUMIs"};
        for (uint32 ii=0; ii<nStats; ii++)
            V[ii]=0;
    };
    
    uint64 numInvalidBarcodes()
    {
        return V[noTooManyWLmatches]+V[noMMtoWLwithoutExact];
    };
    
    uint64 numMappedToTranscriptome()
    {
        return V[yesWLmatch];
    };
    
    /*
    void calcUnique(bool multYes)
    {
        if ( multYes  ) {//multi-genic reads were counted in yesWLmatch
            V[yessubWLmatch_UniqueFeature] = V[yesWLmatch] - V[MultiFeature];
        } else {//no multi-genic reads
            V[yessubWLmatch_UniqueFeature] =  V[yesWLmatch];
        };
    };
    */
    
    uint64 numMappedToTranscriptomeUnique()
    {
        return V[yessubWLmatch_UniqueFeature];
    };
    
    double numSequencingSaturation()
    {
        return 1.0 - double(V[yesUMIs])/V[yessubWLmatch_UniqueFeature]; //yesUMIs is calculated for unqiue-gene reads
    };
};

#endif