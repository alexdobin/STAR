#ifndef H_SoloReadFeatureStats
#define H_SoloReadFeatureStats
#include "IncludeDefine.h"

class SoloReadFeatureStats {
public:
    vector<string> names;
    enum {      nUnmapped,  nNoFeature,  nAmbigFeature,  nAmbigFeatureMultimap,  nTooMany,  nNoExactMatch,  nExactMatch,  nMatch, nMatchUnique, nCellBarcodes,  nUMIs, nStats};
    uint64 V[nStats];    
    SoloReadFeatureStats() 
    {
        names={"nUnmapped","nNoFeature","nAmbigFeature","nAmbigFeatureMultimap","nTooMany","nNoExactMatch","nExactMatch","nMatch","nMatchUnique","nCellBarcodes","nUMIs"};
        for (uint32 ii=0; ii<nStats; ii++)
            V[ii]=0;
    };
    
    uint64 numInvalidBarcodes()
    {
        return V[nTooMany]+V[nNoExactMatch];
    };
    
    uint64 numMappedToTranscriptome()
    {
        return V[nMatch];
    };
    
    /*
    void calcUnique(bool multYes)
    {
        if ( multYes  ) {//multi-genic reads were counted in nMatch
            V[nMatchUnique] = V[nMatch] - V[nAmbigFeature];
        } else {//no multi-genic reads
            V[nMatchUnique] =  V[nMatch];
        };
    };
    */
    
    uint64 numMappedToTranscriptomeUnique()
    {
        return V[nMatchUnique];
    };
    
    double numSequencingSaturation()
    {
        return 1.0 - double(V[nUMIs])/V[nMatchUnique]; //nUMIs is calculated for unqiue-gene reads
    };
};

#endif