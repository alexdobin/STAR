#ifndef H_SoloReadFeatureStats
#define H_SoloReadFeatureStats
#include "IncludeDefine.h"

class SoloReadFeatureStats {
public:
    vector<string> names;
    enum {      nUnmapped,  nNoFeature,  nAmbigFeature,  nAmbigFeatureMultimap,  nTooMany,  nNoExactMatch,  nExactMatch,  nMatch,  nCellBarcodes,  nUMIs, nStats};
    uint64 V[nStats];    
    SoloReadFeatureStats() 
    {
        names={"nUnmapped","nNoFeature","nAmbigFeature","nAmbigFeatureMultimap","nTooMany","nNoExactMatch","nExactMatch","nMatch","nCellBarcodes","nUMIs"};
        for (uint32 ii=0; ii<nStats; ii++)
            V[ii]=0;
    };
    
    uint64 numInvalidBarcodes()
    {
        return V[nTooMany]+V[nNoExactMatch];
    };
    
    uint64 numMappedToTranscriptome()
    {
        return V[nAmbigFeature]+V[nMatch];
    };
    
    uint64 numMappedToTranscriptomeUnique()
    {
        return V[nMatch];
    };
    
    double numSequencingSaturation()
    {
        return 1.0 - double(V[nUMIs])/V[nMatch];
    };
};

#endif