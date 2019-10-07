#ifndef H_SoloReadBarcodeStats
#define H_SoloReadBarcodeStats
#include "IncludeDefine.h"

class SoloReadBarcodeStats {
public:
    vector<string> names;
    enum {      nNoAdapter,  nNoUMI,    nNoCB,   nNinCB,   nNinUMI,   nUMIhomopolymer,  nTooMany,  nNoMatch,   nMismatchesInMultCB,  nExactMatch,    nMismatchOneWL,   nMismatchToMultWL, nStats};
    uint64 V[nStats];    
    SoloReadBarcodeStats() 
    {
        names={"nNoAdapter", "nNoUMI", "nNoCB", "nNinCB", "nNinUMI", "nUMIhomopolymer","nTooMany","nNoMatch", "nMismatchesInMultCB", "nExactMatch", "nMismatchOneWL", "nMismatchToMultWL"};
        for (uint32 ii=0; ii<nStats; ii++)
            V[ii]=0;
    };
    
    uint64 numInvalidBarcodes()
    {
        uint64 n=0;
        for (uint32 ii=0; ii<9; ii++) //first eight number record all invalid barcodes
            n += V[ii];
        
        return n;
    };
};

#endif