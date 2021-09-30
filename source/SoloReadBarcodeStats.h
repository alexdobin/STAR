#ifndef H_SoloReadBarcodeStats
#define H_SoloReadBarcodeStats
#include "IncludeDefine.h"

class SoloReadBarcodeStats {
public:
    vector<string> names;
    enum {      noNoAdapter,  noNoUMI,    noNoCB,   noNinCB,   noNinUMI,   noUMIhomopolymer,  noNoWLmatch,   noTooManyMM,   noTooManyWLmatches,   yesWLmatchExact,   yesOneWLmatchWithMM,   yesMultWLmatchWithMM, nStats};
    uint64 V[nStats];    
    SoloReadBarcodeStats() 
    {
        names={"noNoAdapter", "noNoUMI", "noNoCB", "noNinCB", "noNinUMI", "noUMIhomopolymer","noNoWLmatch", "noTooManyMM", "noTooManyWLmatches", "yesWLmatchExact", "yesOneWLmatchWithMM", "yesMultWLmatchWithMM"};
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