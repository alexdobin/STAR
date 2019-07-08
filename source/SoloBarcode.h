#ifndef CODE_SoloBarcode
#define CODE_SoloBarcode
#include "IncludeDefine.h"

class SoloBarcode {//complex barcodes
public:
    uint32 length;//0=variable length
    uint32 mate;//0 or 1
    int32 posType;//0=start, 1=end
    int32 posAnchor;//0=read start; 1=read end; 2=adapter start; 3=adapter end
    int32 pos;//position with respect to anchor
    uint32 adapterLength;//length of the adapter
    std::vector <uint64> whiteList;//whitelist
    
    bool extractBarcode(string &seqIn, string &qualIn, const uint32 aStart, string &bSeq, string &bQual);
};

#endif