#ifndef CODE_SoloBarcode
#define CODE_SoloBarcode
#include "IncludeDefine.h"

class SoloBarcode {//complex barcodes
public:
    uint32 mate;//0 or 1
    //array: 0=start, 1=end
    int32 anchorType[2];//0=read start; 1=read end; 2=adapter start; 3=adapter end
    int32 anchorDist[2];//distance to anchor
    
    uint32 adapterLength;//length of the adapter
    
    vector<vector<uint64>> wl;//whitelists, one for each length
    uint64 wlFactor;//factor and modulo for converting each whitelist index into global index
    vector<uint32> wlAdd;//additive for each length
    uint32 minLen;//minimum length for this barcode
    uint32 totalSize;//total size of all whitelists
    
    uint64 ilen, icb;//indexes
    
    bool extractBarcode(string &seqIn, string &qualIn, const uint32 aStart, string &bSeq, string &bQual);
    void sortWhiteList();
    void extractPositionsFromString(string &strIn);
};

#endif