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
    
    vector<vector<uint64>> wl;//whitelists, one for each length
    uint64 wlFactor, wlModulo;//factor and modulo for converting each whitelist index into global index
    vector<uint32> wlAdd;//additive for each length
    uint32 minLen;//minimum length for this barcode
    uint32 totalSize;//total size of all whitelists
    
    bool extractBarcode(string &seqIn, string &qualIn, const uint32 aStart, string &bSeq, string &bQual);
    void sortWhiteList();
};

#endif