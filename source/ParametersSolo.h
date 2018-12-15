#ifndef CODE_ParametersSolo
#define CODE_ParametersSolo

#include <unordered_set>
#include <set>

#include "IncludeDefine.h"
class Parameters;

class ParametersSolo {
public:
    string typeStr;
    int type;
    uint32 cbL; //cell barcode length
    uint32 umiL; //umi length
    uint32 bL; //total barcode length
    string soloCBwhitelist;
    
    string strandStr;
    int32 strand;
    
    std::vector <uint32> cbWL;

    uint32 umiMaskLow, umiMaskHigh; //low/high half bit-mask or UMIs 
    
    static const uint64 bufferSize=1048576;
    void initialize(Parameters *pPin);
    
private:
    Parameters *pP;
};

#endif