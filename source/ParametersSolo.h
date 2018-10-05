#ifndef CODE_ParametersSolo
#define CODE_ParametersSolo

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
    
    const uint64 bufferSize=1048576;
    initialize(Parameters *pPin);
    
private:
    Parameters *pP;
};

#endif