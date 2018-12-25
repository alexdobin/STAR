#ifndef CODE_ParametersSolo
#define CODE_ParametersSolo
#include "IncludeDefine.h"

class Parameters;

class ParametersSolo {
public:
    //chemistry, library etc
    string typeStr;
    int type;
    uint32 cbL; //cell barcode length
    uint32 umiL; //umi length
    uint32 bL; //total barcode length
    string soloCBwhitelist;
    std::vector <uint32> cbWL;
    string strandStr;
    int32 strand;    

    //features
    const static vector<string> featureNames;
    vector<string> featureIn;
    bool featureYes[2]; //which features are requested
    //filtering
    char QSbase,QSmax;//quality score base and cutoff
    float cbMinP;//for CBs with non-exact matching to WL, min posterior probability 
    //output
    vector<string> outFileNames;


    
    //constants
    uint32 umiMaskLow, umiMaskHigh; //low/high half bit-mask or UMIs 
    
    //static const uint64 bufferSize=1048576;
    
    void initialize(Parameters *pPin);
    
private:
    Parameters *pP;
};
#endif
