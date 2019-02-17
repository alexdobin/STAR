#ifndef CODE_ParametersSolo
#define CODE_ParametersSolo
#include "IncludeDefine.h"

class Parameters;

class ParametersSolo {
public:
    //chemistry, library etc
    string typeStr;
    int type;
    uint32 cbS,cbL; //cell barcode start,length
    uint32 umiS,umiL; //umi start,length
    uint32 bL; //total barcode length
    string soloCBwhitelist;
    std::vector <uint32> cbWL;
    string strandStr;
    int32 strand;
    //features
    const static vector<string> featureNames;
    vector<string> featureIn;
    vector<uint32> features;
    uint32 nFeatures;
    bool featureYes[2]; //which features are requested
    //filtering
    char QSbase,QSmax;//quality score base and cutoff
    float cbMinP;//for CBs with non-exact matching to WL, min posterior probability
    //algorithms
    vector <string> umiDedup;
    vector <uint32> umiDedupColumns;
    vector <bool> umiDedupYes;
    //output
    vector<string> outFileNames;
    //constants
    uint32 umiMaskLow, umiMaskHigh; //low/high half bit-mask or UMIs 
    
    void initialize(Parameters *pPin);
private:
    Parameters *pP;
};
#endif
