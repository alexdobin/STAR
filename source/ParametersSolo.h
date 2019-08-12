#ifndef CODE_ParametersSolo
#define CODE_ParametersSolo
#include "IncludeDefine.h"
#include "SoloBarcode.h"

class Parameters;

class ParametersSolo {
public:
    //chemistry, library etc
    string typeStr;
    int type;
    bool yes;
    string strandStr;
    int32 strand;   
    
    //simple barcodes
    uint32 cbS,cbL; //cell barcode start,length
    uint32 umiS,umiL; //umi start,length
    uint32 bL; //total barcode length

    vector<string> cbPositionStr;
    string umiPositionStr;
    
    //complex barcodes    
    vector<SoloBarcode> cbV;
    SoloBarcode umiV; //single UMI
    bool adapterYes; //anchor?  
    string adapterSeq; //anchor sequence
    uint32 adapterMismatchesNmax;//max number of mismatches in the anchor
    
    //whitelist - general
    uint64 cbWLsize;
    bool cbWLyes;
    vector<string> soloCBwhitelist;
    vector <uint64> cbWL;    
    vector<string> cbWLstr;
    
    //features
    struct { enum{Gene,GeneFull,SJ,Transcript3p}; } featureTypeInd;
    const static vector<string> featureNames;
    vector<string> featureIn;
    vector<uint32> features, featureInd;
    uint32 nFeatures;
    bool *featureYes; //which features are requested
    
    //filtering
    char QSbase,QSmax;//quality score base and cutoff
    float cbMinP;//for CBs with non-exact matching to WL, min posterior probability
    
    //cell filtering
    struct {
        vector<string> type;
        double cr2maxPercentile;
        double cr2expectedCells;
        double cr2maxMinRatio;
        uint64 cr2maxCellInd;
        uint64 topCells;
    } cellFilter;
    
    //algorithms
    vector <string> umiDedup;
    vector <uint32> umiDedupColumns;
    vector <bool> umiDedupYes;
    int32 CBmatchWLtype;
    
    //output
    vector<string> outFileNames;
    
    //constants
    uint32 umiMaskLow, umiMaskHigh; //low/high half bit-mask or UMIs
    
    bool samAttrYes;//post-processed SAM attributes: error-corrected CB and UMI
    int32 samAttrFeature;//which feature to use for error correction

    void initialize(Parameters *pPin);
    void umiSwapHalves(uint32 &umi);
    void complexWLstrings();
private:
    Parameters *pP;
};
#endif
