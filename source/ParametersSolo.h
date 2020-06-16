#ifndef CODE_ParametersSolo
#define CODE_ParametersSolo

#include <array>

#include "IncludeDefine.h"
#include "SoloBarcode.h"
#include "SoloFeatureTypes.h"

class Parameters;

class ParametersSolo {
public:
    //chemistry, library etc
    string typeStr;
    enum SoloTypes : int32 {None=0, CB_UMI_Simple=1, CB_UMI_Complex=2, CB_samTagOut=3, SmartSeq=4};
    SoloTypes type;
    bool yes;
    string strandStr;
    int32 strand;   
    
    uint32 barcodeRead;//which read is the barcode read = 0,1,2?
    uint32 barcodeStart, barcodeEnd;//start/end of barcode sequence on barcodeRead
    bool barcodeReadYes;
    
    //simple barcodes
    uint32 cbS, cbL; //cell barcode start,length
    uint32 umiS, umiL; //umi start,length
    uint32 bL, cbumiL; //total barcode sequene length, CB+UMI length. Former does may not be equal to the latter

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
    vector<string> featureIn;//string of requested features
    vector<uint32> features;
    uint32 nFeatures;//=features.size(), number of requested features
    
    array<bool,SoloFeatureTypes::N> featureYes; //which features are requested
    array<bool,SoloFeatureTypes::N> readInfoYes;//which features will readInfo (for now only Gene)
    array<int32,SoloFeatureTypes::N> featureInd;//index of each feature - skips unrequested features
    
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
    //vector <bool> umiDedupYes;
    
    struct {
        string type;
        bool mm1; //1 mismatch allowed
        bool mm1_multi; //1 mismatch, multiple matches to WL allowed
        bool oneExact; //CBs require at least one exact match
        bool mm1_multi_pc; //use psedocounts while calculating probabilities of multi-matches
    } CBmatchWL;
    
    struct {
        vector<string> type;
        bool MultiGeneUMI;
    } umiFiltering;
    
    //clusters
    string clusterCBfile;
    
    //output
    vector<string> outFileNames;    
    struct {
    	string featuresGeneField3;
    } outFormat;

    bool samAttrYes;//post-processed SAM attributes: error-corrected CB and UMI
    int32 samAttrFeature;//which feature to use for error correction
    


    //processing
    uint32 redistrReadsNfiles; //numer of files to resditributes reads into
    
    //constants
    uint32 umiMaskLow, umiMaskHigh; //low/high half bit-mask or UMIs

    void initialize(Parameters *pPin);
    void umiSwapHalves(uint32 &umi);
    void complexWLstrings();
private:
    Parameters *pP;
};
#endif
