#ifndef CODE_ParametersSolo
#define CODE_ParametersSolo

#include <array>
#include <unordered_map>
#include <mutex>


#include "IncludeDefine.h"
#include "SoloBarcode.h"
#include "SoloFeatureTypes.h"

class Parameters;
class ParametersSolo;

class UMIdedup {
public:
    const static uint32 tN = 6;
    array<string,tN> typeNames { {"NoDedup", "Exact", "1MM_All", "1MM_Directional", "1MM_CR", "1MM_Directional_UMItools"} };
    enum typeI : int32 { NoDedup=0, Exact=1, All=2, Directional=3, CR=4, Directional_UMItools=5 };
    
    struct {
        uint32_t N;
        array<bool,tN> B;
        bool &NoDedup=B[0], &Exact=B[1], &All=B[2], &Directional=B[3], &CR=B[4], &Directional_UMItools=B[5]; 
    } yes;

    struct {
        //uint32_t N;
        array<uint32_t,tN> I;
        uint32_t &NoDedup=I[0], &Exact=I[1], &All=I[2], &Directional=I[3], &CR=I[4], &Directional_UMItools=I[5];
        uint32_t main; //index for SAM/stats/filtering output
    } countInd; //index in the countCellGennUMI
    
    vector<string> typesIn; //UMIdedup types from user options
    vector<int32> types; //the above converted to typeI numbers
    int32 typeMain; //the type to be used in SAM/stats/filtering output - for now just types[0]
    
    void initialize(ParametersSolo *pS);
    
//protected:
//    int it;
};

class MultiMappers {
public:
    const static uint32 tN = 5;
    array<string,tN> typeNames { {"Unique", "Uniform", "Rescue", "PropUnique", "EM"} };
    enum typeI : int32 { Unique=0, Uniform=1, Rescue=2, PropUnique=3, EM=4 };
    
    struct {
        bool multi; //if multimappers are requested
        uint32_t N;
        array<bool,tN> B;
        bool &Unique=B[0], &Uniform=B[1], &Rescue=B[2], &PropUnique=B[3], &EM=B[4] ;
    } yes;

    struct {
        //uint32_t N;
        array<uint32_t,tN> I;
        uint32_t &Unique=I[0], &Uniform=I[1], &Rescue=I[2], &PropUnique=I[3], &EM=I[4];
        uint32_t main; //index for SAM/stats/filtering output
    } countInd; //index in the countCellGennUMI
    
    vector<string> typesIn; //UMIdedup types from user options
    vector<int32> types; //the above converted to typeI numbers
    int32 typeMain; //the type to be used in SAM/stats/filtering output - for now just types[0]
    
    void initialize(ParametersSolo *pS);
};

class ParametersSolo {
public:
    Parameters *pP;
    bool yes;

    //chemistry, library etc
    string typeStr;
    enum SoloTypes : int32 {None=0, CB_UMI_Simple=1, CB_UMI_Complex=2, CB_samTagOut=3, SmartSeq=4};
    SoloTypes type;
    string strandStr;
    int32 strand;   
    
    uint32 barcodeRead, barcodeReadIn;//which read is the barcode read = 0,1,2?
    uint32 barcodeStart, barcodeEnd;//start/end of barcode sequence on barcodeRead
    bool barcodeReadSeparate;
    
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
    
    //input from SAM files
    vector<string> samAtrrBarcodeSeq, samAtrrBarcodeQual;
    
    //whitelist - general
    uint64 cbWLsize;
    bool cbWLyes;
    vector<string> soloCBwhitelist;
    vector <uint64> cbWL;    
    vector<string> cbWLstr;
    
    MultiMappers multiMap;
    
    //features
    vector<string> featureIn;//string of requested features
    vector<uint32> features;
    uint32 nFeatures;//=features.size(), number of requested features
    
    int32 featureFirst; //which feature is the first on the list
    array<bool,SoloFeatureTypes::N> featureYes; //which features are requested
    array<bool,SoloFeatureTypes::N> readInfoYes;//which features will need readInfo (for now only Gene and GeneFull)
    array<bool,SoloFeatureTypes::N> readIndexYes;//which features will need recording of readIndex (for now only Gene and GeneFull, for multimappers)
    array<bool,SoloFeatureTypes::N> readStatsYes;//which features will need output of read statistics
    array<int32,SoloFeatureTypes::N> featureInd;//index of each feature - skips unrequested features
    
    //filtering
    char QSbase,QSmax;//quality score base and cutoff

    #ifdef MATCH_CellRanger
    double cbMinP;//for CBs with non-exact matching to WL, min posterior probability
    #else
    float cbMinP;//for CBs with non-exact matching to WL, min posterior probability
    #endif
    
    //cell filtering
    struct {
        vector<string> type;
        uint32 topCells;
        
        struct {
            double nExpectedCells;
            double maxPercentile;
            double maxMinRatio;
        } knee;
        
        struct {
            uint32 indMin, indMax; //min/max cell index, sorted by UMI counts,for empty cells
            uint32 umiMin;
            double umiMinFracMedian;
            uint32 candMaxN;
            double FDR;
            uint32 simN;
        } eDcr;//EmptyDrops-CellRanger
        
    } cellFilter;
      
    //CBtype
    struct {
        string typeString;
        int32 type;
        std::unordered_map<string,uint32> strMap;
        std::mutex *strMtx;
    } CBtype;

    //CB match
    struct {
        string type;
        bool mm1; //1 mismatch allowed
        bool mm1_multi; //1 mismatch, multiple matches to WL allowed
        bool oneExact; //CBs require at least one exact match
        bool mm1_multi_pc; //use psedocounts while calculating probabilities of multi-matches
        bool mm1_multi_Nbase; //allow multimatching to WL for CBs with N-bases
        bool EditDist_2; //allow EditDistance <=3
    } CBmatchWL;
    
    //UMIdedup
    UMIdedup umiDedup;
    
    //multi-gene umi
    struct {
        vector<string> type;
        bool MultiGeneUMI       = false;
        bool MultiGeneUMI_All   = false;
        bool yes                = false; //true for non-CR
        bool MultiGeneUMI_CR    = false;
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
    uint32 redistrReadsNfiles; //numer of files to resditribute reads into
    
    //constants
    uint32 umiMaskLow, umiMaskHigh; //low/high half bit-mask or UMIs

    //readStats
    struct {
        string type; //input parameter
        bool yes=false;
    } readStats;

    void initialize(Parameters *pPin);
    void umiSwapHalves(uint32 &umi);
    void complexWLstrings();
    void cellFiltering();

    void init_CBmatchWL();
};
#endif
