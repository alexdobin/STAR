#ifndef H_SoloReadFeature
#define H_SoloReadFeature
#include <set>
#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "SoloReadBarcode.h"

class SoloReadFeature {
public:
    
    uint32 homoPolymer[4];//homopolymer constants
    
    uint32 *cbReadCount;
    
    fstream *strU_0 ,*strU_1, *strU_2; //unique mappers, CB matches whitelist with 0,1>=2 MM

    struct {
        enum {                 nUnmapped,  nNoFeature,  nAmbigFeature,  nAmbigFeatureMultimap,  nTooMany,  nNoExactMatch,  nExactMatch,  nMatch,  nCellBarcodes,  nUMIs, nStats};
        uint64 V[nStats];
        vector<string> names={"nUnmapped","nNoFeature","nAmbigFeature","nAmbigFeatureMultimap","nTooMany","nNoExactMatch","nExactMatch","nMatch","nCellBarcodes","nUMIs",};
    } stats;
        
    string cbSeq, umiSeq, cbQual, umiQual;
    
    SoloReadFeature (int32 feTy, Parameters &Pin, int iChunk);
    void record(SoloReadBarcode &soloBar, uint nTr, set<uint32> &readTrGenes, Transcript *alignOut);
    void addCounts(const SoloReadFeature &soloCBin);
    void addStats(const SoloReadFeature &soloCBin);
    void statsOut(ofstream &streamOut);
    void inputRecords(uint32 **cbP, uint32 *cbReadCountExact);
    
private:
    const int32 featureType;
    
    Parameters &P;
    ParametersSolo &pSolo;
};

#endif
