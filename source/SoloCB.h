#ifndef CODE_SoloCB
#define CODE_SoloCB
#include <set>
#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"

class SoloCB {
public:
    
    uint32 homoPolymer[4];//homopolymer constants
    
    uint32 *cbReadCount, *cbReadCountExact;
    
    fstream *strU_0 ,*strU_1, *strU_2; //unique mappers, CB matches whitelist with 0,1>=2 MM

    struct {
        enum {                 nUnmapped,  nNoFeature,  nAmbigFeature,  nAmbigFeatureMultimap,  nNinBarcode,  nUMIhomopolymer,  nTooMany,  nNoExactMatch,  nNoMatch,  nExactMatch,  nMatch,  nCellBarcodes,  nUMIs, nStats};
        uint64 V[nStats];
        vector<string> names={"nUnmapped","nNoFeature","nAmbigFeature","nAmbigFeatureMultimap","nNinBarcode","nUMIhomopolymer","nTooMany","nNoExactMatch","nNoMatch","nExactMatch","nMatch","nCellBarcodes","nUMIs",};
    } stats;
        
    string cbSeq, umiSeq, cbQual, umiQual;
    
    SoloCB (int32 feTy, Parameters &Pin, int iChunk);
    bool  readCB(uint64 &iReadAll, string &readNameExtra, uint nTr, set<uint32> &readTrGenes, Transcript *alignOut, bool readRecord);
    void addSoloCBcounts(const SoloCB &soloCBin);
    void addSoloCBstats(const SoloCB &soloCBin);
    void statsOut(ofstream &streamOut);
    void inputUMIfeatureCBrecords(uint32 **cbP, uint32 *cbReadCountExact);
    
private:
    const int32 featureType;
    
    Parameters &P;
    ParametersSolo &pSolo;
};

#endif
