#ifndef CODE_SoloCB
#define CODE_SoloCB
#include "IncludeDefine.h"
#include "Parameters.h"

class SoloCB {
public:
    
    uint32 homoPolymer[4];//homopolymer constants
    
    uint32 *cbReadCount, *cbReadCountExact;
    
    fstream *strU_0 ,*strU_1, *strU_2; //uniqe mappers, CB matches whitelist with 0,1>=2 MM

    struct {
        enum {                 nNoGene,  nAmbigGene,  nAmbigGeneMultimap,  nNinBarcode,  nUMIhomopolymer,  nTooMany,  nNoExactMatch,  nNoMatch,  nExactMatch,  nMatch,  nCellBarcodes,  nUMIs, nStats};
        uint64 V[nStats];
        vector<string> names={"nNoGene","nAmbigGene","nAmbigGeneMultimap","nNinBarcode","nUMIhomopolymer","nTooMany","nNoExactMatch","nNoMatch","nExactMatch","nMatch","nCellBarcodes","nUMIs",};
    } stats;
        
    string cbSeq, umiSeq, cbQual, umiQual;
    
    SoloCB (Parameters &Pin, int iChunk);
    void readCB(const uint64 &iReadAll, const string &readNameExtra, const uint nTr, const vector<int32> &readGenes, vector<uint32> &readTranscripts, set<uint32> &readTrGenes);
    void addSoloCBcounts(const SoloCB &soloCBin);
    void addSoloCBstats(const SoloCB &soloCBin);
    void statsOut(ofstream &streamOut);
    void readCBgeneUMIfromFiles(uint32 **cbP, uint32 *cbReadCountExact);
    
private:
    Parameters &P;
    ParametersSolo &pSolo;

};

#endif
