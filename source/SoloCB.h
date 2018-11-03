#ifndef CODE_SoloCB
#define CODE_SoloCB
#include "IncludeDefine.h"
#include "Parameters.h"

class SoloCB {
public:
    uint32* cbReadCount;

    struct {
        enum {nNoGene,nAmbigGene,nAmbigGeneMultimap,nNinBarcode,nTooMany,nNoMatch,nExactMatch,nMatch,nStats};
        uint64 V[nStats];
        vector<string> names={"nNoGene","nAmbigGene","nAmbigGeneMultimap","nNinBarcode","nTooMany","nNoMatch","nExactMatch","nMatch"};
    } stats;
    
    fstream *strU_0; //uniqe mappers, CB matches whitelist
    //ofstream *strU_1; //uniqe mappers, CB with 1MM off whitelist
    
    SoloCB (Parameters &Pin, int iChunk);
    void readCB(const uint64 &iReadAll, const string &readNameExtra, const uint nTr, const vector<int32> &readGenes, vector<uint32> &readTranscripts, set<uint32> &readTrGenes);
    void addSoloCB(const SoloCB &soloCBin);
    void statsOut(ofstream &streamOut);
    void readCBgeneUMIfromFiles(uint32 ** cbP);
    
private:
    Parameters &P;
    ParametersSolo &pSolo;

};

#endif
