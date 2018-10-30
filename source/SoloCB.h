#ifndef CODE_SoloCB
#define CODE_SoloCB
#include "IncludeDefine.h"
#include "Parameters.h"

class SoloCB {
public:
    uint32* cbReadCount;

    struct {
        enum {nNoGene,nAmbigGene,nAmbigGeneMultimap,nNinBarcode,nStats};
        uint64 V[nStats];
//         uint64 nNoGene;
//         uint64 nAmbigGene;
//         uint64 nAmbigGeneMultimap;
//         
//         uint64 nNinBarcode;
    } stats;
    
    ofstream *strU_0; //uniqe mappers, CB matches whitelist
    ofstream *strU_1; //uniqe mappers, CB with 1MM off whitelist
    
    SoloCB (Parameters &Pin, int iChunk);
    void readCB(const string &readNameExtra, const uint nTr, const vector<int32> &readGenes, vector<uint32> &readTranscripts, set<uint32> &readTrGenes);

private:
    Parameters &P;
    ParametersSolo &pSolo;

};

#endif
