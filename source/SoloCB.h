#ifndef CODE_SoloCB
#define CODE_SoloCB
#include "IncludeDefine.h"
#include "Parameters.h"

class SoloCB {
public:
    uint32* cbReadCount;

    ofstream *strU_0; //uniqe mappers, CB matches whitelist
    ofstream *strU_1; //uniqe mappers, CB with 1MM off whitelist
    SoloCB (Parameters &Pin, int iChunk);
    void readCB(const string &readNameExtra, const uint nTr, const vector<int32> &readGene);

private:
    Parameters &P;
    ParametersSolo &pSolo;

};

#endif
