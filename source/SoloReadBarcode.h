#ifndef H_SoloReadBarcode
#define H_SoloReadBarcode
#include <set>
#include "IncludeDefine.h"
#include "Parameters.h"

class SoloReadBarcode {
public:
    uint32 homoPolymer[4];//homopolymer constants       
    string cbSeq, umiSeq, cbQual, umiQual;
    uint32 cbB,umiB;
    int32  cbI;
    int32  cbMatch;//0=exact, 1=1 match with 1MM, 2= >1 matches with 1MM
    string cbMatchString;//CB matches and qualities
    vector<uint32> cbMatchInd;//matches
    uint32 *cbReadCountExact;

    struct {
        enum {                 nNinBarcode,  nUMIhomopolymer,  nTooMany,  nNoMatch,  nStats};
        uint64 V[nStats];
        vector<string> names={"nNinBarcode","nUMIhomopolymer","nTooMany","nNoMatch"};
    } stats;
    
    SoloReadBarcode(Parameters &Pin);
    void getCBandUMI(string &readNameExtra);
    void addCounts(const SoloReadBarcode &rfIn);
    void addStats(const SoloReadBarcode &rfIn);
    void statsOut(ofstream &streamOut);
    
private:  
    Parameters &P;
    ParametersSolo &pSolo;
};

#endif
