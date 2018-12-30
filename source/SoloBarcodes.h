#ifndef CODE_SoloCB
#define CODE_SoloCB
#include <set>
#include "IncludeDefine.h"
#include "Parameters.h"

class SoloBarcodes {
public:
    
    uint32 homoPolymer[4];//homopolymer constants       
    string cbSeq, umiSeq, cbQual, umiQual;
    uint32 cbB,umiB;
    int32  cbI;
    int32  cbMatch;//0=exact, 1=1 match with 1MM, 2= >1 matches with 1MM
    string cbMatchString;//CB matches and qualities
    uint32 *cbReadCountExact;

    struct {
        enum {                 nNinBarcode,  nUMIhomopolymer,  nTooMany,  nNoExactMatch,  nNoMatch,  nExactMatch,  nMatch, nStats};
        uint64 V[nStats];
        vector<string> names={"nNinBarcode","nUMIhomopolymer","nTooMany","nNoExactMatch","nNoMatch","nExactMatch","nMatch"};
    } stats;
    
    
    SoloBarcodes(Parameters &Pin);
    
private:  
    Parameters &P;
    ParametersSolo &pSolo;
};

#endif
