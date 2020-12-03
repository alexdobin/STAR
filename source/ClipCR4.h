#ifndef CODE_ClipCR4
#define CODE_ClipCR4

#include "IncludeDefine.h"
#include "opal/opal.h"

class ClipCR4
{
public:
    int dbN;//number of sequence in the opal "database"
    
    vector<uint8*> storeClip;

    // Results for each sequence in database
    vector<OpalSearchResult> opalRes;
    
    OpalSearchResult** opalResP;
   
    //constructor
    ClipCR4();
    void opalFillOneSeq(uint32 idb, char *seq, uint32 seqL);
    void opalAlign(uint8 *query, uint32 queryLen, int dbN1);
    uint32 polyTail3p(char *seq, uint32 seqLen);
    
private:

    uint32 readLen; //sequence length to align against    
    
    int alphabetLength;
    int gapOpen;
    int gapExt;
    vector<int> scoreMatrix;
    
    // Database
    uint8_t* dbSeqArr;
    uint8_t** dbSeqs;
    int* dbSeqsLen;
    
};

#endif
