#ifndef H_SoloReadBarcode
#define H_SoloReadBarcode
#include <set>
#include "IncludeDefine.h"
#include "Parameters.h"
#include "SoloReadBarcodeStats.h"

class SoloReadBarcode {
private:
    Parameters &P;
    
public:
    ParametersSolo &pSolo;    

    
    uint32 homoPolymer[4];//homopolymer constants
    string cbSeq, umiSeq, cbQual, umiQual, bSeq, bQual;
    vector<string> bStrings; //barcode strings from SAM tags
    string cbSeqCorrected;
    uint64 umiB;
    //int64  cbI;
    int32 cbMatch;//-1: no match, 0: exact, 1: 1 match with 1MM, >1: # of matches with 1MM
    int32 umiCheck;//umi check status
    string cbMatchString;//CB matches and qualities
    vector<uint64> cbMatchInd;//matches
    vector<uint32> cbReadCountExact;
    //map <uint32,uint32> cbReadCountMap;//count read per CB for no WL

    array<uint64,256> qualHist;
    SoloReadBarcodeStats stats;

    SoloReadBarcode(Parameters &Pin);
    void getCBandUMI(char **readSeq, char **readQual, uint64 *readLen, const string &readNameExtraIn, const uint32 &readFilesIndex, const char *readName);
    void addCounts(const SoloReadBarcode &rfIn);
    void addStats(const SoloReadBarcode &rfIn);
    void statsOut(ofstream &streamOut);
    void matchCBtoWL(string &cbSeq1, string &cbQual1, vector<uint64> &cbWL, int32 &cbMatch1, vector<uint64> &cbMatchInd1, string &cbMatchString1);
    bool convertCheckUMI();
    void addStats(const int32 cbMatch1);
    

protected:
    void readNameExtra();
};

#endif
