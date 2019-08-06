#ifndef H_SoloReadFeature
#define H_SoloReadFeature
#include <set>
#include <map>
#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "SoloReadBarcode.h"
#include "SoloCommon.h"
#include "SoloReadFeatureStats.h"

class SoloReadFeature {
public:

    uint32 homoPolymer[4];//homopolymer constants

    uint32 *cbReadCount;
    map <uint32,uint32> cbReadCountMap;
    
    bool readInfoYes;

    fstream *streamReads;

    string cbSeq, umiSeq, cbQual, umiQual;
    
    SoloReadFeatureStats stats;

    SoloReadFeature (int32 feTy, Parameters &Pin, int iChunk);
    void record(SoloReadBarcode &soloBar, uint nTr, set<uint32> &readGene, set<uint32> &readGeneFull, Transcript *alignOut, uint64 iRead, const vector<array<uint32,2>> &readTranscripts);
    void addCounts(const SoloReadFeature &soloCBin);
    void addStats(const SoloReadFeature &soloCBin);
    void statsOut(ofstream &streamOut);
    void inputRecords(uint32 **cbP, uint32 cbPstride, uint32 *cbReadCountExact, ofstream *streamTranscriptsOut, vector<readInfoStruct> &readInfo);

private:
    const int32 featureType;

    Parameters &P;
    ParametersSolo &pSolo;
};

#endif
