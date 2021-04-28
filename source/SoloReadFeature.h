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
#include "ReadAnnotations.h"

class SoloFeature;

class SoloReadFeature {
public:

    uint32 homoPolymer[4];//homopolymer constants

    vector<uint32> cbReadCount;
    map <uintCB,uint32> cbReadCountMap;
    
    vector<uint32> transcriptDistCount;
    
    bool readInfoYes ,readIndexYes;

    fstream *streamReads;

    string cbSeq, umiSeq, cbQual, umiQual;
    
    SoloReadFeatureStats stats;

    SoloReadFeature (int32 feTy, Parameters &Pin, int iChunk);
    void record(SoloReadBarcode &soloBar, uint nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot);
    void addCounts(const SoloReadFeature &soloCBin);
    void addStats(const SoloReadFeature &soloCBin);
    void statsOut(ofstream &streamOut);
    void inputRecords(uint32 **cbP, uint32 cbPstride, vector<uint32> &cbReadCountTotal, vector<readInfoStruct> &readInfo);

private:
    const int32 featureType;

    Parameters &P;
    ParametersSolo &pSolo;
};

#endif
