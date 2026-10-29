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

struct SoloRecordBuf {//one record parsed by inputRecordsBuffered, replayed by replayRecords
    uint64 cb;
    uint32 feature, umi, iread;
    uint32 nReadInc;//0: none, 1: nReadPerCBunique1[cb]++, 2: nReadPerCBmulti1[cb]++
};

class SoloReadFeature {
public:

    uint32 homoPolymer[4];//homopolymer constants

    vector<uint32> cbReadCount;
    map <uintCB,uint32> cbReadCountMap;
    
    vector<uint32> transcriptDistCount;
    
    bool readInfoYes ,readIndexYes;

    fstream *streamReads;

    string cbSeq, umiSeq, cbQual, umiQual;

    SoloReadFlagClass readFlag;
    
    SoloReadFeatureStats stats;

    SoloReadFeature (int32 feTy, Parameters &Pin, int iChunk);
    void record(SoloReadBarcode &soloBar, uint nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot);
    void addCounts(const SoloReadFeature &soloCBin);
    void addStats(const SoloReadFeature &soloCBin);
    void statsOut(ofstream &streamOut);
    void inputRecords(uint32 **cbP, uint32 cbPstride, vector<uint32> &cbReadCountTotal, vector<readInfoStruct> &readInfo, SoloReadFlagClass &readFlagCounts,
                      vector<uint32> &nReadPerCBunique1, vector<uint32> &nReadPerCBmulti1);
    //inputRecords split in two, so that the files can be parsed concurrently: the parse writes only to readInfo (distinct reads per file), stats,
    //readFlagCounts (a per-file copy) and the buffers; replayRecords then makes inputRecords' cbP and nReadPerCB writes in file order
    void inputRecordsBuffered(vector<uint32> &cbReadCountTotal, vector<readInfoStruct> &readInfo, SoloReadFlagClass &readFlagCounts,
                              vector<SoloRecordBuf> &recs, vector<uintCB> &flagCountsOrder);
    void replayRecords(uint32 **cbP, uint32 cbPstride, const vector<SoloRecordBuf> &recs, vector<uint32> &nReadPerCBunique1, vector<uint32> &nReadPerCBmulti1);

private:
    const int32 featureType;

    template <class Sink>
    void inputRecordsT(Sink &sink, vector<uint32> &cbReadCountTotal, vector<readInfoStruct> &readInfo, SoloReadFlagClass &readFlagCounts);

    Parameters &P;
    ParametersSolo &pSolo;
};

#endif
