#include <cmath>
#include "SoloReadFeature.h"
#include "SoloCommon.h"
#include "SoloFeature.h"
#include "soloInputFeatureUMI.h"
#include "serviceFuns.cpp"

namespace {
struct SoloSinkDirect {//writes the shared arrays as each record is parsed
    uint32 **cbP;
    uint32 cbPstride;
    bool readIndexYes;
    vector<uint32> &nReadPerCBunique1, &nReadPerCBmulti1;

    void record(int64 cb, uint32 feature, uint64 umi, uint64 iread)
    {
        cbP[cb][0]=feature;
        cbP[cb][1]=umi;

        if (readIndexYes) {
            cbP[cb][2]=iread;
        };

        cbP[cb]+=cbPstride;
    };
    void countUnique(int64 cb) {nReadPerCBunique1[cb]++;};
    void countMulti(int64 cb) {nReadPerCBmulti1[cb]++;};
    void countsAdd(SoloReadFlagClass &readFlagCounts, int64 cb) {readFlagCounts.countsAdd(cb);};
};

struct SoloSinkBuffer {//stores the records for replayRecords, and the order in which CBs first enter flagCounts
    vector<SoloRecordBuf> &recs;
    vector<uintCB> &flagCountsOrder;

    void record(int64 cb, uint32 feature, uint64 umi, uint64 iread)
    {
        recs.push_back({(uint64)cb, feature, (uint32)umi, (uint32)iread, 0});
    };
    //a count always follows the record() of the same read
    void countUnique(int64) {recs.back().nReadInc=1;};
    void countMulti(int64) {recs.back().nReadInc=2;};
    void countsAdd(SoloReadFlagClass &readFlagCounts, int64 cb)
    {
        auto cbInserted = readFlagCounts.flagCounts.insert({(uintCB)cb, {} });
        if (cbInserted.second)
            flagCountsOrder.push_back(cb);
        for (uint32 ibit=0; ibit<readFlagCounts.nBits; ibit++) {
            (*cbInserted.first).second[ibit] += (uint64) readFlagCounts.checkBit(ibit);
        };
    };
};
};

void SoloReadFeature::inputRecords(uint32 **cbP, uint32 cbPstride, vector<uint32> &cbReadCountTotal, vector<readInfoStruct> &readInfo, SoloReadFlagClass &readFlagCounts,
                                   vector<uint32> &nReadPerCBunique1, vector<uint32> &nReadPerCBmulti1)
{
    SoloSinkDirect sink{cbP, cbPstride, readIndexYes, nReadPerCBunique1, nReadPerCBmulti1};
    inputRecordsT(sink, cbReadCountTotal, readInfo, readFlagCounts);
};

void SoloReadFeature::inputRecordsBuffered(vector<uint32> &cbReadCountTotal, vector<readInfoStruct> &readInfo, SoloReadFlagClass &readFlagCounts,
                                           vector<SoloRecordBuf> &recs, vector<uintCB> &flagCountsOrder)
{
    SoloSinkBuffer sink{recs, flagCountsOrder};
    inputRecordsT(sink, cbReadCountTotal, readInfo, readFlagCounts);
};

void SoloReadFeature::replayRecords(uint32 **cbP, uint32 cbPstride, const vector<SoloRecordBuf> &recs, vector<uint32> &nReadPerCBunique1, vector<uint32> &nReadPerCBmulti1)
{
    SoloSinkDirect sink{cbP, cbPstride, readIndexYes, nReadPerCBunique1, nReadPerCBmulti1};
    for (const auto &r : recs) {
        sink.record(r.cb, r.feature, r.umi, r.iread);
        if (r.nReadInc==1) {
            sink.countUnique(r.cb);
        } else if (r.nReadInc==2) {
            sink.countMulti(r.cb);
        };
    };
};

template <class Sink>
void SoloReadFeature::inputRecordsT(Sink &sink, vector<uint32> &cbReadCountTotal, vector<readInfoStruct> &readInfo, SoloReadFlagClass &readFlagCounts)
{
    streamReads->flush();
    streamReads->seekg(0,std::ios::beg);

    //////////////////////////////////////////// standard features
    uint32 feature;
    uint64 umi, iread, prevIread=(uint64)-1;
    int32 cbmatch;
    int64 cb;
    vector<uint32> trIdDist;
    
    uint64 nReadsIn = 0;

    while (soloInputFeatureUMI(streamReads, featureType, readIndexYes, P.sjAll, iread, cbmatch, feature, umi, trIdDist, readFlagCounts)) {
        if (feature == (uint32)(-1) && !readIndexYes) {//no feature => no record, this can happen for SJs
            streamReads->ignore((uint32)-1, '\n');
            //stats.V[stats.noNoFeature]++; //need separate category for this
            continue;
        };

        bool readIsCounted = false;
        bool featGood = ( feature != (uint32)(-1) );
        bool noMMtoWLwithoutExact = false;
        bool noTooManyWLmatches = false;

        if (cbmatch<=1) {//single match
            *streamReads >> cb;

            if ( pSolo.CBmatchWL.oneExact && cbmatch==1 && cbReadCountTotal[cb]==0 ) {//single 1MM match, no exact matches to this CB
                noMMtoWLwithoutExact = true;

            } else {

                if (!pSolo.cbWLyes) {//if no-WL, the full cbInteger was recorded - now has to be placed in order
                    cb=binarySearchExact<uintCB>(cb, pSolo.cbWL.data(), pSolo.cbWLsize);
                    if (cb+1 == 0)
                        continue; //this cb was not in the tentative WL
                };

                //record feature
                if (featGood) {//good feature, will be counted
                    readIsCounted = true;

                    sink.record(cb, feature, umi, iread);

                } else if (readInfoYes) {//no feature - record readInfo
                    readInfo[iread].cb=cb;
                    readInfo[iread].umi=umi;
                };
            };

        } else {//multiple matches
            #ifdef MATCH_CellRanger
            double ptot=0.0, pmax=0.0, pin;
            #else
            float ptot=0.0, pmax=0.0, pin;
            #endif

            for (uint32 ii=0; ii<(uint32)cbmatch; ii++) {
                uint32 cbin;
                char  qin;
                *streamReads >> cbin >> qin;
                if (cbReadCountTotal[cbin]>0) {//otherwise this cbin does not work
                    qin -= pSolo.QSbase;
                    qin = qin < pSolo.QSmax ? qin : pSolo.QSmax;
                    pin=cbReadCountTotal[cbin]*std::pow(10.0,-qin/10.0);
                    ptot+=pin;
                    if (pin>pmax) {
                        cb=cbin;
                        pmax=pin;
                    };
                };
            };
            if (ptot>0.0 && pmax>=pSolo.cbMinP*ptot) {
                //record feature single-number feature
                if (featGood) {
                    readIsCounted = true;
                    sink.record(cb, feature, umi, iread);
                } else if (readInfoYes) {//no feature - record readInfo
                    readInfo[iread].cb=cb;
                    readInfo[iread].umi=umi;
                };
            } else {
                noTooManyWLmatches = true;
            };
        };

        if ( !readIndexYes || iread != prevIread ) {//no readindex (then only unique-gene reads are recorded) OR only for one align of each read, in case of multimappers
            prevIread = iread;
            if (featGood) {
                if (cbmatch==0) {
                    stats.V[stats.yessubWLmatchExact]++;
                } else if (noMMtoWLwithoutExact) {
                    stats.V[stats.noMMtoWLwithoutExact]++;
                } else if (noTooManyWLmatches) {
                    stats.V[stats.noTooManyWLmatches]++;
                };
            };

            if (readIsCounted) {
                if (feature<geneMultMark) {
                    sink.countUnique(cb);
                } else {
                    sink.countMulti(cb);
                };
            };
            
            if ( pSolo.readStatsYes[featureType] ) {//has to be new iread to avoid muti-counting multi-gene reads
                
                //readIsCounted flag was defined above
                if ( readIsCounted ) {
                    if ( readFlagCounts.checkBit(readFlagCounts.featureU) )
                        readFlagCounts.setBit(readFlagCounts.countedU);
                    if ( readFlagCounts.checkBit(readFlagCounts.featureM) )
                        readFlagCounts.setBit(readFlagCounts.countedM);    
                };

                readFlagCounts.setBit(readFlag.cbMatch); //set this flag for both CB and no-CB, but no-CB will be counted separately below
                //only record this read in readFlagCounts if its CB was defined
                if (cbmatch==0) {//perfect CB match to WL
                    readFlagCounts.setBit(readFlagCounts.cbPerfect);
                    sink.countsAdd(readFlagCounts, cb);
                } else if (cbmatch==1 && !noMMtoWLwithoutExact) {
                    readFlagCounts.setBit(readFlagCounts.cbMMunique);
                    sink.countsAdd(readFlagCounts, cb);
                } else if (cbmatch>1 && !noTooManyWLmatches) {
                    readFlagCounts.setBit(readFlagCounts.cbMMmultiple);
                    sink.countsAdd(readFlagCounts, cb);
                } else {//no CB match
                    readFlagCounts.countsAddNoCB();
                };

                // debug
                /*nReadsIn++;
                uint64 ntot=0;
                for (auto &ii: readFlagCounts.flagCounts)
                    ntot += ii.second[0];

                if (nReadsIn!=ntot)
                    cout << nReadsIn <<" "<< ntot << endl;
                */
                /*//if (featureType==SoloFeatureTypes::Gene && readFlagCounts.flagCounts[cb][readFlagCounts.featureM]+readFlagCounts.flagCounts[cb][readFlagCounts.featureU] != readFlagCounts.flagCounts[cb][readFlagCounts.exonic])
                //    cout << cb <<' '<< iread << endl;
                //if (readFlagCounts.checkBit(readFlagCounts.featureM))
                //    cout << iread <<' '<< readFlagCounts.flagCounts[cb][readFlagCounts.featureM] << endl;
                */
            };
        };
    };
};

