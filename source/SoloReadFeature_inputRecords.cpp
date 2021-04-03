#include <cmath>
#include "SoloReadFeature.h"
#include "SoloCommon.h"
#include "SoloFeature.h"
#include "soloInputFeatureUMI.h"
#include "serviceFuns.cpp"

void SoloReadFeature::inputRecords(uint32 **cbP, uint32 cbPstride, vector<uint32> &cbReadCountTotal, vector<readInfoStruct> &readInfo)
{   
    streamReads->flush();
    streamReads->seekg(0,std::ios::beg);

    //////////////////////////////////////////// standard features
    uint32 feature;
    uint64 umi, iread;
    int32 cbmatch;
    int64 cb;
    vector<uint32> trIdDist;
    while (soloInputFeatureUMI(streamReads, featureType, readIndexYes, P.sjAll, iread, cbmatch, feature, umi, trIdDist)) {
        if (feature == (uint32)(-1) && !readIndexYes) {//no feature => no record, this can happen for SJs
            streamReads->ignore((uint32)-1, '\n');
            //stats.V[stats.nNoFeature]++; //need separate category for this
            continue;
        };

        if (cbmatch<=1) {//single match
            *streamReads >> cb;

            if ( pSolo.CBmatchWL.oneExact && cbmatch==1 && cbReadCountTotal[cb]==0 && feature!=(uint32)(-1) ) {//single 1MM match, no exact matches to this CB
                stats.V[stats.nNoExactMatch]++;
                continue;
            };

            if (!pSolo.cbWLyes) //if no-WL, the full cbInteger was recorded - now has to be placed in order
                cb=binarySearchExact<uintCB>(cb, pSolo.cbWL.data(), pSolo.cbWLsize);

            //record feature ingle-number feature
            if (feature != (uint32)(-1)) {
                cbP[cb][0]=feature;
                cbP[cb][1]=umi;
                if (readIndexYes) {
                    cbP[cb][2]=iread;
                };
                cbP[cb]+=cbPstride;
                if (cbmatch==0)
                    stats.V[stats.nExactMatch]++;
            } else if (readInfoYes) {//no feature - record readInfo
                readInfo[iread].cb=cb;
                readInfo[iread].umi=umi;
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
                if (feature != (uint32)(-1)) {
                    cbP[cb][0]=feature;
                    cbP[cb][1]=umi;
                    if (readIndexYes) {
                        cbP[cb][2]=iread;
                    };
                    cbP[cb]+=cbPstride;
                } else if (readInfoYes) {//no feature - record readInfo
                    readInfo[iread].cb=cb;
                    readInfo[iread].umi=umi;
                };    
            } else if (feature != (uint32)(-1)) {
                stats.V[stats.nTooMany]++;
            };
        };
    };
};

