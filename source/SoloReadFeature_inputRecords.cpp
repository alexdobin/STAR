#include <cmath>
#include "SoloReadFeature.h"
#include "binarySearch2.h"
#include "serviceFuns.cpp"

bool inputFeatureUmi(fstream *strIn, int32 featureType, bool readInfoYes, uint32 &iread, int32 &cbmatch, uint32 &feature, uint32 &umi, array<vector<uint64>,2> &sjAll)
{
    if (!(*strIn >> umi)) //end of file
        return false;

    if (readInfoYes)
        *strIn >> iread;
    
    *strIn >> cbmatch;
    
    if (featureType==0 || featureType==2) {//gene
        *strIn >> feature;
    } else if (featureType==1) {//sj
        uint32 sj[2];
        *strIn >> sj[0] >> sj[1];
        feature=(uint32) binarySearch2(sj[0],sj[1],sjAll[0].data(),sjAll[1].data(),sjAll[0].size());
    };

    return true;
};

void SoloReadFeature::inputRecords(uint32 **cbP, uint32 cbPstride, uint32 *cbReadCountExactTotal)
{   
    {//load exact matches
        streamReads->flush();
        streamReads->seekg(0,ios::beg);
        uint32 feature, umi, iread;
        int32 cbmatch;
        int64 cb;
        while (inputFeatureUmi(streamReads, featureType, readInfoYes, iread, cbmatch, feature, umi, P.sjAll)) {
            if (feature == (uint32)(-1)) {//no feature => no record, this can happen for SJs
                streamReads->ignore((uint32)-1, '\n');
                continue;
            };
            
            if (cbmatch<=1) {//single match
                *streamReads >> cb;
                
                if (cbmatch==1 && cbReadCountExactTotal[cb]==0) {//single 1MM match, no exact matches to this CB
                    stats.V[stats.nNoExactMatch]++;
                    continue;
                };
                
                if (!pSolo.cbWLyes) //if no-WL, the full cbInteger was recorded - now has to be placed in order
                    cb=binarySearchExact<uint64>(cb,pSolo.cbWL.data(),pSolo.cbWL.size());
                
                cbP[cb][0]=feature;
                cbP[cb][1]=umi;
                if (readInfoYes) {
                    cbP[cb][2]=iread;
                };
                cbP[cb]+=cbPstride;
                if (cbmatch==0)
                    stats.V[stats.nExactMatch]++;
            } else {//multiple matches
                float ptot=0.0,pmax=0.0;
                for (uint32 ii=0; ii<(uint32)cbmatch; ii++) {
                    uint32 cbin;
                    char  qin;
                    float pin;
                    *streamReads >> cbin >> qin;
                    if (cbReadCountExactTotal[cbin]>0) {//otherwise this cbin does not work
                        qin -= pSolo.QSbase;
                        qin = qin < pSolo.QSmax ? qin : pSolo.QSmax;
                        pin=cbReadCountExactTotal[cbin]*std::pow(10.0,-qin/10.0);
                        ptot+=pin;
                        if (pin>pmax) {
                            cb=cbin;
                            pmax=pin;
                        };
                    };
                };
                if (ptot>0.0 && pmax>=pSolo.cbMinP*ptot) {
                    cbP[cb][0]=feature;
                    cbP[cb][1]=umi;
                    if (readInfoYes) {
                        cbP[cb][2]=iread;
                    };                
                    cbP[cb]+=cbPstride;
                } else {
                    stats.V[stats.nTooMany]++;
                };
            };
        };
    };

//     if (!pSolo.cbWLyes) //no WL => no mismatch check
//         return;
//     
//     {//1 match
//         strU_1->flush();
//         strU_1->seekg(0,ios::beg);
//         uint32 cb, feature, umi, iread;
//         while (inputFeatureUmi(strU_1,featureType, readInfoYes, iread, feature, umi, P.sjAll)) {
//             *strU_1 >> cb;
//             if (cbReadCountExactTotal[cb]>0) {
//                 if (feature != (uint32)(-1)){
//                     cbP[cb][0]=feature;
//                     cbP[cb][1]=umi;
//                     if (readInfoYes) {
//                         cbP[cb][2]=iread;
//                     };
//                     cbP[cb]+=cbPstride;
//                 };
//             } else {
//                 stats.V[stats.nNoExactMatch]++;
//             };
//         };
//     };
// 
//     {//>1 matches
//         strU_2->flush();
//         strU_2->seekg(0,ios::beg);
//         uint32 cb=0, feature, umi, ncb, iread;
//         while (inputFeatureUmi(strU_2,featureType, readInfoYes, iread, feature, umi, P.sjAll)) {
//             if (feature == (uint32) (-1)) {
//                 strU_2->ignore((uint32) (-1),'\n');//ignore until the end of the line
//                 continue; //nothing to record
//             };
//             *strU_2 >> ncb;
//             float ptot=0.0,pmax=0.0;
//             for (uint32 ii=0; ii<ncb; ii++) {
//                 uint32 cbin;
//                 char  qin;
//                 float pin;
//                 *strU_2 >> cbin >> qin;
//                 if (cbReadCountExactTotal[cbin]>0) {//otherwise this cbin does not work
//                     qin -= pSolo.QSbase;
//                     qin = qin < pSolo.QSmax ? qin : pSolo.QSmax;
//                     pin=cbReadCountExactTotal[cbin]*std::pow(10.0,-qin/10.0);
//                     ptot+=pin;
//                     if (pin>pmax) {
//                         cb=cbin;
//                         pmax=pin;
//                     };
//                 };
//             };
//             if (ptot>0.0 && pmax>=pSolo.cbMinP*ptot) {
//                 cbP[cb][0]=feature;
//                 cbP[cb][1]=umi;
//                 if (readInfoYes) {
//                     cbP[cb][2]=iread;
//                 };                
//                 cbP[cb]+=cbPstride;
//             } else {
//                 stats.V[stats.nTooMany]++;
//             };
//         };
//     };
};
