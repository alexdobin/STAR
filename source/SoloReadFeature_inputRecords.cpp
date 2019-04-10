#include <cmath>
#include "SoloReadFeature.h"
#include "binarySearch2.h"
#include "serviceFuns.cpp"


bool inputFeatureUmi(fstream *strIn, int32 featureType, uint32 &feature, uint32 &umi, array<vector<uint64>,2> &sjAll)
{
    if (!(*strIn >> umi)) //end of file
        return false;

    if (featureType==0 || featureType==2) {//gene
        *strIn >> feature;
    } else if (featureType==1) {//sj
        uint32 sj[2];
        *strIn >> sj[0] >> sj[1];
        feature=(uint32) binarySearch2(sj[0],sj[1],sjAll[0].data(),sjAll[1].data(),sjAll[0].size());
    };

    return true;
};

void SoloReadFeature::inputRecords(uint32 **cbP, uint32 *cbReadCountExactTotal)
{
    {//load exact matches
        strU_0->flush();
        strU_0->seekg(0,ios::beg);
        uint32 feature, umi;
        int64 cb;
        while (inputFeatureUmi(strU_0, featureType, feature, umi, P.sjAll)) {
            *strU_0 >> cb;
            if (!pSolo.cbWLyes) //if no-WL, the full cbInteger was recorded - now has to be placed in order
                cb=binarySearchExact<uint64>(cb,pSolo.cbWL.data(),pSolo.cbWL.size());
            if (feature != (uint32)(-1)){
                cbP[cb][0]=feature;
                cbP[cb][1]=umi;
                cbP[cb]+=2;
                stats.V[stats.nExactMatch]++;
            };
        };
    };

    if (!pSolo.cbWLyes) //no WL => no mismatch check
        return;
    
    {//1 match
        strU_1->flush();
        strU_1->seekg(0,ios::beg);
        uint32 cb, feature, umi;
        while (inputFeatureUmi(strU_1,featureType, feature, umi, P.sjAll)) {
            *strU_1 >> cb;
            if (cbReadCountExactTotal[cb]>0) {
                if (feature != (uint32)(-1)){
                    cbP[cb][0]=feature;
                    cbP[cb][1]=umi;
                    cbP[cb]+=2;
                };
            } else {
                stats.V[stats.nNoExactMatch]++;
            };
        };
    };

    {//2 matches
        strU_2->flush();
        strU_2->seekg(0,ios::beg);
        uint32 cb=0, feature, umi, ncb;
        while (inputFeatureUmi(strU_2,featureType, feature, umi, P.sjAll)) {
            if (feature == (uint32) (-1)) {
                strU_2->ignore((uint32) (-1),'\n');//ignore until the end of the line
                continue; //nothing to record
            };
            *strU_2 >> ncb;
            float ptot=0.0,pmax=0.0;
            for (uint32 ii=0; ii<ncb; ii++) {
                uint32 cbin;
                char  qin;
                float pin;
                *strU_2 >> cbin >> qin;
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
                cbP[cb]+=2;
            } else {
                stats.V[stats.nTooMany]++;
            };
        };
    };
};
