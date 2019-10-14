#include <cmath>
#include "SoloReadFeature.h"
#include "binarySearch2.h"
#include "serviceFuns.cpp"
#include "SoloCommon.h"
#include "SoloFeature.h"

bool inputFeatureUmi(fstream *strIn, int32 featureType, bool readInfoYes, array<vector<uint64>,2> &sjAll, uint32 &iread, int32 &cbmatch, uint32 &feature, uint32 &umi, vector<uint32> &featVecU32)
{
    if (!(*strIn >> umi)) //end of file
        return false;

    if (readInfoYes)
        *strIn >> iread;

    switch (featureType) {
        case SoloFeatureTypes::Gene :
        case SoloFeatureTypes::GeneFull :
            *strIn >> feature;
            break;

        case SoloFeatureTypes::SJ :
            uint32 sj[2];
            *strIn >> sj[0] >> sj[1];
            feature=(uint32) binarySearch2(sj[0],sj[1],sjAll[0].data(),sjAll[1].data(),sjAll[0].size());
            break;

        case SoloFeatureTypes::Transcript3p :
            feature=0;
            uint32 ntr, in1;
            *strIn >> ntr;
            featVecU32.resize(2*ntr);
            for (uint32 ii=0; ii<2*ntr; ii++) {
                *strIn >> in1;
                featVecU32[ii]=in1;
            };
            break;
        };

    *strIn >> cbmatch;

    return true;
};

void SoloReadFeature::inputRecords(uint32 **cbP, uint32 cbPstride, uint32 *cbReadCountTotal, 
                                   ofstream *streamTranscriptsOut, vector<readInfoStruct> &readInfo, SoloFeature **soloFeatAll)
{   
    streamReads->flush();
    streamReads->seekg(0,ios::beg);
    
    switch (featureType) {//for non-standard processing
        
        case SoloFeatureTypes::VelocytoSimple :
        {
            uint32 feature, vtype;
            uint64 iread;
            while (*streamReads >> iread) {//until the end of file
                *streamReads >> feature >> vtype;
                uint64 cb=soloFeatAll[pSolo.featureInd[SoloFeatureTypes::Gene]]->readInfo[iread].cb;
                if (cb+1!=0) {
                    cbP[cb][0]=feature | (vtype << velocytoTypeGeneBitShift ); //encode vType in the top 2 bits;
                    cbP[cb][1]=soloFeatAll[pSolo.featureInd[SoloFeatureTypes::Gene]]->readInfo[iread].umi;
                    cbP[cb]+=cbPstride;
                };
            };
            return;
            break;
        };
    };
    
    //////////////////////////////////////////// standard features
    uint32 feature, umi, iread;
    int32 cbmatch;
    int64 cb;
    vector<uint32> trIdDist;
    while (inputFeatureUmi(streamReads, featureType, readInfoYes, P.sjAll, iread, cbmatch, feature, umi, trIdDist)) {
        if (feature == (uint32)(-1) && !readInfoYes) {//no feature => no record, this can happen for SJs
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
                cb=binarySearchExact<uint64>(cb,pSolo.cbWL.data(),pSolo.cbWLsize);

            //record feature
            if (featureType==SoloFeatureTypes::Transcript3p) {//variable length feature, separate treatment, feature always defined (i.e. !=-1)
                //for now - output all in file
                *streamTranscriptsOut << cb <<' '<< umi <<' '<< trIdDist.size()/2;
                for (auto &tt: trIdDist)
                    *streamTranscriptsOut <<' '<< tt;
                *streamTranscriptsOut << '\n';
                if (cbmatch==0)
                    stats.V[stats.nExactMatch]++;
            } else {//single-number feature
                if (feature != (uint32)(-1)) {
                    cbP[cb][0]=feature;
                    cbP[cb][1]=umi;
                    if (readInfoYes) {
                        cbP[cb][2]=iread;
                    };
                    cbP[cb]+=cbPstride;
                    if (cbmatch==0)
                        stats.V[stats.nExactMatch]++;
                } else {//no feature - record readInfo
                    readInfo[iread].cb=cb;
                    readInfo[iread].umi=umi;
                };
            };
        } else {//multiple matches
            float ptot=0.0,pmax=0.0;
            for (uint32 ii=0; ii<(uint32)cbmatch; ii++) {
                uint32 cbin;
                char  qin;
                float pin;
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
                //record feature
                if (featureType==SoloFeatureTypes::Transcript3p) {//variable length feature
                    //for now - output all in file
                    *streamTranscriptsOut << cb <<' '<< umi <<' '<< trIdDist.size()/2;
                    for (auto &tt: trIdDist)
                        *streamTranscriptsOut <<' '<< tt;
                    *streamTranscriptsOut << '\n';
                } else {//single-number feature
                    if (feature != (uint32)(-1)) {
                        cbP[cb][0]=feature;
                        cbP[cb][1]=umi;
                        if (readInfoYes) {
                            cbP[cb][2]=iread;
                        };
                        cbP[cb]+=cbPstride;
                    } else {//no feature - record readInfo
                        readInfo[iread].cb=cb;
                        readInfo[iread].umi=umi;
                    };    
                };
            } else if (feature != (uint32)(-1)) {
                stats.V[stats.nTooMany]++;
            };
        };
    };
};

