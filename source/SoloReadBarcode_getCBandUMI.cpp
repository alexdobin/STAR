#include "SoloReadBarcode.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"

void SoloReadBarcode::matchCBtoWL(string &cbSeq1, string &cbQual1, uint64 &cbB1, int32 &cbMatch1, vector<uint64> &cbMatchInd1, string &cbMatchString1)
{
    //convert CB and check for Ns
    int64 posN=convertNuclStrToInt64(cbSeq1,cbB1);

    if (!pSolo.cbWLyes) {//no whitelist - no search
        if (posN!=-1) {//Ns are present, discard this read
            stats.V[stats.nNinBarcode]++;
        } else {//no Ns
            //cbI=(int64) cbB1;
            cbMatchInd1.push_back(cbB1);//all possible barcodes are accepted. This will overflow if CB is longer than 31b
            cbMatchString1 = to_string(cbB1);
            cbMatch1=0;
        };
        return;
    };

    if (posN==-2) {//>2 Ns, might already be filtered by Illumina
        stats.V[stats.nNinBarcode]++;
        return;
    } else if (posN==-1) {//no Ns, count only for featureType==gene
        int64 cbI=binarySearchExact<uint64>(cbB1,pSolo.cbWL.data(),pSolo.cbWL.size());
        if (cbI>=0) {//exact match
            cbReadCountExact[cbI]++;//note that this simply counts reads per exact CB, no checks of genes or UMIs
            cbMatchInd1.push_back((uint64) cbI);
            cbMatchString1 = to_string(cbMatchInd1[0]);
            cbMatch1=0;
            return;
        };
    };

    if (posN>=0) {//one N
        int64 cbI=-1;
        uint32 posNshift=2*(pSolo.cbL-1-posN);//shift bits for posN
        for (uint32 jj=0; jj<4; jj++) {
            uint64 cbB11=cbB1^(jj<<posNshift);
            int64 cbI1=binarySearchExact<uint64>(cbB11,pSolo.cbWL.data(),pSolo.cbWL.size());
            if (cbI1>=0) {
                if (cbI>=0) {//had another match already
                    stats.V[stats.nTooMany]++;
                    return;//with N in CB, do not allow matching >1 in WL
                };
                cbI=cbI1;
            };
        };
        if (cbI>=0) {
            cbMatchInd1.push_back((uint64) cbI);
            cbMatchString1 = to_string(cbMatchInd1[0]);
            cbMatch1=1;
            return;
        } else {//no match
            stats.V[stats.nNoMatch]++;
            return;
        };
    };

    //look for 1MM; posN==-1, no Ns
    cbMatch1=0;
    for (uint32 ii=0; ii<pSolo.cbL; ii++) {
        for (uint32 jj=1; jj<4; jj++) {
            int64 cbI1=binarySearchExact<uint64>(cbB1^(jj<<(ii*2)),pSolo.cbWL.data(),pSolo.cbWL.size());
            if (cbI1>=0) {//found match
                //output all
                cbMatchInd1.push_back(cbI1);
                ++cbMatch1;
                cbMatchString1 += ' ' +to_string(cbI1) + ' ' + cbQual1.at(pSolo.cbL-1-ii);
            };
        };
    };
    if (cbMatch1==0) {//no matches
        stats.V[stats.nNoMatch]++;
        cbMatch1=-1;
    } else if (cbMatch1==1) {//1 match, no need to record the quality
        cbMatchString1 = to_string(cbMatchInd1[0]);
    };// else cbMatch contains number of matches, and cbMatchString has CBs and qualities
};   

void SoloReadBarcode::getCBandUMI(string &readNameExtra)
{
    if (pSolo.type==0)
        return;
    //int64 cbI=-999;

    cbMatch=-1;
    cbMatchString="";
    cbMatchInd.clear();
    
    uint32 bLength = readNameExtra.find(' ',pSolo.cbL+pSolo.umiL);
    string bSeq=readNameExtra.substr(0,bLength);
    string bQual=readNameExtra.substr(bLength+1,2*bLength);

    if (pSolo.type==1) {
        cbSeq=bSeq.substr(pSolo.cbS-1,pSolo.cbL);
        umiSeq=bSeq.substr(pSolo.umiS-1,pSolo.umiL);
        cbQual=bQual.substr(pSolo.cbS-1,pSolo.cbL);
        umiQual=bQual.substr(pSolo.umiS-1,pSolo.umiL);
    } else if (pSolo.type==2) {
        uint32 adapterStart=0;
        if (pSolo.adapterYes) {
            if (localAlignHammingDist(bSeq, pSolo.adapterSeq, adapterStart) > pSolo.adapterMismatchesNmax) {
                //TODO: add stats
                return; //no adapter found
            };
        };
        umiSeq="";
        umiQual="";
        if (!pSolo.umiV.extractBarcode(bSeq, bQual, adapterStart, umiSeq, umiQual)) {
            //TODO: add stats
            return;
        };

        cbSeq="";
        cbQual="";
        for (auto &cb : pSolo.cbV) {
            if (!cb.extractBarcode(bSeq, bQual, adapterStart, umiSeq, umiQual)) {
                //TODO: add stats
                return;
            };
        };
    };
    
    //check UMIs, return if bad UMIs
    if (convertNuclStrToInt32(umiSeq,umiB)!=-1) {//convert and check for Ns
        stats.V[stats.nNinBarcode]++;//UMIs are not allowed to have Ns
        return;
    };
    if (umiB==homoPolymer[0] || umiB==homoPolymer[1] || umiB==homoPolymer[2] || umiB==homoPolymer[3]) {
        stats.V[stats.nUMIhomopolymer]++;
        return;
    };

    matchCBtoWL(cbSeq, cbQual, cbB, cbMatch, cbMatchInd, cbMatchString);
};
