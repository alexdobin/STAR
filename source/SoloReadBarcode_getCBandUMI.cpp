#include "SoloReadBarcode.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"

void SoloReadBarcode::getCBandUMI(string &readNameExtra)
{
    if (pSolo.type==0)
        return;
    cbI=-999;

    cbMatch=-1;
    cbMatchString="";

    cbSeq=readNameExtra.substr(pSolo.cbS-1,pSolo.cbL);
    umiSeq=readNameExtra.substr(pSolo.umiS-1,pSolo.umiL);

    uint32 qualStart = readNameExtra.find(' ',pSolo.cbL+pSolo.umiL);
    cbQual=readNameExtra.substr(qualStart+pSolo.cbS,pSolo.cbL);
    umiQual=readNameExtra.substr(qualStart+pSolo.umiS,pSolo.umiL);

    //check UMIs, return if bad UMIs
    if (convertNuclStrToInt32(umiSeq,umiB)!=-1) {//convert and check for Ns
        stats.V[stats.nNinBarcode]++;//UMIs are not allowed to have Ns
        return;
    };
    if (umiB==homoPolymer[0] || umiB==homoPolymer[1] || umiB==homoPolymer[2] || umiB==homoPolymer[3]) {
        stats.V[stats.nUMIhomopolymer]++;
        return;
    };

    //convert CB and check for Ns
    int64 posN=convertNuclStrToInt64(cbSeq,cbB);

    if (!pSolo.cbWLyes) {//no whitelist - no search
        if (posN!=-1) {//Ns are present, discard this read
            stats.V[stats.nNinBarcode]++;
        } else {//no Ns
            cbI=(int64) cbB;//all possible barcodes are accepted. This will overflow if CB is longer than 31b
            cbMatch=0;
        };
        return;
    };

    if (posN==-2) {//>2 Ns, might already be filtered by Illumina
        stats.V[stats.nNinBarcode]++;
        return;
    } else if (posN==-1) {//no Ns, count only for featureType==gene
        cbI=binarySearchExact<uint64>(cbB,pSolo.cbWL.data(),pSolo.cbWL.size());
        if (cbI>=0) {//exact match
            cbReadCountExact[cbI]++;//note that this simply counts reads per exact CB, no checks of genes or UMIs
            cbMatch=0;
            return;
        };
    };

    if (posN>=0) {//one N
        uint32 posNshift=2*(pSolo.cbL-1-posN);//shift bits for posN
        for (uint32 jj=0; jj<4; jj++) {
            uint64 cbB1=cbB^(jj<<posNshift);
            int64 cbI1=binarySearchExact<uint64>(cbB1,pSolo.cbWL.data(),pSolo.cbWL.size());
            if (cbI1>=0) {
                if (cbI>=0) {//had another match already
                    stats.V[stats.nTooMany]++;
                    return;//with N in CB, do not allow matching >1 in WL
                };
                cbI=cbI1;
            };
        };
        if (cbI>=0) {
            cbMatch=1;
            return;
        } else {//no match
            stats.V[stats.nNoMatch]++;
            return;
        };
    };

    //look for 1MM, posN==-1, no Ns
    cbMatch=0;
    cbMatchInd.clear();
    for (uint32 ii=0; ii<pSolo.cbL; ii++) {
        for (uint32 jj=1; jj<4; jj++) {
            int64 cbI1=binarySearchExact<uint64>(cbB^(jj<<(ii*2)),pSolo.cbWL.data(),pSolo.cbWL.size());
            if (cbI1>=0) {//found match
                //output all
                cbI=cbI1;
                cbMatchInd.push_back(cbI1);
                ++cbMatch;
                cbMatchString += ' ' +to_string(cbI1) + ' ' + cbQual.at(pSolo.cbL-1-ii);
            };
        };
    };
    if (cbMatch==0) {//no matches
        stats.V[stats.nNoMatch]++;
        cbMatch=-1;
    };// else cbMatch contains number of matches (1 or >1), and cbMatchString contains matches for >1 case
};
