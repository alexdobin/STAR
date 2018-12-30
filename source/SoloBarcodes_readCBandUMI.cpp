#include "SoloBarcodes.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"

bool SoloBarcodes::readCBandUMI(string &readNameExtra) 
{
    if (pSolo.type==0)
        return false;
    
    cbMatch=-1;
    
    cbSeq=readNameExtra.substr(0,pSolo.cbL);
    umiSeq=readNameExtra.substr(pSolo.cbL,pSolo.umiL);
    cbQual=readNameExtra.substr(pSolo.cbL+pSolo.umiL+1,pSolo.cbL);
    umiQual=readNameExtra.substr(pSolo.cbL+pSolo.umiL+1+pSolo.cbL,pSolo.umiL);
    
    //check UMIs, return if bad UMIs
    if (convertNuclStrToInt32(umiSeq,umiB)!=-1) {//convert and check for Ns
        stats.V[stats.nNinBarcode]++;//UMIs are not allowed to have MMs
        return false;
    };
    if (umiB==homoPolymer[0] || umiB==homoPolymer[1] || umiB==homoPolymer[2] || umiB==homoPolymer[3]) {
        stats.V[stats.nUMIhomopolymer]++;
        return false;
    };       
    
    //convert CB and check for Ns
    int32 posN=convertNuclStrToInt32(cbSeq,cbB);
    cbI=-999;
    if (posN==-2) {//>2 Ns, might already be filtered by Illumina
        stats.V[stats.nNinBarcode]++;
        return false;
    } else if (posN==-1 && featureType==0) {//no Ns, count only for feattureType==gene
        cbI=binarySearchExact(cbB,pSolo.cbWL.data(),pSolo.cbWL.size());
        if (cbI>=0) {//exact match
            cbReadCountExact[cbI]++;//note that this simply counts reads per exact CB, no checks of genes or UMIs
            cbMatch=0;
            return true;
        };
    };
    
    cbI=-999;
    if (posN>=0) {//one N
        uint32 posNshift=2*(pSolo.cbL-1-posN);//shift bits for posN
        for (uint32 jj=0; jj<4; jj++) {
            uint32 cbB1=cbB^(jj<<posNshift);
            int32 cbI1=binarySearchExact(cbB1,pSolo.cbWL.data(),pSolo.cbWL.size());
            if (cbI1>=0) {                        
                if (cbI>=0) {//had another match already
                    stats.V[stats.nTooMany]++;
                    return false;//with N in CB, do not allow matching >1 in WL
                };
                cbI=cbI1;
            };
        };
        if (cbI>=0) {
            cbMatch=1;
            return true;
        } else {//no match
            stats.V[stats.nNoMatch]++;
            return false;
        };
    };
    
    //look for 1MM, posN==-1, no Ns
    cbMatchString("");
    cbI=-999;
    cbMatch=0;
    for (uint32 ii=0; ii<pSolo.cbL; ii++) {
        for (uint32 jj=1; jj<4; jj++) {
            int32 cbI1=binarySearchExact(cbB^(jj<<(ii*2)),pSolo.cbWL.data(),pSolo.cbWL.size());
            if (cbI1>=0) {//found match
                //output all
                cbI=cbI1;
                ++cbMatch;
                cbMatchString += ' ' +to_string(cbI) + ' ' + cbQual.at(pSolo.cbL-1-ii);
                cbReadCount[cbI]++;//this read may be counted in many CBs. This is safe for allocating arrays later. The final number of reads per CB will be calculated later.                   
            };
        };
    };
    if (cbMatch==0) {//no matches    
        stats.V[stats.nNoMatch]++;
        cbMatch=-1;
        return false;
    } else {//one or more matches
        return true;
    };
};
