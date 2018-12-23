#include "SoloCB.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"

void SoloCB::readCB(const uint64 &iReadAll, const string &readNameExtra, const uint nTr, const vector<int32> &readGenes, vector<uint32> &readTranscripts, set<uint32> &readTrGenes) {
    if (pSolo.type==0)
        return;
    
    int32 cbI=-999;
    uint32 cbB, umiB;
    
    cbSeq=readNameExtra.substr(0,pSolo.cbL);
    umiSeq=readNameExtra.substr(pSolo.cbL,pSolo.umiL);
    cbQual=readNameExtra.substr(pSolo.cbL+pSolo.umiL+1,pSolo.cbL);
    umiQual=readNameExtra.substr(pSolo.cbL+pSolo.umiL+1+pSolo.cbL,pSolo.umiL);
    
    int32 posN=convertNuclStrToInt32(cbSeq,cbB);
    if (posN==-2) {//>2 Ns, might already be filtered by Illumina
        stats.V[stats.nNinBarcode]++;
        return;
    } else if (posN==-1) {//no Ns
        cbI=binarySearchExact(cbB,pSolo.cbWL.data(),pSolo.cbWL.size());
        if (cbI>=0) {//exact match
            cbReadCountExact[cbI]++;//note that this simply counts reads per exact CB, no checks of genes or UMIs
        };
    };
    
    //check genes, return if no gene
    if (readTrGenes.size()==0) {
        stats.V[stats.nNoGene]++;
        return;
    } else if (readTrGenes.size()>1) {
        stats.V[stats.nAmbigGene]++;
        if (nTr>1)
            stats.V[stats.nAmbigGeneMultimap]++;
        return;
    };
    
    //check UMIs, return if bad UMIs
    if (convertNuclStrToInt32(umiSeq,umiB)!=-1) {//convert and check for Ns
        stats.V[stats.nNinBarcode]++;//UMIs are not allowed to have MMs
        return;
    };
    if (umiB==homoPolymer[0] || umiB==homoPolymer[1] || umiB==homoPolymer[2] || umiB==homoPolymer[3]) {
        stats.V[stats.nUMIhomopolymer]++;
        return;
    };
    
    if (cbI>=0) {//exact match
        cbReadCount[cbI]++;
        *strU_0 << cbI <<' '<< *readTrGenes.begin() <<' '<< umiB <<'\n';
        return;        
    };
        
    if (posN>=0) {//one N
        uint32 posNshift=2*(pSolo.cbL-1-posN);//shift bits for posN
        for (uint32 jj=0; jj<4; jj++) {
            uint32 cbB1=cbB^(jj<<posNshift);
            int32 cbI1=binarySearchExact(cbB1,pSolo.cbWL.data(),pSolo.cbWL.size());
            if (cbI1>=0) {                        
                if (cbI>=0) {//had another match already
                    stats.V[stats.nTooMany]++;
                    return;//with N in CB, do not allow matching >1 in WL
                };
                cbI=cbI1;
            };
        };
        if (cbI>=0) {
            cbReadCount[cbI]++;
            *strU_1 << cbI <<' '<< *readTrGenes.begin() <<' '<< umiB <<'\n';
        } else {//no match
            stats.V[stats.nNoMatch]++;
        };
        return;
    } else {//look for 1MM, posN==-1, no Ns
        string cbOutString("");
        uint32 ncbOut1=0;
        for (uint32 ii=0; ii<pSolo.cbL; ii++) {
            for (uint32 jj=1; jj<4; jj++) {
                int32 cbI1=binarySearchExact(cbB^(jj<<(ii*2)),pSolo.cbWL.data(),pSolo.cbWL.size());
                if (cbI1>=0) {//found match
                    //simple procedure: accept only if one WL CB exists with 1MM
                    //if (cbI>=0) {//had another match already
                    //    stats.V[stats.nTooMany]++;
                    //    //cout << iReadAll << endl;
                    //    return;    
                    //};
                    //cbI=cbI1;

                    //output all
                    cbI=cbI1;
                    ++ncbOut1;
                    cbOutString += ' ' +to_string(cbI) + ' ' + cbQual.at(pSolo.cbL-1-ii);
                    cbReadCount[cbI]++;//this read may be counted in many CBs. This is safe for allocating arrays later. The final number of reads per CB will be calculated later.                   
                };
            };
        };

        if (ncbOut1==0) {//no matches    
            stats.V[stats.nNoMatch]++;
            return;
        } else if (ncbOut1==1){//only one match
            *strU_1 << cbI <<' '<< *readTrGenes.begin() <<' '<< umiB <<'\n';
        } else {//more than one match
            *strU_2 << *readTrGenes.begin() <<' '<< umiB <<' '<< ncbOut1 <<  cbOutString <<'\n';
        };
    };

};
