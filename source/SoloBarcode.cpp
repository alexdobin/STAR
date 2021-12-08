#include "SoloBarcode.h"
#include "serviceFuns.cpp"
#include "SoloCommon.h"
#include "SequenceFuns.h"
#include "ParametersSolo.h"

void wlAddMismatches(uint32 nMM, uint32 cbLen, vector<uintCB> &wl, vector<uintCB> &wlEd1, vector<uint32> &wlInd1);

void SoloBarcode::sortWhiteList(ParametersSolo *pSolo)
{
    totalSize=0;
    minLen=(uint32)-1;
    wlAdd.resize( wl.size() );
    if (pSolo->CBmatchWL.EditDist_2) {
        wlEd.resize( wl.size() );
        wlEdInd.resize( wl.size() );
    };

    for (uint32 ilen1=1; ilen1 < wl.size(); ilen1++) {//scan through different lengths for this CB
        wlAdd[ilen1]=totalSize;
        if (wl[ilen1].size()>0) {
            if (ilen1<minLen)
                minLen=ilen1;
            std::sort(wl[ilen1].begin(),wl[ilen1].end());//sort
            auto un1=std::unique(wl[ilen1].begin(),wl[ilen1].end());//collapse identical
            wl[ilen1].resize(std::distance(wl[ilen1].begin(),un1));
            totalSize += wl[ilen1].size();

            if (pSolo->CBmatchWL.EditDist_2) {//add mismatches
                wlAddMismatches(2, ilen1, wl[ilen1], wlEd[ilen1], wlEdInd[ilen1]);
            };

        };
    };
};

void SoloBarcode::extractPositionsFromString(string &strIn)
{
    vector<string> p(4);
    splitString(strIn,'_',p);
    anchorType[0] = std::stoi(p[0]);
    anchorType[1] = std::stoi(p[2]);
    anchorDist[0] = std::stoi(p[1]);
    anchorDist[1] = std::stoi(p[3]);
};

void wlAddMismatches(uint32 nMM, uint32 cbLen, vector<uintCB> &wl, vector<uintCB> &wlEd1, vector<uint32> &wlEdInd1)
{
    struct type_cbMMind {
        uintCB cb;
        uint32 ind;
        uint32  mm;
    };

    //max size
    uint64 ntot=wl.size()*(std::pow(cbLen*3,nMM+1)-1)/(cbLen*3-1); //geometric series 1+L+L*L+...

    vector<type_cbMMind> cbMMind;
    cbMMind.reserve(ntot);

    for (uint32 icb=0; icb<wl.size(); icb++) {//fill in 0-MM CBs from wl
        cbMMind.push_back({wl[icb], icb, 0});
    };

    uint64 ind1=0, ind2=wl.size(); //define boundaries of the previous mismatch level
    for (uint32 mm=1; mm<=nMM; mm++) {//main cycle over number of MM
        uint64 ind3=ind2; //current index
        for (uint64 ii=ind1; ii<ind2; ii++) {//add 1 mismatch to previous-mismatch-level sequences
            for (uint32 ll=0; ll<cbLen*2; ll+=2) {//2 bits for each base of CB
                for (uintCB jj=1; jj<4; jj++) {//change to 1 of 3 values. Note that jj=0 corresponds to no change
                    uintCB cbmm = cbMMind[ii].cb^(jj<<ll);
                    //if (cbmm==27 && mm==3)
                    //    cout << ii;
                    cbMMind.push_back({cbmm, cbMMind[ii].ind, mm});
                    ++ind3;
                };
            };
        };

        if (mm==2) {//ins+del only added at mm=ed=2
            for (uint64 ii=0; ii<wl.size(); ii++) {//ins+del added to original barcodes (without mm)
                uintCB cbmm = cbMMind[ii].cb;
                //cout << convertNuclInt64toString(cbmm, cbLen)<<endl;
                for (uint32 ld=0; ld<cbLen*2; ld+=2) {//del
                    uintCB maskd1 = ((uintCB)-1)<<(ld+2);
                    uintCB maskd = (~maskd1)>>2;
                    uintCB cbmmd = (cbmm & maskd) | ((cbmm & maskd1)>>2);
                    //cout << convertNuclInt64toString(cbmmd, cbLen)<<endl;
                    for (uint32 ll=0; ll<cbLen*2; ll+=2) {//ins
                        uintCB cbmm1 = cbmmd<<2;
                        //cout << convertNuclInt64toString(cbmm1, cbLen) <<" "<< ll <<endl;
                        uintCB mask1 = ((uintCB)-1)<<(ll+2);
                        //cout << convertNuclInt64toString(mask1, cbLen)<<endl;
                        uintCB mask = (~mask1)>>2;
                        //cout << convertNuclInt64toString(mask, cbLen)<<endl;                
                        uintCB cbmm2 = (cbmmd & mask) | (cbmm1 & mask1);
                        //cout << convertNuclInt64toString(cbmm2, cbLen)<<endl;
                        
                        for (uintCB jj=0; jj<4; jj++) {
                            uintCB cbmm3 = cbmm2 | (jj<<ll);//cbmm2 has A=00 in ll-position. Change this to one of 4 bases
                            //cout << convertNuclInt64toString(cbmm3, cbLen) <<" "<< jj<<endl;
                            cbMMind.push_back({cbmm3, cbMMind[ii].ind, mm});
                            ++ind3;
                        };
                        //cout << flush;
                    };
                };
            };
        };
        ind1 = ind2;
        ind2 = ind3;
    };

    //select edited SB sequences that have a unique match in 
    std::sort(cbMMind.begin(), cbMMind.end(), [](const type_cbMMind &c1, const type_cbMMind &c2) 
                                                    {return (c1.cb<c2.cb) || (c1.cb==c2.cb && c1.mm<c2.mm) || (c1.cb==c2.cb && c1.mm==c2.mm && c1.ind<c2.ind);
                                                    } );

    uint64 nCBout = 0;
    vector<uint64> nCBmm(nMM+1,0);
    uintCB prevCB = (uintCB)(-1);
    for (uint64 ii=0; ii<cbMMind.size(); ii++) {

        if (ii<cbMMind.size()-1 && cbMMind[ii].cb==cbMMind[ii+1].cb && cbMMind[ii].mm==cbMMind[ii+1].mm && cbMMind[ii].ind==cbMMind[ii+1].ind) {
            cbMMind[ii].mm = nMM+1;
            continue;//skip identical records
        };

        if (    ( ii>0 && cbMMind[ii].cb == prevCB ) 
             || ( ii<cbMMind.size()-1 && cbMMind[ii].cb == cbMMind[ii+1].cb &&  cbMMind[ii].mm == cbMMind[ii+1].mm ) ) {
            cbMMind[ii].mm = nMM+1; //this cb is equal to previous, mark to skip
            //cout << convertNuclInt64toString(cbMMind[ii].cb, cbLen) <<' '<< cbMMind[ii].mm << " 2" << ' '<< convertNuclInt64toString(wl[cbMMind[ii].ind], cbLen)<<' '<<wl[cbMMind[ii].ind]<<'\n';
        } else {
            nCBout++;
            nCBmm[cbMMind[ii].mm]++;
            //cout << convertNuclInt64toString(cbMMind[ii].cb, cbLen) <<' '<< cbMMind[ii].mm << " 1" << ' '<< convertNuclInt64toString(wl[cbMMind[ii].ind], cbLen); //<<' '<< cbMMind[ii].cb;
            //cout <<'\n';
        };
        prevCB = cbMMind[ii].cb;
    };

    wlEd1.resize(nCBout);
    wlEdInd1.resize(nCBout);
    uint32 icb=0;
    for (auto &cb1: cbMMind) {
        if (cb1.mm<=nMM) {
            wlEd1[icb] = cb1.cb;
            wlEdInd1[icb] = cb1.ind;
            ++icb;
        };
    };
};