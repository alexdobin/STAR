#include "SoloBarcode.h"
#include "serviceFuns.cpp"
#include "SoloCommon.h"
#include "SequenceFuns.h"

void wlAddMismatches(uint32 nMM, uint32 cbLen, vector<uintCB> &wl, vector<uint32> &wlInd);

void SoloBarcode::sortWhiteList()
{
    totalSize=0;
    minLen=(uint32)-1;
    wlAdd.resize( wl.size() );
    for (uint32 ilen=1; ilen < wl.size(); ilen++) {//scan through different lengths for this CB
        wlAdd[ilen]=totalSize;
        if (wl[ilen].size()>0) {
            if (ilen<minLen)
                minLen=ilen;
            std::sort(wl[ilen].begin(),wl[ilen].end());//sort
            auto un1=std::unique(wl[ilen].begin(),wl[ilen].end());//collapse identical
            wl[ilen].resize(std::distance(wl[ilen].begin(),un1));
            totalSize += wl[ilen].size();

            //add mismatches
            vector<uint32> wlInd;
            wlAddMismatches(3, ilen, wl[ilen], wlInd);

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

void wlAddMismatches(uint32 nMM, uint32 cbLen, vector<uintCB> &wl, vector<uint32> &wlInd)
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

            if (mm>0)
                continue; //only mm, no indels, for ED=1

            //insertions:
            uintCB maskCB = ~(((uintCB)-1)<<(cbLen*2));
            uintCB cbmm = cbMMind[ii].cb;
            //cout << convertNuclInt64toString(cbmm, cbLen)<<endl;
            for (uint32 shift1=0; shift1<=2; shift1+=2) {
                cbmm >>= shift1;
                for (uint32 ll=0; ll<cbLen*2; ll+=2) {
                    uintCB cbmm1 = cbmm<<2;
                    //cout << convertNuclInt64toString(cbmm1, cbLen) <<" "<< ll <<endl;
                    uintCB mask1 = ((uintCB)-1)<<(ll+2);
                    //cout << convertNuclInt64toString(mask1, cbLen)<<endl;
                    uintCB mask = (~mask1)>>2;
                    //cout << convertNuclInt64toString(mask, cbLen)<<endl;                
                    uintCB cbmm2 = (cbmm & mask) | (cbmm1 & mask1);
                    //cout << convertNuclInt64toString(cbmm2, cbLen)<<endl;
                    
                    for (uintCB jj=0; jj<4; jj++) {
                        uintCB cbmm3 = cbmm2 | (jj<<ll);//cbmm2 has A=00 in ll-position. Change this to one of 4 bases
                        //cout << convertNuclInt64toString(cbmm3, cbLen) <<" "<< jj<<endl;
                        cbMMind.push_back({cbmm3 & maskCB, cbMMind[ii].ind, mm});
                        ++ind3;
                    };
                    //cout << flush;
                };
            };
        };
        ind1 = ind2;
        ind2 = ind3;
    };

    std::sort(cbMMind.begin(), cbMMind.end(), [](const type_cbMMind &c1, const type_cbMMind &c2) 
                                                    {return (c1.cb<c2.cb) || (c1.cb==c2.cb && c1.mm<c2.mm) || (c1.cb==c2.cb && c1.mm==c2.mm && c1.ind<c2.ind);
                                                    } );

    uint64 nCBout = 0;
    vector<uint64> nCBmm(nMM+1,0);
    uintCB prevCB = (uintCB)(-1);
    for (uint64 ii=0; ii<cbMMind.size(); ii++) {

        if (ii<cbMMind.size()-1 && cbMMind[ii].cb==cbMMind[ii+1].cb && cbMMind[ii].mm==cbMMind[ii+1].mm && cbMMind[ii].ind==cbMMind[ii+1].ind) {
            continue;//skip identical records
        };

        if (    ( ii>0 && cbMMind[ii].cb == prevCB ) 
             || ( ii<cbMMind.size()-1 && cbMMind[ii].cb == cbMMind[ii+1].cb &&  cbMMind[ii].mm == cbMMind[ii+1].mm ) ) {
            cbMMind[ii].mm = nMM+1; //this cb is equal to previous, mark to skip
        } else {
            nCBout++;
            nCBmm[cbMMind[ii].mm]++;
            cout << convertNuclInt64toString(cbMMind[ii].cb, cbLen) <<' '<< cbMMind[ii].mm << " 1" << ' '<< convertNuclInt64toString(wl[cbMMind[ii].ind], cbLen); //<<' '<< cbMMind[ii].cb;
            cout <<'\n';
        };
        prevCB = cbMMind[ii].cb;
    };
    //exit(0);
    cout << "aaaa" << endl;
};