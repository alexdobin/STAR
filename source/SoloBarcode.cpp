#include "SoloBarcode.h"
#include "serviceFuns.cpp"
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