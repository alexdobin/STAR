#include "SoloFeature.h"
#include "serviceFuns.cpp"
#include <math.h>
#include <unordered_set>

void SoloFeature::emptyDrops_CR()
{
    
    if (nCB<=pSolo.cellFilter.eDcr.indMin) {
        P.inOut->logMain << "emptyDrops_CR filtering: no empty cells found: nCB=" << nCB <<"   emptyCellMinIndex="<< pSolo.cellFilter.eDcr.indMin << "\n";
        return;
    };
    
    //find genes that were detected in all cells
    unordered_set<uint32> featDet;   
    for (uint32 icb=0; icb<nCB; icb++) {           
        for (uint32 ig=0; ig<nGenePerCB[icb]; ig++) {
            uint32 g1 = countCellGeneUMI[countCellGeneUMIindex[icb]+ig*countMatStride];
            featDet.insert(g1);
        };
    };
    uint32 nFeatDet=featDet.size(); //total number of detected genes - this should have been done already?
    
    //sum gene expression over the collection of empty cells
    typedef struct {uint32 index, count;} IndCount;
    vector<IndCount> indCount(nCB);
    for (uint32 ii=0; ii<nCB; ii++) {
        indCount[ii].index=ii;
        indCount[ii].count=nUMIperCB[ii];
    };
    std::sort(indCount.begin(), indCount.end(), [](const IndCount &ic1, const IndCount &ic2) {
                                                    return ic1.count>ic2.count; //descending order
                                                });
    

    //ambient gene counts
    vector<uint32> ambCount(featuresNumber,0);
    for (uint32 icb=pSolo.cellFilter.eDcr.indMin; icb<min(nCB,pSolo.cellFilter.eDcr.indMax); icb++) {
        for (uint32 ig=0; ig<nGenePerCB[icb]; ig++) {
            uint32 irec = countCellGeneUMIindex[icb]+ig*countMatStride;
            ambCount[countCellGeneUMI[irec+0]] += countCellGeneUMI[irec+1];
        };
    };
    
    //frequencies
    unordered_map<uint32,uint32> ambCountFreq;
    for (auto &ac: ambCount) {
        ambCountFreq[ac]++;
    };
    
     
    return;
};
