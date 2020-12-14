#include "SoloFeature.h"
#include "serviceFuns.cpp"
#include <math.h>

void SoloFeature::cellFiltering()
{    

    if (pSolo.cellFilter.type[0]=="None" ||  nCB<1)
        return;
    if (featureType!=SoloFeatureTypes::Gene && featureType!=SoloFeatureTypes::GeneFull && featureType!=-1) {
        return;
    };

    //simple filtering first
    nUMIperCBsorted=nUMIperCB;
    qsort(nUMIperCBsorted.data(), nCB, sizeof(uint32), funCompareNumbersReverse<uint32>); //sort by gene number

    uint32 nUMImax=0, nUMImin=0;
    if (pSolo.cellFilter.type[0]=="TopCells") {
        nUMImin = nUMIperCBsorted[min(nCB-1,pSolo.cellFilter.topCells)];    
    } else {//other filtering types require simple filtering first
        //find robust max
        uint32 maxind=int(round(pSolo.cellFilter.knee.nExpectedCells*(1.0-pSolo.cellFilter.knee.maxPercentile)));//cell number for robust max
        nUMImax = nUMIperCBsorted[min(nCB-1,maxind)];//robust estimate of the max UMI
        nUMImin = int(round(nUMImax/pSolo.cellFilter.knee.maxMinRatio));
    };
    nUMImin=max(nUMImin,(uint32) 1);//cannot be zero
        
    filteredCells.reset(nCB); //all stats to 0

    for (uint32 icb=0; icb<nCB; icb++) {
        if (nUMIperCB[icb]>=nUMImin) {
            filteredCells.filtVecBool[icb]=true;
            filteredCells.nCellsSimple++;
        };
    };
    
    P.inOut->logMain << "cellFiltering: simple: nUMImax="<< nUMImax <<"; nUMImin="<< nUMImin <<"; nCellsSimple="<< filteredCells.nCellsSimple <<endl;

    if (pSolo.cellFilter.type[0]=="EmptyDrops_CR") {
        emptyDrops_CR();
    };
    
    //filtering is done: filtVecBool=true for kept cells, calculate statistics

    bool *geneDetected = new bool[featuresNumber]; //=true if a gene was detected in at least one cell
    memset((void*) geneDetected, 0, featuresNumber);

    for (uint32 icb=0; icb<nCB; icb++) {
        if (filteredCells.filtVecBool[icb]) {
            
            filteredCells.nCells++;

            filteredCells.nUMIinCells += nUMIperCB[icb]; //nUMIperCB was calculated for umiDedup-main
            
            filteredCells.nReadInCells += nReadPerCB[icb];
            filteredCells.nReadPerCell.push_back(nReadPerCB[icb]);
            
            uint32 ng1 = 0; //number of genes detected in this cell
            for (uint32 ig=0; ig<nGenePerCB[icb]; ig++) {
                uint32 indG1=countCellGeneUMIindex[icb]+ig*countMatStride;
                if (countCellGeneUMI[indG1 + pSolo.umiDedup.countInd.main] > 0) {
                    geneDetected[countCellGeneUMI[indG1]] = true; //gene is present if it's count > 0 for 
                    ng1++;
                };
            };
                        
            filteredCells.nGeneInCells += ng1; //had to recalculate this since some gene counts could be 0
            filteredCells.nGenePerCell.push_back(ng1);
        };
    };   
    
    if (filteredCells.nCells==0) {//all stats were already set to 0
        return;
    };
    
    filteredCells.nGeneDetected=0;
    for (uint32 ii=0; ii<featuresNumber; ii++) {
        if (geneDetected[ii])
            filteredCells.nGeneDetected++;
    };
    
    filteredCells.meanUMIperCell = filteredCells.nUMIinCells / filteredCells.nCells;
    filteredCells.meanReadPerCell = filteredCells.nReadInCells / filteredCells.nCells;
    filteredCells.meanGenePerCell = filteredCells.nGeneInCells / filteredCells.nCells;
    
    sort(filteredCells.nReadPerCell.begin(), filteredCells.nReadPerCell.end());
    sort(filteredCells.nGenePerCell.begin(), filteredCells.nGenePerCell.end());

    filteredCells.medianUMIperCell = nUMIperCBsorted[filteredCells.nCells/2];
    filteredCells.medianGenePerCell = filteredCells.nGenePerCell[filteredCells.nCells/2];
    filteredCells.medianReadPerCell = filteredCells.nReadPerCell[filteredCells.nCells/2];
    
    outputResults(true, outputPrefixFiltered);

    return;
};
