#include "SoloFeature.h"
#include "serviceFuns.cpp"
#include <math.h>

void SoloFeature::cellFiltering()
{    

    if (pSolo.cellFilter.type[0]=="None" ||  nCB<1) {//no filtering, or no cells to filter
        return;
    };
    
    switch(featureType) {
        default:
            return; //filtering not done for other features
            
        case SoloFeatureTypes::Velocyto: {
            filteredCells.reset(nCB);
            
            //Velocyto: use filter from Gene
            SoloFeature &soFeGe= *soloFeatAll[pSolo.featureInd[SoloFeatureTypes::Gene]];
            for (uint32 ic=0; ic<soFeGe.nCB; ic++) {
                if (soFeGe.filteredCells.filtVecBool[ic] && indCBwl[soFeGe.indCB[ic]] != (uint32)-1) {
                    //this cell passde Gene filtering, and is detected in Velocyto
                    filteredCells.filtVecBool[indCBwl[soFeGe.indCB[ic]]] = true;
                };
            };
                        
            nUMIperCBsorted=nUMIperCB;
            std::sort( nUMIperCBsorted.begin(), nUMIperCBsorted.end(), [](const uint32_t &u1, const uint32_t &u2) {return u1>u2;} ); //sort descending
            break;
        };
        case SoloFeatureTypes::Gene:
        case SoloFeatureTypes::GeneFull:
        case SoloFeatureTypes::GeneFull_Ex50pAS:
        case SoloFeatureTypes::GeneFull_ExonOverIntron: 
        case -1: //undefined: matrix loaded from file
        {
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// cell calling (filtering)
            //do the filtering

            //simple filtering first
            nUMIperCBsorted=nUMIperCB;
            std::sort( nUMIperCBsorted.begin(), nUMIperCBsorted.end(), [](const uint32_t &u1, const uint32_t &u2) {return u1>u2;} ); //sort descending

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
        };
    };
    //filtering is done: filtVecBool=true for kept cells
    
    
    //calculate filtered statistics

    std::vector<uint32_t> geneDetected(featuresNumber, 0); //=1 if a gene was detected in at least one cell

    for (uint32 icb=0; icb<nCB; icb++) {
        if (filteredCells.filtVecBool[icb]) {
            
            filteredCells.nCells++;

            filteredCells.nUMIinCells += nUMIperCB[icb]; //nUMIperCB was calculated for umiDedup-main
            
            if (nReadPerCBunique.size()>0) {//for CellFiltering only, read information is not available
                filteredCells.nReadInCells += nReadPerCBtotal[icb];
                filteredCells.nReadInCellsUnique += nReadPerCBunique[icb];            
                filteredCells.nReadPerCellUnique.push_back(nReadPerCBunique[icb]);
            };
            
            uint32 ng1 = 0; //number of genes detected in this cell
            for (uint32 ig=0; ig<nGenePerCB[icb]; ig++) {
                uint32 indG1=countCellGeneUMIindex[icb]+ig*countMatStride;
                if (countCellGeneUMI[indG1 + pSolo.umiDedup.countInd.main] > 0) {
                    geneDetected[countCellGeneUMI[indG1]] = 1; //gene is present if it's count > 0 for 
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
        if (geneDetected[ii]>0)
            filteredCells.nGeneDetected++;
    };
    
    filteredCells.meanUMIperCell = filteredCells.nUMIinCells / filteredCells.nCells;
    filteredCells.meanReadPerCellUnique = filteredCells.nReadInCellsUnique / filteredCells.nCells;
    filteredCells.meanGenePerCell = filteredCells.nGeneInCells / filteredCells.nCells;
    
    std::sort(filteredCells.nReadPerCellUnique.begin(), filteredCells.nReadPerCellUnique.end());
    std::sort(filteredCells.nGenePerCell.begin(), filteredCells.nGenePerCell.end());

    filteredCells.medianUMIperCell = nUMIperCBsorted[filteredCells.nCells/2];
    filteredCells.medianGenePerCell = filteredCells.nGenePerCell[filteredCells.nCells/2];
    filteredCells.medianReadPerCellUnique = filteredCells.nReadPerCellUnique[filteredCells.nCells/2];
    
    //////////////////////////////////////////////////////////////////output filtered matrix
    outputResults(true, outputPrefixFiltered);

    return;
};
