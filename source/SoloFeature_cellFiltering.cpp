#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"
#include "ErrorWarning.h"

void SoloFeature::cellFiltering()
{    
    if (pSolo.cellFilter.type[0]=="None")
        return;
    
    //sort nUperCB
    nUMIperCBsorted=nUMIperCB;
    qsort(nUMIperCBsorted.data(), nCB, sizeof(uint32), funCompareNumbersReverse<uint32>); //sort by gene number

    uint32 nUMImax=0, nUMImin=0;
    if (pSolo.cellFilter.type[0]=="CellRanger2.2") {
        //find robust max
        nUMImax = nUMIperCBsorted[min(nCB-1,pSolo.cellFilter.cr2maxCellInd)];//robust estimate of the max UMI
        nUMImin = int(nUMImax/pSolo.cellFilter.cr2maxMinRatio+0.5);
    } else if (pSolo.cellFilter.type[0]=="TopCells") {
        nUMImin = nUMIperCBsorted[max(nCB-1,pSolo.cellFilter.topCells)];
    };

    cellFilterVec.resize(nCB,false);
    memset(&filteredCells,0,sizeof(filteredCells));

    bool *geneDetected = new bool[Trans.nGe];
    memset((void*) geneDetected, 0, Trans.nGe);

    for (uint32 icb=0; icb<nCB; icb++) {
        if (nUMIperCB[icb]>=nUMImin) {
            cellFilterVec[icb]=true;
            
            filteredCells.nCells++;

            filteredCells.nUMIinCells += nUMIperCB[icb];
            
            filteredCells.nGeneInCells += nGenePerCB[icb];
            filteredCells.nGenePerCell.push_back(nGenePerCB[icb]);
            
            filteredCells.nReadInCells += nReadPerCB[icb];
            filteredCells.nReadPerCell.push_back(nReadPerCB[icb]);
            
            for (uint32 ig=0; ig<nGenePerCB[icb]; ig++) {
                uint32 indG1=countCellGeneUMIindex[icb]+ig*countMatStride;
                geneDetected[countCellGeneUMI[indG1]]=true;            
            };
        };
    };
    
    filteredCells.nGeneDetected=0;
    for (uint32 ii=0; ii<Trans.nGe; ii++) {
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
};
