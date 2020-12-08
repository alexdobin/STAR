#ifndef H_SoloFilteredCells
#define H_SoloFilteredCells

#include "FilteredCells.h"


class SoloFilteredCells {
public:
    array<uint64,12> A;
    
    uint64 &nCells              =A[0];
    uint64 &nReadInCells        =A[1];
    uint64 &medianReadPerCell   =A[2];
    uint64 &meanReadPerCell     =A[3];
    uint64 &nUMIinCells         =A[4];
    uint64 &medianUMIperCell    =A[5];
    uint64 &meanUMIperCell      =A[6];
    uint64 &nGeneInCells        =A[7];
    uint64 &medianGenePerCell   =A[8];
    uint64 &meanGenePerCell     =A[9];
    uint64 &nGeneDetected       =A[10];
    uint64 &nCellsSimple        =A[11];
    
    vector<uint32> nUMIperCell, nReadPerCell, nGenePerCell;
    vector<bool> filterVecBool;
    
    void reset() {
        A={0}; 
        nUMIperCell={};
    };
};

#endif
