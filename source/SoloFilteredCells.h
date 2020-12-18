#ifndef H_SoloFilteredCells
#define H_SoloFilteredCells

class SoloFilteredCells {
public:
    array<uint64,12> C;
    
    uint64 &nCells              =C[0];
    uint64 &nReadInCells        =C[1];
    uint64 &medianReadPerCell   =C[2];
    uint64 &meanReadPerCell     =C[3];
    uint64 &nUMIinCells         =C[4];
    uint64 &medianUMIperCell    =C[5];
    uint64 &meanUMIperCell      =C[6];
    uint64 &nGeneInCells        =C[7];
    uint64 &medianGenePerCell   =C[8];
    uint64 &meanGenePerCell     =C[9];
    uint64 &nGeneDetected       =C[10];
    uint64 &nCellsSimple        =C[11];
    
    vector<uint32> nReadPerCell, nGenePerCell;
    vector<bool> filtVecBool;
    
    void reset(uint64 nCells) {
        C={0}; 
        nReadPerCell.clear(); nReadPerCell.reserve(16384);
        nGenePerCell.clear(); nGenePerCell.reserve(16384);
        filtVecBool.clear();  filtVecBool.resize(nCells, false);
    };
};

#endif
