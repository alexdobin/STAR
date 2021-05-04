#ifndef H_SoloFilteredCells
#define H_SoloFilteredCells

class SoloFilteredCells {
public:
    array<uint64,13> C;
    
    uint64 &nCells              =C[0];
    uint64 &nReadInCells        =C[1];
    uint64 &medianReadPerCellUnique   =C[2];
    uint64 &meanReadPerCellUnique     =C[3];
    uint64 &nUMIinCells         =C[4];
    uint64 &medianUMIperCell    =C[5];
    uint64 &meanUMIperCell      =C[6];
    uint64 &nGeneInCells        =C[7];
    uint64 &medianGenePerCell   =C[8];
    uint64 &meanGenePerCell     =C[9];
    uint64 &nGeneDetected       =C[10];
    uint64 &nCellsSimple        =C[11];
    uint64 &nReadInCellsUnique  =C[12];
    
    vector<uint32> nReadPerCellUnique, nGenePerCell;
    vector<bool> filtVecBool;
    
    void reset(uint64 nCells) {
        C={0}; 
        nReadPerCellUnique.clear(); nReadPerCellUnique.reserve(16384);
        nGenePerCell.clear(); nGenePerCell.reserve(16384);
        filtVecBool.clear();  filtVecBool.resize(nCells, false);
    };
};

#endif
