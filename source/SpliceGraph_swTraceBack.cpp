//#include "SpliceGraph.h"

//void SpliceGraph::swTraceBack(array<typeSeqLen, 2> &alignEnds, array<typeSeqLen, 2> &alignStarts)
//{
    
//     uint32 row = alignEnds[0];    
//     uint32 col = alignEnds[1];
//     uint32 rowT = 0;
//     uint32 colT = 0;
//     while(col > 0 && row > 0 && scoringMatrix[col][row] != 0) {
//         rowT = directionMatrix[col][row].second;
//         colT = directionMatrix[col][row].first;
// 
//         row = rowT;
//         col = colT;
//     };
//     alignStarts[0]=row;
//     alignStarts[1]=col;
//};
