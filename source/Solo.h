#ifndef CODE_Solo
#define CODE_Solo
#include "IncludeDefine.h"
#include "ReadAlignChunk.h"
#include "SoloCB.h"

class Solo {   
public:
//         struct {//uncollapsed read barcodes
//             uint64 nMax;//size of arrays below
//             uint64 N; //number of reads recorded 
//             uint32* cb;
//             uint32* umi;
//             uint32* gene; //gene
//         } reads;
    SoloCB *soloCBall;
        
    Solo(Parameters &Pin);
    soloPostMap(ReadAlignChunk **RAchunk);
       
private:
    Parameters &P;
    ParametersSolo &pSolo;

};

#endif
