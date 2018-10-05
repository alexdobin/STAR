#ifndef CODE_SoloReads
#define CODE_SoloReads
#include "IncludeDefine.h"

class SoloReads {
public:
//         struct {//uncollapsed read barcodes
//             uint64 nMax;//size of arrays below
//             uint64 N; //number of reads recorded 
//             uint32* cb;
//             uint32* umi;
//             uint32* gene; //gene
//         } reads;

    ofstream strUn
    SoloReads (Parameters &Pin);

private:
    Parameters &P;
    ParametersSolo &pSolo;

};

#endif
