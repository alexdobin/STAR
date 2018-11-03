#ifndef CODE_Solo
#define CODE_Solo
#include "IncludeDefine.h"
#include "ReadAlignChunk.h"
#include "SoloCB.h"
#include "Transcriptome.h"

class Solo {   
public:

    SoloCB *soloCBsum, **soloCBall;
        
    uint64 nReadsMapped, nCB; //total number of mapped reads
    
    uint32 *rGeneUMI;//mapped reads sorted by CB
    uint32 *indCB;//index of detected CBs in the whitelist
    uint32 *rCBn;//number of reads for detected CBs in the whitelist
    uint32 **rCBp;//array of pointers to each CB sub-array
    
    Solo(Parameters &Pin, Transcriptome &inTrans);
    void soloPostMap(ReadAlignChunk **RAchunk);
    void collapseUMI(uint32 iCB);
    
private:
    Parameters &P;
    ParametersSolo &pSolo;
    Transcriptome &Trans;
};

#endif
