#ifndef CODE_Solo
#define CODE_Solo
#include "IncludeDefine.h"
#include "ReadAlignChunk.h"
#include "SoloCB.h"
#include "Transcriptome.h"
#include <fstream>

class Solo {   
public:

    SoloCB *soloCBsum, **soloCBall;
        
    uint64 nReadsMapped, nCB; //total number of mapped reads
    
    uint32 *rGeneUMI;//mapped reads sorted by CB
    uint32 *indCB;//index of detected CBs in the whitelist
    uint32 *rCBn;//number of reads for detected CBs in the whitelist
    uint32 **rCBp;//array of pointers to each CB sub-array
    uint32 *nUperCB;//number of UMIs per CB
    uint32 *nGperCB;//number of genes (with >0 UMIs) per CB
    uint32 nCellGeneEntries;//total number of non-zero cell/gene combinations (entries in the output matrix)
    
    ofstream *soloStatsStream;
    
    Solo(Parameters &Pin, Transcriptome &inTrans);
    void soloPostMap(ReadAlignChunk **RAchunk);
    void collapseUMI(uint32 *rGU, uint32 rN, uint32 &nGenes, uint32 &nUtot, uint32 *umiArray);
    void outputNumUMIperGeneCB();    


private:
    static const uint32 umiArrayStride=3;
    
    Parameters &P;
    ParametersSolo &pSolo;
    Transcriptome &Trans;
};

#endif
