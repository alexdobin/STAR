#ifndef H_SoloFeature
#define H_SoloFeature
#include "IncludeDefine.h"
#include "ReadAlignChunk.h"
#include "Transcriptome.h"
#include <fstream>

#include "SoloRead.h"

class SoloFeature {   
public:

    SoloReadFeature *readFeatSum, **readFeatAll;
    SoloReadBarcode *readBarSum, **readBarAll;

    uint64 nReadsMapped, nCB; //total number of mapped reads
    
    uint32 *rGeneUMI;//mapped reads sorted by CB
    uint32 *indCB;//index of detected CBs in the whitelist
    uint32 *rCBn;//number of reads for detected CBs in the whitelist
    uint32 **rCBp;//array of pointers to each CB sub-array
    uint32 *nUperCB;//number of UMIs per CB
    uint32 *nGperCB;//number of genes (with >0 UMIs) per CB
    uint32 nCellGeneEntries;//total number of non-zero cell/gene combinations (entries in the output matrix)
    
    ofstream *statsStream;
    
    array<vector<uint64>,2> sjAll;
    
    SoloFeature(int feTy, Parameters &Pin, Transcriptome &inTrans);
    void processRecords(ReadAlignChunk **RAchunk);
    void collapseUMI(uint32 *rGU, uint32 rN, uint32 &nGenes, uint32 &nUtot, uint32 *umiArray);
    void outputResults();    

private:   
    const int32 featureType;
    
    Parameters &P;
    ParametersSolo &pSolo;
    Transcriptome &Trans;
    
    static const uint32 umiArrayStride=3;    
};

#endif
