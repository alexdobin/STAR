#ifndef H_SoloFeature
#define H_SoloFeature
#include "IncludeDefine.h"
#include "ReadAlignChunk.h"
#include "Transcriptome.h"
#include <fstream>

#include "SoloCommon.h"
#include "SoloRead.h"

class SoloFeature {
private:
    const int32 featureType;

    Parameters &P;
    Transcriptome &Trans;

    static const uint32 umiArrayStride=3;
    enum {rguG, rguU, rguR};
    uint32 rguStride;
    
    SoloFeature **soloFeatAll;
    
public:
    ParametersSolo &pSolo;

    SoloReadFeature *readFeatSum, **readFeatAll;
    SoloReadBarcode *readBarSum;

    uint64 nReadsMapped, nCB, nReadsInput; //total number of mapped reads

    uint32 *rGeneUMI;//mapped reads sorted by CB
    uint32 *indCB;//index of detected CBs in the whitelist
    uint32 *rCBn;//number of reads for detected CBs in the whitelist
    uint32 **rCBp;//array of pointers to each CB sub-array

    vector<uint32> nUMIperCB, nUMIperCBsorted;//number of UMIs per CB, and the same sorted (descendant)
    vector<uint32> nGenePerCB;//number of genes (with >0 UMIs) per CB
    vector<uint32> nReadPerCB;//number of reads per CB
    
    vector<uint32> countCellGeneUMI;//sparsified matrix for the counts, each entry is: geneID count1 count2 ... countNcounts
    vector<uint32> countCellGeneUMIindex;//index of CBs in the count matrix
    uint32 countMatStride; //number of counts per entry in the count matrix
    
    //vector<vector<array<uint32,4>>> cellFeatureCounts; //another way to collect the cell/feature counts
    
    vector<bool> cellFilterVec;
    struct {
        uint64 nCells, nReadInCells, medianReadPerCell, meanReadPerCell, nUMIinCells, medianUMIperCell, meanUMIperCell, nGeneInCells, medianGenePerCell, meanGenePerCell, nGeneDetected;
        vector<uint32> nUMIperCell, nReadPerCell, nGenePerCell;
    } filteredCells;
    
    string outputPrefix;
    ofstream *streamTranscriptsOut;
    
    array<vector<uint64>,2> sjAll;
    
    vector<readInfoStruct> readInfo; //corrected CB/UMI information for each read

    SoloFeature(int32 feTy, Parameters &Pin, Transcriptome &inTrans, SoloReadBarcode *readBarSumIn, SoloFeature **soloFeatAll);
    void processRecords(ReadAlignChunk **RAchunk);
    void sumThreads(ReadAlignChunk **RAchunk);
    void countCBgeneUMI();
    void countVelocyto();
    void collapseUMI(uint32 iCB, uint32 *umiArray);
    void outputResults(bool cellFilterYes);
    void addBAMtags(char *&bam0, uint32 &size0, char* bam1);
    void statsOutput();
    void cellFiltering();
};

#endif
