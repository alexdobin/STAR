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
    
public:
    ParametersSolo &pSolo;

    SoloReadFeature *readFeatSum, **readFeatAll;
    SoloReadBarcode *readBarSum, **readBarAll;

    uint64 nReadsMapped, nCB; //total number of mapped reads

    uint32 *rGeneUMI;//mapped reads sorted by CB
    uint32 *indCB;//index of detected CBs in the whitelist
    uint32 *rCBn;//number of reads for detected CBs in the whitelist
    uint32 **rCBp;//array of pointers to each CB sub-array

    vector<uint32> nUMIperCB;//number of UMIs per CB
    vector<uint32> nGenePerCB;//number of genes (with >0 UMIs) per CB
    vector<uint64> nReadPerCB;//number of reads per CB
    
    vector<bool> cellFilterVec;
    struct {
        uint64 nCells, nReadInCells, medianReadPerCell, meanReadPerCell, nUMIinCells, medianUMIperCell, meanUMIperCell, nGeneInCells, medianGenePerCell, meanGenePerCell, nGeneDetected;
        vector<uint32> nUMIperCell, nReadPerCell, nGenePerCell;
    } filteredCells;
    
    string outputPrefix;
    ofstream *streamTranscriptsOut;
    
    array<vector<uint64>,2> sjAll;
    
    vector<readInfoStruct> readInfo; //corrected CB/UMI information for each read

    SoloFeature(int feTy, Parameters &Pin, Transcriptome &inTrans);
    void processRecords(ReadAlignChunk **RAchunk);
    void collapseUMI(uint32 *rGU, uint32 rN, uint32 &nGenes, uint32 &nUtot, uint32 *umiArray, uint64 cellBarcode);
    void outputResults(bool cellFilterYes);
    void addBAMtags(char *&bam0, uint32 &size0, char* bam1);
    void statsOutput();
    void cellFiltering();
};

#endif
