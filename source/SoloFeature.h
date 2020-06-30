#ifndef H_SoloFeature
#define H_SoloFeature

#include <fstream>
#include <unordered_map>
#include <unordered_set>

#include "IncludeDefine.h"
#include "ReadAlignChunk.h"
#include "Transcriptome.h"

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
    uint32 featuresNumber; //number of features (i.e. genes, SJs, etc)

    uint32 *rGeneUMI;//mapped reads sorted by CB
    uint32 *indCB;//index of detected CBs in the whitelist
    vector<uint32> indCBwl; //reverse of indCB: index of WL CBs in detected CB list
    uint32 *rCBn;//number of reads for detected CBs in the whitelist
    uint32 **rCBp;//array of pointers to each CB sub-array

    vector<uint32> nUMIperCB, nUMIperCBsorted;//number of UMIs per CB, and the same sorted (descendant)
    vector<uint32> nGenePerCB;//number of genes (with >0 UMIs) per CB
    vector<uint32> nReadPerCB;//number of reads per CB
    
    vector<uint32> countCellGeneUMI;//sparsified matrix for the counts, each entry is: geneID count1 count2 ... countNcounts
    vector<uint32> countCellGeneUMIindex;//index of CBs in the count matrix
    uint32 countMatStride; //number of counts per entry in the count matrix
    
    vector<unordered_map<uint32, unordered_set<uint64>>> cbFeatureUMImap; //for SmartSeq counting
    
    vector<bool> cellFilterVec;
    struct {
        uint64 nCells, nReadInCells, medianReadPerCell, meanReadPerCell, nUMIinCells, medianUMIperCell, meanUMIperCell, nGeneInCells, medianGenePerCell, meanGenePerCell, nGeneDetected;
        vector<uint32> nUMIperCell, nReadPerCell, nGenePerCell;
    } filteredCells;
    
    string outputPrefix;
    
    array<vector<uint64>,2> sjAll;
    
    vector<readInfoStruct> readInfo; //corrected CB/UMI information for each read
    
    vector<uint32> redistrFilesCBindex, redistrFilesCBfirst; //redistr file for each CB, CB boundaries in redistributed files
    vector<uint64> redistrFilesNreads; //number of reads in each file
    vector <fstream*> redistrFilesStreams;

    SoloFeature(int32 feTy, Parameters &Pin, Transcriptome &inTrans, SoloReadBarcode *readBarSumIn, SoloFeature **soloFeatAll);
    void processRecords(ReadAlignChunk **RAchunk);
    void sumThreads(ReadAlignChunk **RAchunk);
    void countSmartSeq();
    void countCBgeneUMI();
    void countVelocyto();
    void quantTranscript();
    void collapseUMI(uint32 iCB, uint32 *umiArray);
    void outputResults(bool cellFilterYes);
    void addBAMtags(char *&bam0, uint32 &size0, char* bam1);
    void statsOutput();
    void cellFiltering();
    void redistributeReadsByCB();
};

#endif
