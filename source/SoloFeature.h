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
#include "ReadAlignChunk.h"

#include "SoloFilteredCells.h"

class SoloFeature {
private:
    Parameters &P;
    ReadAlignChunk **RAchunk;    
    Transcriptome &Trans;
    
    SoloFeature **soloFeatAll;
    
    static const uint32 umiArrayStride=3;
    enum {rguG, rguU, rguR};
    uint32 rguStride;
    
public:
    ParametersSolo &pSolo;

    SoloReadFeature *readFeatSum, **readFeatAll;
    SoloReadBarcode *readBarSum;

    const int32 featureType;   

    uint64 nReadsMapped, nReadsInput; //total number of mapped reads
    uint32 nCB;
    uint32 featuresNumber; //number of features (i.e. genes, SJs, etc)

    uint32 *rGeneUMI;//mapped reads sorted by CB
    uint32 *rCBn;//number of reads for detected CBs in the whitelist
    uint32 **rCBp;//array of pointers to each CB sub-array

    vector<uint32> indCB;//index of detected CBs in the whitelist
    vector<uint32> indCBwl; //reverse of indCB: index of WL CBs in detected CB list
    vector<uint32> nUMIperCB, nUMIperCBsorted;//number of UMIs per CB, and the same sorted (descendant)
    vector<uint32> nGenePerCB, nGenePerCBmulti;//number of genes (with >0 UMIs) per CB
    vector<uint32> nReadPerCB;//number of reads per CB. With multimappers: all aligns per CB
    vector<uint32> nReadPerCBunique, nReadPerCBtotal; //number of unique and multiple reads per CB
    uint32 nReadPerCBmax;

    vector<double> nUMIperCBmulti;

    vector<uint32> countCellGeneUMI;//sparsified matrix for the counts, each entry is: geneID count1 count2 ... countNcounts
    vector<uint32> countCellGeneUMIindex;//index of CBs in the count matrix
    uint32 countMatStride; //number of counts per entry in the count matrix
    
    struct {
        vector<double> m;
        vector<uint32> i;
        uint32 s;
    } countMatMult;
    
    vector<unordered_map<uint32, unordered_set<uint64>>> cbFeatureUMImap; //for SmartSeq counting
       
    string outputPrefix, outputPrefixFiltered;
    
    SoloFilteredCells filteredCells;
    
    array<vector<uint64>,2> sjAll;
    
    vector<readInfoStruct> readInfo; //corrected CB/UMI information for each read
    SoloReadFlagClass readFlagCounts;

    
    vector<uint32> redistrFilesCBindex, redistrFilesCBfirst; //redistr file for each CB, CB boundaries in redistributed files
    vector<uint64> redistrFilesNreads; //number of reads in each file
    vector <fstream*> redistrFilesStreams;

    SoloFeature(Parameters &Pin, ReadAlignChunk **RAchunk, Transcriptome &inTrans, int32 feTy, SoloReadBarcode *readBarSumIn, SoloFeature **soloFeatAll);
    void clearLarge(); //clear large vectors
    void processRecords();
    void sumThreads();
    void countSmartSeq();
    void countCBgeneUMI();
    void countVelocyto();
    void quantTranscript();
    
    void collapseUMI(uint32 iCB, uint32 *umiArray);
    void collapseUMI_CR(uint32 iCB, uint32 *umiArray);
    void collapseUMIall();
    void collapseUMIperCB(uint32 iCB, vector<uint32> &umiArray, vector<uint32> &gID,  vector<uint32> &gReadS);

    uint32 umiArrayCorrect_CR         (const uint32 nU0, uintUMI *umiArr, const bool readInfoRec, const bool nUMIyes, unordered_map <uintUMI,uintUMI> &umiCorr);
    uint32 umiArrayCorrect_Directional(const uint32 nU0, uintUMI *umiArr, const bool readInfoRec, const bool nUMIyes, unordered_map <uintUMI,uintUMI> &umiCorr, const int32 dirCountAdd);
    uint32 umiArrayCorrect_Graph      (const uint32 nU0, uintUMI *umiArr, const bool readInfoRec, const bool nUMIyes, unordered_map <uintUMI,uintUMI> &umiCorr);

    void outputResults(bool cellFilterYes, string outputPrefixMat);
    void addBAMtags(char *&bam0, uint32 &size0, char* bam1);
    void statsOutput();
    void redistributeReadsByCB();
    
    void cellFiltering();
    void emptyDrops_CR();
    void loadRawMatrix();
};

#endif
