#ifndef H_SoloCommon
#define H_SoloCommon

#include <array>
#include <unordered_map>

typedef struct{
    uint64 cb; 
    uint32 umi;
} readInfoStruct;

typedef struct{
    uint32 tr; 
    uint8  type;
} trTypeStruct;

typedef uint32 uintUMI;
typedef uint64 uintCB;
typedef uint32 uintRead;

#define uintUMIbits 32
#define velocytoTypeGeneBits 4
#define velocytoTypeGeneBitShift 28
#define geneMultMark (((uint32)1)<<31)

class SoloReadFlagClass
{
public:
    bool yes=false;
    typedef uint32 typeFlag;
    typeFlag flag=0;
    enum: uint32 {cbMatch, cbPerfect, cbMMunique, cbMMmultiple, genomeU, genomeM, featureU, featureM, exonic, intronic, exonicAS, intronicAS, mito, countedU, countedM, nBits};
    const vector<string> statNames={"cbMatch", "cbPerfect", "cbMMunique", "cbMMmultiple", "genomeU", "genomeM", "featureU", "featureM", "exonic", "intronic", "exonicAS", "intronicAS", "mito", "countedU", "countedM"};

    /* lookup table, probably not efficient
        typeFlag bitMask[nBits], bitInt[nBits];
        
        void SoloReadFlagClass() {
            for (uint32 ibit=0; ibit< nBits; ibit++) {
                bitInt[ii] = ((typeFlag)1) << ibit;
                bitMask[ii] = ~ bitInt[ii];
            };

        };
    */

    unordered_map < uintCB, array<uint64,nBits> > flagCounts;
    array<uint64,nBits> flagCountsNoCB={};

    void setBit(uint32 ibit) {
        flag |= ((typeFlag)1) << ibit;
    };

    typeFlag checkBit(uint32 ibit) {
        return (flag>>ibit) & ((typeFlag)1);
    };

    void countsAdd(uintCB cb) 
    {//adds flag bits to the count for a given cb
        auto cbInserted = flagCounts.insert({cb, {} });
        for (uint32 ibit=0; ibit<nBits; ibit++) {
            (*cbInserted.first).second[ibit] += (uint64) checkBit(ibit);
        };
    };

    void countsAddNoCB()
    {//adds flag bits to the count for no-CB reads
        for (uint32 ibit=0; ibit<nBits; ibit++)
            flagCountsNoCB[ibit] += (uint64) checkBit(ibit);
    };

    void countsAddNoCBarray(array<uint64,nBits> &arrIn)
    {
        for (uint32 ibit=0; ibit<nBits; ibit++)
            flagCountsNoCB[ibit] += arrIn[ibit];
    };

};

#endif