#include "IncludeDefine.h"
#include "Genome.h"

bool sjAlignSplit(uint a1,uint aLength, const Genome &mapGen, uint &a1D, uint &aLengthD, uint &a1A, uint &aLengthA, uint &isj) {
    uint sj1=(a1-mapGen.sjGstart)%mapGen.sjdbLength;
    if (sj1<mapGen.sjdbOverhang && sj1+aLength>mapGen.sjdbOverhang) {//align crosses the junctions
        isj=(a1-mapGen.sjGstart)/mapGen.sjdbLength;
        aLengthD=mapGen.sjdbOverhang-sj1;
        aLengthA=aLength-aLengthD;
        a1D=mapGen.sjDstart[isj]+sj1;
        a1A=mapGen.sjAstart[isj];
        return true;
    } else {
        return false;
    };
};
