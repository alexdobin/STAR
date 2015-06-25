#include "IncludeDefine.h"                                                                                                                                                     
#include "Parameters.h"

bool sjAlignSplit(uint a1,uint aLength,Parameters* P, uint &a1D, uint &aLengthD, uint &a1A, uint &aLengthA, uint &isj) {
    uint sj1=(a1-P->sjGstart)%P->sjdbLength;
    if (sj1<P->sjdbOverhang && sj1+aLength>P->sjdbOverhang) {//align crosses the junctions
        isj=(a1-P->sjGstart)/P->sjdbLength;
        aLengthD=P->sjdbOverhang-sj1;
        aLengthA=aLength-aLengthD;
        a1D=P->sjDstart[isj]+sj1;
        a1A=P->sjAstart[isj];
        return true;
    } else {
        return false;
    };
};
