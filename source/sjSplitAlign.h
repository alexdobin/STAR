#ifndef SJSPLITALIGN_H
#define SJSPLITALIGN_H

#include "IncludeDefine.h"                                                                                                                                                     
#include "Parameters.h"

#ifdef __cplusplus
extern "C" {
#endif

bool sjAlignSplit(uint a1, uint aLength, Parameters* P, uint &a1D, uint &aLengthD, uint &a1A, uint &aLengthA, uint &isj); 

#ifdef __cplusplus
}
#endif

#endif SJSPLITALIGN_H