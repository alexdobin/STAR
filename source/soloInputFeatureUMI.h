#ifndef H_soloInputFeatureUMI
#define H_soloInputFeatureUMI

#include <fstream>
#include <array>
#include <vector>
#include "IncludeDefine.h"
#include "SoloCommon.h"

bool soloInputFeatureUMI(fstream *strIn, int32 featureType, bool readInfoYes, array<vector<uint64>,2> &sjAll, uint64 &iread, 
                            int32 &cbmatch, uint32 &feature, uint64 &umi, vector<uint32> &featVecU32, SoloReadFlagClass &readFlagCounts);

#endif
