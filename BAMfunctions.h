#ifndef DEF_BAMfunctions
#define DEF_BAMfunctions

#include "IncludeDefine.h"
#include SAMTOOLS_BGZF_H
        void outBAMwriteHeader (BGZF* fp, const string &samh, const vector <string> &chrn, const vector <uint> &chrl);
        
#endif