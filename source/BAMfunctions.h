#ifndef DEF_BAMfunctions
#define DEF_BAMfunctions

#include "IncludeDefine.h"

#include <htslib/bgzf.h>
#include <htslib/sam.h>

void outBAMwriteHeader (BGZF* fp, const string &samh, const vector <string> &chrn, const vector <uint> &chrl);
int bam_read1_fromArray(char *bamChar, bam1_t *b);
string bam_cigarString (bam1_t *b);

#endif
