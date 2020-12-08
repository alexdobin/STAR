#ifndef DEF_BAMfunctions
#define DEF_BAMfunctions

#include "IncludeDefine.h"
#include SAMTOOLS_BGZF_H
#include SAMTOOLS_SAM_H
#include "ErrorWarning.h"


void outBAMwriteHeader (BGZF* fp, const string &samh, const vector <string> &chrn, const vector <uint> &chrl);
int bam_read1_fromArray(char *bamChar, bam1_t *b);
string bam_cigarString (bam1_t *b);
        
int reg2bin(int beg, int end);
int bamAttrArrayWrite(int32 attr, const char* tagName, char* attrArray );
int bamAttrArrayWrite(float attr, const char* tagName, char* attrArray );
int bamAttrArrayWrite(char attr, const char* tagName, char* attrArray );
int bamAttrArrayWrite(string &attr, const char* tagName, char* attrArray );
int bamAttrArrayWrite(const vector<char> &attr, const char* tagName, char* attrArray );
int bamAttrArrayWrite(const vector<int32> &attr, const char* tagName, char* attrArray );
int bamAttrArrayWriteSAMtags(string &attrStr, char *attrArray, Parameters &P);

template <class TintType>
TintType bamAttributeInt(const char *bamAux, const char *attrName) {//not tested!!!
    const char *attrStart=strstr(bamAux,attrName);
    if (attrStart==NULL) return (TintType) -1;
    switch (attrStart[2]) {
        case ('c'):
            return (TintType) *(int8_t*)(attrStart+3);
        case ('s'):
            return (TintType) *(int16_t*)(attrStart+3);
        case ('i'):
            return (TintType) *(int32_t*)(attrStart+3);
        case ('C'):
            return (TintType) *(uint8_t*)(attrStart+3);
        case ('S'):
            return (TintType) *(uint16_t*)(attrStart+3);
        case ('I'):
            return (TintType) *(uint32_t*)(attrStart+3);
    };
};

template <typename intType>
int bamAttrArrayWriteInt(intType xIn, const char* tagName, char* attrArray, Parameters &P) {//adapted from samtools
    attrArray[0]=tagName[0];attrArray[1]=tagName[1];
    #define ATTR_RECORD_INT(_intChar,_intType,_intValue) attrArray[2] = _intChar; *(_intType*)(attrArray+3) = (_intType) _intValue; return 3+sizeof(_intType)
    int64 x = (int64) xIn;
    if (x < 0) {
        if (x >= -127) {
            ATTR_RECORD_INT('c',int8_t,x);
        } else if (x >= -32767) {
            ATTR_RECORD_INT('s',int16_t,x);
        } else {
            ATTR_RECORD_INT('i',int32_t,x);
            if (!(x>=-2147483647)) {
                ostringstream errOut;
                errOut <<"EXITING because of FATAL BUG: integer out of range for BAM conversion: "<< x <<"\n";
                errOut <<"SOLUTION: contact Alex Dobin at dobin@cshl.edu\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);
            };
        };
    } else {
        if (x <= 255) {
            ATTR_RECORD_INT('C',uint8_t,x);
        } else if (x <= 65535) {
            ATTR_RECORD_INT('S',uint16_t,x);
        } else {
            ATTR_RECORD_INT('I',uint32_t,x);
            if (!(x<=4294967295)) {
                ostringstream errOut;
                errOut <<"EXITING because of FATAL BUG: integer out of range for BAM conversion: "<< x <<"\n";
                errOut <<"SOLUTION: contact Alex Dobin at dobin@cshl.edu\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);
            };
        };
    };
};
        
#endif