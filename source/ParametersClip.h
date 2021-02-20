#ifndef CODE_ParametersClip
#define CODE_ParametersClip

#include "IncludeDefine.h"
#include "ClipMate.h"
#include "ClipCR4.h"
#include <array>
#include <vector>

class Parameters;

class ReadClipInput
{
public:
    vector <uint32> N;
    vector <uint32> NafterAd;
    vector <string> adSeq;
    vector <double> adMMp;
};

class ParametersClip
{//
public:
    //bool yes; //trimming is performed

    vector<string> adapterType;
    
    array<ReadClipInput,2> in;

    void initialize(Parameters *pPin);
    void initializeClipMates(vector<vector<ClipMate>> &clipMates);
        
private:
    Parameters *pP;
};

#endif
