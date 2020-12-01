#ifndef CODE_ParametersClip
#define CODE_ParametersClip

#include "IncludeDefine.h"
#include "ClipCR4.h"

class Parameters;

class ClipMate
{
public:
    //clip parameters
	int type; //standard sequence clip: 0=5p, 1=3p, -1=no clip; 10/11 = 10X CR4 5/3p clip
	uint32 N;
	uint32 NafterAd;
	string adSeq;
	vector<char> adSeqNum;
	double adMMp;
    
    //clipped infor from clipChunk
    char clippedInfo;
    
    //clip results
    uint32 clippedAdN;  //adapter bases clipped
    uint32 clippedAdMM; //adapter mismatches
    uint32 clippedN; //total number of bases clipped
    
    ClipCR4 *cr4; //CR4 clipping structure

	void initialize(uint32 Nin, const string &adSeqIn, uint32 afterAdNin, double adMMpIn);
	uint32 clip(uint &Lread, char *SeqNum);
    void clipChunk(char *chArr, uint64 chSize);

};


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
	bool yes; //trimming is performed

    vector<string> adapterType;
    
	array<ReadClipInput,2> in;

	vector<vector<ClipMate>> clipMates; //5p,3p; two mates
	
	void initialize(Parameters *pPin);
        
private:
	Parameters *pP;
};

#endif
