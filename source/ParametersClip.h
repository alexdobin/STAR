#ifndef CODE_ParametersClip
#define CODE_ParametersClip

#include "IncludeDefine.h"

class Parameters;

class ClipMate
{
public:
	int type; //0=5p, 1=3p, -1=no clip
	uint32 N;
	uint32 NafterAd;
	string adSeq;
	vector<char> adSeqNum;
	double adMMp;

	void initialize(uint32 Nin, const string &adSeqIn, uint32 afterAdNin, double adMMpIn);
	void clip(uint &Lread, char *SeqNum, uint32 &clippedN);

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

	array<ReadClipInput,2> in;

	vector<vector<ClipMate>> clipMates; //5p,3p; two mates

	void initialize(Parameters *pPin);
        
private:
	Parameters *pP;
};

#endif
