#include "ParametersClip.h"
#include "Parameters.h"
#include "SequenceFuns.h"

void ClipMate::initialize(uint32 Nin, const string &adSeqIn, uint32 NafterAdin, double adMMpIn)
{
	N=Nin;

	adSeq=adSeqIn;
	if (adSeq=="-") {
		adSeq="";
	} else {
        if ( adSeq=="polyA") {
            adSeqNum.clear(); //it should be empty, but just in case...
            adSeqNum.resize(DEF_readSeqLengthMax, 0); //fill with A=0
        } else {
            adSeqNum.resize(adSeq.size(),0);
            convertNucleotidesToNumbers(adSeq.data(), adSeqNum.data(), adSeqNum.size());
        };
	};

	if (N==0 && adSeq=="")
		type=-1;

	NafterAd=NafterAdin;
	adMMp=adMMpIn;
};

void ParametersClip::initialize(Parameters *pPin)
{
	pP=pPin;
	yes=false;
    
    clipMates.resize(min(2LLU, pP->readNmates));//3rd read can only be barcode - do not trim it (for now).

	if ( clipMates.size()>1 ) {//duplicate values for the 2nd mate.
		for (int ip=0; ip<2; ip++) {
			if ( in[ip].N.size()==1 )
				in[ip].N.push_back(in[ip].N[0]);
			if ( in[ip].NafterAd.size()==1 )
				in[ip].NafterAd.push_back(in[ip].NafterAd[0]);
			if ( in[ip].adSeq.size()==1 )
				in[ip].adSeq.push_back(in[ip].adSeq[0]);
			if ( in[ip].adMMp.size()==1 )
				in[ip].adMMp.push_back(in[ip].adMMp[0]);
		};
	};

	for (uint32 im=0; im<clipMates.size(); im++) {
		for (int ip=0; ip<2; ip++) {//fill in the ip
			clipMates[im].resize(2);
			clipMates[im][ip].type=ip;
			clipMates[im][ip].initialize(in[ip].N[im], in[ip].adSeq[im], in[ip].NafterAd[im], in[ip].adMMp[im]);
		};
	};

};
