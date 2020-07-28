#include "ParametersClip.h"
#include "Parameters.h"
#include "SequenceFuns.h"

void ClipMate::clip(uint &Lread, char *SeqNum, uint32 &clippedN)
{
	clippedN=0;

	if (type<0)
		return; //no clip for this mate

	uint LreadOld=Lread;

	if (N>0) {//clip N bases
		if (Lread>N) {
			Lread -= N; // for 3p this is all
			clippedN += N;
			if (type==0) {//5p
				memmove(SeqNum, SeqNum+N, Lread);
			};
		} else {
			Lread=0;
			clippedN=LreadOld;
			return;
		}
	};

	uint32 clippedAdN=0;
	if (adSeq.length()>0) {//clip adapter
		if (type==0) {//5p

		} else {//3p
			clippedAdN = Lread-localSearch(SeqNum, Lread, adSeqNum.data(), adSeqNum.size(), adMMp);
		};

		Lread -= clippedAdN;
		clippedN += clippedAdN;
	};

    if (NafterAd>0) {
    	if (Lread > NafterAd) {
    		Lread -= NafterAd;
    		clippedN += NafterAd;
			if (type==0) {//5p
				memmove(SeqNum, SeqNum+NafterAd, Lread);
			};
    	} else {
    		Lread=0;
			clippedN=LreadOld;
			return;
    	};
    };
};
