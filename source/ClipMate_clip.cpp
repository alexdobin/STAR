#include "ParametersClip.h"
#include "Parameters.h"
#include "SequenceFuns.h"

uint32 ClipMate::clip(uint &Lread, char *SeqNum)
{
	clippedN=0;

	if (type<0)
		return 0; //no clip for this mate

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
		}
	};

	if (adSeq.length()>0) {//clip adapter
		if (type==0) {//5p
            //not implemented yet
		} else {//3p            
            //clippedAdN = Lread-localSearch(SeqNum, Lread, adSeqNum.data(), adSeqNum.size(), adMMp);
			uint64 clippedAdN1 = Lread-localSearch(SeqNum, Lread, adSeqNum.data(), adSeqNum.size(), adMMp);
            
            //clippedAdN = localSearchGeneral(SeqNum, Lread, adSeqNum, 0, Lread, adMMp,   clippedAdMM);
            vector<uint32> vecMM({20, 23, 26, 30, 40, 50, 60, 70, 80, 90});
            clippedAdN = localSearchGeneral(SeqNum, Lread, adSeqNum, 0, Lread, adMMp, vecMM,   clippedAdMM);
            
            Lread=Lread;
		};

		Lread -= clippedAdN;
		clippedN += clippedAdN;
	};

    if (NafterAd>0) {
    	if (Lread > NafterAd) {
    		Lread -= NafterAd;
    		clippedN += NafterAd;
			if (type==0) {//5p. For 3p, no need to move sequence
				memmove(SeqNum, SeqNum+NafterAd, Lread);
			};
    	} else {//0-length after clipping
    		Lread=0;
			clippedN=LreadOld;
    	};
    };
    
    return clippedN ;
};
