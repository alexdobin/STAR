#include "ParametersClip.h"
#include "Parameters.h"
#include "SequenceFuns.h"

uint32 ClipMate::clip(uint &Lread, char *seqNum)
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
				memmove(seqNum, seqNum+N, Lread);
			};
		} else {
			Lread=0;
			clippedN=LreadOld;
		}
	};

	if (adSeq.length()>0) {//clip adapter
        switch (type) {
            /* not implemented yet
            case 0: {//5p - not tested yet
                vector<uint32> vecMM({20, 22, 24, 30, 40, 50, 60, 70, 80, 90});
                clippedAdN = localSearchGeneral(seqNum, Lread, adSeqNum, -(int32)adSeqNum.size()+1, (int32)Lread-(int32)adSeqNum.size(), adMMp, vecMM,   clippedAdMM);     
                memmove(seqNum, seqNum+clippedAdN, Lread-clippedAdN);
                break;
            };
            */
            case 1: {//3p            
                clippedAdN = Lread-localSearch(seqNum, Lread, adSeqNum.data(), adSeqNum.size(), adMMp);
                break;
                /* new way, not tested properly
                uint64 clippedAdN1 = Lread-localSearch(seqNum, Lread, adSeqNum.data(), adSeqNum.size(), adMMp);
                
                //clippedAdN = localSearchGeneral(seqNum, Lread, adSeqNum, 0, Lread, adMMp,   clippedAdMM);
                vector<uint32> vecMM({20, 23, 26, 30, 40, 50, 60, 70, 80, 90});
                clippedAdN = localSearchGeneral(seqNum, Lread, adSeqNum, Lread-1, -1, adMMp, vecMM,   clippedAdMM);
                
                Lread=Lread;
                */
            };
            case 10: {//5p: CR4
                clippedAdN = min( (uint32)clippedInfo, (uint32)Lread );
                memmove(seqNum, seqNum+clippedAdN, Lread-clippedAdN);
                break;
            };
            case 11: {//3p: CR4, polyA
                clippedAdN = cr4->polyTail3p(seqNum, Lread);
            };
		};

		Lread -= clippedAdN;
		clippedN += clippedAdN;
	};

    if (NafterAd>0) {
    	if (Lread > NafterAd) {
    		Lread -= NafterAd;
    		clippedN += NafterAd;
			if (type==0) {//5p. For 3p, no need to move sequence
				memmove(seqNum, seqNum+NafterAd, Lread);
			};
    	} else {//0-length after clipping
    		Lread=0;
			clippedN=LreadOld;
    	};
    };
    
    return clippedN;
};
