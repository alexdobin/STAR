#include "ClipMate.h"
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
    
    if (type==10)
        cr4 = new ClipCR4;

    NafterAd=NafterAdin;
    adMMp=adMMpIn;
};


