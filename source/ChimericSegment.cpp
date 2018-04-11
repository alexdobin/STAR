#include "ChimericSegment.h"

ChimericSegment::ChimericSegment(Parameters &Pin, Transcript &alignIn) : P(Pin), pCh(Pin.pCh), align(alignIn)
{
    if ( (align.intronMotifs[1]==0 && align.intronMotifs[2]==0) || (align.intronMotifs[1]>0 && align.intronMotifs[2]>0)) {//strand is undefined
        str=0;
    } else if ( (align.Str==0) == (align.intronMotifs[1]>0)) {//strand the same as RNA. 
        //This assumes that the aligns have consistent strands, i.e. only intronMotifs[1]>0 OR intronMotifs[2]>0
        str=1;
    } else {//strand opposite to RNA
        str=2;
    };
    roS=align.Str==0 ? align.exons[0][EX_R] : align.Lread - align.exons[align.nExons-1][EX_R] - align.exons[align.nExons-1][EX_L];
    roE=align.Str==0 ? align.exons[align.nExons-1][EX_R] + align.exons[align.nExons-1][EX_L] - 1 : align.Lread - align.exons[0][EX_R] - 1;
    if (roS>align.readLength[0]) roS--;
    if (roE>align.readLength[0]) roE--;
};

bool ChimericSegment::segmentCheck()
{
    bool segGood = true;
    segGood = segGood && align.rLength >= pCh.segmentMin; //mapped length >= chim segmentMin
    segGood = segGood && align.intronMotifs[0]==0; //no non-canonical unannotated juncions. 
    return segGood;
    
        //this is already tested for each align with default --outFilterIntronStrands RemoveInconsistentStrands
        //segGood = segGood && (align.intronMotifs[1]==0 || align.intronMotifs[2]==0); //consistent intron motifs. 
        //this is not requiered since seg2 is tested for length
        //   segGood = segGood && (align.exons[align.nExons-1][EX_R] + align.exons[align.nExons-1][EX_L] + P.pCh.segmentMin <= Lread
        //             || align.exons[0][EX_R] >= P.pCh.segmentMin); //uncovered by seg1 read length is <= segmentMin

};