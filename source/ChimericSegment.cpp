#include "ChimericSegment.h"

ChimericSegment::ChimericSegment(Parameters &Pin, Transcript &alignIn) : P(Pin), align(alignIn)
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