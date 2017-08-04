#include "ChimericSegment.h"

ChimericSegment::ChimericSegment(Parameters &Pin, Transcript &alignIn, uint Lread, uint *readLengthIn) : P(Pin), align(alignIn), readLength(readLengthIn)
{
    if (align.intronMotifs[1]==0 && align.intronMotifs[2]==0) {//strand is undefined
        str=0;
    } else if ( (align.Str==0) == (align.intronMotifs[1]>0)) {//strand the same as RNA
        str=1;
    } else {//strand opposite to RNA
        str=2;
    };
    uint roS=align.Str==0 ? align.exons[0][EX_R] : Lread - align.exons[align.nExons-1][EX_R] - align.exons[align.nExons-1][EX_L];
    uint roE=align.Str==0 ? align.exons[align.nExons-1][EX_R] + align.exons[align.nExons-1][EX_L] - 1 : Lread - align.exons[0][EX_R] - 1;
    if (roS>readLength[0]) roS--;
    if (roE>readLength[0]) roE--;
};