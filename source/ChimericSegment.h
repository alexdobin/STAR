#ifndef CODE_ChimericSegment
#define CODE_ChimericSegment

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"

class ChimericSegment
{//
    public:
        Transcript  &align;     //alignment
        uint roS,roE,str;     //start/end/strand in original read coordinates
        Parameters &P;
        ChimericSegment(Parameters &Pin, Transcript &alignIn); //allocate
    private:
};

#endif