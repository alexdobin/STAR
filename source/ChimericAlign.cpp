#include "ChimericAlign.h"

ChimericAlign::ChimericAlign(ChimericSegment &seg1in, ChimericSegment &seg2in, int chimScoreIn) : seg1(seg1in), seg2(seg2in), P(seg1in.P), chimScore(chimScoreIn)
{
    stitchingDone=false;
};