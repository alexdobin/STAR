#ifndef CODE_extendAlign
#define CODE_extendAlign

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"

struct ExtendAlignResult {//the fields extendAlign sets in a Transcript
    uint extendL;
    intScore maxScore;
    uint nMatch;
    uint nMM;
};

//res gets the values a freshly reset Transcript would hold after the Transcript* overload; returns its result for such a Transcript
bool extendAlign( char* R, char* G, uint rStart, uint gStart, int dR, int dG, uint L, uint Lprev, uint nMMprev, uint nMMmax, double pMMmax, bool extendToEnd, ExtendAlignResult &res );

bool extendAlign( char* R, char* G, uint rStart, uint gStart, int dR, int dG, uint L, uint Lprev, uint nMMprev, uint nMMmax, double pMMmax, bool extendToEnd, Transcript* trA );

inline void addExtension(Transcript &trA, const ExtendAlignResult &ext) {//trA.add() of a reset Transcript holding ext
    trA.maxScore+=ext.maxScore;
    trA.nMatch+=ext.nMatch;
    trA.nMM+=ext.nMM;
};

#endif
