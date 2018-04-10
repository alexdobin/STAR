#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "extendAlign.h"

bool extendAlign( char* R, char* G, uint rStart, uint gStart, int dR, int dG, uint L, uint Lprev, uint nMMprev, uint nMMmax, double pMMmax, bool extendToEnd, Transcript* trA ) {

// find the maximum score

int iS,iG;

int Score=0, nMatch=0, nMM=0;
trA->maxScore=0;

R=R+rStart;
G=G+gStart;

if (extendToEnd) {//end to end extension

    int iExt;
    for (iExt=0;iExt<(int) L;iExt++) {
        iS=dR*iExt;
        iG=dG*iExt;

        if ((gStart+iG)==(uint)(-1) || G[iG]==5) {//prohibit extension through chr boundary
            trA->extendL=0;
            trA->maxScore=-999999999;
            trA->nMatch=0;
            trA->nMM=nMMmax+1;
            return true;
//             return false;
        };
        if (R[iS]==MARK_FRAG_SPACER_BASE) break; //no extension through the spacer between fragments

        if (R[iS]>3 || G[iG]>3) continue;//no penalties for Ns in reads or genome

        if (G[iG]==R[iS]) {//Match
            nMatch++;
            Score += scoreMatch;
        } else {
            nMM++;
            Score -= scoreMatch;
        };
    };

    if (iExt>0) {
        trA->extendL=iExt;
        trA->maxScore=Score;
        trA->nMatch=nMatch;
        trA->nMM=nMM;
        return true;
    } else {
        return false;
    };

};


for (int i=0;i<(int) L;i++) {
    iS=dR*i;
    iG=dG*i;

    if ((gStart+iG)==(uint)(-1) || G[iG]==5 || R[iS]==MARK_FRAG_SPACER_BASE) break; //no extension through chr boundary, or through the spacer between fragments
    if (R[iS]>3 || G[iG]>3) continue;//no penalties for Ns in reads or genome

    if (G[iG]==R[iS]) {//Match
        nMatch++;
        Score += scoreMatch;
        if (Score>trA->maxScore) {//record new maximum
            if (nMM+nMMprev <= min(pMMmax*double(Lprev+i+1), double(nMMmax)) ) {//check nMM, if too many mismatches - do not record this maximum. Do not break - there might be still hope to make a long extension
                trA->extendL=i+1;
                trA->maxScore=Score;
                trA->nMatch=nMatch;
                trA->nMM=nMM;
            };
        };
    } else {//MM
        if (nMM+nMMprev >= min(pMMmax*double(Lprev+L), double(nMMmax)) ) {//there is no hope to extend it further, break
            break;
        };

        nMM++;
        Score -= scoreMatch;
    };
};

// decide of the extension worked
bool extDone =  trA->extendL > 0;

return extDone;

};

