#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "extendAlign.h"
#include <climits>

namespace {

//local extension with steps DR, DG of +-1: nIter already stops where the chromosome-boundary test would break, and the
//"no hope" test breaks at the mismatch where nMM>=nMMstop. The record test is the original expression, evaluated only until
//it holds: with pMMmax>=0 and no wrap in Lprev+i+1 it is monotone in i, so it stays true until nMM changes
template <int DR, int DG>
inline bool extendLocalUnit(const char* R, const char* G, int nIter, uint Lprev, uint nMMprev, uint nMMmax, double pMMmax, int nMMstop, ExtendAlignResult &res)
{
    int Score=0, nMatch=0, nMM=0;
    intScore maxScore=0;
    int extendL=0, bestMatch=0, bestMM=0;
    bool recordOK=false;

    for (int i=0;i<nIter;i++) {
        const char g=G[DG*i];
        if (g==5) break;
        const char r=R[DR*i];
        if (r==MARK_FRAG_SPACER_BASE) break;
        if (r>3 || g>3) continue;

        if (g==r) {//Match
            nMatch++;
            Score += scoreMatch;
            if (Score>maxScore) {
                if (!recordOK)
                    recordOK = nMM+nMMprev <= min(pMMmax*double(Lprev+i+1), double(nMMmax));
                if (recordOK) {
                    extendL=i+1;
                    maxScore=Score;
                    bestMatch=nMatch;
                    bestMM=nMM;
                };
            };
        } else {//MM
            if (nMM>=nMMstop)
                break;
            nMM++;
            Score -= scoreMatch;
            recordOK=false;
        };
    };

    if (extendL==0)
        return false;
    res.extendL=extendL;
    res.maxScore=maxScore;
    res.nMatch=bestMatch;
    res.nMM=bestMM;
    return true;
};

};

bool extendAlign( char* R, char* G, uint rStart, uint gStart, int dR, int dG, uint L, uint Lprev, uint nMMprev, uint nMMmax, double pMMmax, bool extendToEnd, ExtendAlignResult &res ) {

// find the maximum score

res.extendL=0;
res.maxScore=0;
res.nMatch=0;
res.nMM=0;

int iS,iG;

int Score=0, nMatch=0, nMM=0;

R=R+rStart;
G=G+gStart;

if (extendToEnd) {//end to end extension

    int iExt;
    for (iExt=0;iExt<(int) L;iExt++) {
        iS=dR*iExt;
        iG=dG*iExt;

        if ((gStart+iG)==(uint)(-1) || G[iG]==5) {//prohibit extension through chr boundary
            res.extendL=0;
            res.maxScore=-999999999;
            res.nMatch=0;
            res.nMM=nMMmax+1;
            return true;
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
        res.extendL=iExt;
        res.maxScore=Score;
        res.nMatch=nMatch;
        res.nMM=nMM;
        return true;
    } else {
        return false;
    };

};

const int nL=(int) L;
if (nL<=0)
    return false;

const uint kExact=(uint)1<<52; //integers up to here are exact doubles
if ( (dR==1 || dR==-1) && (dG==1 || dG==-1) && pMMmax>=0.0 && Lprev<kExact && nMMprev<kExact ) {
    //nMMstop: smallest nMM with nMM+nMMprev >= noHope, the original "no hope" test, which is monotone in nMM
    const double noHope=min(pMMmax*double(Lprev+L), double(nMMmax));
    int nMMstop=-1;
    if (double(nMMprev)>=noHope) {
        nMMstop=0;
    } else if (noHope<=double(kExact)) {//nMMprev < noHope <= 2^52, both exact: the test is nMM+nMMprev >= ceil(noHope)
        uint c=(uint) noHope;
        if (double(c)<noHope)
            c++;
        c-=nMMprev;
        nMMstop = c<(uint)INT_MAX ? (int)c : INT_MAX;
    };

    if (nMMstop>=0) {//else NaN or a limit beyond 2^52: the general loop below
        const uint iB = dG==1 ? ~gStart : gStart+1; //i at which gStart+dG*i is (uint)(-1), the chromosome-boundary break
        const int nIter = iB<(uint)nL ? (int)iB : nL;
        if (dR==1) {
            if (dG==1)
                return extendLocalUnit<1,1>(R, G, nIter, Lprev, nMMprev, nMMmax, pMMmax, nMMstop, res);
            else
                return extendLocalUnit<1,-1>(R, G, nIter, Lprev, nMMprev, nMMmax, pMMmax, nMMstop, res);
        } else {
            if (dG==1)
                return extendLocalUnit<-1,1>(R, G, nIter, Lprev, nMMprev, nMMmax, pMMmax, nMMstop, res);
            else
                return extendLocalUnit<-1,-1>(R, G, nIter, Lprev, nMMprev, nMMmax, pMMmax, nMMstop, res);
        };
    };
};

int extendL=0, bestMatch=0, bestMM=0;
intScore maxScore=0;

for (int i=0;i<nL;i++) {
    iS=dR*i;
    iG=dG*i;

    if ((gStart+iG)==(uint)(-1) || G[iG]==5 || R[iS]==MARK_FRAG_SPACER_BASE) break; //no extension through chr boundary, or through the spacer between fragments
    if (R[iS]>3 || G[iG]>3) continue;//no penalties for Ns in reads or genome

    if (G[iG]==R[iS]) {//Match
        nMatch++;
        Score += scoreMatch;
        if (Score>maxScore) {//record new maximum
            if (nMM+nMMprev <= min(pMMmax*double(Lprev+i+1), double(nMMmax)) ) {//check nMM, if too many mismatches - do not record this maximum. Do not break - there might be still hope to make a long extension
                extendL=i+1;
                maxScore=Score;
                bestMatch=nMatch;
                bestMM=nMM;
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
if (extendL==0)
    return false;
res.extendL=extendL;
res.maxScore=maxScore;
res.nMatch=bestMatch;
res.nMM=bestMM;
return true;

};

bool extendAlign( char* R, char* G, uint rStart, uint gStart, int dR, int dG, uint L, uint Lprev, uint nMMprev, uint nMMmax, double pMMmax, bool extendToEnd, Transcript* trA ) {
    //the fields of trA are written only where the extension sets them, and the local extension returns trA->extendL>0, as before
    ExtendAlignResult res;
    if (extendAlign(R, G, rStart, gStart, dR, dG, L, Lprev, nMMprev, nMMmax, pMMmax, extendToEnd, res)) {
        trA->extendL=res.extendL;
        trA->maxScore=res.maxScore;
        trA->nMatch=res.nMatch;
        trA->nMM=res.nMM;
        return true;
    };
    trA->maxScore=0;
    return !extendToEnd && trA->extendL>0;
};
