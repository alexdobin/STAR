#include "SuffixArrayFuns.h"
#include "PackedArray.h"

inline uint medianUint2(uint a, uint b) {
    // returns (a+b)/2
    return a/2 + b/2 + (a%2 + b%2)/2;
};

uint compareSeqToGenome(Genome &mapGen, char** s2, uint S, uint N, uint L, uint iSA, bool dirR, bool& compRes) {
    // compare s to g, find the maximum identity length
    // s2[0] read sequence; s2[1] complementary sequence
    // S position to start search from in s2[0],s2[1]
    //dirR forward or reverse direction search on read sequence

    register int64 ii;

    uint SAstr=mapGen.SA[iSA];
    bool dirG = (SAstr>>mapGen.GstrandBit) == 0; //forward or reverse strand of the genome
    SAstr &= mapGen.GstrandMask;

    char *g=mapGen.G;

    if (dirR && dirG) {//forward on read, forward on genome
        char* s  = s2[0] + S + L;
        g += SAstr + L;
        for (ii=0;(uint) ii < N-L; ii++)
        {
            if (s[ii]!=g[ii])
            {
                if (s[ii]>g[ii])
                {
                    compRes=true;
                    return ii+L;
                } else
                {
                    compRes=false;
                    return ii+L;
                };
            };
        };
//         if (s[ii]>g[ii]) {compRes=true;} else {compRes=false;};
        return N; //exact match
    } else if (dirR && !dirG) {
        char* s  = s2[1] + S + L;
        g += mapGen.nGenome-1-SAstr - L;
        for (ii=0; (uint) ii < N-L; ii++)
        {
            if (s[ii]!=g[-ii])
            {
                if (s[ii]>g[-ii] || g[-ii]>3)
                {
                    compRes=false;
                    return ii+L;
                } else
                {
                    compRes=true;
                    return ii+L;
                };
            };
        };
        return N;
    } else if (!dirR && dirG) {
        char* s  = s2[1] + S - L;
        g += SAstr + L;
        for (ii=0; (uint) ii < N-L; ii++)
        {
            if (s[-ii]!=g[ii])
            {
                if (s[-ii]>g[ii]) {
                    compRes=true;
                    return ii+L;

                } else
                {
                    compRes=false;
                    return ii+L;
                };
            };
        };
        return N;
    } else {//if (!dirR && !dirG)
        char* s  = s2[0] + S - L;
        g += mapGen.nGenome-1-SAstr - L;
        for (ii=0; (uint) ii < N-L; ii++)
        {
            if (s[-ii]!=g[-ii])
            {
                if (s[-ii]>g[-ii] || g[-ii]>3)
                {
                    compRes=false;
                    return ii+L;
                } else
                {
                    compRes=true;
                    return ii+L;
                };
            };
        };
        return N;
    };
};

uint findMultRange(Genome &mapGen, uint i3, uint L3, uint i1, uint L1, uint i1a, uint L1a, uint i1b, uint L1b, char** s, bool dirR, uint S) {
    // given SA index i3 and identity length L3, return the index of the farthest element with the same length, starting from i1,L1 or i1a,L1a, or i1b,L1b

    bool compRes;

    if (L1<L3) { //search between i1 and i3
        L1b=L1; i1b=i1; i1a=i3;
    }
    else {
        if (L1a<L1) {//search between i1a and i1b, else: search bewtween i1a and i1b
            L1b=L1a; i1b=i1a; i1a=i1;
        };
    };
    while ( (i1b+1<i1a)|(i1b>i1a+1) ) { //L1a is the target length, i1a...i1b is the initial range, i1c,L1c is the value in the middle
        uint i1c=medianUint2(i1a,i1b);
        //uint L1c=identityLength(&g[mapGen.SA[i3]+L1b],&g[mapGen.SA[i1c]+L1b],L3-L1b)+L1b;
        uint L1c=compareSeqToGenome(mapGen,s,S,L3,L1b,i1c,dirR,compRes);
        if (L1c==L3) {
            i1a=i1c;
        }
        else { //L1c<L3, move i1c
            i1b=i1c;L1b=L1c;
        };
    };
    return i1a;
};

uint maxMappableLength(Genome &mapGen, char** s, uint S, uint N, uint i1, uint i2, bool dirR, uint& L, uint* indStartEnd) {
    // find minimum mappable length of sequence s to the genome g with suffix array SA; length(s)=N; [i1 i2] is
    // returns number of mappings (1=unique);range indStartEnd; min mapped length = L
    // binary search in SA space

    bool compRes;

    uint L1,L2,i3,L3,L1a,L1b,L2a,L2b,i1a,i1b,i2a,i2b;

    L1=compareSeqToGenome(mapGen,s,S,N,L,i1,dirR,compRes);
    L2=compareSeqToGenome(mapGen,s,S,N,L,i2,dirR,compRes);

//     L1=identityLength(&s[L],&g[mapGen.SA[i1]]);
//     L2=identityLength(&s[L],&g[mapGen.SA[i2]]);
    L= min(L1,L2);

    L1a=L1;L1b=L1;i1a=i1;i1b=i1;
    L2a=L2;L2b=L2;i2a=i2;i2b=i2;

    i3=i1;L3=L1; //in case i1+1>=i2 an not iteration of the loope below is ever made
    while (i1+1<i2) {//main binary search loop
        i3=medianUint2(i1,i2);
        L3=compareSeqToGenome(mapGen,s,S,N,L,i3,dirR,compRes);

        if (L3==N) break; //found exact match, exit the binary search

        if (compRes) { //move 1 to 3
            if (L3>L1) {
               L1b=L1a; L1a=L1; i1b=i1a; i1a=i1;
            };
            i1=i3;L1=L3;
        }
        else {
            if (L3>L2) { //move 2 to 3
               L2b=L2a; L2a=L2; i2b=i2a; i2a=i2;
            };
            i2=i3;L2=L3;
        };
        L= min(L1,L2);

    };

    if (L3<N) {//choose longest alignment length between L1 and L2
        if (L1>L2) {
            i3=i1;L3=L1;
        } else {
            i3=i2;L3=L2;
        };
    };
    // now i3,L3 is the "best" alignment, i.e. longest length

    // find the range of SA indices in which the identiyLength is the same
    i1=findMultRange(mapGen,i3,L3,i1,L1,i1a,L1a,i1b,L1b,s,dirR,S);
    i2=findMultRange(mapGen,i3,L3,i2,L2,i2a,L2a,i2b,L2b,s,dirR,S);

    L=L3; //output
    indStartEnd[0]=i1;
    indStartEnd[1]=i2;

    return i2-i1+1;
};


int compareRefEnds (Genome &mapGen, uint64 SAstr,  uint64 gInsert, bool strG, bool strR)
{
    if ( strG)
    {// + strand g
       return strR ? (SAstr < gInsert ? 1:-1) : 1;
    } else
    {// - strand g
       return strR ? -1 : ( gInsert==-1LLU ? -1 : ( SAstr < mapGen.nGenome-gInsert ? 1:-1) );
    };
};

uint compareSeqToGenome1(Genome &mapGen, char** s2, uint S, uint N, uint L, uint iSA, bool dirR, uint64 gInsert, int & compRes) {
    // compare s to g, find the maximum identity length
    // s2[0] read sequence; s2[1] complementary sequence
    // S position to start search from in s2[0],s2[1]
    // dirR: strand of the s

    //different treatment of 5 (spacer) in the sequence and genome
    // 5 is allowed in the sequence
    // 5 in the genome is < than 5 in the sequence

    //TODO no need for complementary sequence

    register int64 ii;

    uint SAstr=mapGen.SA[iSA];
    bool dirG = (SAstr>>mapGen.GstrandBit) == 0; //forward or reverse strand of the genome
    SAstr &= mapGen.GstrandMask;
    char *g=mapGen.G;

    if (dirG) {//forward on read, forward on genome
        char* s  = s2[0] + S + L;
        g += SAstr + L;
        for (ii=0;(uint) ii < N-L; ii++)
        {
            if (s[ii]!=g[ii])
            {
                if (s[ii]>g[ii])
                {
                    compRes=1;
                    return ii+L;
                } else
                {
                    compRes=-1;
                    return ii+L;
                };
            } else if (s[ii]==GENOME_spacingChar)
            {//this already implies the s[ii]==g[ii]
                compRes=compareRefEnds (mapGen, SAstr, gInsert, dirG, dirR);
                return ii+L;
            };
        };
//         if (s[ii]>g[ii]) {compRes=true;} else {compRes=false;};
        return N; //exact match
    }
    else
    {
        char* s  = s2[1] + S + L;
        g += mapGen.nGenome-1-SAstr - L;
        for (ii=0; (uint) ii < N-L; ii++)
        {
            if (s[ii]!=g[-ii])
            {
                char s1=s[ii],g1=g[-ii];
                if (s1<4) s1=3-s1;
                if (g1<4) g1=3-g1;
                if (s1>g1) {
                    compRes=1;
                    return ii+L;
                } else
                {
                    compRes=-1;
                    return ii+L;
                };
                break;
            } else if (s[ii]==GENOME_spacingChar)
            {//this already implies the s[ii]==g[ii]
                compRes=compareRefEnds (mapGen, SAstr, gInsert, dirG, dirR);
                return ii+L;
            };
        };
        return N;
    };
};


uint suffixArraySearch1(Genome &mapGen, char** s, uint S, uint N, uint64 gInsert, bool strR, uint i1, uint i2, uint L)
{
    // binary search in SA space
    // s[0],s[1] - query sequence, complementary sequence
    // S - start offset
    // N - sequence length
    // g - genome sequence
    // gInsert - index where the sequence insertion happened
    // SA - suffix array
    // strR - strand of the query sequence
    // i1,i2 = starting indices in SA
    // L - starting length
    // output: first SA index > searched string, i.e. g[SA[index-1]]<s<g[SA[index]]

    int compRes;

    uint L1=compareSeqToGenome1(mapGen,s,S,N,L,i1,strR,gInsert,compRes);
    if (compRes<0)
    {// the sequence is smaller than the first index of the SA, cannot proceed
        L=L1;
        return 0;
    };

    uint L2=compareSeqToGenome1(mapGen, s,S,N,L,i2,strR,gInsert,compRes);
    if (compRes>0)
    {//the sequence is bigger than the last SA index, return a huge number
        L=L2;
        return -2llu;
    };

    L=min(L1,L2);

    uint i3=i1,L3=L1; //in case i1+1>=i2 an not iteration of the loope below is ever made
    while (i1+1<i2) {//main binary search loop
        i3=medianUint2(i1,i2);
        L3=compareSeqToGenome1(mapGen,s,S,N,L,i3,strR,gInsert,compRes);//cannot do this because these sj sequences contains spacers=5
        if (L3==N) {//this should not really happen
            L=N;
            return i3;
//             cerr << "Bug L3==N"<<endl;
//             exit(-1); //found exact match of the whole read length, exit the binary search
        };

        if (compRes>0)
        { //move 1 to 3
            i1=i3;L1=L3;
        } else if (compRes<0)
        {//move 2 to 3
            i2=i3;L2=L3;
        }
        L= min(L1,L2);
    };
    return i2; //index at i2 is always bigger than the sequence
};

uint funCalcSAiFromSA(char* gSeq, PackedArray& gSA, Genome &mapGen, uint iSA, int L, int & iL4)
{
    uint SAstr=gSA[iSA];
    bool dirG = (SAstr>>mapGen.GstrandBit) == 0; //forward or reverse strand of the genome
    SAstr &= mapGen.GstrandMask;
    iL4=-1;
    register uint saind=0;
    if (dirG)
    {
        register uint128 g1=*( (uint128*) (gSeq+SAstr) );
        for (int ii=0; ii<L; ii++)
        {
            register char g2=(char) g1;
            if (g2>3)
            {
                iL4=ii;
                saind <<= 2*(L-ii);
                return saind;
            };
            saind=saind<<2;
            saind+=g2;
            g1=g1>>8;
        };
        return saind;
    } else
    {
        register uint128 g1=*( (uint128*) (gSeq+mapGen.nGenome-SAstr-16) );
        for (int ii=0; ii<L; ii++)
        {
            register char g2=(char) (g1>>(8*(15-ii)));
            if (g2>3)
            {
                iL4=ii;
                saind <<= 2*(L-ii);
                return saind;
            };
            saind=saind<<2;
            saind+=3-g2;
        };
        return saind;
    };

};

int64 funCalcSAi(char *G, uint iL) {
    int64 ind1=0;
    for (uint iL1=0;iL1<=iL;iL1++) {
        uint g=(uint) G[iL1];
        if (g>3) {
            return -ind1;
        } else {
            ind1 <<= 2;
            ind1 += g;
        };
   };
   return ind1;
};
