#include "SequenceFuns.h"
#include <assert.h>  

void complementSeqNumbers(char* ReadsIn, char* ReadsOut, uint Lread) {//complement the numeric sequences
    for (uint jj=0;jj<Lread;jj++) {//0-3 map to 3-0; every other byte, negative chars included, maps to itself
        unsigned char c=(unsigned char)ReadsIn[jj];
        ReadsOut[jj]=char(c<4 ? 3-c : c);
    };
};

//per-byte tables indexed by the unsigned input byte (a switch per base compiles to a mispredicted indirect jump)
//convertNucleotidesToNumbers: A/a 0, C/c 1, G/g 2, T/t 3, every other byte 4
static const char convNuclTable[256]={
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
    4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
};
//revComplementNucleotides: the IUPAC complement of A C G T N R Y K M S W B D V H and their lower case; every other byte maps to itself
static const unsigned char revComplNuclTable[256]={
    0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
    16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
    32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,
    48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,
    64,84,86,71,72,69,70,67,68,73,74,77,76,75,78,79,
    80,81,89,83,65,85,66,87,88,82,90,91,92,93,94,95,
    96,116,118,103,104,101,102,99,100,105,106,109,108,107,110,111,
    112,113,121,115,97,117,98,119,120,114,122,123,124,125,126,127,
    128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,
    144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,
    160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,
    176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,
    192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,
    208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,
    224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,
    240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,
};

void revComplementNucleotides(char* ReadsIn, char* ReadsOut, uint Lread) {//complement the numeric sequences
    for (uint jj=0;jj<Lread;jj++) {
        ReadsOut[jj]=char(revComplNuclTable[(unsigned char)ReadsIn[Lread-1-jj]]);
    };
};

void revComplementNucleotides(string &seq) {//complement the numeric sequences
    string seq1(seq);
    for (uint jj=0;jj<seq.size();jj++) {
        seq[jj]=char(revComplNuclTable[(unsigned char)seq1[seq.size()-1-jj]]);
    };
};



char  nuclToNumBAM(char cc){
    switch (cc) {//=ACMGRSVTWYHKDBN
        case ('='): cc=0;break;
        case ('A'): case ('a'): cc=1;break;
        case ('C'): case ('c'): cc=2;break;
        case ('M'): case ('m'): cc=3;break;
        case ('G'): case ('g'): cc=4;break;
        case ('R'): case ('r'): cc=5;break;
        case ('S'): case ('s'): cc=6;break;
        case ('V'): case ('v'): cc=7;break;
        case ('T'): case ('t'): cc=8;break;
        case ('W'): case ('w'): cc=9;break;
        case ('Y'): case ('y'): cc=10;break;
        case ('H'): case ('h'): cc=11;break;
        case ('K'): case ('k'): cc=12;break;
        case ('D'): case ('d'): cc=13;break;
        case ('B'): case ('b'): cc=14;break;
        case ('N'): case ('n'): cc=15;break;
        default: cc=15;
    };
    return cc;
};

void nuclPackBAM(char* ReadsIn, char* ReadsOut, uint Lread) {//pack nucleotides for BAM
    for (uint jj=0;jj<Lread/2;jj++) {
        ReadsOut[jj]=nuclToNumBAM(ReadsIn[2*jj])<<4 | nuclToNumBAM(ReadsIn[2*jj+1]);
    };
    if (Lread%2==1) {
        ReadsOut[Lread/2]=nuclToNumBAM(ReadsIn[Lread-1])<<4;
    };
};

void convertNucleotidesToNumbers(const char* R0, char* R1, const uint Lread) {//transform sequence  from ACGT into 0-1-2-3 code
    for (uint jj=0;jj<Lread;jj++) {
        R1[jj]=convNuclTable[(unsigned char)R0[jj]];//anything else is converted to N
    };
};

void convertCapitalBasesToNum(uint8_t *rS, uint64_t N)
{//only capital bases are allowed
    for (uint64_t ib=0; ib<N; ib++) {
        switch (rS[ib]) {
            case 'A':
                rS[ib]=0;
                break;
            case 'C':
                rS[ib]=1;
                break;        
            case 'G':
                rS[ib]=2;
                break;   
            case 'T':
                rS[ib]=3;
                break;   
            default:
                rS[ib]=4;
        };
    };
};

uint convertNucleotidesToNumbersRemoveControls(const char* R0, char* R1, const uint Lread) {//transform sequence  from ACGT into 0-1-2-3 code
    uint iR1=0;
    for (uint jj=0;jj<Lread;jj++) {
        switch (int(R0[jj])){
            case (65): case(97):
                R1[jj]=char(0);break;//A
            case (67): case(99):
                R1[jj]=char(1);break;//C
            case (71): case(103):
                R1[jj]=char(2);break;//G
            case (84): case(116):
                R1[jj]=char(3);break;//T
            default:
                if (int(R0[jj]) < 32) {//control characters are skipped
                    continue;
                } else {//all non-control non-ACGT characters are convreted to N
                    R1[jj]=char(4);//anything else
                };
        };
        ++iR1;
    };
    return iR1;
};


char convertNt01234(const char R0) {//transform sequence  from ACGT into 0-1-2-3 code
    switch(R0)
    {
        case('a'):
        case('A'):
            return 0;
            break;
        case('c'):
        case('C'):
            return 1;
            break;
        case('g'):
        case('G'):
            return 2;
            break;
        case('t'):
        case('T'):
            return 3;
            break;
        default:
            return 4;
    };
};

int32 convertNuclStrToInt32(const string S, uint32 &intOut) {
    intOut=0;
    int32 posN=-1;
    for (uint32 ii=0; ii<S.size(); ii++) {
        uint32 nt = (uint32) convertNt01234(S.at(ii));
        if (nt>3) {//N
            if (posN>=0)
                return -2; //two Ns
            posN=ii;
            nt=0;
        };
        intOut = intOut << 2;
        intOut +=nt;
        //intOut += nt<<(2*ii);
    };
    return posN;
};

string convertNuclInt32toString(uint32 nuclNum, const uint32 L) {
    string nuclOut(L,'N');
    string nuclChar="ACGT";

    for (uint32 ii=1; ii<=L; ii++) {
        nuclOut[L-ii] = nuclChar[nuclNum & 3];
        nuclNum = nuclNum >> 2;
    };

    return nuclOut;
};

int64 convertNuclStrToInt64(const string S, uint64 &intOut) {
    intOut=0;
    int64 posN=-1;
    for (uint64 ii=0; ii<S.size(); ii++) {
        uint64 nt = (uint64) convertNt01234(S[ii]);
        if (nt>3) {//N
            if (posN>=0)
                return -2; //two Ns
            posN=ii;
            nt=0;
        };
        intOut = intOut << 2;
        intOut +=nt;
        //intOut += nt<<(2*ii);
    };
    return posN;
};

string convertNuclInt64toString(uint64 nuclNum, const uint32 L) {
    string nuclOut(L,'N');
    string nuclChar="ACGT";

    for (uint64 ii=1; ii<=L; ii++) {
        nuclOut[L-ii] = nuclChar[nuclNum & 3];
        nuclNum = nuclNum >> 2;
    };

    return nuclOut;
};


uint chrFind(uint Start, uint i2, uint* chrStart) {// find chromosome from global locus
    uint i1=0, i3;
    while (i1+1<i2) {
        i3=(i1+i2)/2;
        if ( chrStart[i3] > Start ) {
            i2=i3;
        } else {
            i1=i3;
        };
    };
    return i1;
};

uint localSearch(const char *x, uint nx, const char *y, uint ny, double pMM){
    //find the best alignment of two short sequences x and y
    //pMM is the maximum percentage of mismatches
    uint nMatch=0, nMM=0, nMatchBest=0, nMMbest=0, ixBest=nx;
    for (uint ix=0;ix<nx;ix++) {
        nMatch=0; nMM=0;
        for (uint iy=0;iy<min(ny,nx-ix);iy++) {
            if (x[ix+iy]>3) continue;
            if (x[ix+iy]==y[iy]) {
                nMatch++;
            } else {
                nMM++;
            };
        };

        if ( ( nMatch>nMatchBest || (nMatch==nMatchBest && nMM<nMMbest) ) && double(nMM)/double(nMatch)<=pMM) {
            ixBest=ix;
            nMatchBest=nMatch;
            nMMbest=nMM;
        };
    };
    return ixBest;
};

uint localSearchNisMM(const char *x, uint nx, const char *y, uint ny, double pMM){
    //find the best alignment of two short sequences x and y
    //pMM is the maximum percentage of mismatches
    //Ns in x OR y are considered mismatches
    uint nMatch=0, nMM=0, nMatchBest=0, nMMbest=0, ixBest=nx;
    for (uint ix=0;ix<nx;ix++) {
        nMatch=0; nMM=0;
        for (uint iy=0;iy<min(ny,nx-ix);iy++) {
            if (x[ix+iy]==y[iy] && y[iy]<4) {
                nMatch++;
            } else {
                nMM++;
            };
        };

        if ( ( nMatch>nMatchBest || (nMatch==nMatchBest && nMM<nMMbest) ) && double(nMM)/double(nMatch)<=pMM) {
            ixBest=ix;
            nMatchBest=nMatch;
            nMMbest=nMM;
        };
    };
    return ixBest;
};

uint32 localAlignHammingDist(const string &text, const string &query, uint32 &pos)
{
    uint32 distBest=query.size();
    if (text.size()<query.size()) {//query is longer than text, no match
    	return text.size()+1;
    };
    for (uint32 ii=0; ii<text.size()-query.size()+1; ii++) {
        uint32 dist1=0;
        for (uint32 jj=0; jj<query.size(); jj++) {
            if (query[jj]!='N' && text[jj+ii]!=query[jj]) {//N in query does not count as mismatch
                ++dist1;
            };
        };
        if (dist1<distBest) {
            distBest=dist1;
            pos=ii;
        };
    };
    return distBest;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
uint32 localSearchGeneral(const char *text, const uint32 textLen, const vector<char> &query, const int32 textStart, const int32 textEnd, double pMM, vector <uint32> vecMM, uint32 &nMM)
{
    assert(textEnd <= (int32)textLen);
    assert(textStart + (int32)query.size() >= 0);

    nMM=0;
    
    uint32 nMatchBest=0;
    int32 posBest=textLen;
    uint32 clippedL = 0;

    
    int32 dirSearch = (textStart <= textEnd ? 1 : -1); //search direction
    
    for (int32 pos=textStart; pos!=textEnd; pos+=dirSearch) {
        int32 qs = max(0, -pos);
        int32 qe = min((uint32)query.size(), (uint32)(textLen-pos) );
                       
        uint32 nMatch1=0, nMM1=0;

        for (uint32 iq=qs; iq<qe; iq++) {
            if (text[pos+iq]>3) 
                continue; //Ns in the text are not counted as matches or mismatches
            if (text[pos+iq]==query[iq]) {
                nMatch1++;
            } else {
                nMM1++;
                if ( nMM1 >= vecMM.size() ) {
                    nMatch1=0;
                    break; //too many mismatches
                };
            };
        };
        
        //if ( (nMatch1>nMatchBest || (nMatch1==nMatchBest && nMM1<nMM)) && double(nMM1)<=double(nMatch1)*pMM ) {
        if ( (nMatch1>nMatchBest || (nMatch1==nMatchBest && nMM1<nMM)) && nMM1<vecMM.size() && (qe-qs)>=vecMM[nMM1]) {
            posBest=pos;
            nMatchBest=nMatch1;
            nMM=nMM1;
            clippedL = (uint32)(textStart <= textEnd ? posBest+(int32)query.size(): -posBest+(int32)textLen );
        };        
    };
        
    return clippedL;
};
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
uint qualitySplit(char* r,uint L, uint maxNsplit, uint  minLsplit, uint** splitR) {
    //splits the read r[L] by quality scores q[L], outputs in splitR - split coordinate/length - per base
    //returns number of good split regions
    uint iR=0,iS=0,iR1,LgoodMin=0, iFrag=0;
    while ( (iR<L) & (iS<maxNsplit) ) { //main cycle
        //find next good base
        while ( iR<L && r[iR]>3 ) {
            if (r[iR]==MARK_FRAG_SPACER_BASE) 
                iFrag++; //count read fragments
            iR++;
        };

        if (iR==L) break; //exit when reached end of read

        iR1=iR;

        //find the next bad base
        while ( iR<L && r[iR]<=3 ) {
            iR++;
        };

        if ( (iR-iR1)>LgoodMin ) LgoodMin=iR-iR1;
        if ( (iR-iR1)<minLsplit ) continue; //too short for a good region

        splitR[0][iS]=iR1;      //good region start
        splitR[1][iS]=iR-iR1;   //good region length
        splitR[2][iS]=iFrag;    //good region fragment
        iS++;
    };

    if (iS==0) splitR[1][0]=LgoodMin; //output min good piece length

    return iS;
};

