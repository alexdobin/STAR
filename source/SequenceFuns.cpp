#include "SequenceFuns.h"
#include <assert.h>  

void complementSeqNumbers(char* ReadsIn, char* ReadsOut, uint Lread) {//complement the numeric sequences
    for (uint jj=0;jj<Lread;jj++) {
        switch (int(ReadsIn[jj])){
            case (3): ReadsOut[jj]=char(0);break;
            case (2): ReadsOut[jj]=char(1);break;
            case (1): ReadsOut[jj]=char(2);break;
            case (0): ReadsOut[jj]=char(3);break;
            default:  ReadsOut[jj]=ReadsIn[jj];
        };
    };
};

void revComplementNucleotides(char* ReadsIn, char* ReadsOut, uint Lread) {//complement the numeric sequences
    for (uint jj=0;jj<Lread;jj++) {
        switch (ReadsIn[Lread-1-jj]){
            case ('A'): ReadsOut[jj]='T';break;
            case ('C'): ReadsOut[jj]='G';break;
            case ('G'): ReadsOut[jj]='C';break;
            case ('T'): ReadsOut[jj]='A';break;
            case ('N'): ReadsOut[jj]='N';break;
            case ('R'): ReadsOut[jj]='Y';break;
            case ('Y'): ReadsOut[jj]='R';break;
            case ('K'): ReadsOut[jj]='M';break;
            case ('M'): ReadsOut[jj]='K';break;
            case ('S'): ReadsOut[jj]='S';break;
            case ('W'): ReadsOut[jj]='W';break;
            case ('B'): ReadsOut[jj]='V';break;
            case ('D'): ReadsOut[jj]='H';break;
            case ('V'): ReadsOut[jj]='B';break;
            case ('H'): ReadsOut[jj]='D';break;

            case ('a'): ReadsOut[jj]='t';break;
            case ('c'): ReadsOut[jj]='g';break;
            case ('g'): ReadsOut[jj]='c';break;
            case ('t'): ReadsOut[jj]='a';break;
            case ('n'): ReadsOut[jj]='n';break;
            case ('r'): ReadsOut[jj]='y';break;
            case ('y'): ReadsOut[jj]='r';break;
            case ('k'): ReadsOut[jj]='m';break;
            case ('m'): ReadsOut[jj]='k';break;
            case ('s'): ReadsOut[jj]='s';break;
            case ('w'): ReadsOut[jj]='w';break;
            case ('b'): ReadsOut[jj]='v';break;
            case ('d'): ReadsOut[jj]='h';break;
            case ('v'): ReadsOut[jj]='b';break;
            case ('h'): ReadsOut[jj]='d';break;

            default:   ReadsOut[jj]=ReadsIn[Lread-1-jj];
        };
    };
};

void revComplementNucleotides(string &seq) {//complement the numeric sequences
    string seq1(seq);
    for (uint jj=0;jj<seq.size();jj++) {
        switch (seq1[seq.size()-1-jj]){
            case ('A'): seq[jj]='T';break;
            case ('C'): seq[jj]='G';break;
            case ('G'): seq[jj]='C';break;
            case ('T'): seq[jj]='A';break;
            case ('N'): seq[jj]='N';break;
            case ('R'): seq[jj]='Y';break;
            case ('Y'): seq[jj]='R';break;
            case ('K'): seq[jj]='M';break;
            case ('M'): seq[jj]='K';break;
            case ('S'): seq[jj]='S';break;
            case ('W'): seq[jj]='W';break;
            case ('B'): seq[jj]='V';break;
            case ('D'): seq[jj]='H';break;
            case ('V'): seq[jj]='B';break;
            case ('H'): seq[jj]='D';break;

            case ('a'): seq[jj]='t';break;
            case ('c'): seq[jj]='g';break;
            case ('g'): seq[jj]='c';break;
            case ('t'): seq[jj]='a';break;
            case ('n'): seq[jj]='n';break;
            case ('r'): seq[jj]='y';break;
            case ('y'): seq[jj]='r';break;
            case ('k'): seq[jj]='m';break;
            case ('m'): seq[jj]='k';break;
            case ('s'): seq[jj]='s';break;
            case ('w'): seq[jj]='w';break;
            case ('b'): seq[jj]='v';break;
            case ('d'): seq[jj]='h';break;
            case ('v'): seq[jj]='b';break;
            case ('h'): seq[jj]='d';break;

            default:   seq[jj]=seq1[seq.size()-1-jj];
        };
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
                            R1[jj]=char(4);//anything else is converted to N
                    };
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

