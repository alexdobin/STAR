#include "SequenceFuns.h"

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

void convertNucleotidesToNumbers(const char* R0, char* R1, uint Lread) {//transform sequence  from ACGT into 0-1-2-3 code    
    for (uint jj=0;jj<Lread;jj++) {
                    switch (int(R0[jj])){
                        case (65): case(97):  R1[jj]=char(0);break;//A
                        case (67): case(99):  R1[jj]=char(1);break;//C           
                        case (71): case(103): R1[jj]=char(2);break;//G                       
                        case (84): case(116): R1[jj]=char(3);break;//T                                
//                         case (78): R1[jj]=char(9);break;//N
                        default:   R1[jj]=char(9);//anything else
                    };
                };
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

uint qualitySplit(char* r, char* q, uint L, char Qsplit, uint maxNsplit, uint  minLsplit, uint** splitR) {
    //splits the read r[L] by quality scores q[L], outputs in splitR - split coordinate/length - per base
    //returns number of good split regions
    uint iR=0,iS=0,iR1,LgoodMin=0, iFrag=0;
    while ( (iR<L) & (iS<maxNsplit) ) { //main cycle
        //find next good base
        while ( (iR<L) && ( (q[iR]<Qsplit) || (r[iR]>3) ) ) {
            if (r[iR]==MARK_FRAG_SPACER_BASE) iFrag++; //count read fragments
            iR++;
        };
        
        if (iR==L) break; //exit when reached end of read
        
        iR1=iR;
        
        //find the next bad base
        while ( iR<L && q[iR]>=Qsplit && r[iR]<=3 ) {
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

