#include "OutSJ.h"
#include "ErrorWarning.h"

OutSJ::OutSJ (uint nSJmax, Parameters &Pin, Genome &genomeIn) : oneSJ(genomeIn), P(Pin), mapGen(genomeIn)  {//do I need P?

    data = new char [oneSJ.dataSize*nSJmax]; //allocate big array of SJ loci and properties
    memset(data,0,oneSJ.dataSize*nSJmax);
    N=0;//initialize the counter
};


int compareSJ(const void* i1, const void* i2) {//compare SJs from the data structure
    uint s1=*( (uint*)i1 );
    uint s2=*( (uint*)i2 );

    if (s1>s2) {
        return 1;
    } else if (s1<s2) {
        return -1;
    } else {
        uint32 g1=*( (uint32*)( (char*)i1 + sizeof(uint)) );
        uint32 g2=*( (uint32*)( (char*)i2 + sizeof(uint)) );
        if (g1>g2) {
            return 1;
        } else if (g1<g2) {
            return -1;
        } else {
            return 0;
        };
    };
};

void OutSJ::collapseSJ() {//collapse junctions. Simple version now: re-sort everything
                                //TODO: sort only unsorted portion, and merge two sorted lists
                                //TODO: if there are many collapsed junctions, add the new junction to them w/o sorting
                                //TODO: stranded version
    //sort by start and gap
    if (N==0) return;
    qsort((void*) data, N, oneSJ.dataSize, compareSJ);
    //collapse
    uint isj1=0; //last collapsed junction
    char* isj1P=data;
    for (uint isj=1;isj<N;isj++) {//cycle through all non-collapsed junctions
        char* isjP=data+isj*oneSJ.dataSize;
        if ( compareSJ( (void*) isjP, (void*) isj1P ) == 0 ) {
            oneSJ.collapseOneSJ(isj1P, isjP, P);
        } else {//originate new junction: copy from isj to isj1+1
            isj1++;
            isj1P=data+isj1*oneSJ.dataSize;
            if (isj!=isj1)
            {
                memcpy(isj1P,isjP,oneSJ.dataSize);
            };
        };
    };
    N=isj1+1;
};

Junction::Junction(Genome &genomeIn) : mapGen(genomeIn) {
};

//////////////////////////////////////////////////// oneJunctionWrite
void Junction::junctionPointer(char* sjPoint, uint isj) {//
    char* d1=sjPoint+isj*dataSize;
    start=(uint*) (d1+startP);
    gap=(uint32*) (d1+gapP);
    strand=d1+strandP;
    motif=d1+motifP;
    annot=d1+annotP;
    countUnique=(uint32*) (d1+countUniqueP);
    countMultiple=(uint32*) (d1+countMultipleP);
    overhangLeft=(uint16*) (d1+overhangLeftP);
    overhangRight=(uint16*) (d1+overhangRightP);
};

void Junction::outputStream(ostream &outStream) {
    uint sjChr=mapGen.chrBin[*start >> mapGen.pGe.gChrBinNbits];
    outStream << mapGen.chrName.at(sjChr) <<"\t"<< *start + 1 - mapGen.chrStart[sjChr] <<"\t"<<*start + *gap - mapGen.chrStart[sjChr] \
            <<"\t"<< int(*strand) <<"\t"<< int(*motif) <<"\t"<< int (*annot) <<"\t"<< *countUnique <<"\t"<< *countMultiple \
            <<"\t"<< *overhangLeft << endl;
};

void Junction::collapseOneSJ(char* isj1P, char* isjP, Parameters& P) {//collapse isj junction into isj1: increase counts in isj1. choose max overhangs, motif, annot
    *(uint32*)(isj1P+countUniqueP)   += *(uint32*)(isjP+countUniqueP);
    *(uint32*)(isj1P+countMultipleP) += *(uint32*)(isjP+countMultipleP);

    if (*(uint16*)(isj1P+overhangLeftP) < *(uint16*)(isjP+overhangLeftP) ) {
        *(uint16*)(isj1P+overhangLeftP) = *(uint16*)(isjP+overhangLeftP) ;
    };
    if (*(uint16*)(isj1P+overhangRightP) < *(uint16*)(isjP+overhangRightP) ) {
        *(uint16*)(isj1P+overhangRightP) = *(uint16*)(isjP+overhangRightP) ;
    };

    if (*(isj1P+motifP) != *(isjP+motifP) ) {
            uint s1=*(uint*)(isj1P+startP);
            uint c1=mapGen.chrBin[ s1 >> mapGen.pGe.gChrBinNbits];
            
            stringstream errOut;
            errOut <<"EXITING becaues of BUG: different motifs for the same junction while collapsing junctions\n" \
                   << mapGen.chrName[c1] <<" "<< s1-mapGen.chrStart[c1]+1 <<" "<<s1-mapGen.chrStart[c1]+1 + *(uint32*)(isj1P+gapP) <<" "<<int(*(char*)(isj1P+motifP)) <<" "<<int(*(char*)(isjP+motifP)) \
                   <<" "<<int(*(char*)(isj1P+annotP)) <<" "<<int(*(char*)(isjP+annotP))<<"\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);\
//         *(isj1P+motifP) = *(isjP+motifP) ;
    };
    if (*(isj1P+annotP) < *(isjP+annotP) ) {
            stringstream errOut;
            errOut <<"EXITING becaues of BUG: different annotation status for the same junction while collapsing junctions:"\
                   <<*(uint*)(isj1P+startP) <<" "<<*(uint32*)(isj1P+gapP) <<" "<<int(*(char*)(isj1P+annotP)) <<" "<<int(*(char*)(isjP+annotP))<<"\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);\

//         *(isj1P+annotP) = *(isjP+annotP) ;
    };

}
