#include "Transcriptome.h"
#include "ReadAlign.h"
#include "serviceFuns.cpp"
#include "Quantifications.h"

void Transcriptome::geneCounts(Transcript &a, Quantifications &q) {

    //search end of alignment among the starts of exons
    uint64 g1=a.exons[a.nExons][EX_G];//end of the alignment
    uint64 ex1=binarySearch1<uint64>(g1, exG.s, exG.nEx);


    
};
