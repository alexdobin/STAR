#include "ChimericAlign.h"

ChimericAlign::ChimericAlign(ChimericSegment &seg1in, ChimericSegment &seg2in, int chimScoreIn, Genome &genomeIn, ReadAlign *RAin) 
                              : seg1(seg1in), seg2(seg2in),chimScore(chimScoreIn), P(seg1in.P), pCh(P.pCh), mapGen(genomeIn), RA(RAin) {
    stitchingDone=false;
    
    al1=&seg1.align;
    al2=&seg2.align;
    
    if (al1->roStart > al2->roStart) 
        swap (al1,al2);
            
    ex1 = al1->Str==1 ? 0 : al1->nExons-1;
    ex2 = al2->Str==0 ? 0 : al2->nExons-1;
};

bool ChimericAlign::chimericCheck() {
    bool chimGood=true;
    
    chimGood = chimGood && al1->exons[ex1][EX_iFrag] <= al2->exons[ex2][EX_iFrag];//otherwise - strange configuration, both segments contain two mates
        //if ( trChim[0].exons[e0][EX_iFrag] > trChim[1].exons[e1][EX_iFrag] ) {//strange configuration, rare, similar to the next one
        //    chimN=0;//reject such chimeras
            //good test example:
            //CTTAGCTAGCAGCGTCTTCCCAGTGCCTGGAGGGCCAGTGAGAATGGCACCCTCTGGGATTTTTGCTCCTAGGTCT
            //TTGAGGTGAAGTTCAAAGATGTGGCTGGCTGTGAGGAGGCCGAGCTAGAGATCATGGAATTTGTGAATTTCTTGAA
        //} else 
    
    //junction overhangs too short for chimerically spliced mates
    chimGood = chimGood && (al1->exons[ex1][EX_iFrag] < al2->exons[ex2][EX_iFrag] || (al1->exons[ex1][EX_L] >= pCh.junctionOverhangMin &&  al2->exons[ex2][EX_L] >= pCh.junctionOverhangMin) );

    return chimGood;
};
