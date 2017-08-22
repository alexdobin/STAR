#include "ChimericAlign.h"

void ChimericAlign::chimericStitching(char *genSeq, char *readSeq) {
 
    if (stitchingDone)
        return;
    
    stitchingDone=true;
    
    chimStr = max(seg1.str,seg2.str); //segment strands are either equal, or one is zero - select the non-zero strand

    chimRepeat1=0; chimRepeat2=0; chimJ1=0; chimJ2=0; chimMotif=0;

    if ( al1.exons[ex1][EX_iFrag] < al2.exons[ex2][EX_iFrag] ) {//mates bracket the chimeric junction
        chimMotif=-1;
        if (al1.Str==1) {//negative strand
            chimJ1=al1.exons[ex1][EX_G]-1;
        } else {
            chimJ1=al1.exons[ex1][EX_G]+al1.exons[ex1][EX_L];
        };
        if (al2.Str==0) {//positive strand
            chimJ2=al2.exons[ex2][EX_G]-1;
        } else {
            chimJ2=al2.exons[ex2][EX_G]+al2.exons[ex2][EX_L];
        };
    } else {//chimeric junctions is within one of the mates, check and shift chimeric junction if necessary
            uint roStart0 = al1.Str==0 ? al1.exons[ex1][EX_R] : al1.Lread - al1.exons[ex1][EX_R] - al1.exons[ex1][EX_L];
            uint roStart1 = al2.Str==0 ? al2.exons[ex2][EX_R] : al1.Lread - al2.exons[ex2][EX_R] - al2.exons[ex2][EX_L];

            uint jR, jRbest=0;
            int jScore=0,jMotif=0,jScoreBest=-999999,jScoreJ=0;
            for (jR=0; jR<roStart1+al2.exons[ex2][EX_L]-roStart0-1; jR++) {//scan through the exons to find a canonical junction, and check for mismatches

                if (jR==al1.readLength[0]) jR++; //skip the inter-mate base

                char bR=readSeq[roStart0+jR];

                char b0,b1;
                if (al1.Str==0) {
                    b0=genSeq[al1.exons[ex1][EX_G]+jR];
                } else {
                    b0=genSeq[al1.exons[ex1][EX_G]+al1.exons[ex1][EX_L]-1-jR];
                    if (b0<4) b0=3-b0;
                };

                if (al2.Str==0) {
                    b1=genSeq[al2.exons[ex2][EX_G]-roStart1+roStart0+jR];
                } else {
                    b1=genSeq[al2.exons[ex2][EX_G]+al2.exons[ex2][EX_L]-1+roStart1-roStart0-jR];
                    if (b1<4) b1=3-b1;
                };

                if ( ( P.pCh.filter.genomicN && (b0>3 || b1>3) ) || bR>3) {//chimera is not called if there are Ns in the genome or in the read
                    chimScore=0;
                    return;
                };

                char b01,b02,b11,b12;
                if (al1.Str==0) {
                    b01=genSeq[al1.exons[ex1][EX_G]+jR+1];
                    b02=genSeq[al1.exons[ex1][EX_G]+jR+2];
                } else {
                    b01=genSeq[al1.exons[ex1][EX_G]+al1.exons[ex1][EX_L]-1-jR-1];
                    if (b01<4) b01=3-b01;
                    b02=genSeq[al1.exons[ex1][EX_G]+al1.exons[ex1][EX_L]-1-jR-2];
                    if (b02<4) b02=3-b02;
                };
                if (al2.Str==0) {
                    b11=genSeq[al2.exons[ex2][EX_G]-roStart1+roStart0+jR-1];
                    b12=genSeq[al2.exons[ex2][EX_G]-roStart1+roStart0+jR];
                } else {
                    b11=genSeq[al2.exons[ex2][EX_G]+al2.exons[ex2][EX_L]-1+roStart1-roStart0-jR+1];
                    if (b11<4) b11=3-b11;
                    b12=genSeq[al2.exons[ex2][EX_G]+al2.exons[ex2][EX_L]-1+roStart1-roStart0-jR];
                    if (b12<4) b12=3-b12;
                };

                jMotif=0;
                if (b01==2 && b02==3 && b11==0 && b12==2) {//GTAG
                    if (chimStr!=2) {
                        jMotif=1;
                    };
                } else if(b01==1 && b02==3 && b11==0 && b12==1) {//CTAC
                    if (chimStr!=1) {
                        jMotif=2;
                    };
                };

                if (bR==b0 && bR!=b1) {
                    jScore++;
                } else if (bR!=b0 && bR==b1) {
                    jScore--;
                };

                jScoreJ =jMotif==0 ? jScore +  P.pCh.scoreJunctionNonGTAG : jScore ;

                if ( jScoreJ > jScoreBest || (jScoreJ == jScoreBest && jMotif>0) ) {
                    chimMotif=jMotif;
                    jRbest=jR;
                    jScoreBest=jScoreJ;
                };
            };//jR cycle
            

            //shift junction in trChim
            if (al1.Str==1) {
                al1.exons[ex1][EX_R] +=al1.exons[ex1][EX_L]-jRbest-1;
                al1.exons[ex1][EX_G] +=al1.exons[ex1][EX_L]-jRbest-1;
                al1.exons[ex1][EX_L]=jRbest+1;
                chimJ1=al1.exons[ex1][EX_G]-1;
            } else {
                al1.exons[ex1][EX_L]=jRbest+1;
                chimJ1=al1.exons[ex1][EX_G]+al1.exons[ex1][EX_L];
            };

            if (al2.Str==0) {
                al2.exons[ex2][EX_R] +=roStart0+jRbest+1-roStart1;
                al2.exons[ex2][EX_G] +=roStart0+jRbest+1-roStart1;
                al2.exons[ex2][EX_L]=roStart1+al2.exons[ex2][EX_L]-roStart0-jRbest-1;
                chimJ2=al2.exons[ex2][EX_G]-1;
            } else {
                al2.exons[ex2][EX_L]=roStart1+al2.exons[ex2][EX_L]-roStart0-jRbest-1;
                chimJ2=al2.exons[ex2][EX_G]+al2.exons[ex2][EX_L];
            };
            //find repeats
            char b0,b1;
            for (jR=0;jR<100;jR++) {//forward check
                if (al1.Str==0) {
                    b0=genSeq[chimJ1+jR];
                } else {
                    b0=genSeq[chimJ1-jR];
                    if (b0<4) b0=3-b0;
                };

                if (al2.Str==0) {
                    b1=genSeq[chimJ2+1+jR];
                } else {
                    b1=genSeq[chimJ2-1-jR];
                    if (b1<4) b1=3-b1;
                };
                if (b0!=b1) break;
            };
            chimRepeat2=jR;
            for (jR=0;jR<100;jR++) {//reverse check
                if (al1.Str==0) {
                    b0=genSeq[chimJ1-1-jR];
                } else {
                    b0=genSeq[chimJ1+1+jR];
                    if (b0<4) b0=3-b0;
                };

                if (al2.Str==0) {
                    b1=genSeq[chimJ2-jR];
                } else {
                    b1=genSeq[chimJ2+jR];
                    if (b1<4) b1=3-b1;
                };
                if (b0!=b1) break;
            };
            chimRepeat1=jR;
    };//chimeric junction is within a mate

    if (chimMotif>=0 && (al1.exons[ex1][EX_L]<P.pCh.junctionOverhangMin+chimRepeat1 || al2.exons[ex2][EX_L]<P.pCh.junctionOverhangMin+chimRepeat2) ) {
        //filter out cases where linear junctions that are very close to chimeric junction
        chimScore=0;
        return;
    };
};
        