#include "blocksOverlap.h"

uint blocksOverlap(Transcript &t1, Transcript &t2) {//calculate overlap between blocks of two transcripts

    uint i1=0, i2=0, nOverlap=0;
    while (i1<t1.nExons && i2<t2.nExons) {//scan through all exons
        uint rs1=t1.exons[i1][EX_R];
        uint rs2=t2.exons[i2][EX_R];
        uint re1=t1.exons[i1][EX_R]+t1.exons[i1][EX_L];//1st base after the end
        uint re2=t2.exons[i2][EX_R]+t2.exons[i2][EX_L];
        uint gs1=t1.exons[i1][EX_G];
        uint gs2=t2.exons[i2][EX_G];

        if (rs1>=re2) {//t1 block is on the right to t2, no hope of overlap
            i2++;
        } else if (rs2>=re1) {//t2 block is on the right to t1, no hope of overlap
            i1++;
        } else if (gs1-rs1 != gs2-rs2) {//no overlap
            if (re1>=re2) i2++;//1 is on the right of 2
            if (re2>=re1) i1++;//2 is on the right of 1
        } else {//overlap
            nOverlap += min(re1,re2) - max(rs1,rs2);
            if (re1>=re2) i2++;//1 is on the right of 2
            if (re2>=re1) i1++;//2 is on the right of 1
        };
    };

    //debug
//     uint nO1=0;
//     for (uint ir=0;ir<t1.Lread;ir++) {
//         if (t1.gMap[ir]==t2.gMap[ir] && t1.gMap[ir]>0) nO1++;
//     };
//
//     if (nOverlap!=nO1) {
//         exit(255);
//     };
//

    return nOverlap;
};

