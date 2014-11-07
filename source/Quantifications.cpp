#include "Quantifications.h"

Quantifications::Quantifications (uint32 nGeIn) {

    geneCounts.nGe=nGeIn;
    geneCounts.uStr = new uintQ* [2];

    for (int istr=0; istr<2; istr++) {
        geneCounts.uStr[istr] = new uintQ [geneCounts.nGe];
        for (uint32 ii=0; ii<geneCounts.nGe; ii++) {
            geneCounts.uStr[istr][ii]=0;
        };
    };
};
