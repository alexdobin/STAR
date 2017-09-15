#include "IncludeDefine.h"
#include "Parameters.h"
//#include "Transcript.h"
#include "ReadAlign.h"
#include "SequenceFuns.h"

int ReadAlign::createExtendWindowsWithAlign(uint a1, uint aStr) {


    uint aBin = (a1 >> P.winBinNbits); //align's bin
    uint iBinLeft=aBin, iBinRight=aBin;
    uintWinBin* wB=winBin[aStr];
    uint iBin=-1, iWin=-1, iWinRight=-1;


    if (wB[aBin]==uintWinBinMax) {//proceed if there is no window at this bin
        //check neighboring bins

        bool flagMergeLeft=false;
        if (aBin>0) {//merge left only if there are bins on the left
            for (iBin=aBin-1;  iBin >= ( aBin>P.winAnchorDistNbins ? aBin-P.winAnchorDistNbins : 0 );  --iBin) {//go left, find windows in Anchor range
                if (wB[iBin]<uintWinBinMax) {
                    flagMergeLeft=true;
                    break;
                };
                if (iBin==0) break;
            };
            flagMergeLeft = flagMergeLeft && (mapGen.chrBin[iBin>>P.winBinChrNbits]==mapGen.chrBin[aBin>>P.winBinChrNbits]);
            if (flagMergeLeft) {//this align can be merged into the existing window
                iWin=wB[iBin];
                iBinLeft=WC[iWin][WC_gStart];
                for (uint ii=iBin+1; ii<=aBin; ii++) {//mark al bins with the existing windows ID
                    wB[ii]=iWin;
                };
            };
        };

        bool flagMergeRight=false;
        if (aBin+1<P.winBinN) {//merge left only if there are bins on the right
            for (iBin=aBin+1;  iBin<min(aBin+P.winAnchorDistNbins+1,P.winBinN);  ++iBin) {//go right, find windows in Anchor range
                if (wB[iBin]<uintWinBinMax) {
                    flagMergeRight=true;
                    break;
                };
            };

            flagMergeRight = flagMergeRight && (mapGen.chrBin[iBin>>P.winBinChrNbits]==mapGen.chrBin[aBin>>P.winBinChrNbits]);
            if (flagMergeRight) {//this align can be merged into the existing window
                while (wB[iBin]==wB[iBin+1]) ++iBin; //extend through all bins of the right window
                iBinRight=iBin;
                iWinRight=wB[iBin];
//                 if (iBin!=WC[iWinRight][WC_gEnd]) {//debug, switch off!!!
//                     cerr <<"BUG in createWindows"<<endl<<flush;
//                     exit(0);
//                 };
                if (!flagMergeLeft) iWin=wB[iBin];//left window overwrites right
                for (uint ii=aBin; ii<=iBin; ii++) {//mark al bins with the existing windows ID
                    wB[ii]=iWin;
                };
            };
        };


        if (!flagMergeLeft && !flagMergeRight) {//no merging, a new window was added
            wB[aBin]=iWin=nW; //add new window ID for now, may change it later
            WC[iWin][WC_Chr]=mapGen.chrBin[aBin >> P.winBinChrNbits];
            WC[iWin][WC_Str]=aStr;
            WC[iWin][WC_gEnd]=WC[iWin][WC_gStart]=aBin;
            ++nW;
            if (nW>=P.alignWindowsPerReadNmax) {
                nW=P.alignWindowsPerReadNmax-1;
                return EXIT_createExtendWindowsWithAlign_TOO_MANY_WINDOWS; //too many windows, do not record TODO: record a marker
            };
        } else {//record windows after merging
            WC[iWin][WC_gStart]=iBinLeft;
            WC[iWin][WC_gEnd]=iBinRight;
            if (flagMergeLeft && flagMergeRight) {//kill right window, it was merged with the left one
                WC[iWinRight][WC_gStart]=1;
                WC[iWinRight][WC_gEnd]=0;
            };
        };
    };
    return 0;
};

