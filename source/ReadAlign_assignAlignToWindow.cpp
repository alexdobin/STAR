#include "IncludeDefine.h"
#include "Parameters.h"
#include "ReadAlign.h"
#include "ErrorWarning.h"

void ReadAlign::assignAlignToWindow(uint a1, uint aLength, uint aStr, uint aNrep, uint aFrag, uint aRstart, bool aAnchor, uint sjA) {

    uint iW=winBin[aStr][a1>>P->winBinNbits];

    if (iW==uintWinBinMax || (!aAnchor && aLength < WALrec[iW]) ) return; //alignment does not belong to any window, or it's shorter than rec-length

    
    //check if this alignment overlaps with any other alignment in the window, record the longest of the two
    {//do not check for overlap if this is an sj-align
        uint iA;
        for (iA=0; iA<nWA[iW]; iA++) {
            if (aFrag==WA[iW][iA][WA_iFrag] && WA[iW][iA][WA_sjA]==sjA \
                && a1+WA[iW][iA][WA_rStart]==WA[iW][iA][WA_gStart]+aRstart \
                && ( (aRstart>=WA[iW][iA][WA_rStart] && aRstart<WA[iW][iA][WA_rStart]+WA[iW][iA][WA_Length]) \
                  || (aRstart+aLength>=WA[iW][iA][WA_rStart] && aRstart+aLength<WA[iW][iA][WA_rStart]+WA[iW][iA][WA_Length]) ) ) {//this piece overlaps with iA
                break;
            };
        };
        if (iA<nWA[iW]) {//found overlap
            if (aLength>WA[iW][iA][WA_Length]) {//replace 
                
                uint iA0;//iA0 is where the align has to be inserted
                for (iA0=0;iA0<nWA[iW];iA0++) 
                {//find the insertion point TODO binary search
                    if (iA0!=iA && aRstart<WA[iW][iA0][WA_rStart]) 
                    {//do not compare with the piece to be removed
                        break;
                    };
                };
                
                if (iA0>iA)
                {//true insertion place since iA will be removed
                    --iA0;
                };
                
                if (iA0<iA) {//shift aligns down
                    for (uint iA1=iA;iA1>iA0;iA1--) {//shift aligns to free up insertion point
                        for (uint ii=0;ii<WA_SIZE;ii++) {
                                WA[iW][iA1][ii]=WA[iW][iA1-1][ii];
                        };
                    };
                } else if (iA0>iA) {//shift aligns up
                    for (uint iA1=iA;iA1<iA0;iA1++) {//shift aligns to free up insertion point
                        for (uint ii=0;ii<WA_SIZE;ii++) {
                                WA[iW][iA1][ii]=WA[iW][iA1+1][ii];
                        };
                    };
                };
                
                
                WA[iW][iA0][WA_rStart]=aRstart;
                WA[iW][iA0][WA_Length]=aLength;
                WA[iW][iA0][WA_gStart]=a1;
                WA[iW][iA0][WA_Nrep]=aNrep;                
                WA[iW][iA0][WA_Anchor]=int(aAnchor);//=0 if not, =1 if yes
                WA[iW][iA0][WA_iFrag]=aFrag;
                WA[iW][iA0][WA_sjA]=sjA;
                
            };
            return; //do not record new align
        };
    };          

    if (nWA[iW]==P->seedPerWindowNmax) {//too many aligns per window,  re-calculate min-length, remove the shortest one,

        WALrec[iW]=Lread+1; 
        for (uint iA=0; iA<nWA[iW]; iA++) {//find the new min-length    
            if (WA[iW][iA][WA_Anchor]!=1) WALrec[iW]=min(WALrec[iW],WA[iW][iA][WA_Length]); //protect the anchors - they are not counted for min-length
        };


        if (WALrec[iW]==Lread+1) {//this could happen if there are too many anchors
            mapMarker=MARKER_TOO_MANY_ANCHORS_PER_WINDOW;                    
            nW=0;
            return;
        };


        if (!aAnchor && aLength < WALrec[iW]) return; //alignment is shorter than min-length, do not record - unless it's an anchor                

        uint iA1=0;
        for (uint iA=0; iA<nWA[iW]; iA++) {//remove the shortest aligns
            if ( WA[iW][iA][WA_Anchor]==1 || WA[iW][iA][WA_Length] > WALrec[iW] ) {//re-record the anchors and long aligns
                for (uint ii=0; ii<WA_SIZE; ii++) WA[iW][iA1][ii]=WA[iW][iA][ii]; //re-record the iA-th alignment into iA1-th place
                iA1++;                     
            };
        };

        nWA[iW]=iA1;

        if (!aAnchor && aLength <= WALrec[iW]) {//current align was removed, zero out its nWAP
            nWAP[iW]=0;
        };

    };

    if ( aAnchor || aLength > WALrec[iW] ) {
        if (nWA[iW]>=P->seedPerWindowNmax) {
            exitWithError("BUG: iA>=P->seedPerWindowNmax in stitchPieces, exiting",std::cerr, P->inOut->logMain, EXIT_CODE_BUG, *P);            
        };
     
        uint iA;                                      
        for (iA=0; iA<nWA[iW]; iA++) {//find the insertion point in case aligns are not sorted by aRstart
                                    //TODO binary search
            if (aRstart<WA[iW][iA][WA_rStart]) break;
        };
        for (uint iA1=nWA[iW];iA1>iA;iA1--) {//shift aligns for free up insertion point
            for (uint ii=0;ii<WA_SIZE;ii++) {
                    WA[iW][iA1][ii]=WA[iW][iA1-1][ii];
            };
        };
        
        // now iW is the window to which this align belongs, record it

        // This is the piece that breaks
//         bt full
//#0  0x000000000044e64f in ReadAlign::assignAlignToWindow (this=0x61c494e0, a1=5768917221, aLength=16, aStr=0, aNrep=232, aFrag=0, aRstart=3906, aAnchor=false, sjA=18446744073709551615)
//    at ReadAlign_assignAlignToWindow.cpp:118
//        iA = 0
//        iW = 10000
//#1  0x000000000043e926 in ReadAlign::stitchPieces (this=0x61c494e0, R=0x647016d0, Q=0x6486fa80, 
//    G=0x2aaaad6810d8 "\003\001\003\003\002\003\001\001\003\003\003\002\003\002\003\003\002\001\003\003\001\002\003\003\001\003", SA=..., Lread=8204) at ReadAlign_stitchPieces.cpp:180
//        a1 = 5768917221
//        aStr = 0
//        aRstart = 3906
//        iSA = 2474311944
//        aFrag = 0
//        aLength = 16
//        aDir = 1
//        aNrep = 232
//        aAnchor = false
//        iP = 2030
//        swWinCovMax = 2
//        iW1 = 2
//        trNtotal = 2
//#2  0x00000000004415af in ReadAlign::mapOneRead (this=0x61c494e0) at ReadAlign_mapOneRead.cpp:102
//        seedSearchStartLmax = 50
//#3  0x000000000044f278 in ReadAlign::oneRead (this=0x61c494e0) at ReadAlign_oneRead.cpp:70
//        readStatus = {1, 0}
//#4  0x00000000004453f4 in ReadAlignChunk::mapChunk (this=0x61c49230) at ReadAlignChunk_mapChunk.cpp:25
//        readStatus = 0
//#5  0x0000000000444b44 in ReadAlignChunk::processChunks (this=0x61c49230) at ReadAlignChunk_processChunks.cpp:144
//        __FUNCTION__ = "processChunks"
//#6  0x000000000047475f in ThreadControl::threadRAprocessChunks (RAchunk=0x61c49230) at ThreadControl.h:22
//No locals.
//#7  0x00002aaaabacc851 in start_thread () from /lib64/libpthread.so.0
//No symbol table info available.
//#8  0x00002aaaabdca90d in clone () from /lib64/libc.so.6
//No symbol table info available.
        if (iW >=P->seedPerWindowNmax ) {
            exitWithError("BUG: iW>=P->seedPerWindowNmax in stitchPieces, exiting",std::cerr, P->inOut->logMain, EXIT_CODE_BUG, *P);
        }
        exitWithError("BUG: iW >=P->seedPerWindowNmax in assignAlignToWindow, exiting",std::cerr, P->inOut->logMain, EXIT_CODE_BUG, *P);
        WA[iW][iA][WA_rStart]=aRstart;
        WA[iW][iA][WA_Length]=aLength;
        WA[iW][iA][WA_gStart]=a1;
        WA[iW][iA][WA_Nrep]=aNrep;                
        WA[iW][iA][WA_Anchor]=int(aAnchor);//=0 if not, =1 if yes
        WA[iW][iA][WA_iFrag]=aFrag;
        WA[iW][iA][WA_sjA]=sjA;

        nWA[iW]++;
        nWAP[iW]++;
        if (aAnchor && WlastAnchor[iW]<iA) {
            WlastAnchor[iW]=iA; //record the index of the last anchor
        }
    };
};
