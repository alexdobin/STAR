#include "ReadAlign.h"
#include "SuffixArrayFuns.h"
#include "ErrorWarning.h"

uint ReadAlign::maxMappableLength2strands(uint pieceStartIn, uint pieceLengthIn, uint iDir, uint iSA1, uint iSA2, uint& maxLbest, uint iFrag) {
    //returns number of mappings, maxMappedLength=mapped length
    uint Nrep=0, indStartEnd[2], maxL;

    uint NrepAll[P->genomeSAsparseD], indStartEndAll[P->genomeSAsparseD][2], maxLall[P->genomeSAsparseD];
    maxLbest=0;

    bool dirR = iDir==0;

    for (uint iDist=0; iDist<min(pieceLengthIn,P->genomeSAsparseD); iDist++) {//cycle through different distances
        uint pieceStart;
        uint pieceLength=pieceLengthIn-iDist;

        //calculate full index
        uint Lmax=min(P->genomeSAindexNbases,pieceLength);
        uint ind1=0;
        if (dirR) {//forward search
            pieceStart=pieceStartIn+iDist;
            for (uint ii=0;ii<Lmax;ii++) {//calculate index TODO: make the index calculation once for the whole read and store it
                ind1 <<=2LLU;
                ind1 += ((uint) Read1[0][pieceStart+ii]);
            };
        } else {//reverse search
            pieceStart=pieceStartIn-iDist;
            for (uint ii=0;ii<Lmax;ii++) {//calculate index TODO: make the index calculation once for the whole read and store it
                ind1 <<=2LLU;
                ind1 += ( 3-((uint) Read1[0][pieceStart-ii]) );
            };
        };

        //find SA boundaries
        uint Lind=Lmax;
        while (Lind>0) {//check the precense of the prefix for Lind
            iSA1=SAi[P->genomeSAindexStart[Lind-1]+ind1];
            if ((iSA1 & P->SAiMarkAbsentMaskC) == 0) {//prefix exists
                break;
            } else {//this prefix does not exist, reduce Lind
                --Lind;
                ind1 = ind1 >> 2;
            };
        };

        if (P->genomeSAindexStart[Lind-1]+ind1+1 < P->genomeSAindexStart[Lind]) {//we are not at the end of the SA
            iSA2=((SAi[P->genomeSAindexStart[Lind-1]+ind1+1] & P->SAiMarkNmask) & P->SAiMarkAbsentMask) - 1;
        } else {
            iSA2=P->nSA-1;
        };


    //#define SA_SEARCH_FULL

    #ifdef SA_SEARCH_FULL
        //full search of the array even if the index search gave maxL
        maxL=0;
        Nrep = maxMappableLength(Read1, pieceStart, pieceLength, G, SA, iSA1 & P->SAiMarkNmask, iSA2, dirR, maxL, indStartEnd, P);
    #else
        if (Lind < P->genomeSAindexNbases && (iSA1 & P->SAiMarkNmaskC)==0 ) {//no need for SA search
            indStartEnd[0]=iSA1;
            indStartEnd[1]=iSA2;
            Nrep=indStartEnd[1]-indStartEnd[0]+1;
            maxL=Lind;
        } else if (iSA1==iSA2) {//unique align already, just find maxL
            if ((iSA1 & P->SAiMarkNmaskC)!=0) {
                ostringstream errOut;
                errOut  << "BUG: in ReadAlign::maxMappableLength2strands";
                exitWithError(errOut.str(), std::cerr, P->inOut->logMain, EXIT_CODE_BUG, *P);
            };
            indStartEnd[0]=indStartEnd[1]=iSA1;
            Nrep=1;
            bool comparRes;
            maxL=compareSeqToGenome(Read1,pieceStart, pieceLength, Lind, G, SA, iSA1, dirR, comparRes, P);
        } else {//SA search, pieceLength>maxL
        if ( (iSA1 & P->SAiMarkNmaskC)==0 ) {//no N in the prefix
                maxL=Lind;
            } else {
                maxL=0;
            };
            Nrep = maxMappableLength(Read1, pieceStart, pieceLength, G, SA, iSA1 & P->SAiMarkNmask, iSA2, dirR, maxL, indStartEnd, P);
        };
    #endif

        if (maxL+iDist > maxLbest) {//this idist is better
            maxLbest=maxL+iDist;
        };
        NrepAll[iDist]=Nrep;
        indStartEndAll[iDist][0]=indStartEnd[0];
        indStartEndAll[iDist][1]=indStartEnd[1];
        maxLall[iDist]=maxL;
    };

    for (uint iDist=0; iDist<min(pieceLengthIn,P->genomeSAsparseD); iDist++) {//cycle through different distances, store the ones with largest maxL
        if ( (maxLall[iDist]+iDist) == maxLbest) {
            storeAligns(iDir, (dirR ? pieceStartIn+iDist : pieceStartIn-iDist), NrepAll[iDist], maxLall[iDist], indStartEndAll[iDist], iFrag);
        };
    };
    return Nrep;
};
