#include "IncludeDefine.h"
#include "Transcript.h"
// local alignment with Smith-Waterman algorithm
intSWscore alignSmithWaterman(char *R, uint rL, char *G, uint gL, \
        intSWscore pMatch, intSWscore pMismatch, intSWscore pGapOpen, intSWscore pGapExtend, \
        char* T, uint Tsize, Transcript &trA) {

    intSWscore *H=new intSWscore[rL+1];

    uint rL1=rL+1;
    if (rL1*(gL+1)>Tsize) return (intSWscore) 0;

    intSWscore *E=new intSWscore[rL1];

    memset(H,0,sizeof(H[0])*(rL1));
    memset(E,0,sizeof(E[0])*(rL1));



    intSWscore maxH=0;
    uint maxHr=0, maxHg=0;

    for (uint ig=1;ig<=gL;ig++) {//cycle over colums
        intSWscore F=(intSWscore) 0;
        intSWscore HdiagPrev=0;
        for (uint ir=1;ir<=rL;ir++) {//cycle over rows

            E[ir]=max( E[ir]-pGapExtend, H[ir]-pGapOpen );
            E[ir]=max( E[ir], (intSWscore) 0 );

            F = max( F-pGapExtend, H[ir-1]-pGapOpen );
            F = max( F, (intSWscore) 0);

            intSWscore Hdiag = G[ig-1]==R[ir-1] ? HdiagPrev+pMatch : HdiagPrev-pMismatch;

//             if (H[ir]>E[ir] & H[ir]>F)

            HdiagPrev=H[ir];

            if (F>Hdiag && F>E[ir]) {//insertion (gap in read)
                H[ir]=F;
                T[ir+ig*rL1]=1;
            } else if (Hdiag>F && Hdiag>E[ir]) {//match-mismatch
                H[ir]=Hdiag;
                T[ir+ig*rL1]=2;
            } else {//deletion (gap in genome)
                H[ir]=E[ir];
                T[ir+ig*rL1]=3;
            };

            if (H[ir]<0) {
                H[ir]=0;
            };

            if (H[ir]==0) {
                T[ir+ig*rL1]=0;
            };

//             Hdiag=max(Hdiag,E[ir]);
//             Hdiag=max(Hdiag,F);
//             H[ir]=max(Hdiag,(intSWscore) 0);

            if (H[ir]>maxH) {
                maxH=H[ir];
                maxHr=ir;
                maxHg=ig;
            };
            #ifdef DEBUG_SW
                stdOut << setw(2)<<H[ir]<<"  ";
                if (ir==rL) stdOut<<endl;
            #endif
        };
    };

    //traceback the alignment
    //last operation should always be match
    if (T[maxHr+maxHg*rL1]!=2) {
        cerr <<"BUG in SW, last operation is mot match, EXITING"<<endl<<flush;
        exit(1);
    };
    trA.nExons=0;
    uint ig=maxHg, ir=maxHr;
    char prevOper=0;

    trA.nExons=-1;
    while (T[ir+ig*rL1]>0 && ir>0 && ig>0) {
        if (T[ir+ig*rL1]==2) {
            if (prevOper==2) {//increase length
                trA.exons[trA.nExons][EX_L]++;
            } else {//new exon
                ++trA.nExons;
                trA.exons[trA.nExons][EX_L]=1;
                trA.exons[trA.nExons][EX_R]=ir;
                trA.exons[trA.nExons][EX_G]=ig;
                prevOper=2;
            };
            --ir;
            --ig;
        } else if (T[ir+ig*rL1]==1) {//gap in read
            --ir;
            prevOper=1;
        } else {//gap in genome
            --ig;
            prevOper=3;
        };
    };

    ++trA.nExons;
    for (uint ii=0;ii<trA.nExons;ii++) {//subtract length
        trA.exons[ii][EX_R] -= trA.exons[ii][EX_L]; //note that exon loci have extra +1 because T matrix is shifted by +1
        trA.exons[ii][EX_G] -= trA.exons[ii][EX_L];
    };

    for (uint ii=0;ii<trA.nExons/2;ii++) {//reverse order
        for (uint jj=0;jj<EX_SIZE;jj++) {
            swap(trA.exons[ii][jj],trA.exons[trA.nExons-1-ii][jj]);
        };
    };

//     for (uint ii=1;ii<trA.nExons;ii++) {//subtract EX_G of the first exon
//         trA.exons[ii][EX_G] -= trA.exons[0][EX_G];
//     };



    #ifdef DEBUG_SW
            for (uint jj=0;jj<3;jj++) {
                for (uint ii=0;ii<trA.nExons;ii++) {
                    stdOut << trA.exons[jj][trA.nExons-1-ii]<<",";
                };
                stdOut <<"\t";
            };
            stdOut <<endl;
            stdOut << maxHr <<"\t"<< maxHg << endl;
            for (int ir=0;ir<=rL;ir++){
                for (int ig=0;ig<=gL;ig++){
                    stdOut << setw(2) << int(T[ir+ig*rL1]) <<"  ";
                };
                stdOut <<endl;
            };
    #endif
//     uint ig=maxHg, ir=maxHr;
//     uint* rMap=new uint [rL];
//     int H1, F1, E1;
//     while (ir>=1 && ig>=1) {
//
//         rMap[ir-1]=ig;
//
//     };

    return maxH;
};
