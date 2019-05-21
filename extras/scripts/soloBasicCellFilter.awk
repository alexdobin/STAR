BEGIN {
    maxCells=3000;
    maxPercentile=0.99;
    maxMinRatio=10;
    nHeaderLines=3;
    fOutCB=ARGV[1] ".filtered";
    fOutMat=ARGV[2] ".filtered";
}

{
    if (ARGIND==1) {# read barcodes
        CB[NR]=$1;
    } else if (FNR<=nHeaderLines) {
         A[FNR]=$0;
         if (FNR==nHeaderLines) {
             nGenes=$1;
         };
    } else {
         cellG[FNR]=$1;
         cellI[FNR]=$2;
         cellN[FNR]=$3;
         cellTot[$2]+=$3;
         nLines=FNR;
    };
}
END {
    asort(cellTot,cellTotSort);
    nMax=cellTotSort[-int((1-maxPercentile)*maxCells)+length(cellTot)];
    nMin=nMax/maxMinRatio;

    nCell=0;
    for (ii=1; ii<=length(CB); ii++) {
        if (cellTot[ii]>=nMin) {
            print CB[ii] > fOutCB;
            nCell++;
            cellInew[ii]=nCell;
            print nCell,cellTot[ii] > fOutCB ".counts";
        };
    };

    print cellTotSort[length(cellTotSort)]+0,nMax,nMin,nCell;

    nMat=0;
    for (ii=nHeaderLines+1; ii<=nLines; ii++) {
        if (cellTot[cellI[ii]]>=nMin) {
            nMat++;
        };
    };
    
    for (ii=1;ii<nHeaderLines;ii++) {
        print A[ii] > fOutMat;
    }; 

    print nGenes,nCell+0,nMat+0 > fOutMat;

    for (ii=nHeaderLines+1; ii<=nLines; ii++) {
        if (cellTot[cellI[ii]]>=nMin) {
            print cellG[ii],cellInew[cellI[ii]],cellN[ii] > fOutMat;
        };
    };

};
