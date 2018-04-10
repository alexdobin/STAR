#include "ReadAlign.h"
#include "SequenceFuns.h"
#include "Stats.h"
#include "serviceFuns.cpp"

int ReadAlign::mapOneRead() {

    #ifdef OFF_BEFORE_SEEDING
        #warning OFF_BEFORE_SEEDING
        nW=0;
        return 0;
    #endif

    revertStrand = false; //the 2nd read is awlays on opposite strand. 1st and 2nd reads have already been reversed.

    if (Lread>0) {
        Nsplit=qualitySplit(Read1[0], Qual1[0], Lread, P.Qsplit, P.maxNsplit, P.minLsplit, splitR);
    } else {
        Nsplit=0;
    };

    resetN(); //reset aligns counters to 0

    //reset/initialize a transcript
    trInit->reset();
    trInit->Chr=0;    trInit->Str=0; trInit->roStr=0;    trInit->cStart=0;     trInit->gLength=0; //to generate nice output of 0 for non-mapped reads
    trInit->iRead=iRead;
    trInit->Lread=Lread;
    trInit->nExons=0;
    trInit->readLengthOriginal=readLengthOriginal;
    trInit->readLengthPairOriginal=readLengthPairOriginal;
    trInit->readLength=readLength;
    trInit->readNmates=readNmates;
    trInit->readName=readName;
    
    trBest=trInit;

    uint seedSearchStartLmax=min(P.seedSearchStartLmax,(uint) (P.seedSearchStartLmaxOverLread*(Lread-1)));
    // align all good pieces
    for (uint ip=0; ip<Nsplit; ip++) {

        uint Nstart = P.seedSearchStartLmax>0 && seedSearchStartLmax<splitR[1][ip] ? splitR[1][ip]/seedSearchStartLmax+1 : 1;
        uint Lstart = splitR[1][ip]/Nstart;
        bool flagDirMap=true;
        for (uint iDir=0; iDir<2; iDir++) {//loop over two directions

            uint Lmapped, L;

            for (uint istart=0; istart<Nstart; istart++) {

//               #ifdef COMPILE_FOR_LONG_READS
//
//               #else
                if (flagDirMap || istart>0) {//check if the 1st piece in reveree direction does not need to be remapped
                    Lmapped=0;
                    while ( istart*Lstart + Lmapped + P.minLmap < splitR[1][ip] ) {//map until unmapped portion is <=minLmap

                        uint Shift = iDir==0 ? ( splitR[0][ip] + istart*Lstart + Lmapped ) : \
                                   ( splitR[0][ip] + splitR[1][ip] - istart*Lstart-1-Lmapped); //choose Shift for forward or reverse

                        //uint seedLength=min(splitR[1][ip] - Lmapped - istart*Lstart, P.seedSearchLmax);
                        uint seedLength=splitR[1][ip] - Lmapped - istart*Lstart;
                        maxMappableLength2strands(Shift, seedLength, iDir, 0, mapGen.nSA-1, L, splitR[2][ip]);//L=max mappable length, unique or multiple
                        if (iDir==0 && istart==0 && Lmapped==0 && Shift+L == splitR[1][ip] ) {//this piece maps full length and does not need to be mapped from the opposite direction
                            flagDirMap=false;
                        };
                        Lmapped+=L;
                    };//while ( istart*Lstart + Lmapped + P.minLmap < splitR[1][ip] )
                };//if (flagDirMap || istart>0)

                if (P.seedSearchLmax>0) {//search fixed length. Not very efficient, need to improve
                    uint Shift = iDir==0 ? ( splitR[0][ip] + istart*Lstart ) : \
                                   ( splitR[0][ip] + splitR[1][ip] - istart*Lstart-1); //choose Shift for forward or reverse
                    uint seedLength = min(P.seedSearchLmax, iDir==0 ? (splitR[0][ip] + splitR[1][ip]-Shift):(Shift+1) );
                    maxMappableLength2strands(Shift, seedLength, iDir, 0, mapGen.nSA-1, L, splitR[2][ip]);//L=max mappable length, unique or multiple
                };


//               #endif
            };//for (uint istart=0; istart<Nstart; istart++)
        };
    };

    #ifdef OFF_AFTER_SEEDING
        #warning OFF_AFTER_SEEDING
        return 0;
    #endif

    if (Lread<P.outFilterMatchNmin) {//read is too short (trimmed too much?)
        mapMarker=MARKER_READ_TOO_SHORT;
        trBest->rLength=0; //min good piece length
        nW=0;
    } else if (Nsplit==0) {//no good pieces
        mapMarker=MARKER_NO_GOOD_PIECES;
        trBest->rLength=splitR[1][0]; //min good piece length
        nW=0;
    } else if (Nsplit>0 && nA==0) {
        mapMarker=MARKER_ALL_PIECES_EXCEED_seedMultimapNmax;
        trBest->rLength=multNminL;
        nW=0;
    } else if (Nsplit>0 && nA>0) {//otherwise there are no good pieces, or all pieces map too many times: read cannot be mapped
//         qsort((void*) PC, nP, sizeof(uint)*PC_SIZE, funCompareUint2);//sort PC by rStart and length
        stitchPieces(Read1, Lread);
    };
    
    return 0;
};
