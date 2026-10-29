#include "ReadAlign.h"
#include "SequenceFuns.h"
#include "Stats.h"
#include "serviceFuns.cpp"
#include "SuffixArrayFuns.h"

//The first seed search of each chain in mapOneRead's loops does not depend on the others, and each is three dependent
//misses (SAi, SA, genome byte); loading one level for all of them before the loops overlaps the misses.
//Prefetches and reads only: the searches are unchanged.
#if defined(__GNUC__)
//always_inline: GCC 13 deletes calls to an out-of-line function whose only effects are prefetches
__attribute__((always_inline))
#endif
static inline void prefetchChainFirstSeeds(Genome &mapGen, Parameters &P, char *R, uint **splitR, uint Nsplit, uint seedSearchStartLmax)
{
    const uint nMax=32;
    uint saiInd[nMax], Lseed[nMax], iSAa[nMax], iSAb[nMax], iSAc[nMax];
    uint n=0;
    const uint Lsai=P.pGe.gSAindexNbases;

    //the chains' first seeds, as the loops compute them with Lmapped=0 (iDir=1, istart=0 is included though it runs only when
    //the forward search does not map the whole piece)
    for (uint ip=0; ip<Nsplit && n<nMax; ip++) {
        uint Nstart = P.seedSearchStartLmax>0 && seedSearchStartLmax<splitR[1][ip] ? splitR[1][ip]/seedSearchStartLmax+1 : 1;
        uint Lstart = splitR[1][ip]/Nstart;
        for (uint iDir=0; iDir<2; iDir++) {
            for (uint istart=0; istart<Nstart && n<nMax && istart*Lstart + P.seedMapMin < splitR[1][ip]; istart++) {
                uint Shift = iDir==0 ? ( splitR[0][ip] + istart*Lstart ) : ( splitR[0][ip] + splitR[1][ip] - istart*Lstart-1);
                uint Lmax=min(Lsai, splitR[1][ip] - istart*Lstart);
                uint ind1=0;
                if (iDir==0) {
                    for (uint ii=0;ii<Lmax;ii++)
                        ind1 = (ind1<<2) + ((uint) R[Shift+ii]);
                } else {
                    for (uint ii=0;ii<Lmax;ii++)
                        ind1 = (ind1<<2) + ( 3-((uint) R[Shift-ii]) );
                };
                saiInd[n]=mapGen.genomeSAindexStart[Lmax-1]+ind1;
                Lseed[n]=Lmax;
                prefetchPackedWords(mapGen.SAi, saiInd[n], saiInd[n]+1);
                ++n;
            };
        };
    };

    //prefetch targets for a present, N-free full-length prefix (others are skipped): its SA range, which ends before the next prefix's start
    const uint saiEnd=mapGen.genomeSAindexStart[Lsai];
    uint m=0;
    for (uint k=0; k<n; k++) {
        if (Lseed[k]!=Lsai || saiInd[k]+1>=saiEnd)
            continue;
        uint i1=mapGen.SAi[saiInd[k]];
        if ( (i1 & (mapGen.SAiMarkAbsentMaskC | mapGen.SAiMarkNmaskC)) != 0 )
            continue;
        uint i2=mapGen.SAi[saiInd[k]+1];
        if ( (i2 & mapGen.SAiMarkAbsentMaskC) == 0 ) {
            i2 = (i2 & mapGen.SAiMarkNmask) - 1;
        } else if ( (mapGen.SAi[saiEnd-1] & mapGen.SAiMarkAbsentMaskC) == 0 ) {//an absent next prefix holds the next present prefix's start, trusted only when the level ends present; maxMappableLength2strands searches to nSA-1 here, so this range only aims the prefetches
            i2 = (i2 & mapGen.SAiMarkNmask & mapGen.SAiMarkAbsentMask) - 1;
        } else {
            continue;
        };
        prefetchPackedWords(mapGen.SA, i1, i1);
        if (i2!=i1)
            prefetchPackedWords(mapGen.SA, i2, i2);
        uint i3=i1;
        if (i1+1<i2) {
            i3=i1/2 + i2/2 + (i1%2 + i2%2)/2;
            prefetchPackedWords(mapGen.SA, i3, i3);
        };
        iSAa[m]=i1;
        iSAb[m]=i2;
        iSAc[m]=i3;
        ++m;
    };

    //the genome byte where each entry compare (both range ends) and the first midpoint's compare start; the midpoint's
    //starts at min(L1,L2)>=Lsai, usually in the same cache line
    for (uint k=0; k<m; k++) {
        prefetchSuffixStart(mapGen, iSAa[k], Lsai);
        if (iSAb[k]!=iSAa[k])
            prefetchSuffixStart(mapGen, iSAb[k], Lsai);
        if (iSAc[k]!=iSAa[k])
            prefetchSuffixStart(mapGen, iSAc[k], Lsai);
    };
};

int ReadAlign::mapOneRead() {

    #ifdef OFF_BEFORE_SEEDING
        #warning OFF_BEFORE_SEEDING
        nW=0;
        return 0;
    #endif

    revertStrand = false; //the 2nd read is awlays on opposite strand. 1st and 2nd reads have already been reversed.

    if (Lread>0) {
        Nsplit=qualitySplit(Read1[0], Lread, P.maxNsplit, P.seedSplitMin, splitR);
        // splitR[0][fragnum] => good region start position   (from SequenceFuns.cpp)
        // splitR[1][fragnum] => good reagion length
        // splitR[2][fragnum] => fragnum ?
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
    trInit->readNmates=readNmates; //not readNends: this is alignment
    trInit->readName=readName;

    trBest=trInit;

    uint seedSearchStartLmax=min(P.seedSearchStartLmax, // 50
                                  (uint) (P.seedSearchStartLmaxOverLread*(Lread-1))); // read length
    prefetchChainFirstSeeds(mapGen, P, Read1[0], splitR, Nsplit, seedSearchStartLmax);
    // align all good pieces
    for (uint ip=0; ip<Nsplit; ip++) {


        // if the good piece is long, then try multiple times to map it
        uint Nstart = P.seedSearchStartLmax>0 && seedSearchStartLmax<splitR[1][ip] ? splitR[1][ip]/seedSearchStartLmax+1 : 1;
        uint Lstart = splitR[1][ip]/Nstart;  // length of segment to map.
        bool flagDirMap=true;

        for (uint iDir=0; iDir<2; iDir++) {//loop over two directions (leftToRight and RightToLeft)

            uint Lmapped, L;

            for (uint istart=0; istart<Nstart; istart++) {  // each try to map a segment.

//               #ifdef COMPILE_FOR_LONG_READS
//
//               #else
                if (flagDirMap || istart>0) {//check if the 1st piece in reveree direction does not need to be remapped
                    Lmapped=0;  // length of segment mapped so far.

                    // begin mapping starting from segment start position (istart*Lstart)
                    while ( istart*Lstart + Lmapped + P.seedMapMin < splitR[1][ip] ) {//map until unmapped portion is <=minLmap (default: 5)

                        // Shift is the position in the read to begin mapping from.
                        uint Shift = iDir==0 ? ( splitR[0][ip] + istart*Lstart + Lmapped ) : \
                                   ( splitR[0][ip] + splitR[1][ip] - istart*Lstart-1-Lmapped); //choose Shift for forward or reverse

                        //uint seedLength=min(splitR[1][ip] - Lmapped - istart*Lstart, P.seedSearchLmax);
                        uint seedLength=splitR[1][ip] - Lmapped - istart*Lstart; // what's left of the read to align.
                        maxMappableLength2strands(Shift, seedLength, iDir, 0, mapGen.nSA-1, L, splitR[2][ip]);//L=max mappable length, unique or multiple
                        if (iDir==0 && istart==0 && Lmapped==0 && Shift+L == splitR[1][ip] ) {//this piece maps full length and does not need to be mapped from the opposite direction
                            flagDirMap=false;
                        };
                        Lmapped+=L;
                    };//while ( istart*Lstart + Lmapped + P.minLmap < splitR[1][ip] )
                };//if (flagDirMap || istart>0)

                if (P.seedSearchLmax>0) {//search fixed length. Not very efficient, need to improve
                    // off by default.
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
