#include "sjdbBuildIndex.h"
// #include "sjdbLoadFromStream.h"
// #include "sjdbPrepare.h"
#include "ErrorWarning.h"
#include "SuffixArrayFuns.h"
#include "SequenceFuns.h"
#include "serviceFuns.cpp"
#include "IncludeDefine.h"
#include "streamFuns.h"
#include "binarySearch2.h"
#include "ErrorWarning.h"
#include <cmath>

#include "funCompareUintAndSuffixes.h"

void sjdbBuildIndex (Parameters &P, char *Gsj, char *G, PackedArray &SA, PackedArray &SA2, PackedArray &SAi, Genome &mapGen, Genome &mapGen1) {

    #define SPACER_CHAR GENOME_spacingChar

    if (mapGen.sjdbN==0)
    {//no junctions to insert
        return;
    };

    time_t rawtime;
    time ( &rawtime );
    P.inOut->logMain   << timeMonthDayTime(rawtime) << " ..... inserting junctions into the genome indices" <<endl;
    *P.inOut->logStdOut  << timeMonthDayTime(rawtime) << " ..... inserting junctions into the genome indices" <<endl;

    uint nGsj=mapGen.sjdbLength*mapGen.sjdbN;
    for (uint ii=1; ii<=mapGen.sjdbN; ii++)
    {
        Gsj[ii*mapGen.sjdbLength-1]=SPACER_CHAR; //to make sure this is > than any genome char
    };
    Gsj[nGsj*2]=SPACER_CHAR;//mark the end of the text

    for (uint ii=0; ii<nGsj; ii++) {//reverse complement junction sequences
        Gsj[nGsj*2-1-ii]=Gsj[ii]<4 ? 3-Gsj[ii] : Gsj[ii]; //reverse complement
    };

    char* G1c=new char[nGsj*2+1];
    complementSeqNumbers(Gsj, G1c, nGsj*2+1);

    uint32* oldSJind=new uint32[mapGen1.sjdbN];

//     uint nIndicesSJ1=mapGen.sjdbOverhang;
    uint   nIndicesSJ1=mapGen.sjdbLength;//keep all indices - this is pre-2.4.1 of generating the genome

    uint64* indArray=new uint64[2*mapGen.sjdbN*(nIndicesSJ1+1)*2];//8+4 bytes for SA index and index in the genome * nJunction * nIndices per junction * 2 for reverse compl
    uint64 sjNew=0;
    #pragma omp parallel num_threads(P.runThreadN)
    #pragma omp for schedule (dynamic,1000) reduction(+:sjNew)
    for (uint isj=0; isj<2*mapGen.sjdbN; isj++) {//find insertion points for each of the sequences

        char** seq1=new char*[2];
        seq1[0]=Gsj+isj*mapGen.sjdbLength;
        seq1[1]=G1c+isj*mapGen.sjdbLength;

        uint isj1=isj<mapGen.sjdbN ? isj : 2*mapGen.sjdbN-1-isj;
        int sjdbInd = mapGen1.sjdbN==0 ? -1 : binarySearch2(mapGen.sjdbStart[isj1],mapGen.sjdbEnd[isj1],mapGen1.sjdbStart,mapGen1.sjdbEnd,mapGen1.sjdbN);
        if (sjdbInd<0)
        {//count new junctions
            ++sjNew;
        } else
        {//record new index of the old junctions
            oldSJind[sjdbInd]=isj1;
        };

        for (uint istart1=0; istart1<nIndicesSJ1;istart1++) {

            uint istart=istart1;
//             uint istart=isj<mapGen.sjdbN ? istart1 : istart1+1; //for rev-compl junction, shift by one base to start with the 1st non-spacer base
            uint ind1=2*(isj*nIndicesSJ1+istart1);
            if (sjdbInd>=0 || seq1[0][istart]>3)
            {//no index for already included junctions, or suffices starting with N
                indArray[ind1]=-1;
            } else
            {
                //indArray[ind1] =  suffixArraySearch(seq1, istart, mapGen.sjdbLength-istart1, G, SA, true, 0, mapGen.nSA-1, 0, P) ;
                indArray[ind1] =  suffixArraySearch1(mapGen, seq1, istart, 10000, -1LLU, true, 0, mapGen.nSA-1, 0) ;
                //-1LLU results in suffixes for the new junctions to be always included in SA *after* the suffixes of the old junctions
                //for identical suffixes, this may result in unstable ordering
                indArray[ind1+1] = isj*mapGen.sjdbLength+istart;
            };
        };
    };
    sjNew = sjNew/2;//novel junctions were double counted on two strands

    time ( &rawtime );
    P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished SA search: number of new junctions=" << sjNew <<", old junctions="<<mapGen.sjdbN-sjNew<<endl;

    uint nInd=0;//true number of new indices
    for (uint ii=0; ii<2*mapGen.sjdbN*nIndicesSJ1; ii++) {//remove entries that cannot be inserted, this cannot be done in the parallel cycle above
        if (indArray[ii*2]!= (uint) -1) {
            indArray[nInd*2]=indArray[ii*2];
            indArray[nInd*2+1]=indArray[ii*2+1];
            ++nInd;
        };
    };

    g_funCompareUintAndSuffixes_G=Gsj;
    qsort((void*) indArray, nInd, 2*sizeof(uint64), funCompareUintAndSuffixes);
    time ( &rawtime );
    P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished sorting SA indicesL nInd="<<nInd <<endl;

    indArray[2*nInd]=-999; //mark the last junction
    indArray[2*nInd+1]=-999; //mark the last junction

    mapGen.nGenome=mapGen.chrStart[mapGen.nChrReal]+nGsj;
    mapGen.nSA+=nInd;

    uint GstrandBit1 = (uint) floor(log(mapGen.nGenome)/log(2))+1;
    if (GstrandBit1<32) GstrandBit1=32; //TODO: use simple access function for SA
    if ( GstrandBit1 != mapGen.GstrandBit)
    {//too many junctions were added - GstrandBit changed
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: cannot insert junctions on the fly because of strand GstrandBit problem\n";
        errOut << "SOLUTION: please contact STAR author at https://groups.google.com/forum/#!forum/rna-star\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };

    SA2.defineBits(mapGen.GstrandBit+1,mapGen.nSA);
    uint nGsjNew=sjNew*mapGen.sjdbLength; //this is the actual number of bytes added to the genome, while nGsj is the total size of all junctions

    uint N2bit= 1LLU << mapGen.GstrandBit;
    uint strandMask=~N2bit;

    /*testing
    PackedArray SAo;
    SAo.defineBits(mapGen.GstrandBit+1,mapGen.nSA);
    SAo.allocateArray();
    ifstream oldSAin("./DirTrue/SA");
    oldSAin.read(SAo.charArray,SAo.lengthByte);
    oldSAin.close();
    */

    uint isj=0, isa2=0;
    for (uint isa=0;isa<mapGen1.nSA;isa++) {
        while (isa==indArray[isj*2]) {//insert sj index before the existing index
            uint ind1=indArray[isj*2+1];
            if (ind1<nGsj) {
                ind1+=mapGen.chrStart[mapGen.nChrReal];
            } else {//reverse strand
                ind1=(ind1-nGsj) | N2bit;
            };
            SA2.writePacked(isa2,ind1);
            /*testing
            if (SA2[isa2]!=SAo[isa2]) {
               cout <<isa2 <<" "<< SA2[isa2]<<" "<<SAo[isa2]<<endl;
               //sleep(100);
            };
            */
            ++isa2; ++isj;
        };

        uint ind1=SA[isa];
        if ( (ind1 & N2bit)>0 )
        {//- strand
            uint ind1s = mapGen1.nGenome - (ind1 & strandMask);
            if (ind1s>mapGen.chrStart[mapGen.nChrReal])
            {//this index was an old sj, may need to shift it
                uint sj1 = (ind1s-mapGen.chrStart[mapGen.nChrReal])/mapGen.sjdbLength;//old junction index
                ind1s += (oldSJind[sj1]-sj1)*mapGen.sjdbLength;
                ind1 = (mapGen.nGenome - ind1s) | N2bit;
            } else
            {
                ind1+=nGsjNew; //reverse complementary indices are all shifted by the length of junctions
            };
        } else
        {//+ strand
            if (ind1>mapGen.chrStart[mapGen.nChrReal])
            {//this index was an old sj, may need to shift it
                uint sj1 = (ind1-mapGen.chrStart[mapGen.nChrReal])/mapGen.sjdbLength;//old junction index
                ind1 += (oldSJind[sj1]-sj1)*mapGen.sjdbLength;
            };
        };

        SA2.writePacked(isa2,ind1);
            /*testing
            if (SA2[isa2]!=SAo[isa2]) {
               cout <<isa2 <<" "<< SA2[isa2]<<" "<<SAo[isa2]<<endl;
               //sleep(100);
            };
            */
        ++isa2;
    };

    for (;isj<nInd;isj++) {//insert last new indices after the last old index
        uint ind1=indArray[isj*2+1];
        if (ind1<nGsj) {
            ind1+=mapGen.chrStart[mapGen.nChrReal];
        } else {//reverse strand
            ind1=(ind1-nGsj) | N2bit;
        };
        SA2.writePacked(isa2,ind1);
            /*testing
            if (SA2[isa2]!=SAo[isa2]) {
               cout <<isa2 <<" "<< SA2[isa2]<<" "<<SAo[isa2]<<endl;
               //sleep(100);
            };
            */
        ++isa2;
    };

    time ( &rawtime );
    P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished inserting junction indices" <<endl;

    //SAi insertions
    for (uint iL=0; iL < mapGen.pGe.gSAindexNbases; iL++) {
        uint iSJ=0;
        uint ind0=mapGen.genomeSAindexStart[iL]-1;//last index that was present in the old genome
        for (uint ii=mapGen.genomeSAindexStart[iL];ii<mapGen.genomeSAindexStart[iL+1]; ii++) {//scan through the longest index
            uint iSA1=SAi[ii];
            uint iSA2=iSA1 & mapGen.SAiMarkNmask & mapGen.SAiMarkAbsentMask;

            if ( iSJ<nInd && (iSA1 &  mapGen.SAiMarkAbsentMaskC)>0 )
            {//index missing from the old genome
                uint iSJ1=iSJ;
                int64 ind1=funCalcSAi(Gsj+indArray[2*iSJ+1],iL);
                while (ind1 < (int64)(ii-mapGen.genomeSAindexStart[iL]) && indArray[2*iSJ]-1<iSA2) {
                    ++iSJ;
                    ind1=funCalcSAi(Gsj+indArray[2*iSJ+1],iL);
                };
                if (ind1 == (int64)(ii-mapGen.genomeSAindexStart[iL]) ) {
                    SAi.writePacked(ii,indArray[2*iSJ]-1+iSJ+1);
                    for (uint ii0=ind0+1; ii0<ii; ii0++) {//fill all the absent indices with this value
                        SAi.writePacked(ii0,(indArray[2*iSJ]-1+iSJ+1) | mapGen.SAiMarkAbsentMaskC);
                    };
                    ++iSJ;
                    ind0=ii;
                } else {
                    iSJ=iSJ1;
                };
            } else
            {//index was present in the old genome
                while (iSJ<nInd && indArray[2*iSJ]-1+1<iSA2) {//for this index insert "smaller" junctions
                    ++iSJ;
                };

                while (iSJ<nInd && indArray[2*iSJ]-1+1==iSA2) {//special case, the index falls right behind SAi
                    if (funCalcSAi(Gsj+indArray[2*iSJ+1],iL) >= (int64) (ii-mapGen.genomeSAindexStart[iL]) ) {//this belongs to the next index
                        break;
                    };
                    ++iSJ;
                };

                SAi.writePacked(ii,iSA1+iSJ);

                for (uint ii0=ind0+1; ii0<ii; ii0++) {//fill all the absent indices with this value
                    SAi.writePacked(ii0,(iSA2+iSJ) | mapGen.SAiMarkAbsentMaskC);
                };
                ind0=ii;
            };
        };

    };

    for (uint isj=0;isj<nInd;isj++) {
        int64 ind1=0;
        for (uint iL=0; iL < mapGen.pGe.gSAindexNbases; iL++) {
            uint g=(uint) Gsj[indArray[2*isj+1]+iL];
            ind1 <<= 2;
            if (g>3) {//this iSA contains N, need to mark the previous
                for (uint iL1=iL; iL1 < mapGen.pGe.gSAindexNbases; iL1++) {
                    ind1+=3;
                    int64 ind2=mapGen.genomeSAindexStart[iL1]+ind1;
                    for (; ind2>=0; ind2--) {//find previous index that is not absent
                        if ( (SAi[ind2] & mapGen.SAiMarkAbsentMaskC)==0 ) {
                            break;
                        };
                    };
                    SAi.writePacked(ind2,SAi[ind2] | mapGen.SAiMarkNmaskC);
                    ind1 <<= 2;
                };
                break;
            } else {
                ind1 += g;
            };
        };
    };
    time ( &rawtime );
    P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished SAi" <<endl;

    //change parameters, most parameters are already re-defined in sjdbPrepare.cpp
    SA.defineBits(mapGen.GstrandBit+1,mapGen.nSA);//same as SA2
    SA.pointArray(SA2.charArray);
    mapGen.nSAbyte=SA.lengthByte;
    mapGen.sjGstart=mapGen.chrStart[mapGen.nChrReal];
    memcpy(G+mapGen.chrStart[mapGen.nChrReal],Gsj, nGsj);

    /* testing
    PackedArray SAio=SAi;
    SAio.allocateArray();
    ifstream oldSAiin("./DirTrue/SAindex");
    oldSAiin.read(SAio.charArray,8*(P.pGe.gSAindexNbases+2));//skip first bytes
    oldSAiin.read(SAio.charArray,SAio.lengthByte);
    oldSAiin.close();


    for (uint ii=0;ii<mapGen.nSA;ii++) {
        if (SA2[ii]!=SAo[ii]) {
            cout <<ii <<" "<< SA2[ii]<<" "<<SAo[ii]<<endl;
        };
    };


    for (uint iL=0; iL < P.pGe.gSAindexNbases; iL++) {
        for (uint ii=mapGen.genomeSAindexStart[iL];ii<mapGen.genomeSAindexStart[iL+1]; ii++) {//scan through the longets index
                if ( SAio[ii]!=SAi[ii] ) {
                    cout <<ii<<" "<<SAio[ii]<<" "<<SAi[ii]<<endl;
                };
        };
    };
    */

    /*
    ofstream genomeOut("/home/dobin/Genome");
    fstreamWriteBig(genomeOut,G,mapGen.nGenome+nGsj,"777","777",P);
    genomeOut.close();
    genomeOut.open("/home/dobin/SA");
    fstreamWriteBig(genomeOut,SA2.charArray,SA2.lengthByte,"777","777",P);
    genomeOut.close();
    */

    delete [] indArray;
    delete [] G1c;
    delete [] oldSJind;

};
