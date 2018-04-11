#include "genomeSAindex.h"
#include "TimeFunctions.h"
#include "SuffixArrayFuns.h"
#include "ErrorWarning.h"

void genomeSAindex(char * G, PackedArray & SA, Parameters & P, PackedArray & SAi, Genome &mapGen)
 {
    mapGen.genomeSAindexStart = new uint [mapGen.pGe.gSAindexNbases+1];
    mapGen.genomeSAindexStart[0]=0;
    for (uint ii=1;ii<=mapGen.pGe.gSAindexNbases;ii++) {//L-mer indices starts
        mapGen.genomeSAindexStart[ii] = mapGen.genomeSAindexStart[ii-1] + ( 1LLU<<(2*ii) );
    };
    mapGen.nSAi = mapGen.genomeSAindexStart[mapGen.pGe.gSAindexNbases];

    /* testing
    //     uint* SAi1=new uint[mapGen.nSAi];

    PackedArray SAio;
    SAio.defineBits(mapGen.GstrandBit+3,mapGen.nSAi);
    SAio.allocateArray();
    ifstream oldSAiin("./DirTrue/SAindex");
    oldSAiin.read(SAio.charArray,8*(mapGen.pGe.gSAindexNbases+2));//skip first bytes
    oldSAiin.read(SAio.charArray,SAio.lengthByte);
    oldSAiin.close();
    */


    mapGen.SAiMarkNbit=mapGen.GstrandBit+1;
    mapGen.SAiMarkAbsentBit=mapGen.GstrandBit+2;

    mapGen.SAiMarkNmaskC=1LLU << mapGen.SAiMarkNbit;
    mapGen.SAiMarkNmask=~mapGen.SAiMarkNmaskC;
    mapGen.SAiMarkAbsentMaskC=1LLU << mapGen.SAiMarkAbsentBit;
    mapGen.SAiMarkAbsentMask=~mapGen.SAiMarkAbsentMaskC;


    SAi.defineBits(mapGen.GstrandBit+3,mapGen.nSAi);
    SAi.allocateArray();

    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain    << timeMonthDayTime(rawTime) <<" ... generating Suffix Array index\n" <<flush;
    *P.inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... generating Suffix Array index\n" <<flush;

    /*testing
    PackedArray SA1=SA;
    uint* ind0=new uint[mapGen.pGe.gSAindexNbases];

    for (uint ii=0; ii<mapGen.pGe.gSAindexNbases; ii++) {
        ind0[ii]=-1;//this is needed in case "AAA...AAA",i.e. indPref=0 is not present in the genome for some lengths
    };
    uint* SAi1=new uint[mapGen.nSAi];

    for (uint isa=0; isa<mapGen.nSA; isa++) {//for all suffixes
        if (isa%100000000==0) P.inOut->logMain  << isa*100/mapGen.nSA << "% " << flush;

        uint SAstr=SA1[isa];
        bool dirG = (SAstr>>mapGen.GstrandBit) == 0; //forward or reverse strand of the genome
        SAstr &= mapGen.GstrandMask;
        if (!dirG) SAstr=mapGen.nGenome-1-SAstr;

        uint indPref=0;
        for (uint iL=0; iL < mapGen.pGe.gSAindexNbases; iL++) {//calculate index

            indPref <<= 2;

            uint g1= (uint) G[dirG ? SAstr+iL : SAstr-iL]; //reverese if (-) strand

            if (g1>3) {//if N, this suffix does not belong in SAi
                for (uint iL1=iL; iL1 < mapGen.pGe.gSAindexNbases; iL1++) {
                    SAi1[mapGen.genomeSAindexStart[iL1]+ind0[iL1]] |= mapGen.SAiMarkNmaskC;
                };
                break;
            };

            if (!dirG) g1=3-g1; //complement if (-) strand

            indPref += (uint) g1;

            if ( indPref > ind0[iL] || isa==0 ) {//new && good index, record it
                SAi1[mapGen.genomeSAindexStart[iL]+indPref]=isa;
                for (uint ii=ind0[iL]+1; ii<indPref; ii++) {//index is not present, record to the last present suffix
                    SAi1[mapGen.genomeSAindexStart[iL]+ii] = isa | mapGen.SAiMarkAbsentMaskC;
                };
                ind0[iL]=indPref;
            } else if ( indPref < ind0[iL] ) {
                ostringstream errOut;
                errOut << "BUG: next index is smaller than previous, EXITING\n" <<flush;
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            };
        };
    };//for (uint isa=0; isa<mapGen.nSA; isa++)
    */

    genomeSAindexChunk(G, SA, P, SAi, 0, SA.length-1, mapGen);

    time(&rawTime);
    P.inOut->logMain    << timeMonthDayTime(rawTime) <<" ... completed Suffix Array index\n" <<flush;
    *P.inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... completed Suffix Array index\n" <<flush;

//     for (uint ii=1;ii<=mapGen.pGe.gSAindexNbases-1;ii++) {//L-mer indices starts
//         cout <<ii<<endl;
//         for (uint jj=mapGen.genomeSAindexStart[ii-1]; jj<mapGen.genomeSAindexStart[ii]; jj++)
//         {
//             if (SAi[jj]!=SAio[jj])
//             {
//                 cout <<ii <<" "<< jj<<" "<<jj-mapGen.genomeSAindexStart[ii-1]<<" "<<SAi[jj]<<" "<<SAio[jj]<<" "<<endl;
//                 sleep(100);
//             };
//         };
//     };


 };


void genomeSAindexChunk(char * G, PackedArray & SA, Parameters & P, PackedArray & SAi, uint iSA1, uint iSA2, Genome &mapGen)
{
    uint* ind0=new uint[mapGen.pGe.gSAindexNbases];
    uint* ind0a=new uint[mapGen.pGe.gSAindexNbases];

    for (uint ii=0; ii<mapGen.pGe.gSAindexNbases; ii++) {
        ind0[ii]=-1;//this is needed in case "AAA...AAA",i.e. indPref=0 is not present in the genome for some lengths
        ind0a[ii]=-1;//this is needed in case "AAA...AAA",i.e. indPref=0 is not present in the genome for some lengths
    };

    PackedArray SAi1;
    SAi1=SAi;
    SAi1.allocateArray();

    uint isaStep=mapGen.nSA/(1llu<<(2*mapGen.pGe.gSAindexNbases))+1;
//     isaStep=8;

    uint isa=iSA1;
    int iL4;
    uint indFull=funCalcSAiFromSA(G,SA,mapGen,isa,mapGen.pGe.gSAindexNbases,iL4);
    while (isa<=iSA2) {//for all suffixes

        /* testing
        uint SAstr=SA[isa];
        bool dirG = (SAstr>>mapGen.GstrandBit) == 0; //forward or reverse strand of the genome
        SAstr &= mapGen.GstrandMask;
        if (!dirG) SAstr=mapGen.nGenome-1-SAstr;
        uint indPref1=0;
        */

        for (uint iL=0; iL < mapGen.pGe.gSAindexNbases; iL++) {//calculate index
            /*{//testing: old way
                indPref1 <<= 2;

                uint g1= (uint) G[dirG ? SAstr+iL : SAstr-iL]; //reverese if (-) strand

                if (g1>3) {//if N, this suffix does not belong in SAi
                    for (uint iL1=iL; iL1 < mapGen.pGe.gSAindexNbases; iL1++) {
                        SAi1.writePacked(mapGen.genomeSAindexStart[iL1]+ind0[iL1],SAi[mapGen.genomeSAindexStart[iL1]+ind0[iL1]] | mapGen.SAiMarkNmaskC);
                    };
                } else //relying on the true code to break iL cycle
                {
                    if (!dirG) g1=3-g1; //complement if (-) strand

                    indPref1 += (uint) g1;

                    if ( indPref1 > ind0a[iL] || isa==0 ) {//new && good index, record it
                        SAi1.writePacked(mapGen.genomeSAindexStart[iL]+indPref1, isa);
                        for (uint ii=ind0a[iL]+1; ii<indPref1; ii++) {//index is not present, record to the last present suffix
                            SAi1.writePacked(mapGen.genomeSAindexStart[iL]+ii, isa | mapGen.SAiMarkAbsentMaskC);
                        };
                        ind0a[iL]=indPref1;
                    } else if ( indPref1 < ind0a[iL] ) {
                        ostringstream errOut;
                        errOut << "BUG: next index is smaller than previous, EXITING\n" <<flush;
                        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                    };
                };
            };
            */

            uint indPref = indFull >> (2*(mapGen.pGe.gSAindexNbases-1-iL));
//             if (indPref!=indPref1)
//                 cout<< iL <<" "<< isa <<" "<< indPref <<" "<<indPref1<<endl;


            if ( (int)iL==iL4 ) {//this suffix contains N and does not belong in SAi
                for (uint iL1=iL; iL1 < mapGen.pGe.gSAindexNbases; iL1++) {
                    SAi.writePacked(mapGen.genomeSAindexStart[iL1]+ind0[iL1],SAi[mapGen.genomeSAindexStart[iL1]+ind0[iL1]] | mapGen.SAiMarkNmaskC);
//                     if (SAi[mapGen.genomeSAindexStart[iL]+ind0[iL1]] != SAi1[mapGen.genomeSAindexStart[iL]+ind0[iL1]])
//                         cout<< iL <<" "<< isa <<" "<< indPref <<" "<<indPref1<<endl;

                };
                break;//break the iL cycle
            };

            if ( indPref > ind0[iL] || isa==0 ) {//new && good index, record it
                //testing
//                 if (funCalcSAiFromSA(G,SA,isa,iL+1,P)!=indPref)
//                     cout<< iL <<" "<< isa <<" "<< indPref <<" "<<funCalcSAiFromSA(G,SA,isa,iL+1,P)<<endl;

                SAi.writePacked(mapGen.genomeSAindexStart[iL]+indPref, isa);
//                 if (SAi[mapGen.genomeSAindexStart[iL]+indPref] != SAi1[mapGen.genomeSAindexStart[iL]+indPref])
//                     cout<< iL <<" "<< isa <<" "<< indPref <<" "<<indPref1<<endl;

                for (uint ii=ind0[iL]+1; ii<indPref; ii++) {//index is not present, record to the last present suffix
                    SAi.writePacked(mapGen.genomeSAindexStart[iL]+ii, isa | mapGen.SAiMarkAbsentMaskC);
                };
                ind0[iL]=indPref;

            } else if ( indPref < ind0[iL] ) {
                ostringstream errOut;
                errOut << "BUG: next index is smaller than previous, EXITING\n" <<flush;
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            };


        };

        //find next index not equal to the current one
        funSAiFindNextIndex(G, SA, isaStep, isa, indFull, iL4, mapGen);//indFull and iL4 have been already defined at the previous step
//         isa++;
//         indFull=funCalcSAiFromSA(G,SA,isa,mapGen.pGe.gSAindexNbases,P,iL4);
    };//isa cycle
    delete [] ind0;

 };

void funSAiFindNextIndex(char * G, PackedArray & SA, uint isaStep, uint & isa, uint & indFull, int & iL4, Genome &mapGen)
 {
    uint indFullPrev=indFull;
    int iL4prev=iL4;
    isa+=isaStep;
    while (isa<mapGen.nSA && (indFull=funCalcSAiFromSA(G,SA,mapGen,isa,mapGen.pGe.gSAindexNbases,iL4))==indFullPrev && iL4==iL4prev)
    {//make large step in isa while the indFull/iL4 are still the same
        isa+=isaStep;
    };
    if (isa>=mapGen.nSA)
    {//reached the end of the SA
        indFull=funCalcSAiFromSA(G,SA,mapGen,mapGen.nSA-1,mapGen.pGe.gSAindexNbases,iL4);
        if (indFull==indFullPrev && iL4==iL4prev)
        {
            isa=mapGen.nSA;//no more indices, the last one is equal to the previous
            return;
        };
    };

    {//binary search
        uint i1=isa-isaStep;
        uint i2=min(isa,mapGen.nSA-1);
        while (i1+1<i2)
        {
            isa=i1/2 + i2/2 + (i1%2 + i2%2)/2;
            if ((indFull=funCalcSAiFromSA(G,SA,mapGen,isa,mapGen.pGe.gSAindexNbases,iL4))==indFullPrev && iL4==iL4prev)
            {
                i1=isa;
            } else
            {
                i2=isa;
            };
        };
        if (isa==i1)
        {
            isa=i2;
            indFull=funCalcSAiFromSA(G,SA,mapGen,isa,mapGen.pGe.gSAindexNbases,iL4);
        };
    };
};
