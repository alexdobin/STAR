#include "genomeSAindex.h" 
#include "TimeFunctions.h"
#include "SuffixArrayFuns.h"
#include "ErrorWarning.h"

void genomeSAindex(char * G, PackedArray & SA, Parameters * P, PackedArray & SAi)
 {
    P->genomeSAindexStart = new uint [P->genomeSAindexNbases+1];
    P->genomeSAindexStart[0]=0;
    for (uint ii=1;ii<=P->genomeSAindexNbases;ii++) {//L-mer indices starts
        P->genomeSAindexStart[ii] = P->genomeSAindexStart[ii-1] + ( 1LLU<<(2*ii) );
    };
    P->nSAi = P->genomeSAindexStart[P->genomeSAindexNbases];
    
    ///* testing
//     uint* SAi1=new uint[P->nSAi];

    PackedArray SAio;
    SAio.defineBits(P->GstrandBit+3,P->nSAi);
    SAio.allocateArray();
    ifstream oldSAiin("./DirTrue/SAindex");
    oldSAiin.read(SAio.charArray,8*(P->genomeSAindexNbases+2));//skip first bytes
    oldSAiin.read(SAio.charArray,SAio.lengthByte);
    oldSAiin.close();      
    //*/
    

    P->SAiMarkNbit=P->GstrandBit+1;
    P->SAiMarkAbsentBit=P->GstrandBit+2;
    
    P->SAiMarkNmaskC=1LLU << P->SAiMarkNbit;
    P->SAiMarkNmask=~P->SAiMarkNmaskC;
    P->SAiMarkAbsentMaskC=1LLU << P->SAiMarkAbsentBit;
    P->SAiMarkAbsentMask=~P->SAiMarkAbsentMaskC;

//     uint* SAi=new uint[P->nSAi];
    
    SAi.defineBits(P->GstrandBit+3,P->nSAi);
    SAi.allocateArray();   
    
    time_t rawTime;
    time(&rawTime);    
    P->inOut->logMain    << timeMonthDayTime(rawTime) <<" ... Generating Suffix Array index\n" <<flush;   
    *P->inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... Generating Suffix Array index\n" <<flush; 
    
    genomeSAindexChunk(G, SA, P, SAi, 0, P->nSA-1);

    
    time(&rawTime);    
    P->inOut->logMain    << timeMonthDayTime(rawTime) <<" ... Packing Suffix Array index\n" <<flush;   
    *P->inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... Packing Suffix Array index\n" <<flush; 
      
    for (uint ii=1;ii<=P->genomeSAindexNbases-1;ii++) {//L-mer indices starts
        cout <<ii<<endl;
        for (uint jj=P->genomeSAindexStart[ii-1]; jj<P->genomeSAindexStart[ii]; jj++)
        {
            if (SAi[jj]!=SAio[jj])
                cout <<ii <<" "<< jj<<SAi[jj]<<" "<<SAio[jj]<<" "<<endl;
        };
    };
//     //pack SAi
//     SAip.defineBits(P->GstrandBit+3,P->nSAi);//SAi uses an extra bit compared to SA because it needs to store values > nSA
//     SAip.pointArray((char*) SAi);
//     for (uint ii=0;ii<SAip.length;ii++) {
//         SAip.writePacked(ii,SAi[ii]);
//     };
    
 };
 
 
 void genomeSAindexChunk(char * G, PackedArray & SA, Parameters * P, PackedArray & SAi, uint iSA1, uint iSA2)
 {
         uint* ind0=new uint[P->genomeSAindexNbases];
//     uint* ind0a=new uint[P->genomeSAindexNbases];

    for (uint ii=0; ii<P->genomeSAindexNbases; ii++) {
        ind0[ii]=-1;//this is needed in case "AAA...AAA",i.e. indPref=0 is not present in the genome for some lengths
//         ind0a[ii]=-1;//this is needed in case "AAA...AAA",i.e. indPref=0 is not present in the genome for some lengths
    };    
    uint isaStep=P->nSA/(1llu<<(2*P->genomeSAindexNbases));
//     uint isareport=0;
    for (uint isa=iSA1; isa<=iSA2; isa++) {//for all suffixes

//         if (isareport+100000000<isa) 
//         {
//             P->inOut->logMain  << isa*100/P->nSA << "% " << flush;         
//             cout << isa*100/P->nSA << "% " << flush; 
//         };
        uint indFull, indFullPrev;
        int iL4, iL4prev;

        if (isaStep>4 && isa>0)
        {//skip indices with the same prefix

            isa+=isaStep-1;//-1 since +1 was added in the cycle
            while (isa<P->nSA && funCalcSAiFromSA(G,SA,isa,P->genomeSAindexNbases,P,iL4)==indFullPrev && iL4==iL4prev)
            {
                isa+=isaStep;
            };
            if (isa>=P->nSA)
            {
                isa=P->nSA-1;
                indFull=funCalcSAiFromSA(G,SA,isa,P->genomeSAindexNbases,P,iL4);
            } else
            {//binary search
                uint i1=isa-isaStep;
                uint i2=isa;
                while (i1+1<i2)
                {
                    isa=i1/2 + i2/2 + (i1%2 + i2%2)/2;
                    if ((indFull=funCalcSAiFromSA(G,SA,isa,P->genomeSAindexNbases,P,iL4))==indFullPrev && iL4==iL4prev)
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
                    indFull=funCalcSAiFromSA(G,SA,isa,P->genomeSAindexNbases,P,iL4);
                };
            };
        } else
        {
            indFull=funCalcSAiFromSA(G,SA,isa,P->genomeSAindexNbases,P,iL4);
        };
        
        indFullPrev=indFull;
        iL4prev=iL4;
        
        /* testing
        uint SAstr=SA[isa];
        bool dirG = (SAstr>>P->GstrandBit) == 0; //forward or reverse strand of the genome
        SAstr &= P->GstrandMask;
        if (!dirG) SAstr=P->nGenome-1-SAstr;
        uint indPref1=0;
         */
      
        for (uint iL=0; iL < P->genomeSAindexNbases; iL++) {//calculate index
//             {//testing: old way
//                 indPref1 <<= 2;
// 
//                 uint g1= (uint) G[dirG ? SAstr+iL : SAstr-iL]; //reverese if (-) strand
// 
//                 if (g1>3) {//if N, this suffix does not belong in SAi
//                     for (uint iL1=iL; iL1 < P->genomeSAindexNbases; iL1++) {
//                         SAi1[P->genomeSAindexStart[iL1]+ind0a[iL1]] |= P->SAiMarkNmaskC;
//                     };
//                 } else //relying on the true code to break iL cycle
//                 {
//                     if (!dirG) g1=3-g1; //complement if (-) strand
// 
//                     indPref1 += (uint) g1;
// 
//                     if ( indPref1 > ind0a[iL] || isa==0 ) {//new && good index, record it
//                         SAi1[P->genomeSAindexStart[iL]+indPref1]=isa;
//                         for (uint ii=ind0a[iL]+1; ii<indPref1; ii++) {//index is not present, record to the last present suffix
//                             SAi1[P->genomeSAindexStart[iL]+ii] = isa | P->SAiMarkAbsentMaskC;
//                         };
//                         ind0a[iL]=indPref1;
//                     } else if ( indPref1 < ind0a[iL] ) {
//                         ostringstream errOut;
//                         errOut << "BUG: next index is smaller than previous, EXITING\n" <<flush;
//                         exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
//                     };
//                 };
//             };

            uint indPref = indFull >> (2*(P->genomeSAindexNbases-1-iL));
//             if (indPref!=indPref1)
//                 cout<< iL <<" "<< isa <<" "<< indPref <<" "<<indPref1<<endl;
            

            if ( iL==iL4 ) {//this suffix contains N and does not belong in SAi
                for (uint iL1=iL; iL1 < P->genomeSAindexNbases; iL1++) {
                    SAi.writePacked(P->genomeSAindexStart[iL1]+ind0[iL1],SAi[P->genomeSAindexStart[iL1]+ind0[iL1]] | P->SAiMarkNmaskC);
//                     if (SAi[P->genomeSAindexStart[iL]+ind0[iL1]] != SAio[P->genomeSAindexStart[iL]+ind0[iL1]])
//                         cout<< iL <<" "<< isa <<" "<< indPref <<" "<<indPref1<<endl;
                    
                };
                break;//break the iL cycle
            };

            if ( indPref > ind0[iL] || isa==0 ) {//new && good index, record it
                //testing
//                 if (funCalcSAiFromSA(G,SA,isa,iL+1,P)!=indPref)
//                     cout<< iL <<" "<< isa <<" "<< indPref <<" "<<funCalcSAiFromSA(G,SA,isa,iL+1,P)<<endl;
                        
                SAi.writePacked(P->genomeSAindexStart[iL]+indPref, isa);
//                 if (SAi[P->genomeSAindexStart[iL]+indPref] != SAio[P->genomeSAindexStart[iL]+indPref])
//                     cout<< iL <<" "<< isa <<" "<< indPref <<" "<<indPref1<<endl;
                
                for (uint ii=ind0[iL]+1; ii<indPref; ii++) {//index is not present, record to the last present suffix
                    SAi.writePacked(P->genomeSAindexStart[iL]+ii, isa | P->SAiMarkAbsentMaskC); 
                };
                ind0[iL]=indPref;

            } else if ( indPref < ind0[iL] ) {
                ostringstream errOut;
                errOut << "BUG: next index is smaller than previous, EXITING\n" <<flush;
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
            };
            
            
        };
    };//for (uint isa=0; isa<P->nSA; isa++)
    delete [] ind0;

 };