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

void sjdbBuildIndex (Parameters *P, Parameters *P1, char *Gsj, char *G, PackedArray &SA, PackedArray &SA2, PackedArray &SAi) {
    
    #define SPACER_CHAR GENOME_spacingChar

    if (P->sjdbN==0)
    {//no junctions to insert
        return;
    };
    
    time_t rawtime;
    time ( &rawtime );
    P->inOut->logMain   << timeMonthDayTime(rawtime) << " ..... Inserting junctions into the genome indices" <<endl;    
    *P->inOut->logStdOut  << timeMonthDayTime(rawtime) << " ..... Inserting junctions into the genome indices" <<endl;
    
    uint nGsj=P->sjdbLength*P->sjdbN;
    for (uint ii=1; ii<=P->sjdbN; ii++) 
    {
        Gsj[ii*P->sjdbLength-1]=SPACER_CHAR; //to make sure this is > than any genome char
    };
    Gsj[nGsj*2]=SPACER_CHAR;//mark the end of the text

    for (uint ii=0; ii<nGsj; ii++) {//reverse complement junction sequences
        Gsj[nGsj*2-1-ii]=Gsj[ii]<4 ? 3-Gsj[ii] : Gsj[ii]; //reverse complement
    };

    char* G1c=new char[nGsj*2+1];
    complementSeqNumbers(Gsj, G1c, nGsj*2+1);

    uint32* oldSJind=new uint32[P1->sjdbN];
    
//     uint nIndicesSJ1=P->sjdbOverhang;
    uint   nIndicesSJ1=P->sjdbLength;//keep all indices - this is pre-2.4.1 of generating the genome
    
    uint64* indArray=new uint64[2*P->sjdbN*(nIndicesSJ1+1)*2];//8+4 bytes for SA index and index in the genome * nJunction * nIndices per junction * 2 for reverse compl
    uint64 sjNew=0;
    #pragma omp parallel num_threads(P->runThreadN)
    #pragma omp for schedule (dynamic,1000) reduction(+:sjNew)
    for (uint isj=0; isj<2*P->sjdbN; isj++) {//find insertion points for each of the sequences

        char** seq1=new char*[2];
        seq1[0]=Gsj+isj*P->sjdbLength;
        seq1[1]=G1c+isj*P->sjdbLength;
        
        uint isj1=isj<P->sjdbN ? isj : 2*P->sjdbN-1-isj;
        int sjdbInd = P1->sjdbN==0 ? -1 : binarySearch2(P->sjdbStart[isj1],P->sjdbEnd[isj1],P1->sjdbStart,P1->sjdbEnd,P1->sjdbN);
        if (sjdbInd<0) 
        {//count new junctions
            ++sjNew;
        } else 
        {//record new index of the old junctions
            oldSJind[sjdbInd]=isj1;
        };
        
        for (uint istart1=0; istart1<nIndicesSJ1;istart1++) {
            
            uint istart=istart1;
//             uint istart=isj<P->sjdbN ? istart1 : istart1+1; //for rev-compl junction, shift by one base to start with the 1st non-spacer base
            uint ind1=2*(isj*nIndicesSJ1+istart1);
            if (sjdbInd>=0 || seq1[0][istart]>3) 
            {//no index for already included junctions, or suffices starting with N
                indArray[ind1]=-1;
            } else 
            {
                //indArray[ind1] =  suffixArraySearch(seq1, istart, P->sjdbLength-istart1, G, SA, true, 0, P->nSA-1, 0, P) ;
                indArray[ind1] =  suffixArraySearch1(seq1, istart, 10000, G, -1LLU, SA, true, 0, P->nSA-1, 0, P) ;
                indArray[ind1+1] = isj*P->sjdbLength+istart;
            };
        };
    };
    sjNew = sjNew/2;//novel junctions were double counted on two strands
    
    time ( &rawtime );
    P->inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished SA search: number of new junctions=" << sjNew <<", old junctions="<<P->sjdbN-sjNew<<endl;
    
    uint nInd=0;//true number of new indices
    for (uint ii=0; ii<2*P->sjdbN*nIndicesSJ1; ii++) {//remove entries that cannot be inserted, this cannot be done in the parallel cycle above
        if (indArray[ii*2]!= (uint) -1) {
            indArray[nInd*2]=indArray[ii*2];
            indArray[nInd*2+1]=indArray[ii*2+1];
            ++nInd;
        };
    };

    globalGenomeArray=Gsj;
    qsort((void*) indArray, nInd, 2*sizeof(uint64), funCompareUintAndSuffixes);
    time ( &rawtime );
    P->inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished sorting SA indicesL nInd="<<nInd <<endl;

    indArray[2*nInd]=-999; //mark the last junction
    indArray[2*nInd+1]=-999; //mark the last junction
    
    P->nGenome=P->chrStart[P->nChrReal]+nGsj;    
    P->nSA+=nInd;
    
    uint GstrandBit1 = (uint) floor(log(P->nGenome)/log(2))+1;
    if (GstrandBit1<32) GstrandBit1=32; //TODO: use simple access function for SA
    if ( GstrandBit1 != P->GstrandBit) 
    {//too many junctions were added - GstrandBit changed
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: cannot insert junctions on the fly because of strand GstrandBit problem\n";
        errOut << "SOLUTION: please contact STAR author at https://groups.google.com/forum/#!forum/rna-star\n";
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_GENOME_FILES, *P);
    };
    
    SA2.defineBits(P->GstrandBit+1,P->nSA);
    uint nGsjNew=sjNew*P->sjdbLength; //this is the actual number of bytes added to the genome, while nGsj is the total size of all junctions
    
    uint N2bit= 1LLU << P->GstrandBit;
    uint strandMask=~N2bit;
    
    /*testing
    PackedArray SAo;
    SAo.defineBits(P->GstrandBit+1,P->nSA);
    SAo.allocateArray();
    ifstream oldSAin("./DirTrue/SA");
    oldSAin.read(SAo.charArray,SAo.lengthByte);
    oldSAin.close();
    */
    
    uint isj=0, isa2=0;
    for (uint isa=0;isa<P1->nSA;isa++) {   
        while (isa==indArray[isj*2]) {//insert sj index before the existing index
            uint ind1=indArray[isj*2+1];
            if (ind1<nGsj) {
                ind1+=P->chrStart[P->nChrReal];
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
            uint ind1s = P1->nGenome - (ind1 & strandMask);
            if (ind1s>P->chrStart[P->nChrReal])
            {//this index was an old sj, may need to shift it
                uint sj1 = (ind1s-P->chrStart[P->nChrReal])/P->sjdbLength;//old junction index
                ind1s += (oldSJind[sj1]-sj1)*P->sjdbLength;
                ind1 = (P->nGenome - ind1s) | N2bit;
            } else
            {
                ind1+=nGsjNew; //reverse complementary indices are all shifted by the length of junctions
            };
        } else
        {//+ strand
            if (ind1>P->chrStart[P->nChrReal])
            {//this index was an old sj, may need to shift it
                uint sj1 = (ind1-P->chrStart[P->nChrReal])/P->sjdbLength;//old junction index
                ind1 += (oldSJind[sj1]-sj1)*P->sjdbLength;
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
            ind1+=P->chrStart[P->nChrReal];
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
    P->inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished inserting junction indices" <<endl;
    
    //SAi insertions
    for (uint iL=0; iL < P->genomeSAindexNbases; iL++) {
        uint iSJ=0;
        uint ind0=P->genomeSAindexStart[iL]-1;//last index that was present in the old genome
        for (uint ii=P->genomeSAindexStart[iL];ii<P->genomeSAindexStart[iL+1]; ii++) {//scan through the longest index
            uint iSA1=SAi[ii];
            uint iSA2=iSA1 & P->SAiMarkNmask & P->SAiMarkAbsentMask;
            
            if ( iSJ<nInd && (iSA1 &  P->SAiMarkAbsentMaskC)>0 ) 
            {//index missing from the old genome
                uint iSJ1=iSJ;
                int64 ind1=funCalcSAi(Gsj+indArray[2*iSJ+1],iL);
                while (ind1 < (int64)(ii-P->genomeSAindexStart[iL]) && indArray[2*iSJ]-1<iSA2) {
                    ++iSJ;
                    ind1=funCalcSAi(Gsj+indArray[2*iSJ+1],iL);
                };
                if (ind1 == (int64)(ii-P->genomeSAindexStart[iL]) ) {
                    SAi.writePacked(ii,indArray[2*iSJ]-1+iSJ+1);
                    for (uint ii0=ind0+1; ii0<ii; ii0++) {//fill all the absent indices with this value
                        SAi.writePacked(ii0,(indArray[2*iSJ]-1+iSJ+1) | P->SAiMarkAbsentMaskC);
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
                    if (funCalcSAi(Gsj+indArray[2*iSJ+1],iL) >= (int64) (ii-P->genomeSAindexStart[iL]) ) {//this belongs to the next index
                        break;
                    };
                    ++iSJ;
                };   
                
                SAi.writePacked(ii,iSA1+iSJ);
                
                for (uint ii0=ind0+1; ii0<ii; ii0++) {//fill all the absent indices with this value
                    SAi.writePacked(ii0,(iSA2+iSJ) | P->SAiMarkAbsentMaskC);
                };
                ind0=ii;
            };
        };

    };

    for (uint isj=0;isj<nInd;isj++) {
        int64 ind1=0;
        for (uint iL=0; iL < P->genomeSAindexNbases; iL++) {
            uint g=(uint) Gsj[indArray[2*isj+1]+iL];
            ind1 <<= 2;
            if (g>3) {//this iSA contains N, need to mark the previous
                for (uint iL1=iL; iL1 < P->genomeSAindexNbases; iL1++) {
                    ind1+=3;
                    int64 ind2=P->genomeSAindexStart[iL1]+ind1;
                    for (; ind2>=0; ind2--) {//find previous index that is not absent
                        if ( (SAi[ind2] & P->SAiMarkAbsentMaskC)==0 ) {
                            break;
                        };
                    };
                    SAi.writePacked(ind2,SAi[ind2] | P->SAiMarkNmaskC);
                    ind1 <<= 2;
                };
                break;
            } else {
                ind1 += g;
            };
        };
    };
    time ( &rawtime );
    P->inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished SAi" <<endl;
    
    //change parameters, most parameters are already re-defined in sjdbPrepare.cpp
    SA.defineBits(P->GstrandBit+1,P->nSA);//same as SA2
    SA.pointArray(SA2.charArray);
    P->nSAbyte=SA.lengthByte;
    P->sjGstart=P->chrStart[P->nChrReal];
    memcpy(G+P->chrStart[P->nChrReal],Gsj, nGsj);
    
    /* testing
    PackedArray SAio=SAi;
    SAio.allocateArray();
    ifstream oldSAiin("./DirTrue/SAindex");
    oldSAiin.read(SAio.charArray,8*(P->genomeSAindexNbases+2));//skip first bytes
    oldSAiin.read(SAio.charArray,SAio.lengthByte);
    oldSAiin.close();  
    

    for (uint ii=0;ii<P->nSA;ii++) {
        if (SA2[ii]!=SAo[ii]) {
            cout <<ii <<" "<< SA2[ii]<<" "<<SAo[ii]<<endl;
        };
    };


    for (uint iL=0; iL < P->genomeSAindexNbases; iL++) {
        for (uint ii=P->genomeSAindexStart[iL];ii<P->genomeSAindexStart[iL+1]; ii++) {//scan through the longets index
                if ( SAio[ii]!=SAi[ii] ) {
                    cout <<ii<<" "<<SAio[ii]<<" "<<SAi[ii]<<endl;
                };
        };
    };    
    */
    
    /*
    ofstream genomeOut("/home/dobin/Genome");
    fstreamWriteBig(genomeOut,G,P->nGenome+nGsj,"777","777",P);
    genomeOut.close(); 
    genomeOut.open("/home/dobin/SA");
    fstreamWriteBig(genomeOut,SA2.charArray,SA2.lengthByte,"777","777",P);
    genomeOut.close();
    */
    
    delete [] indArray;
    delete [] G1c;
    delete [] oldSJind;       
    
};
