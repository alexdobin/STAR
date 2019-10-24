#include "Transcript.h"
bool Transcript::transformGenome(Genome &genOut, Transcript & A)
{
    uint32 nB=0;//number of out blocks
    auto &coBl=genOut.genomeOut.convBlocks;
    
    //uint64 icb=0;
    
    auto cBit=coBl.begin();
    for (uint32 iA=0; iA<nExons; iA++) {//loop over blocks
        
        uint64 r1=exons[iA][EX_R];//r start of a-block
        uint64 len=exons[iA][EX_L];        
        uint64 g1=exons[iA][EX_G];//g start of a-block
        uint64 g2=g1+len-1;//end of a-block
        
        
        //find c-block whose start is on the left of a-block
        array<uint64,3> g1array =  {g1,0,0};
        cBit = std::upper_bound(cBit, coBl.end(), g1array,
                               [](const array<uint64,3> &t1, const array<uint64,3> &t2)
                               {
                                  return t1[0] < t2[0];
                               });
        
        --cBit;//upperbound finds element > value, need to step back
        //icb=binarySearchStride(coBl.data(), coBl.size(), exons[iA][EX_G], icb, 3);

               
        uint64 b1=(*cBit)[0];
        uint64 b2=(*cBit)[0]+(*cBit)[1]-1;
        uint64 b1o=(*cBit)[2];
        
        //    **********---********----********
        //    b1      b2
        //        ++++++++++++++++
        //        g1            g2
        
        if (g1 <= b2 ) {//fill the first o-block
            A.exons[nB][EX_iFrag]=exons[iA][EX_iFrag];
            A.exons[nB][EX_R] = r1;
            A.exons[nB][EX_G] = b1o + g1 - b1;
            
            if (g2<=b2) {//g2 inside this c-block
                A.exons[nB][EX_L]=len;
            } else {//g2 is past the end of this c-block
                A.exons[nB][EX_L]=b2-g1+1;
            };
//             if (genOut.genomeOut.gapsAreJunctions && g2>=b2) {//add 0-lenght block to mark splice junction
//                     ++nB;
//                     A.exons[nB][EX_R] = A.exons[nB-1][EX_R]+A.exons[nB-1][EX_L];
//                     A.exons[nB][EX_G] = (*cBit)[2]+(*cBit)[1]; //end of this c-block in new genome
//                     A.exons[nB][EX_L] = 0;
//                     A.canonSJ[nB]=1;
//                     A.sjAnnot[nB]=1;  
//                     ++nB;
//                     A.exons[nB][EX_R] = A.exons[nB-1][EX_R];
//                     A.exons[nB][EX_G] = (*(cBit+1))[2]; //start of the next c-block in new genome
//                     A.exons[nB][EX_L] = 0;
//                     A.canonSJ[nB]=-1;
//                     A.sjAnnot[nB]=0;  
//             };
            ++nB;
        };
        
        ++cBit;
        while (g2 >= (*cBit)[0]) {//until c-block start is to the right of g2 
            A.exons[nB][EX_iFrag]=exons[iA][EX_iFrag];            
            A.exons[nB][EX_G]=(*cBit)[2];
            A.exons[nB][EX_R]=r1+(*cBit)[0]-g1;
            
            if (g2 < (*cBit)[0]+(*cBit)[1]) {//g2 inside this c-block
                A.exons[nB][EX_L]=g2-(*cBit)[0]+1;
            } else {//g2 is past the end of this c-block
                A.exons[nB][EX_L]=(*cBit)[1];//full block length
            };
//             if (genOut.genomeOut.gapsAreJunctions && g2 >= (*cBit)[0]+(*cBit)[1]-1) {//add 0-lenght block to mark splice junction
//                 ++nB;
//                 A.exons[nB][EX_R] = A.exons[nB-1][EX_R]+A.exons[nB-1][EX_L];
//                 A.exons[nB][EX_G] = (*cBit)[2]+(*cBit)[1]; //end of this c-block in new genome
//                 A.exons[nB][EX_L] = 0;
//                 A.canonSJ[nB]=1;
//                 A.sjAnnot[nB]=1;  
//                 ++nB;
//                 A.exons[nB][EX_R] = A.exons[nB-1][EX_R];
//                 A.exons[nB][EX_G] = (*(cBit+1))[2]; //start of the next c-block in new genome
//                 A.exons[nB][EX_L] = 0;
//                 A.canonSJ[nB]=-1;
//                 A.sjAnnot[nB]=0;  
//             };
            ++nB;
            ++cBit;
        };
        --cBit;
        
    };
    
    A.nExons=nB;
    if (nB>0) {
        A.Str = Str;

//         if (A.exons[0][EX_G]>=genOut.nGenome) {//convert to +strand
//             for (uint32 ib=0; ib<nB; ib++) {
//                 A.exons[ib][EX_G]=2*genOut.nGenome-A.exons[ib][EX_G]-A.exons[ib][EX_L];
//                 A.exons[ib][EX_R]=Lread-1-A.exons[ib][EX_R]-A.exons[ib][EX_L];
//             };
//             for (uint32 ib=0; ib<nB/2; ib++) {
//                 swap(A.exons[ib],A.exons[nB-ib-1]);
//             };
//         };
//         
        
        A.Chr = genOut.chrBin[A.exons[0][EX_G] >> genOut.pGe.gChrBinNbits];
    };
    
    {//recalculate canonSJ, sjAnnot
        A.sjAnnot[0]=0;
        A.canonSJ[0]=-1;
        for (uint64 ia=0; ia<A.nExons-1; ia++) {
            A.sjAnnot[ia]=0;//TODO: check annotations in the new genome coordinates
            if (A.exons[ia+1][EX_iFrag]!=A.exons[ia][EX_iFrag]) {
                A.canonSJ[ia]=-3;
                continue;
            };
            uint64 gapG=A.exons[ia+1][EX_G]-(A.exons[ia][EX_G]+A.exons[ia][EX_L]);
            uint64 gapR=A.exons[ia+1][EX_R]-A.exons[ia][EX_R]-A.exons[ia][EX_L];
            if (gapR>0) {//insertion
                A.canonSJ[ia]=-2;
            } else if (gapG>=genOut.P.alignIntronMin) {//junction
                A.canonSJ[ia]=0;
            } else {//deletion
                A.canonSJ[ia]=-1;
            };
        };
    };
    
    return (nB>0);
};

