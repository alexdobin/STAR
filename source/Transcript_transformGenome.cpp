#include "Transcript.h"
#include "binarySearch2.h"

bool Transcript::transformGenome(Genome &genOut, Transcript & A)
{
    uint32 nB=0;//number of out blocks
    auto &coBl=genOut.genomeOut.convBlocks;
    
    //uint64 icb=0;
    
    //auto cBit=coBl.begin();
    for (uint32 iA=0; iA<nExons; iA++) {//loop over blocks
        
        uint64 r1=exons[iA][EX_R];//r start of a-block
        uint64 len=exons[iA][EX_L];        
        uint64 g1=exons[iA][EX_G];//g start of a-block
        uint64 g2=g1+len-1;//end of a-block
        
        
        //find c-block whose start is on the left of a-block
        array<uint64,3> g1array =  {g1,0,0};
        auto cBit=coBl.begin(); //TODO might reuse previous value of cBit, but need to take care of overlapping mates
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
        
        /*    **********---********----********
         *    b1      b2
         *        ++++++++++++++++
         *        g1            g2
        */
        
        if (g1 <= b2 ) {//fill the first o-block
            A.exons[nB][EX_iFrag]=exons[iA][EX_iFrag];
            A.exons[nB][EX_R] = r1;
            A.exons[nB][EX_G] = b1o + g1 - b1;
            
            if (g2<=b2) {//g2 inside this c-block
                A.exons[nB][EX_L]=len;
            } else {//g2 is past the end of this c-block
                A.exons[nB][EX_L]=b2-g1+1;
            };
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
            ++nB;
            ++cBit;
        };
        --cBit;
    };

    if (nB==0)
        return false; //transformation did not work

    {//merge exons w/o R/G gaps
        uint32 nB1 = 1;
        for ( uint32 ib=1; ib<nB; ib++) {
            
            if (nB1!=ib) {//copy ib into nB1
                for (uint32 ii=0; ii<EX_SIZE; ii++)
                    A.exons[nB1][ii] = A.exons[ib][ii];
            };
            
            if (A.exons[nB1][EX_iFrag]!=A.exons[nB1-1][EX_iFrag]) {//gap between mates, do not adjust
                ++nB1;
                continue;
            };

            uint64 gapR = A.exons[nB1][EX_R] - A.exons[nB1-1][EX_R] - A.exons[nB1-1][EX_L];
            uint64 gapG = A.exons[nB1][EX_G] - A.exons[nB1-1][EX_G] - A.exons[nB1-1][EX_L]; 
            
            if ( gapR == gapG) {//eliminate both gaps by increasing previous exon length
                A.exons[nB1-1][EX_L] += A.exons[nB1][EX_L] + gapR;

            } else {//new exon
                uint64 minGap = min(gapR, gapG);
                if (minGap>0) {//flush the gap to the left
                    A.exons[nB1][EX_L] += minGap;
                    A.exons[nB1][EX_G] -= minGap;
                    A.exons[nB1][EX_R] -= minGap;
                };
                ++nB1;
            };
        };
        nB = nB1;
    };
    
    A.nExons=nB;   
    A.Str = Str;        
    A.Chr = genOut.chrBin[A.exons[0][EX_G] >> genOut.pGe.gChrBinNbits];
    
    {//recalculate canonSJ, sjAnnot
        for (uint64 ia=0; ia<A.nExons-1; ia++) {
            A.sjAnnot[ia]=0;
            
            if (A.exons[ia+1][EX_iFrag]!=A.exons[ia][EX_iFrag]) {//mate gap
                A.canonSJ[ia]=-3;
                continue;
            };
            
            uint64 jS=A.exons[ia][EX_G]+A.exons[ia][EX_L];
            uint64 jE=A.exons[ia+1][EX_G]-1;//intron start/end
            int sjdbInd=binarySearch2(jS, jE, genOut.sjdbStart, genOut.sjdbEnd, genOut.sjdbN);
            if (sjdbInd>=0) {//annotated
                A.sjAnnot[ia] = 1;
                A.canonSJ[ia] = genOut.sjdbMotif[sjdbInd];
                if (genOut.sjdbMotif[sjdbInd]==0) {//shift to match annotations
                    if (A.exons[ia][EX_L] <= genOut.sjdbShiftLeft[sjdbInd]) {
                        return false; //this align is not allowed
                    };
                    A.exons[ia][EX_L] -= genOut.sjdbShiftLeft[sjdbInd];
                    A.exons[ia+1][EX_G] -= genOut.sjdbShiftLeft[sjdbInd];
                };

            } else {//unannotated
                uint64 gapG=jE-jS+1;
                uint64 gapR=A.exons[ia+1][EX_R]-A.exons[ia][EX_R]-A.exons[ia][EX_L];
                if (gapR>0) {//insertion
                    A.canonSJ[ia]=-2;
                    
                } else if (gapG>=genOut.P.alignIntronMin) {//junction
                    A.canonSJ[ia]=0;
                    if        ( genOut.G[jS]==2 && genOut.G[jS+1]==3 && genOut.G[jE-1]==0 && genOut.G[jE]==2 ) {//GTAG
                        A.canonSJ[ia]=1;
                    } else if ( genOut.G[jS]==1 && genOut.G[jS+1]==3 && genOut.G[jE-1]==0 && genOut.G[jE]==1 ) {//CTAC
                        A.canonSJ[ia]=2;
                    } else if ( genOut.G[jS]==2 && genOut.G[jS+1]==1 && genOut.G[jE-1]==0 && genOut.G[jE]==2 ) {//GCAG
                        A.canonSJ[ia]=3;
                    } else if ( genOut.G[jS]==1 && genOut.G[jS+1]==3 && genOut.G[jE-1]==2 && genOut.G[jE]==1 ) {//CTGC
                        A.canonSJ[ia]=4;
                    } else if ( genOut.G[jS]==0 && genOut.G[jS+1]==3 && genOut.G[jE-1]==0 && genOut.G[jE]==1 ) {//ATAC
                        A.canonSJ[ia]=5;
                    } else if ( genOut.G[jS]==2 && genOut.G[jS+1]==3 && genOut.G[jE-1]==0 && genOut.G[jE]==3 ) {//GTAT
                        A.canonSJ[ia]=6;
                    };                
                    
                } else {//deletion
                    A.canonSJ[ia]=-1;
                };
            };
        };
    };
    
    return true;
};

