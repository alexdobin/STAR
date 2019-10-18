#include "Transcript.h"
bool Transcript::convertGenomeCigar(Genome &genOut, Transcript & A)
{
    auto &coBl=genOut.genomeOut.convBlocks;    
    auto cBit=coBl.begin();
    cBit = std::upper_bound(cBit, coBl.end(), array<uint64,3> {gStart,0,0},
                               [](const array<uint64,3> &t1, const array<uint64,3> &t2)
                               {
                                  return t1[0] < t2[0];
                               });
    --cBit;
    
    uint64 icB = cBit-coBl.begin();
    
    A.gStart=gStart-coBl[icB][0]+coBl[icB][2];
    
    A.cigar.clear();
    A.cigar.reserve(cigar.size()+100); //add 100 for junctions
    
    uint64 gEnd=gStart;//gEnd is end+1
    uint32 l1=0, l2=0;
    for (auto &cc : cigar) {
        if (cc[0]==BAM_CIGAR_M)
            l1+=cc[1];
        if (cc[0]==BAM_CIGAR_M || cc[0]==BAM_CIGAR_D || cc[0]==BAM_CIGAR_N) {//M,D,N: blocks that have to be split if they overlap blocks boundaries
            uint64 gStart1=gEnd;
            gEnd += cc[1];
            uint64 bE;
            while (( bE=coBl[icB][0]+coBl[icB][1] ) <= gEnd) {//MDN-block ends at the c-block end, or goes past
                uint32 len1 = (uint32)( bE-gStart1 );//first part of the a-block
                
                if (cc[0]==BAM_CIGAR_M)
                    l2+=len1;
                
                if (cc[0]==A.cigar.back()[0]) {//this can only happen for N
                    A.cigar.back()[1] += len1;//add to previous gap
                } else {
                    A.cigar.push_back({cc[0], len1});
                };
                
                len1 = (uint32)( coBl[icB+1][2]-coBl[icB][2]-coBl[icB][1] );//insert junction gap
                
                if (A.cigar.back()[0]==BAM_CIGAR_N) {
                    A.cigar.back()[1] += len1;//add to previous gap
                } else {
                    A.cigar.push_back({BAM_CIGAR_N, len1});
                };
                
                gStart1=bE;
                ++icB;
            };
            
            if (gEnd>gStart1)
                A.cigar.push_back({ cc[0], (uint32) (gEnd-gStart1) });
              
            if(cc[0]==BAM_CIGAR_M)
                l2 += gEnd-gStart1;
            if (l1!=l2)
                cout << l1 <<" "<< l2 <<endl;                
            
            if (gEnd<coBl[icB][0]) //step back
                --icB;
                
        } else {//S,I
            l1+=cc[1];
            l2+=cc[1];
            A.cigar.push_back(cc);
        };
    };

    if (A.cigar.back()[0] == BAM_CIGAR_N)
        A.cigar.pop_back();
    
    A.Str = Str;

    if (A.gStart >= genOut.nGenome) {//convert to +strand
        A.Str = 1-A.Str;
        A.gStart = 2*genOut.nGenome - (A.gStart+gEnd-gStart);
        std::reverse(A.cigar.begin(), A.cigar.end());
    };
        
    A.Chr = genOut.chrBin[A.gStart >> genOut.pGe.gChrBinNbits];
    
    return (true); //conversion is always possible
};
