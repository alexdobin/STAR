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
    uint32 gGapNew=0;
    for (uint32 iC=0; iC<cigar.size(); iC++) {
        auto cc=cigar[iC];        
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
                gGapNew += len1;
                
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
            if (l1!=l2) {
                cout << l1 <<" "<< l2 <<endl;                
            };
            
            if (gEnd<coBl[icB][0]) //step back
                --icB;
                
        } else {//S,I
            l1+=cc[1];
            l2+=cc[1];
            A.cigar.push_back(cc);
        };        
    };
    if (l1!=l2) {
        cout << l1 <<" "<< l2 <<endl;                
    };
    
    //remove last junction, or the junction before last S - it could appear if the last cigar M lands at the end of the c-block
    if (A.cigar.back()[0] == BAM_CIGAR_N) {
        gGapNew-=A.cigar.back()[1];
        A.cigar.pop_back(); 
    };
    if (A.cigar.back()[0] == BAM_CIGAR_S && (*(A.cigar.end()-2))[0] == BAM_CIGAR_N) {
        gGapNew-=(*(A.cigar.end()-2))[1];
        A.cigar.erase(A.cigar.end()-2);
    };
    
    A.gLength=gEnd-gStart+gGapNew;
    A.Str = Str;

    if (A.gStart >= genOut.genomeOut.nMinusStrandOffset) {//convert to +strand
        A.Str = 1-A.Str;
        A.gStart = 2*genOut.genomeOut.nMinusStrandOffset - (A.gStart+A.gLength);
        std::reverse(A.cigar.begin(), A.cigar.end());
    };
        
    A.Chr = genOut.chrBin[A.gStart >> genOut.pGe.gChrBinNbits];
    
    return (true); //conversion is always possible
};
