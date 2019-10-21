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
    //uint32 l1=0, l2=0;
    //uint32 gGapNew=0;
    for (uint32 iC=0; iC<cigar.size(); iC++) {
        auto cc=cigar[iC];        
        //if (cc[0]==BAM_CIGAR_M)
        //    l1+=cc[1];
        if (cc[0]==BAM_CIGAR_M || cc[0]==BAM_CIGAR_D || cc[0]==BAM_CIGAR_N) {//M,D,N: blocks that have to be split if they overlap blocks boundaries
            uint64 gStart1=gEnd;
            gEnd += cc[1];
            uint64 bE;
            while (( bE=coBl[icB][0]+coBl[icB][1] ) <= gEnd) {//MDN-block ends at the c-block end, or goes past
                uint32 len1 = (uint32)( bE-gStart1 );//first part of the a-block
                
                //if (cc[0]==BAM_CIGAR_M)
                //    l2+=len1;

                A.cigar.push_back({cc[0], len1});
                
                len1 = (uint32)( coBl[icB+1][2]-coBl[icB][2]-coBl[icB][1] );//insert junction gap
                //gGapNew += len1;
                
                A.cigar.push_back({BAM_CIGAR_N, len1});
                
                gStart1=bE;
                ++icB;
            };
            
            if (gEnd>gStart1)
                A.cigar.push_back({ cc[0], (uint32) (gEnd-gStart1) });
              
            //if(cc[0]==BAM_CIGAR_M)
            //    l2 += gEnd-gStart1;
            //if (l1!=l2) {
            //    cout << l1 <<" "<< l2 <<endl;                
            //};
            
            if (gEnd<coBl[icB][0]) //step back
                --icB;
                
        } else {//S,I
            //l1+=cc[1];
            //l2+=cc[1];
            A.cigar.push_back(cc);
        };        
    };
    //if (l1!=l2) {
    //    cout << l1 <<" "<< l2 <<endl;                
    //};
    
    {//merge consecutive identical operations
        auto c0=A.cigar.begin();
        for (auto c1=A.cigar.begin()+1; c1!=A.cigar.end(); c1++) {
            if ( c1[0][0]==c0[0][0] ) {//add length to the prev element
                c0[0][1] += c1[0][1];
            } else if ( c1[0][0]==BAM_CIGAR_N && c0[0][0]==BAM_CIGAR_I && c0[-1][0]==BAM_CIGAR_N){
                //NIN : add junction to the junction before I, since I placement does not matter
                //relying on I operation not being the first one, so c0 is at least cigar.begin+1
                c0[-1][1] += c1[0][1];
            } else{
                c0++;
                *c0=*c1;
            };
        };
        A.cigar.resize(c0+1-A.cigar.begin());
    };
     
    {//remove last junction, or the junction before last S - it could appear if the last cigar M lands at the end of the c-block
        if (A.cigar.back()[0] == BAM_CIGAR_N) {
            //gGapNew-=A.cigar.back()[1];
            A.cigar.pop_back(); 
        };
        if (A.cigar.back()[0] == BAM_CIGAR_S && (*(A.cigar.end()-2))[0] == BAM_CIGAR_N) {
            //gGapNew-=(*(A.cigar.end()-2))[1];
            A.cigar.erase(A.cigar.end()-2);
        };
    };
    
    {//remove start junctions with too short overhangs
        //TODO: loop over until all such junctions are removed
        array<uint32,5> nOp={0,0,0,0,0};
        for (auto c1=A.cigar.begin(); c1!=A.cigar.end(); c1++) {
            nOp[(*c1)[0]] += (*c1)[1];
            if ( (*c1)[0] == BAM_CIGAR_M ) {
                if (nOp[BAM_CIGAR_M] >= genOut.P.alignSJDBoverhangMin) //TODO Transcript class should have P
                    break; //overhang is long enough
            } else if ( (*c1)[0]==BAM_CIGAR_N ) {//short overhang for this junction
                //find the next M
                c1++;
                while ( (*c1)[0]!=BAM_CIGAR_M ) {
                    nOp[(*c1)[0]] += (*c1)[1];
                    c1++;
                };
                //remove junction
                A.cigar.erase(A.cigar.begin(), c1-1);//keep one element at the beginning to convert it to S
                A.cigar.at(0)={BAM_CIGAR_S, nOp[BAM_CIGAR_M]+nOp[BAM_CIGAR_S]+nOp[BAM_CIGAR_I]};
                A.gStart += nOp[BAM_CIGAR_M]+nOp[BAM_CIGAR_D]+nOp[BAM_CIGAR_N];
                break;
            };
        };
    };
    
    {//remove terminal junctions with too short overhangs
        //TODO: loop over until all such junctions are removed
        array<uint32,5> nOp={0,0,0,0,0};
        for (auto c1=A.cigar.end()-1; c1!=A.cigar.begin(); c1--) {
            nOp[(*c1)[0]] += (*c1)[1];
            if ( (*c1)[0] == BAM_CIGAR_M ) {
                if (nOp[BAM_CIGAR_M] >= genOut.P.alignSJDBoverhangMin) //TODO Transcript class should have P
                    break; //overhang is long enough
            } else if ( (*c1)[0]==BAM_CIGAR_N ) {//short overhang for this junction
                //find the next M
                c1--;
                while ( (*c1)[0]!=BAM_CIGAR_M ) {
                    nOp[(*c1)[0]] += (*c1)[1];
                    c1--;
                };
                //remove junction
                A.cigar.erase(c1+1, A.cigar.end());
                A.cigar.push_back({BAM_CIGAR_S, nOp[BAM_CIGAR_M]+nOp[BAM_CIGAR_S]+nOp[BAM_CIGAR_I]});
                break;
            };
        };
    };
    
    //TODO: also need to recalculate nMM
    A.gLength=0;
    for (const auto &cc : A.cigar)
        if (cc[0]==BAM_CIGAR_M || cc[0]==BAM_CIGAR_D || cc[0]==BAM_CIGAR_N)
            A.gLength += cc[1];
        
    A.Str = Str;
    if (A.gStart >= genOut.genomeOut.nMinusStrandOffset) {//convert to +strand
        A.Str = 1-A.Str;
        A.gStart = 2*genOut.genomeOut.nMinusStrandOffset - (A.gStart+A.gLength);
        std::reverse(A.cigar.begin(), A.cigar.end());
    };
        
    A.Chr = genOut.chrBin[A.gStart >> genOut.pGe.gChrBinNbits];
    
    return (true); //conversion is always possible
};
