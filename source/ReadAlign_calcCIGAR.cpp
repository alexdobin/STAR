#include "ReadAlign.h"

void ReadAlign::calcCIGAR(Transcript const &trOut, uint nMates, uint iExMate, uint leftMate) {
    matesCIGAR.clear();
    for (uint imate=0;imate<nMates;imate++) {

        uint iEx1 = (imate==0 ? 0 : iExMate+1);
        uint iEx2 = (imate==0 ? iExMate : trOut.nExons-1);
        uint Mate=trOut.exons[iEx1][EX_iFrag];
        uint Str= trOut.Str;
        
        samStreamCIGAR.str(std::string());

        uint trimL;
        if (Str==0 && Mate==0) {
            trimL=clip5pNtotal[Mate];
        } else if (Str==0 && Mate==1) {
            trimL=clip3pNtotal[Mate];
        } else if (Str==1 && Mate==0) {
            trimL=clip3pNtotal[Mate];
        } else {
            trimL=clip5pNtotal[Mate];
        };

        uint trimL1 = trimL + trOut.exons[iEx1][EX_R] - (trOut.exons[iEx1][EX_R]<readLength[leftMate] ? 0 : readLength[leftMate]+1);
        if (trimL1>0) {
            samStreamCIGAR << trimL1 << "S"; //initial trimming
        };

        for (uint ii=iEx1;ii<=iEx2;ii++) {
            if (ii>iEx1) {//record gaps
                uint gapG=trOut.exons[ii][EX_G]-(trOut.exons[ii-1][EX_G]+trOut.exons[ii-1][EX_L]);
                uint gapR=trOut.exons[ii][EX_R]-trOut.exons[ii-1][EX_R]-trOut.exons[ii-1][EX_L];
                //it's possible to have a D or N and I at the same time
                if (gapR>0){
                    samStreamCIGAR << gapR;
                    samStreamCIGAR << "I";
                };
                if (trOut.canonSJ[ii-1]>=0 || trOut.sjAnnot[ii-1]==1) {//junction: N
                    samStreamCIGAR << gapG;
                    samStreamCIGAR << "N";
                } else if (gapG>0) {//deletion: N
                    samStreamCIGAR << gapG;
                    samStreamCIGAR << "D";
                };
            };
            samStreamCIGAR << trOut.exons[ii][EX_L] << "M";
        };

        uint trimR1=(trOut.exons[iEx1][EX_R]<readLength[leftMate] ? \
            readLengthOriginal[leftMate] : readLength[leftMate]+1+readLengthOriginal[Mate]) \
            - trOut.exons[iEx2][EX_R]-trOut.exons[iEx2][EX_L] - trimL;
        if ( trimR1 > 0 ) {
            samStreamCIGAR << trimR1 << "S"; //final trimming
        };
        matesCIGAR.push_back(samStreamCIGAR.str());
    };
};