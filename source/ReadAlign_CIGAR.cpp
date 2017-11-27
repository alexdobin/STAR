#include "ReadAlign.h"

uint ReadAlign::alignCIGAR
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
                    samStreamSJmotif <<','<< trOut.canonSJ[ii-1] + (trOut.sjAnnot[ii-1]==0 ? 0 : SJ_SAM_AnnotatedMotifShift); //record junction type
//                     samStreamSJannot <<','<< (int) trOut.sjAnnot[ii-1]; //record annotation type
                    samStreamSJintron <<','<< trOut.exons[ii-1][EX_G] + trOut.exons[ii-1][EX_L] + 1 - P->chrStart[trOut.Chr] <<','\
                                   << trOut.exons[ii][EX_G] - P->chrStart[trOut.Chr]; //record intron loci
                } else if (gapG>0) {//deletion: N
                    samStreamCIGAR << gapG;
                    samStreamCIGAR << "D";
                };
            };
            samStreamCIGAR << trOut.exons[ii][EX_L] << "M";
        };

        string SJmotif = samStreamSJmotif.str();
        string SJintron = samStreamSJintron.str();
//         string SJannot = samStreamSJannot.str();

        if (SJmotif.length()==0) {//no junctions recorded, mark with -1
            SJmotif=",-1";
            SJintron=",-1";
//             SJannot=",-1";
        };

        uint trimR1=(trOut.exons[iEx1][EX_R]<readLength[leftMate] ? \
            readLengthOriginal[leftMate] : readLength[leftMate]+1+readLengthOriginal[Mate]) \
            - trOut.exons[iEx2][EX_R]-trOut.exons[iEx2][EX_L] - trimL;
        if ( trimR1 > 0 ) {
            samStreamCIGAR << trimR1 << "S"; //final trimming
        };
        CIGAR=samStreamCIGAR.str();