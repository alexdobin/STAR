#include "Transcript.h"
#include "SequenceFuns.h"

string Transcript::generateCigarP() 
{//generates CIGARp string for the transcript with "p" operation for PE reads
 //p is a special CIGAR operation to encode gap between mates. This gap is negative for overlapping mates.

    string CIGAR;
    stringstream samStreamCIGAR;

    uint leftMate=0;
    if (readNmates>1) leftMate=Str;

    uint trimL=exons[0][EX_R] - (exons[0][EX_R]<readLengthOriginal[leftMate] ? 0 : readLengthOriginal[leftMate]+1);
    if (trimL>0) {
        samStreamCIGAR << trimL << "S"; //initial trimming
    };

    for (uint ii=0;ii<nExons;ii++) {//cycle over all exons, record CIGAR
        if (ii>0) {//record gaps
            uint gapG=exons[ii][EX_G]-(exons[ii-1][EX_G]+exons[ii-1][EX_L]);

            if (exons[ii][EX_G] >= (exons[ii-1][EX_G]+exons[ii-1][EX_L]) ) {//
                if (canonSJ[ii-1]==-3) {//gap between mates
                    //soft clipping of the second mate
                    uint s1=readLengthOriginal[leftMate]-(exons[ii-1][EX_R]+exons[ii-1][EX_L]);
                    uint s2=exons[ii][EX_R]-(readLengthOriginal[leftMate]+1);
                    if (s1>0){
                        samStreamCIGAR << s1 << "S";
                    };
                    samStreamCIGAR << gapG << "p";
                    if (s2>0){
                        samStreamCIGAR << s2 << "S";
                    };

                } else {
                    //it's possible to have a D or N and I for at the same time
                    uint gapR=exons[ii][EX_R]-exons[ii-1][EX_R]-exons[ii-1][EX_L]; //gapR>0 always
                    if (gapR>0){
                        samStreamCIGAR << gapR << "I";
                    };
                    if (canonSJ[ii-1]>=0 || sjAnnot[ii-1]==1) {//junction: N
                        samStreamCIGAR << gapG << "N";
                    } else if (gapG>0) {//deletion
                        samStreamCIGAR << gapG << "D";
                    };
                };
            } else {//mates overlap
                samStreamCIGAR << "-" << (exons[ii-1][EX_G]+exons[ii-1][EX_L]) - exons[ii][EX_G] << "p";
            };
        };
        samStreamCIGAR << exons[ii][EX_L] << "M";
    };


    trimL=(exons[nExons-1][EX_R]<readLengthOriginal[leftMate] ? readLengthOriginal[leftMate] : readLengthPairOriginal) - exons[nExons-1][EX_R]-exons[nExons-1][EX_L];
    if ( trimL > 0 ) {
        samStreamCIGAR << trimL << "S"; //final trimming
    };
    CIGAR=samStreamCIGAR.str();

    return CIGAR;
};
