#include "ReadAlign.h"
#include "SequenceFuns.h"

string ReadAlign::outputTranscriptCIGARp(Transcript const &trOut) {//generates CIGARp string for the transcript
                                                                   //p is a special CIGAR operation to encode gap between mates. This gap is negative for overlapping mates

    string CIGAR;
    samStreamCIGAR.str(std::string());

    uint leftMate=0;
    if (P.readFilesIn.size()>1) leftMate=trOut.Str;

    uint trimL=trOut.exons[0][EX_R] - (trOut.exons[0][EX_R]<readLengthOriginal[leftMate] ? 0 : readLengthOriginal[leftMate]+1);
    if (trimL>0) {
        samStreamCIGAR << trimL << "S"; //initial trimming
    };

    for (uint ii=0;ii<trOut.nExons;ii++) {//cycle over all exons, record CIGAR
        if (ii>0) {//record gaps
            uint gapG=trOut.exons[ii][EX_G]-(trOut.exons[ii-1][EX_G]+trOut.exons[ii-1][EX_L]);

            if (trOut.exons[ii][EX_G] >= (trOut.exons[ii-1][EX_G]+trOut.exons[ii-1][EX_L]) ) {//
                if (trOut.canonSJ[ii-1]==-3) {//gap between mates
                    //soft clipping of the second mate
                    uint s1=readLengthOriginal[leftMate]-(trOut.exons[ii-1][EX_R]+trOut.exons[ii-1][EX_L]);
                    uint s2=trOut.exons[ii][EX_R]-(readLengthOriginal[leftMate]+1);
                    if (s1>0){
                        samStreamCIGAR << s1 << "S";
                    };
                    samStreamCIGAR << gapG << "p";
                    if (s2>0){
                        samStreamCIGAR << s2 << "S";
                    };

                } else {
                    //it's possible to have a D or N and I for at the same time
                    uint gapR=trOut.exons[ii][EX_R]-trOut.exons[ii-1][EX_R]-trOut.exons[ii-1][EX_L]; //gapR>0 always
                    if (gapR>0){
                        samStreamCIGAR << gapR << "I";
                    };
                    if (trOut.canonSJ[ii-1]>=0 || trOut.sjAnnot[ii-1]==1) {//junction: N
                        samStreamCIGAR << gapG << "N";
                    } else if (gapG>0) {//deletion
                        samStreamCIGAR << gapG << "D";
                    };
                };
            } else {//mates overlap
                samStreamCIGAR << "-" << (trOut.exons[ii-1][EX_G]+trOut.exons[ii-1][EX_L]) - trOut.exons[ii][EX_G] << "p";
            };
        };
        samStreamCIGAR << trOut.exons[ii][EX_L] << "M";
    };


    trimL=(trOut.exons[trOut.nExons-1][EX_R]<readLengthOriginal[leftMate] ? readLengthOriginal[leftMate] : readLengthPairOriginal) - trOut.exons[trOut.nExons-1][EX_R]-trOut.exons[trOut.nExons-1][EX_L];
    if ( trimL > 0 ) {
        samStreamCIGAR << trimL << "S"; //final trimming
    };
    CIGAR=samStreamCIGAR.str();

    return CIGAR;



};
