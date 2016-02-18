#include "Transcriptome.h"
#include "ReadAlign.h"
#include "serviceFuns.cpp"


int alignToTranscript(Transcript &aG, uint trS1, uint8 trStr1, uint32 *exSE1, uint32 *exLenCum1, uint16 exN1, Transcript &aT) {

    //find exon that overlaps beginning of the read
    uint32 g1=aG.exons[0][EX_G]-trS1;//start of the transcript
    uint32 ex1=binarySearch1<uint32>(g1, exSE1, 2*exN1);
    if (ex1>=2*exN1) return 0; //align start is to the right of all exons

    if (ex1%2==1) {//beginning of the read >=end of an exon
        if (exSE1[ex1]==g1) {//first base of the read is exactly the last base of the exon
            --ex1;
        } else {
            return 0;//beginning of the read is past the end of an exon, align does not belong to this transcript
        };
    };
    ex1=ex1/2; //this is the first exon of the alignment

    aT.nExons=0;
    aT.primaryFlag=false;

    aG.canonSJ[aG.nExons-1]=-999; //marks the last exons
    for (uint32 iab=0; iab<aG.nExons; iab++) {//scan through all blocks of the align

//         g1+=aG.exons[iab][EX_L]-1;//last base of the block
        if (aG.exons[iab][EX_G]+aG.exons[iab][EX_L]>exSE1[2*ex1+1]+trS1+1) {//block extends past exon end
            return 0;
        };

        if (iab==0 || aG.canonSJ[iab-1]<0) {
            aT.exons[aT.nExons][EX_R]=aG.exons[iab][EX_R];
            aT.exons[aT.nExons][EX_G]=aG.exons[iab][EX_G]-trS1-exSE1[2*ex1]+exLenCum1[ex1];
            aT.exons[aT.nExons][EX_L]=aG.exons[iab][EX_L];
            aT.exons[aT.nExons][EX_iFrag]=aG.exons[iab][EX_iFrag];
            if (aT.nExons>0) aT.canonSJ[aT.nExons-1]=aG.canonSJ[iab-1];
            ++aT.nExons;
        } else {
            aT.exons[aT.nExons-1][EX_L]+=aG.exons[iab][EX_L];
        };
        switch (aG.canonSJ[iab]) {
            case -999: //last exon
                if (trStr1==2) {//convert align coordinates if on the -strand
                    uint32 trlength=exLenCum1[exN1-1]+exSE1[2*exN1-1]-exSE1[2*exN1-2]+1; //transcript length
                    for (uint32 iex=0; iex<aT.nExons; iex++)
                    {
                        aT.exons[iex][EX_R]=aG.Lread-(aT.exons[iex][EX_R]+aT.exons[iex][EX_L]);
                        aT.exons[iex][EX_G]=trlength-(aT.exons[iex][EX_G]+aT.exons[iex][EX_L]);
                    };
                    for (uint32 iex=0; iex<aT.nExons/2; iex++)
                    {
                        swap(aT.exons[iex][EX_R],aT.exons[aT.nExons-1-iex][EX_R]);
                        swap(aT.exons[iex][EX_G],aT.exons[aT.nExons-1-iex][EX_G]);
                        swap(aT.exons[iex][EX_L],aT.exons[aT.nExons-1-iex][EX_L]);
                        swap(aT.exons[iex][EX_iFrag],aT.exons[aT.nExons-1-iex][EX_iFrag]);
                    };
                    for (uint32 iex=0; iex<(aT.nExons-1)/2; iex++)
                    {
                        swap(aT.canonSJ[iex],aT.canonSJ[aT.nExons-2-iex]);
                    };

                };
                for (uint32 iex=0; iex<aT.nExons; iex++)
                {//no junctions in the transcritomic coordinates
                    aT.sjAnnot[iex]=0;
                    aT.shiftSJ[iex][0]=0;
                    aT.shiftSJ[iex][1]=0;
                    aT.sjStr[iex]=0;
                };

                return 1; //reached the end of blocks, align is consistend with this transcript
                break;
            case -3: //mate connection
                ex1=binarySearch1<uint32>(aG.exons[iab+1][EX_G]-trS1, exSE1, 2*exN1);
                if (ex1%2==1) {//beginning of the mext mate in the middle of the exon?
                    return 0; //align does not belong to this transcript
                } else {
                    ex1=ex1/2; //this is the first exon of the second mate
                };
                break;
            case -2: //insertion
                break;
            case -1: //deletion
                break;
            default://junctions
                if ( aG.exons[iab][EX_G]+aG.exons[iab][EX_L]==exSE1[2*ex1+1]+trS1+1 && aG.exons[iab+1][EX_G]==exSE1[2*(ex1+1)]+trS1 ) {
                    //junction matches transcript junction
                    ++ex1;
                } else {
                    return 0;
                };
        };
    };
    return 0; //this should not happen
};

uint32 Transcriptome::quantAlign (Transcript &aG, Transcript *aTall) {
    uint32 nAtr=0; //number of alignments to the transcriptome

    //binary search through transcript starts
    uint32 tr1=binarySearch1a<uint>(aG.exons[0][EX_G], trS, nTr);
    if (tr1==(uint32) -1) return 0; //alignment outside of range of all transcripts

    uint aGend=aG.exons[aG.nExons-1][EX_G];

    ++tr1;
    do {//cycle back through all the transcripts
        --tr1;
        if (aGend<=trE[tr1]) {//this transcript contains the read
                int aStatus=alignToTranscript(aG, trS[tr1], trStr[tr1], exSE+2*trExI[tr1], exLenCum+trExI[tr1], trExN[tr1], aTall[nAtr]);
                if (aStatus==1) {//align conforms with the transcript
                    aTall[nAtr].Chr = tr1;
                    aTall[nAtr].Str = trStr[tr1]==1 ? aG.Str : 1-aG.Str; //TODO strandedness
                    ++nAtr;
                };
        };
    } while (trEmax[tr1]>=aGend && tr1>0);

    return nAtr;
};
