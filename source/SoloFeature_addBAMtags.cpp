#include "SoloFeature.h"
#include "BAMfunctions.h"
#include "SequenceFuns.h"

void SoloFeature::addBAMtags(char *&bam0, uint32 &size0, char *bam1)
{//add extra tags to the BAM record
    
    uint64 iread = * ((uint64*) (bam0+size0));
    iread = iread >> 32; //iRead was encoded in the upper 32 bitsls
    //cout << iread <<" "<< convertNuclInt64toString(readInfo[iread].cb,  pSolo.cbL) <<" "<< convertNuclInt64toString(readInfo[iread].umi, pSolo.umiL)<<endl;

    if (readInfo[iread].cb + 1 != 0) {//otherwise do not add CB/UB tags            
        string cb  = pSolo.cbWLstr[readInfo[iread].cb];
        string umi = convertNuclInt64toString(readInfo[iread].umi, pSolo.umiL);
        memcpy(bam1, bam0, size0);

        size0 += bamAttrArrayWrite(cb,  "CB", bam1+size0);
        size0 += bamAttrArrayWrite(umi, "UB", bam1+size0);
        uint32 *bam1i = (uint32*) bam1;
        bam1i[0] = size0-sizeof(uint32);
        bam0=bam1;
    };
};