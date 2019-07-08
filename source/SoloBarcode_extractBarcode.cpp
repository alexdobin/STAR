#include "SoloBarcode.h"

bool SoloBarcode::extractBarcode(string &seqIn, string &qualIn, const uint32 aStart, string &bSeq, string &bQual)
{//input: sequence seqIn, anchor start aStart
 //output: start position of the barcode
    int32 bStart=0, bEnd=0, aPos=0;
    switch (posAnchor) {
        case 0: //read start
            aPos=0;
            break;
        case 1: //read end
            aPos=seqIn.size()-1;
            break;
        case 2: //adapter start
            aPos=aStart;
            break;
        case 3: //adapter end
            aPos=(int32)aStart+adapterLength-1;
            break;
    };
    
    switch (posType) {
        case 0: //start
            bStart=aPos+pos;
            if (length>0) {
                bEnd=bStart+length-1;
            } else {//extend to the end of the sequence
                bEnd=seqIn.size()-1;
            };
            break;
        case 1: //end
            bEnd=aPos+pos;
            if (length>0) {
                bStart=bEnd+1-length;
            } else {//extend to the start of the sequence
                bStart=0;
            };
            break;
    };
                
    if (bStart<0 || bEnd>(int32)seqIn.size()) //something went wrong
        return false;

    bSeq +=seqIn.substr(bStart,bEnd-bStart+1);
    bQual+=qualIn.substr(bStart,bEnd-bStart+1);
    return true;
};