#include "SoloBarcode.h"
#include <array>

bool SoloBarcode::extractBarcode(string &seqIn, string &qualIn, const uint32 adapterStart, string &bSeq, string &bQual)
{//input: sequence seqIn, adapter start adapterStart
 //output: start position of the barcode   
    array<int32,2> pos={0,0};
    for (uint32 ii=0; ii<2; ii++) {
        switch (anchorType[ii]) {//this calculates the position of the anchor base
            case 0: //read start
                pos[ii]=0;
                break;
            case 1: //read end
                pos[ii]=seqIn.size()-1;
                break;
            case 2: //adapter start
                pos[ii]=(int32)adapterStart;
                break;
            case 3: //adapter end
                pos[ii]=(int32)adapterStart+adapterLength-1;
                break;
        };
        pos[ii]+=anchorDist[ii];
    };
                
    bSeq="";
    bQual="";
    if (pos[0]<0 || pos[1]>(int32)seqIn.size() || pos[0]>pos[1]) //something went wrong
        return false;

    bSeq =seqIn.substr(pos[0],pos[1]-pos[0]+1);
    bQual=qualIn.substr(pos[0],pos[1]-pos[0]+1);
    return true;
};
