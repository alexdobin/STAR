#include "BAMbinSortByCoordinate.h"
#include "ErrorWarning.h"
#include "serviceFuns.cpp"
#include "BAMfunctions.h"

void BAMbinSortByCoordinate(uint32 iBin, uint binN, uint binS, uint nThreads, string dirBAMsort, Parameters &P, Genome &mapGen) {

    if (binS==0) return; //nothing to do for empty bins
    //allocate arrays
    char *bamIn=new char[binS+1];
    uint *startPos=new uint[binN*3];

    uint bamInBytes=0;
    //load all aligns
    for (uint it=0; it<nThreads; it++) {
        string bamInFile=dirBAMsort+to_string(it)+"/"+to_string((uint) iBin);
        ifstream bamInStream;
        bamInStream.open(bamInFile.c_str(),std::ios::binary | std::ios::ate);//open at the end to get file size
        uint s1=bamInStream.tellg();
        if (s1>0)
        {
            bamInStream.seekg(std::ios::beg);
            bamInStream.read(bamIn+bamInBytes,s1);//read the whole file
        };
        bamInBytes += bamInStream.gcount();
        bamInStream.close();
        remove(bamInFile.c_str());
    };
    if (bamInBytes!=binS) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: number of bytes expected from the BAM bin does not agree with the actual size on disk: ";
        errOut << binS <<"   "<< bamInBytes <<"   "<< iBin <<"\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, 1, P);
    };

    //extract coordinates

    for (uint ib=0,ia=0;ia<binN;ia++) {
        uint32 *bamIn32=(uint32*) (bamIn+ib);
        startPos[ia*3]  =( ((uint) bamIn32[1]) << 32) | ( (uint)bamIn32[2] );
        startPos[ia*3+2]=ib;
        ib+=bamIn32[0]+sizeof(uint32);//note that size of the BAM record does not include the size record itself
        startPos[ia*3+1]=*( (uint*) (bamIn+ib) ); //read order
        ib+=sizeof(uint);
    };

    //sort
    qsort((void*) startPos, binN, sizeof(uint)*3, funCompareArrays<uint,3>);

    BGZF *bgzfBin;
    bgzfBin=bgzf_open((dirBAMsort+"/b"+to_string((uint) iBin)).c_str(),("w"+to_string((long long) P.outBAMcompression)).c_str());
    outBAMwriteHeader(bgzfBin,P.samHeaderSortedCoord,mapGen.chrNameAll,mapGen.chrLengthAll);
    //send ordered aligns to bgzf one-by-one
    for (uint ia=0;ia<binN;ia++) {
        char* ib=bamIn+startPos[ia*3+2];
        bgzf_write(bgzfBin,ib, *((uint32*) ib)+sizeof(uint32) );
    };

    bgzf_flush(bgzfBin);
    bgzf_close(bgzfBin);
    //release memory
    delete [] bamIn;
    delete [] startPos;
};