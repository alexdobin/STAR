#include "BAMbinSortUnmapped.h"
#include "ErrorWarning.h"
#include "BAMfunctions.h"

void BAMbinSortUnmapped(uint32 iBin, uint nThreads, string dirBAMsort, Parameters &P, Genome &mapGen) {

    BGZF *bgzfBin;
    bgzfBin=bgzf_open((dirBAMsort+"/b"+to_string((uint) iBin)).c_str(),("w"+to_string((long long) P.outBAMcompression)).c_str());
    outBAMwriteHeader(bgzfBin,P.samHeaderSortedCoord,mapGen.chrNameAll,mapGen.chrLengthAll);


    vector<string> bamInFile;
    std::map <uint,uint> startPos;

    for (uint it=0; it<nThreads; it++) {//files from all threads, and BySJout
        bamInFile.push_back(dirBAMsort+to_string(it)+"/"+to_string((uint) iBin));
        bamInFile.push_back(dirBAMsort+to_string(it)+"/"+to_string((uint) iBin)+".BySJout");
    };
    vector<uint32> bamSize(bamInFile.size(),0);//record sizes

    //allocate arrays
    char **bamIn=new char* [bamInFile.size()];
    ifstream *bamInStream = new ifstream [bamInFile.size()];

    for (uint it=0; it<bamInFile.size(); it++) {//initialize
        bamIn[it] = new char [BAMoutput_oneAlignMaxBytes];

        bamInStream[it].open(bamInFile.at(it).c_str());//opean all files

        bamInStream[it].read(bamIn[it],sizeof(int32));//read BAM record size
        if (bamInStream[it].good()) {
            bamSize[it]=((*(uint32*)bamIn[it])+sizeof(int32));//true record size +=4 (4 bytes for uint-iRead)
            bamInStream[it].read(bamIn[it]+sizeof(int32),bamSize.at(it)-sizeof(int32)+sizeof(uint));//read the rest of the record, including last uint = iRead
            startPos[*(uint*)(bamIn[it]+bamSize.at(it))]=it;//startPos[iRead]=it : record the order of the files to output
        } else {//nothing to do here, file is empty, do not record it
        };
    };

    //send ordered aligns to bgzf one-by-one
    while (startPos.size()>0) {
        uint it=startPos.begin()->second;
        uint startNext=startPos.size()>1 ? (++startPos.begin())->first : (uint) -1;

        while (true) {
            bgzf_write(bgzfBin, bamIn[it], bamSize.at(it));
            bamInStream[it].read(bamIn[it],sizeof(int32));//read record size
            if (bamInStream[it].good()) {
                 bamSize[it]=((*(uint32*)bamIn[it])+sizeof(int32));
                 bamInStream[it].read(bamIn[it]+sizeof(int32),bamSize.at(it)-sizeof(int32)+sizeof(uint));//read the rest of the record, including la$
                 uint iRead=*(uint*)(bamIn[it]+bamSize.at(it));
                 if (iRead>startNext) {//this read from this chunk is > than a read from another chunk
                     startPos[iRead]=it;
                     break;
                 };
            } else {//nothing to do here, reached the end of the file
                 break;
            };
        };
        startPos.erase(startPos.begin());
    };

    bgzf_flush(bgzfBin);
    bgzf_close(bgzfBin);


    for (uint it=0; it<bamInFile.size(); it++) {//destroy at the end
        bamInStream[it].close();
        remove(bamInFile.at(it).c_str());
        delete [] bamIn[it];
    };
    delete [] bamIn;
    delete [] bamInStream;
};
