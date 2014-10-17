#ifndef CODE_BAMoutput
#define CODE_BAMoutput

#include "IncludeDefine.h"
#include SAMTOOLS_BGZF_H


class BAMoutput {//
public:
    //sorted output
    BAMoutput (uint64 bamArraySizeIn, uint32 nBinsIn, uint64 genomeSize, int iChunk, string tmpDir);
    void coordOneAlign (char *bamIn, uint bamSize, uint chrStart, uint iRead);
    void coordFlush ();
    //unsorted output
    BAMoutput (uint64 bamArraySizeIn, BGZF *bgzfBAMin);
    void unsortedOneAlign (char *bamIn, uint bamSize, uint bamSize2);
    void unsortedFlush ();
    
    uint32 nBins; //number of bins to split genome into
    uint* binTotalN; //total number of aligns in each bin
    uint* binTotalBytes;//total size of aligns in each bin
private:
    uint64 bamArraySize; //this size will be allocated
    char* bamArray; //large array to store the bam alignments, pre-sorted
    uint64 binSize;//storage size of each bin
    uint64 binGlen;//bin genomic length
    char **binStart; //pointers to starts of the bins
    uint64 *binBytes, binBytes1;//number of bytes currently written to each bin
    ofstream **binStream;//output streams for each bin
    BGZF *bgzfBAM;
};

#endif