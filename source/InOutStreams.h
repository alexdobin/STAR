#ifndef INOUTSTREAMS_DEF
#define INOUTSTREAMS_DEF

#include "IncludeDefine.h"
#include SAMTOOLS_BGZF_H

class InOutStreams {
    public:
    ostream *logStdOut, *outSAM;
    ofstream logStdOutFile, outSAMfile;
    BGZF *outBAMfileUnsorted, *outBAMfileCoord, *outQuantBAMfile;

    ofstream outChimSAM, outChimJunction, logMain, logProgress, logFinal, outUnmappedReadsStream[MAX_N_MATES];
    ifstream readIn[MAX_N_MATES];

    //compilation-optional streams
    ofstream outLocalChains;

    InOutStreams();
    ~InOutStreams();
};

#endif
