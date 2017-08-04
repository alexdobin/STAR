#ifndef CODE_ChimericTranscript
#define CODE_ChimericTranscript

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"

class ChimericTranscript
{//
    public:
        Transcript  **chTrs; //all chimeric transcripts
        uint        nCh;     //number of recorded (best) chimeric transcripts
        uint        nChSize;  //size of the chTrs array, will be increased if nCh > nChSize
        
        ChimericTranscript(Parameters &Pin); //allocate
    private:
};

#endif