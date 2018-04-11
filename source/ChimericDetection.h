#ifndef CODE_ChimericDetection
#define CODE_ChimericDetection

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ChimericAlign.h"
#include "Genome.h"

class ReadAlign;

class ChimericDetection {
    private:
        Parameters &P;
        ReadAlign *RA;
        Transcript ***trAll;
        uint nW, *nWinTr;
        char** Read1;
        Genome &outGen;        
        uint *readLength;
        
    public:
        uint chimN;
        vector <ChimericAlign> chimAligns;
        bool chimRecord;
        int chimScoreBest;
                
        ChimericDetection(Parameters &Pin, Transcript ***trAll, uint *nWinTr, char** Read1in, Genome &genomeIn, fstream *ostreamChimJunctionIn, ReadAlign *RA);
        bool chimericDetectionMult(uint nWin, uint *readLengthIn);
        fstream *ostreamChimJunction;
};

#endif