#ifndef CODE_ChimericDetection
#define CODE_ChimericDetection

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ChimericAlign.h"
#include "Genome.h"

class ChimericDetection {
    public:
        uint chimN;
        vector <ChimericAlign> chimAligns;
        bool chimRecord;
        int chimScoreBest;
        
        vector <Transcript> alignsModified; //separate storage for modified aligns
        
        ChimericDetection(Parameters &Pin, Transcript ***trAll, uint *nWinTr, char** Read1in, Genome &mapGenIn, fstream &ostreamChimJunctionIn);
        bool chimericDetectionMult(uint nWin, uint *readLengthIn);
        
    private:
        Parameters &P;
        Transcript ***trAll;
        uint nW, *nWinTr;
        char** Read1;
        uint *readLength;
        Genome &outGen;
        fstream &ostreamChimJunction;
};

#endif