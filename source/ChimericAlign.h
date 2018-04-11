#ifndef CODE_ChimericAlign
#define CODE_ChimericAlign

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ChimericSegment.h"
#include "Genome.h"

#include <memory>

class ReadAlign;

class ChimericAlign
{//
    public:
        ChimericSegment  seg1, seg2;     //two chimeric segments
//         std::unique_ptr <Transcript> al1, al2; //two chimeric alignments - modified by chimeric switching
        Transcript *al1, *al2;
        uint ex1, ex2;
        
        uint chimJ1, chimJ2, chimRepeat1, chimRepeat2;
        int chimMotif, chimStr, chimScore;
        
        ChimericAlign(ChimericSegment &seg1in, ChimericSegment &seg2in, int chimScoreIn, Genome &genomeIn, ReadAlign *RAin); //allocate
        void chimericJunctionOutput(fstream &outStream, uint chimN);
        void chimericStitching(char *genSeq, char *readSeq);
        bool chimericCheck();
        
        bool stitchingDone;
        
    private:
        Parameters &P;
        ParametersChimeric &pCh;     
        Genome &mapGen;        
        ReadAlign *RA;

};

#endif