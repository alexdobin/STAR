#ifndef H_SuperTranscript
#define H_SuperTranscript

#include "IncludeDefine.h"
#include "Parameters.h"

struct sjInfo {
    uint32 start;
    uint32 end;
    uint32 tr;
    uint32 super;
};

class SuperTranscript {
    private:
    Parameters &P;
public:
    vector<uint8> seqConcat;//concatenated sequences of supertranscripts, a.k.a. Condensed Genome 
    vector<vector<uint8>> seq;//sequences of supertranscripts
    vector<pair<uint64, uint64>> startEndInFullGenome;//superTr start end in normal genome
    vector<uint64> trIndex;//superTr's index this tr belongs to
    vector<pair<uint64, uint64>> trStartEnd;//tr start/end in the superTr it belongs to
    
    vector<uint8*> seqp;//pointers to sequence for each superTr
    vector<uint64> &length;//superTr lengths == Genome.chrLength
    vector<sjInfo> sj; //all splice junctions
    vector<vector<array<uint32,2>>> sjC; //collapsed splice junctions
    
    SuperTranscript(Parameters &P, vector<uint64> &chrLength) : P(P), length(chrLength) {};
    void sjCollapse();
    void load(char *G, vector<uint64> &chrStart);
};
#endif
