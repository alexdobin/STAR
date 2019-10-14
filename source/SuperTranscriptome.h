#ifndef H_SuperTranscriptome
#define H_SuperTranscriptome

#include "IncludeDefine.h"
#include "Parameters.h"

struct sjInfo {
    uint32 start;
    uint32 end;
    uint32 tr;
    uint32 super;
};

class SuperTranscript {//one supertranscript
public:
    uint8 *seqP;//pointer to sequence
    uint32 length;
    vector<array<uint32,3>> sjC;//collapsed junctions
    vector<uint32> sjDonor;//SJ donor coordinates, sorted
};


class SuperTranscriptome {
    private:
    Parameters &P;
public:   
    vector<uint8> seqConcat;//concatenated sequences of supertranscripts, a.k.a. Condensed Genome 
    vector<vector<uint8>> seq;//sequences of supertranscripts
    vector<uint64> trIndex;//superTr's index this tr belongs to
    vector<array<uint64,2>> trStartEnd;//tr start/end in the superTr it belongs to
    vector<sjInfo> sj;//all junctions
    
    vector<SuperTranscript> superTrs;
    uint32 sjNmax, sjDonorNmax;//max number of SJs per superTr, SJ donors
    uint32 N; //number of superTr
    
    
    SuperTranscriptome(Parameters &P) : P(P) {};
    void sjCollapse();
    void load(char *G, vector<uint64> &chrStart, vector<uint64> &chrLength);
};
#endif
