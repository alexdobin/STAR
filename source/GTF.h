#ifndef H_GTF
#define H_GTF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "SjdbClass.h"
#include "Genome.h"
#include "SuperTranscriptome.h"

class GTF {
private:
    Genome &genome;
    Parameters &P;
    string dirOut;
    SjdbClass &sjdbLoci;
public:
    bool gtfYes; 

    enum {exT,exS,exE,exG,exL}; //indexes in the exonLoci array
    uint64 exonN;
    vector<array<uint64,exL>> exonLoci;
    vector<uint32> transcriptStrand;
    vector <string> transcriptID, geneID;    
    vector<array<string,2>> geneAttr;
    std::map <string,uint64> transcriptIDnumber, geneIDnumber;

    vector<vector<uint8>> transcriptSeq;//sequences of normal transcripts
    vector<array<uint64,2>> transcriptStartEnd;//normal transcripts start/end in the normal genome

    SuperTranscriptome superTr;
    
    GTF(Genome &genomeIn, Parameters &Pin, string dirOutIn, SjdbClass &sjdbLociIn);
    uint64 transcriptGeneSJ();
    void superTranscript();
};

#endif
