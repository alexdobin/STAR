#ifndef H_GTF
#define H_GTF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "SjdbClass.h"
#include "Genome.h"

class GTF {
public:
    bool gtfYes; 

    enum {exT,exS,exE,exG,exL}; //indexes in the exonLoci array
    uint64 exonN;
    vector<array<uint64,exL>> exonLoci;
    vector<uint32> transcriptStrand;
    vector <string> transcriptID, geneID;    
    
    vector<vector<uint8>> transcriptSeq;//sequences of normal transcripts
    vector<vector<uint8>> superTrSeq;//sequences of supertranscripts
    vector<pair<uint64, uint64>> superTrStartEnd;//superTr start end in normal genome
    vector<uint8> superTrGenomeSeq;//condensed genome made of supertranscripts
    vector<pair<uint64, uint64>> transcriptStartEnd;//normal transcripts start/end in the normal genome
    vector<uint64> transcriptSuperTrIndex;//superTr's index this tr belongs to
    vector<pair<uint64, uint64>> transcriptSuperTrStartEnd;//tr start/end in the superTr it belongs to
    
    enum {sjS,sjE,sjT,sjSu,sjL};//start end transcript superTranscript
    vector<array<uint32,sjL>> spliceJunctions; //splice junctions
    
    vector<array<string,2>> geneAttr;
    std::map <string,uint64> transcriptIDnumber, geneIDnumber;
    

    GTF(Genome &genomeIn, Parameters &Pin, string dirOutIn, SjdbClass &sjdbLociIn);
    uint64 transcriptGeneSJ();
    void superTranscript();
private:
    Genome &genome;
    Parameters &P;
    string dirOut;
    SjdbClass &sjdbLoci;
};

#endif
