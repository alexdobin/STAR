#ifndef H_GTF
#define H_GTF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "SjdbClass.h"
#include "Genome.h"

class GTF{
public:
    bool gtfYes; 

    enum {exT,exS,exE,exG,exL}; //indexes in the exonLoci array
    uint64 exonN;
    vector<array<uint64,exL>> exonLoci;
    vector<uint32> transcriptStrand;
    vector <string> transcriptID, geneID;    
    
    vector<vector<uint8>> sequenceOfNormalTranscripts;
    vector<vector<uint8>> sequenceOfSuperTranscripts;
    vector<pair<uint64, uint64>> superTranscriptIntervals;
    vector<uint8> sequenceOfCondensedGenome;
    vector<pair<uint64, uint64>> transcriptIntervals;
    vector<uint64> normalTranscriptSuperTindex;
    vector<pair<uint64, uint64>> normalTranscriptIntervalsInST;
    
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
