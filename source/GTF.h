#ifndef H_GTF
#define H_GTF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "SjdbClass.h"
#include "Genome.h"
#include "SuperTranscriptome.h"
#include <thread>

class GTFearly;

class GTF {
private:
    Genome &genome;
    Parameters &P;
    string dirOut;
    SjdbClass &sjdbLoci;

    GTF(Genome &genomeIn, Parameters &Pin, const string &dirOut, SjdbClass &sjdbLoci, GTFearly &early);
    void loadExons(ostream &logOut, ostream &errStream, string *fatalError);

    std::thread writeThread[2];//transcriptGeneSJ's table writes: writeTables, writeSjdbList
    void writeTables(ofstream *exgeOutP, ofstream *geOutP, const string dirOut);
    void writeSjdbList(const string dirOut, uint64 sjdbN1, uint64 sjdbN2);
public:
    bool gtfYes; 

    enum {exT,exS,exE,exG,exL}; //indexes in the exonLoci array
    uint64 exonN;
    vector<array<uint64,exL>> exonLoci;
    vector<uint32> transcriptStrand;
    vector <string> transcriptID, geneID;    
    vector<array<string,2>> geneAttr;

    vector<vector<uint8>> transcriptSeq;//sequences of normal transcripts
    vector<array<uint64,2>> transcriptStartEnd;//normal transcripts start/end in the normal genome

    SuperTranscriptome superTrome;
    
    GTF(Genome &genomeIn, Parameters &Pin, const string &dirOut, SjdbClass &sjdbLoci);
    ~GTF();
    static void loadEarly(Genome &genome, Parameters &P, const string &dirOut);
    uint64 transcriptGeneSJ(const string &dirOut, bool writeLater=false);//with writeLater, the tables it writes are complete after writeJoin()
    void writeJoin();
    void superTranscript();
};

#endif
