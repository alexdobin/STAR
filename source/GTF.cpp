#include "GTF.h"

#include "ErrorWarning.h"
#include "streamFuns.h"
#include "TimeFunctions.h"

#include <thread>
#include <system_error>
#include <unordered_map>
#include <cstring>

namespace {
class GtfAttrParser {//the value of each attribute name in an exon line's last field
//the original search: ';', '=', '\t' and '"' become spaces, and a name's value is the space-delimited word after the first " name "
public:
    GtfAttrParser(const vector<vector<string>> &slotNamesIn);
    void parse(const string &attrField, vector<string> &exAttr);//exAttr[ii]: the value of the last of slot ii's names that has one, else ""
private:
    vector<vector<string>> slotNames;
    bool byWord;//every name is non-empty and has no space, so " name " can only match one whole word
    vector<string> names;//distinct names
    vector<vector<uint32>> slotIndex;//slot -> indexes in names, in search order
    vector<uint8> state;//per name: 0 not seen, 1 seen and waiting for its value, 2 has a value
    vector<size_t> valStart, valLen;
    bool sep[256];
};

GtfAttrParser::GtfAttrParser(const vector<vector<string>> &slotNamesIn) : slotNames(slotNamesIn)
{
    byWord=true;
    slotIndex.resize(slotNames.size());
    for (uint32 ii=0; ii<slotNames.size(); ii++) {
        for (auto &name : slotNames[ii]) {
            if (name.empty() || name.find(' ')!=string::npos)
                byWord=false;
            uint32 k=0;
            while (k<names.size() && names[k]!=name)
                k++;
            if (k==names.size())
                names.push_back(name);
            slotIndex[ii].push_back(k);
        };
    };
    state.resize(names.size());
    valStart.resize(names.size());
    valLen.resize(names.size());
    for (int c=0; c<256; c++)
        sep[c]=false;
    sep[(uint8)' ']=sep[(uint8)';']=sep[(uint8)'=']=sep[(uint8)'\t']=sep[(uint8)'"']=true;
};

void GtfAttrParser::parse(const string &attrField, vector<string> &exAttr)
{
    for (auto &a : exAttr)
        a.clear();

    if (!byWord) {//a name may span words or match between two spaces: search as the original did
        string oneLine1(attrField);
        replace(oneLine1.begin(),oneLine1.end(),';',' ');
        replace(oneLine1.begin(),oneLine1.end(),'=',' ');
        replace(oneLine1.begin(),oneLine1.end(),'\t',' ');
        replace(oneLine1.begin(),oneLine1.end(),'\"',' ');
        for (uint32 ii=0; ii<slotNames.size(); ii++) {
            for (auto &attr1 : slotNames[ii]) {
                size_t pos1=oneLine1.find(" " + attr1 + " ");
                if (pos1!=string::npos)
                    pos1=oneLine1.find_first_not_of(" ", pos1+attr1.size()+1);
                if (pos1!=string::npos)
                    exAttr[ii]=oneLine1.substr(pos1, oneLine1.find_first_of(" ",pos1)-pos1);
            };
        };
        return;
    };

    //" name " matches a word equal to name with a separator before it (one after it only matters when a value follows); its value is the next word
    const char *s=attrField.data();
    size_t n=attrField.size();
    uint32 nOpen=names.size();//names still not seen or waiting for a value
    std::fill(state.begin(), state.end(), 0);
    size_t i=0;
    while (nOpen>0) {
        while (i<n && sep[(uint8)s[i]])
            i++;
        if (i==n)
            break;
        size_t b=i;
        while (i<n && !sep[(uint8)s[i]])
            i++;
        size_t len=i-b;
        for (uint32 k=0; k<names.size(); k++) {
            if (state[k]==1) {
                state[k]=2;
                valStart[k]=b;
                valLen[k]=len;
                nOpen--;
            } else if (state[k]==0 && b>0 && names[k].size()==len && memcmp(names[k].data(), s+b, len)==0) {
                state[k]=1;
            };
        };
    };
    for (uint32 ii=0; ii<slotIndex.size(); ii++) {
        for (uint32 k : slotIndex[ii]) {
            if (state[k]==2)
                exAttr[ii].assign(s+valStart[k], valLen[k]);
        };
    };
};
}

//whitespace as istream >> sees it in the classic locale
static inline bool gtfIsSpace(char c)
{
    return c==' ' || (c>='\t' && c<='\r');
};

//next token at or after pos, as >> string reads it: [b,e), pos left at e; false when only whitespace is left
static inline bool gtfNextToken(const string &s, size_t &pos, size_t &b, size_t &e)
{
    size_t n=s.size();
    while (pos<n && gtfIsSpace(s[pos]))
        ++pos;
    if (pos==n)
        return false;
    b=pos;
    while (pos<n && !gtfIsSpace(s[pos]))
        ++pos;
    e=pos;
    return true;
};

//next token if it is 1 to 19 decimal digits: >> uint64 reads it whole and without overflow
static inline bool gtfNextUint(const string &s, size_t &pos, uint64 &v)
{
    size_t b, e;
    if (!gtfNextToken(s,pos,b,e) || e-b>19)
        return false;
    uint64 v1=0;
    for (size_t ii=b; ii<e; ii++) {
        if (s[ii]<'0' || s[ii]>'9')
            return false;
        v1=v1*10+(uint64)(s[ii]-'0');
    };
    v=v1;
    return true;
};

//next non-whitespace character, as >> char reads it
static inline bool gtfNextChar(const string &s, size_t &pos, char &c)
{
    size_t n=s.size();
    while (pos<n && gtfIsSpace(s[pos]))
        ++pos;
    if (pos==n)
        return false;
    c=s[pos++];
    return true;
};

class GTFearly {//GTF exon lines parsed on a thread while genomeLoad reads the genome arrays
public:
    SjdbClass sjdbLociUnused;//the parse does not touch sjdbLoci
    ostringstream logOut, stdOut, errStream;//its Log.out, stdout and stderr lines, written by the GTF that takes the parse
    string fatalError;//its fatal input error, raised there
    GTF *gtf;
    std::thread thread;
    GTFearly() : gtf(NULL) {};
    ~GTFearly() {
        if (thread.joinable())
            thread.join();
        delete gtf;
    };
};

void GTF::loadEarly(Genome &genome, Parameters &P, const string &dirOut)
{//start parsing the exon lines on a thread; reads only chrName, chrStart, chrLength and chrNameIndex of the genome, which must not change until the parse is taken
    std::shared_ptr<GTFearly> early=std::make_shared<GTFearly>();
    GTFearly *e=early.get();

    time_t rawTime;//localtime is not thread-safe: stamp the line here
    time(&rawTime);
    e->logOut << timeMonthDayTime(rawTime) <<" ..... processing annotations GTF\n";
    e->stdOut << timeMonthDayTime(rawTime) <<" ..... processing annotations GTF\n";

    try {
        e->thread=std::thread([e, &genome, &P, dirOut]() {e->gtf=new GTF(genome, P, dirOut, e->sjdbLociUnused, *e);});
    } catch (const std::system_error &) {//no thread: the GTF constructor parses as usual
        return;
    };
    genome.gtfEarly=early;
};

GTF::GTF(Genome &genome, Parameters &P, const string &dirOut, SjdbClass &sjdbLoci, GTFearly &early)
                 : genome(genome), P(P), dirOut(dirOut), sjdbLoci(sjdbLoci), superTrome(P)
{//runs on the thread started by loadEarly
    gtfYes=true;
    loadExons(early.logOut, early.errStream, &early.fatalError);
};

GTF::GTF(Genome &genome, Parameters &P, const string &dirOut, SjdbClass &sjdbLoci) 
                 : genome(genome), P(P), dirOut(dirOut), sjdbLoci(sjdbLoci), superTrome(P)
 {//initialize; load gtf file; returns number of added junctions
                     
    if (genome.sjdbOverhang==0 || genome.pGe.sjdbGTFfile=="-") {//no GTF
        gtfYes=false;
        return;
    };
    gtfYes=true;
    mkdir(dirOut.c_str(), P.runDirPerm);

    if (genome.gtfEarly) {//genomeLoad started the parse on a thread: wait for it and take its result
        std::shared_ptr<GTFearly> early=genome.gtfEarly;
        genome.gtfEarly.reset();
        early->thread.join();
        P.inOut->logMain << early->logOut.str() <<flush;
        *P.inOut->logStdOut << early->stdOut.str() <<flush;
        std::cerr << early->errStream.str() <<flush;
        if (early->fatalError!="")
            exitWithError(early->fatalError,std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        exonN=early->gtf->exonN;
        exonLoci=std::move(early->gtf->exonLoci);
        transcriptStrand=std::move(early->gtf->transcriptStrand);
        transcriptID=std::move(early->gtf->transcriptID);
        geneID=std::move(early->gtf->geneID);
        geneAttr=std::move(early->gtf->geneAttr);
        return;
    };
    
    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ..... processing annotations GTF\n" <<flush;
    *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ..... processing annotations GTF\n" <<flush;

    loadExons(P.inOut->logMain, std::cerr, NULL);
};

void GTF::loadExons(ostream &logOut, ostream &errStream, string *fatalError)
{//parse the exon lines; with fatalError, a fatal input error is stored there instead of exiting
    auto fatal=[&](const string &message) {
        if (fatalError==NULL)
            exitWithError(message,std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        *fatalError=message;
    };

    std::unordered_map <string,uint64> transcriptIDnumber, geneIDnumber;//numbered by first appearance
    string lastTrID, lastGeneID;//exon lines come grouped by transcript: the previous line's IDs and numbers skip the lookups
    uint64 lastTrN=0, lastGeneN=0;
    bool haveLastTr=false, haveLastGene=false;
    string lastChr;
    uint64 lastChrI=0;
    bool haveLastChr=false;

    GtfAttrParser attrParser({ {genome.pGe.sjdbGTFtagExonParentTranscript}, {genome.pGe.sjdbGTFtagExonParentGene}, genome.pGe.sjdbGTFtagExonParentGeneName, genome.pGe.sjdbGTFtagExonParentGeneType }); //trID, gID, gName, gBiotype
    vector<string> exAttr(4); //trID, gID, gName, gBiotype

    ifstream sjdbStreamIn ( genome.pGe.sjdbGTFfile.c_str() );
    if (sjdbStreamIn.fail()) {
        ostringstream errOut;
        errOut << "FATAL error, could not open file pGe.sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<"\n";
        fatal(errOut.str());
        return;
    };

    if (genome.chrNameIndex.size()==0) {
        for (uint64 ii=0;ii<genome.nChrReal;ii++) {
            genome.chrNameIndex[genome.chrName[ii]]=ii;
        };
    };

    //one pass; the fields before the attributes are split at the whitespace istream >> skips, and a line whose
    //coordinates are not plain digits, or which ends before the strand, is read with istringstream
    exonN=0;
    uint64 exonLinesN=0; //lines with the exon feature, whether or not they are used
    exonLoci.clear();
    string oneLine, chr1, ddd2, oneLine1;
    const string &featureExon=genome.pGe.sjdbGTFfeatureExon;
    while (sjdbStreamIn.good()) {

        getline(sjdbStreamIn,oneLine);
        size_t pos=0, chrB=0, chrE=0, ftB=0, ftE=0, b, e;
        if (gtfNextToken(oneLine,pos,chrB,chrE) && gtfNextToken(oneLine,pos,b,e))
            gtfNextToken(oneLine,pos,ftB,ftE);
        if ((chrE==chrB || oneLine[chrB]!='#') && oneLine.compare(ftB,ftE-ftB,featureExon)==0) {//exonic line, process
            ++exonLinesN;
            chr1.assign(oneLine,chrB,chrE-chrB);

            if (genome.pGe.sjdbGTFchrPrefix!="-") chr1=genome.pGe.sjdbGTFchrPrefix + chr1;

            if (!haveLastChr || chr1!=lastChr) {
                auto chrIt=genome.chrNameIndex.find(chr1);
                if (chrIt==genome.chrNameIndex.end()) {//chr not in Genome
                    logOut << "WARNING: while processing sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<": chromosome '"<<chr1<<"' not found in Genome fasta files for line:\n";
                    logOut << oneLine <<"\n"<<flush;
                    continue; //do not process exons/transcripts on missing chromosomes
                };
                lastChr=chr1;
                lastChrI=chrIt->second;
                haveLastChr=true;
            };
            uint64 chrI=lastChrI;

            uint64 ex1=0,ex2=0;
            char str1=0;
            if (gtfNextUint(oneLine,pos,ex1) && gtfNextUint(oneLine,pos,ex2) && gtfNextToken(oneLine,pos,b,e) && gtfNextChar(oneLine,pos,str1)) {
                if (gtfNextToken(oneLine,pos,b,e)) {
                    oneLine1.assign(oneLine,e,string::npos);//the last field, from the separator after the frame
                } else {
                    oneLine1.clear();//no frame: >> fails and getline reads nothing
                };
            } else {
                istringstream oneLineStream (oneLine);
                oneLineStream >> ddd2 >> ddd2 >> ddd2 >> ex1 >> ex2 >> ddd2 >> str1 >> ddd2; //read all fields except the last
                oneLine1.clear();
                getline(oneLineStream, oneLine1);//get the last field
            };
            if ( ex2 > genome.chrLength[chrI] ) {
            	warningMessage("while processing sjdbGTFfile=" + genome.pGe.sjdbGTFfile + ", line:\n" + oneLine + "\n exon end = " + to_string(ex2) +  \
            			       " is larger than the chromosome " + chr1 + " length = " + to_string(genome.chrLength[chrI] ) + " , will skip this exon\n", \
							   errStream, logOut, P);
            	continue;
            };
            
            attrParser.parse(oneLine1, exAttr);

            if (exAttr[0]=="") {//no transcript ID
                logOut << "WARNING: while processing pGe.sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<": no transcript_id for line:\n";
                logOut << oneLine <<"\n"<<flush;
                exAttr[0]="tr_" + chr1 +"_"+ to_string(ex1) +"_"+ to_string(ex2) +"_"+ to_string(exonN); //unique name for the transcript
            };
            
            if (exAttr[1]=="") {//no gene ID
                logOut << "WARNING: while processing pGe.sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<": no gene_id for line:\n";
                logOut << oneLine <<"\n"<<flush;
                exAttr[1]="MissingGeneID";
            };
            
            if (exAttr[2]=="") {//no gene name
                exAttr[2]=exAttr[1];
            };
            
            if (exAttr[3]=="") {//no gene name
                exAttr[3]="MissingGeneType";
            };
            
            if (!haveLastTr || exAttr[0]!=lastTrID) {
                auto trIt=transcriptIDnumber.insert(std::pair <string,uint64> (exAttr[0],(uint64) transcriptIDnumber.size())).first;//insert new element if necessary with a new numeric value
                if (transcriptID.size() < transcriptIDnumber.size()) {//new transcript
                    transcriptID.push_back(exAttr[0]);
                    if (str1=='+') {
                       transcriptStrand.push_back(1);
                    } else if (str1=='-') {
                       transcriptStrand.push_back(2);
                    } else {
                       transcriptStrand.push_back(0);
                    };
                };
                lastTrID=exAttr[0];
                lastTrN=trIt->second;
                haveLastTr=true;
            };

            if (!haveLastGene || exAttr[1]!=lastGeneID) {
                auto gIt=geneIDnumber.insert(std::pair <string,uint64> (exAttr[1],(uint64) geneIDnumber.size())).first;//insert new element if necessary with a new numeric value
                if (geneID.size() < geneIDnumber.size()) {//new gene is added
                    geneID.push_back(exAttr[1]);
                    geneAttr.push_back({exAttr[2],exAttr[3]});
                };
                lastGeneID=exAttr[1];
                lastGeneN=gIt->second;
                haveLastGene=true;
            };

            array<uint64,exL> exon1;
            exon1[exT]=lastTrN;
            exon1[exS]=ex1+genome.chrStart[chrI]-1;
            exon1[exE]=ex2+genome.chrStart[chrI]-1;
            exon1[exG]=lastGeneN;
            exonLoci.push_back(exon1);
            ++exonN;
        };//if (chr1.substr(0,1)!="#" && featureType=="exon")
    };//

    if (exonLinesN==0)         {
        ostringstream errOut;
        errOut << "Fatal INPUT FILE error, no ""exon"" lines in the GTF file: " << genome.pGe.sjdbGTFfile <<"\n";
        errOut << "Solution: check the formatting of the GTF file, it must contain some lines with ""exon"" in the 3rd column.\n";
        errOut << "          Make sure the GTF file is unzipped.\n";
        errOut << "          If exons are marked with a different word, use --sjdbGTFfeatureExon .\n";
        fatal(errOut.str());
        return;
    };

    if (exonN==0) {
        ostringstream errOut;
        errOut << "Fatal INPUT FILE error, no valid ""exon"" lines in the GTF file: " << genome.pGe.sjdbGTFfile <<"\n";
        errOut << "Solution: check the formatting of the GTF file. One likely cause is the difference in chromosome naming between GTF and FASTA file.\n";
        fatal(errOut.str());
        return;
    };
    
    return;
};
