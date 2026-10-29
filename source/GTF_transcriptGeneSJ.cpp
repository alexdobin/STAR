#include "GTF.h"
#include "serviceFuns.cpp"
#include "streamFuns.h"
#include <algorithm>
#include <cstdlib>
#include <tuple>

struct GTFexgeRow {//exon start/end/strand/gene/transcript, sorted by all five
    uint64 exStart, exEnd, exStrand, geID, trID;
    bool operator<(const GTFexgeRow &b) const {
        return std::tie(exStart,exEnd,exStrand,geID,trID) < std::tie(b.exStart,b.exEnd,b.exStrand,b.geID,b.trID);
    };
};
struct GTFextrRow {//transcript start/end/ID, exon start/end and gene, sorted by the first five
    uint64 trStart, trEnd, trID, exStart, exEnd, geID;
    bool operator<(const GTFextrRow &b) const {
        return std::tie(trStart,trEnd,trID,exStart,exEnd) < std::tie(b.trStart,b.trEnd,b.trID,b.exStart,b.exEnd);
    };
};
struct GTFsjRow {//junction start/end/strand/gene, sorted by start/end
    uint64 start, end, str, gene;
    bool operator<(const GTFsjRow &b) const {
        return start<b.start || (start==b.start && end<b.end);
    };
};

static GTF *gtfWritePending=NULL;//its table writes are finished before exit()
static void gtfWriteJoinAtExit()
{
    if (gtfWritePending!=NULL)
        gtfWritePending->writeJoin();
};

uint64 GTF::transcriptGeneSJ(const string &dirOut, bool writeLater)
{//sort exonLoci by transcript ID and exon coordinates
 //fills sjdbLoci from GTF junctions
 //the tables are written on threads: complete on return or, with writeLater, after writeJoin()
    if (!gtfYes)
        return 0;

    exonN=exonLoci.size();
    std::stable_sort(exonLoci.begin(), exonLoci.end(), [](const array<uint64,exL> &a, const array<uint64,exL> &b) {//the order of glibc's qsort (a stable merge sort) with funCompareUint2
        return a[exT]<b[exT] || (a[exT]==b[exT] && a[exS]<b[exS]);
    });

    writeJoin();
    ofstream *exgeOut=&ofstrOpen(dirOut+"/exonGeTrInfo.tab",ERROR_OUT,P);

    static bool atExitSet=false;
    if (!atExitSet) {
        atexit(gtfWriteJoinAtExit);
        atExitSet=true;
    };

    ofstream *geOut=&ofstrOpen(dirOut+"/geneInfo.tab",ERROR_OUT,P);

    gtfWritePending=this;
    try {
        writeThread[0]=std::thread(&GTF::writeTables, this, exgeOut, geOut, dirOut);
    } catch (const std::system_error &) {//no thread: write here
        writeTables(exgeOut, geOut, dirOut);
    };

    //make junctions
    vector<GTFsjRow> sjLoci;
    sjLoci.reserve(exonN);
    uint64 trIDn=exonLoci[0][exT];
    for (uint64 iex=1; iex<exonN; iex++) {
        if (trIDn==exonLoci[iex][exT]) {
                uint64 chr1=genome.chrBin[exonLoci[iex][exS] >> genome.pGe.gChrBinNbits];
            if ( exonLoci[iex][exS]<=exonLoci[iex-1][exE]+1 ) {//touching - nothing to do
            } else if ( exonLoci[iex][exS]<=exonLoci[iex-1][exE] ) {//overlapping
                P.inOut->logMain << "WARNING: while processing pGe.sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<": overlapping exons:\n";
                P.inOut->logMain << genome.chrName[chr1] <<"\t"<< exonLoci[iex-1][exS]+1-genome.chrStart[chr1] << "\t"<< exonLoci[iex-1][exE]+1-genome.chrStart[chr1]  <<"\n";
                P.inOut->logMain << genome.chrName[chr1] <<"\t"<< exonLoci[iex][exS]  +1-genome.chrStart[chr1] << "\t"<< exonLoci[iex][exE]  +1-genome.chrStart[chr1]  <<"\n";
            } else {
                sjLoci.push_back({exonLoci[iex-1][exE]+1, exonLoci[iex][exS]-1, (uint64) transcriptStrand[trIDn], exonLoci[iex][exG]+1});//genes are numbered from 1
            };
        } else {
            trIDn=exonLoci[iex][exT];
        };
    };

    std::stable_sort(sjLoci.begin(), sjLoci.end());//the order of glibc's qsort with funCompareUint2

    char strandChar[3]={'.','+','-'};
    uint64 sjdbN1=sjdbLoci.chr.size();
    sjdbLoci.gene.resize(sjdbN1); //need to resize in case sjdbLoci was loaded from files without gene attribute. TODO make sure gene is always present
    for (uint64 ii=0;ii<sjLoci.size();ii++) {
        if ( ii==0 || (sjLoci[ii].start!=sjLoci[ii-1].start)
                   || (sjLoci[ii].end!=sjLoci[ii-1].end)
                   || (sjLoci[ii].str!=sjLoci[ii-1].str) ) {
            uint64 chr1=genome.chrBin[sjLoci[ii].start >> genome.pGe.gChrBinNbits];
            sjdbLoci.chr.push_back(genome.chrName[chr1]);
            sjdbLoci.start.push_back(sjLoci[ii].start+1-genome.chrStart[chr1]);
            sjdbLoci.end.push_back(sjLoci[ii].end+1-genome.chrStart[chr1]);
            sjdbLoci.str.push_back(strandChar[sjLoci[ii].str]);
            sjdbLoci.gene.push_back({sjLoci[ii].gene});
        } else {
            sjdbLoci.gene.back().insert(sjLoci[ii].gene);
        };
    };

    try {//until writeJoin, sjdbLoci's chr, start, end, str and gene are only read
        writeThread[1]=std::thread(&GTF::writeSjdbList, this, dirOut, sjdbN1, (uint64) sjdbLoci.chr.size());
    } catch (const std::system_error &) {
        writeSjdbList(dirOut, sjdbN1, sjdbLoci.chr.size());
    };

    sjdbLoci.priority.resize(sjdbLoci.chr.size(),20);

    P.inOut->logMain << "Processing pGe.sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<", found:\n";
    P.inOut->logMain << "\t\t"  << transcriptID.size() <<" transcripts\n" << "\t\t"  << exonN << " exons (non-collapsed)\n" << "\t\t"  << sjdbLoci.chr.size()-sjdbN1 << " collapsed junctions\n";
    P.inOut->logMain << "Total junctions: " <<sjdbLoci.chr.size()<<"\n";
    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ..... finished GTF processing\n\n" <<flush;

    if (!writeLater)
        writeJoin();

    return sjdbLoci.chr.size()-sjdbN1;
};

void GTF::writeJoin()
{
    for (auto &t : writeThread) {
        if (t.joinable())
            t.join();
    };
    if (gtfWritePending==this)
        gtfWritePending=NULL;
};

GTF::~GTF()
{
    writeJoin();
};

void GTF::writeTables(ofstream *exgeOutP, ofstream *geOutP, const string dirOut)
{//exonGeTrInfo.tab, geneInfo.tab, transcriptInfo.tab and exonInfo.tab; reads only exonN, the sorted exonLoci, transcriptStrand, transcriptID, geneID and geneAttr
    {//exon-gene data structures: exon start/end/strand/gene/transcript
        //re-sort exons by exons loci
        vector<GTFexgeRow> exgeLoci(exonN);
        for (uint64 iex=0; iex<exonN; iex++) {
            exgeLoci[iex]={exonLoci[iex][exS], exonLoci[iex][exE], transcriptStrand[exonLoci[iex][exT]], exonLoci[iex][exG], exonLoci[iex][exT]};
        };

        std::stable_sort(exgeLoci.begin(), exgeLoci.end());//the order of glibc's qsort with funCompareArrays<uint64,5>

        ofstream &exgeOut=*exgeOutP;
        exgeOut<<exonN<<"\n";
        for (uint64 iex=0; iex<exonN; iex++) {
             exgeOut<<exgeLoci[iex].exStart <<"\t"<<  exgeLoci[iex].exEnd <<"\t"<< exgeLoci[iex].exStrand \
              <<"\t"<< exgeLoci[iex].geID <<"\t"<< exgeLoci[iex].trID <<"\n"; //the last value, transript-number, is worng here since tranascripts are re-sorted later
        };
        exgeOut.close();
        delete exgeOutP;

        ofstream &geOut=*geOutP;
        geOut << geneID.size() << "\n";
        for (uint64 ig=0; ig<geneID.size(); ig++) {//just geneID for now
            geOut << geneID[ig] <<"\t"<< geneAttr[ig][0] <<"\t"<< geneAttr[ig][1] <<"\n";
        };
        geOut.close();
        delete geOutP;
    };

    {//exon-transcript data structures
        //re-sort transcripts by transcript start/end
        vector<GTFextrRow> extrLoci(exonN); //this also contains transcripts start and end

        uint64 trex1=0;
        for (uint64 iex=0; iex<=exonN; iex++) {
            if (iex==exonN || exonLoci[iex][exT] != exonLoci[trex1][exT]) {
                for (uint64 iex1=trex1; iex1<iex; iex1++) {//go back and fill the trend
                    extrLoci[iex1].trEnd=exonLoci[iex-1][exE];
                };
                if (iex==exonN) break;
                trex1=iex;
            };
            extrLoci[iex].trStart=exonLoci[trex1][exS];
            extrLoci[iex].trID=exonLoci[iex][exT];
            extrLoci[iex].exStart=exonLoci[iex][exS];
            extrLoci[iex].exEnd=exonLoci[iex][exE];
            extrLoci[iex].geID=exonLoci[iex][exG];
        };

        std::stable_sort(extrLoci.begin(), extrLoci.end());//the order of glibc's qsort with funCompareArrays<uint64,5>

        ofstream trOut ((dirOut+"/transcriptInfo.tab").c_str());
        trOut<<transcriptID.size() << "\n";
        ofstream exOut ((dirOut+"/exonInfo.tab").c_str());
        exOut<<exonN<<"\n";

        uint64 trid=extrLoci[0].trID;
        uint64 trex=0;
        uint64 trstart=extrLoci[0].trStart;
        uint64 trend=extrLoci[0].trEnd;
        uint64 exlen=0;
        for (uint64 iex=0;iex<=exonN; iex++) {
            if (iex==exonN || extrLoci[iex].trID != trid) {//start of the new transcript
                //write out previous transcript
                trOut << transcriptID.at(trid) <<"\t"<< extrLoci[iex-1].trStart<<"\t"<< extrLoci[iex-1].trEnd \
                       <<"\t"<< trend << "\t"<< (uint64) transcriptStrand[trid]  <<"\t"<< iex-trex <<"\t"<<trex<<"\t"<<extrLoci[iex-1].geID<<"\n";
                if (iex==exonN) break;
                trid=extrLoci[iex].trID;
                trstart=extrLoci[iex].trStart;
                trex=iex;
                trend=max(trend,extrLoci[iex-1].trEnd);
                exlen=0;
            };
            exOut << extrLoci[iex].exStart-trstart <<"\t"<< extrLoci[iex].exEnd-trstart <<"\t"<< exlen <<"\n";
            exlen+=extrLoci[iex].exEnd-extrLoci[iex].exStart+1;
        };
        trOut.close();
        exOut.close();
    };
};

void GTF::writeSjdbList(const string dirOut, uint64 sjdbN1, uint64 sjdbN2)
{//sjdbList.fromGTF.out.tab: the junctions sjdbN1..sjdbN2-1 transcriptGeneSJ added to sjdbLoci
    ofstream sjdbList ((dirOut+"/sjdbList.fromGTF.out.tab").c_str());
    for (uint64 ii=sjdbN1;ii<sjdbN2; ii++) {
        sjdbList << sjdbLoci.chr.at(ii)<<"\t"<< sjdbLoci.start.at(ii) << "\t"<< sjdbLoci.end.at(ii)  <<"\t"<< sjdbLoci.str.at(ii);

        auto gg=sjdbLoci.gene[ii].cbegin();//iterator for genes
        sjdbList <<"\t"<< *gg;
        ++gg;
        for (; gg!=sjdbLoci.gene[ii].cend(); gg++)
            sjdbList <<","<< *gg;
        sjdbList<<"\n";
    };
    sjdbList.close();
};
