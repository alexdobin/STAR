#include "Transcriptome.h"
#include "streamFuns.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"
#include "serviceFuns.cpp"

Transcriptome::Transcriptome (Parameters &Pin, TranscriptomeEarly *earlyIn) : P(Pin), early(earlyIn) {

    if (!P.quant.yes)
        return;

    if (!P.pGe.transform.outQuant) {//standard
        trInfoDir = P.pGe.sjdbGTFfile=="-" ? P.pGe.gDir : P.sjdbInsert.outDir; //if GTF file is given at the mapping stage, it's always used for transcript info
    } else {//transformed genome
        trInfoDir = P.pGeOut.gDir;
    };

    ifstream &geStream = trOpen(trInfoDir+"/geneInfo.tab", ERROR_OUT, "SOLUTION: utilize --sjdbGTFfile /path/to/annotations.gtf option at the genome generation step or mapping step");
    geStream >> nGe;
    geID.resize(nGe);
    geName.resize(nGe);
    geBiotype.resize(nGe);
    geStream.ignore(999,'\n');
    for (uint ii=0;ii<nGe;ii++) {
        string line1;
        getline(geStream,line1);
        istringstream stream1(line1);
        stream1 >> geID[ii] >> geName[ii] >> geBiotype[ii];
    };
    geStream.close();

    if ( P.quant.trSAM.yes || P.quant.gene.yes || P.quant.geneFull_Ex50pAS.yes ) {//load exon-transcript structures
        //load tr and ex info
        ifstream & trinfo = trOpen(trInfoDir+"/transcriptInfo.tab", ERROR_OUT, "SOLUTION: utilize --sjdbGTFfile /path/to/annotantions.gtf option at the genome generation step or mapping step");
        trinfo >> nTr;
        trS=new uint [nTr];
        trE=new uint [nTr];
        trEmax=new uint [nTr];
        trExI=new uint32 [nTr];
        trExN=new uint16 [nTr];
        trStr=new uint8 [nTr];
        trID.resize(nTr);
        trGene=new uint32 [nTr];
        trLen=new uint32 [nTr];
        
        for (uint32 itr=0; itr<nTr; itr++) {
            uint16 str1;
            trinfo >> trID[itr] >> trS[itr] >> trE[itr] >> trEmax[itr] >> str1 >> trExN[itr] >> trExI[itr] >> trGene[itr];
            trStr[itr]=str1;

            if (!trinfo.good()) {
                ostringstream errOut;
                errOut <<"EXITING because of FATAL GENOME INDEX FILE error: transcriptInfo.tab is corrupt, or is incompatible with the current STAR version\n";
                errOut <<"SOLUTION: re-generate genome index";
                fatal(errOut.str(), EXIT_CODE_GENOME_FILES);
            };

        };
        logOut() << "Loaded transcript database, nTr="<<nTr<<endl;
        trinfo.close();

        ifstream & exinfo = trOpen(trInfoDir+"/exonInfo.tab", ERROR_OUT, "SOLUTION: utilize --sjdbGTFfile /path/to/annotantions.gtf option at the genome generation step or mapping step");
        exinfo >> nEx;
        exSE = new uint32 [2*nEx];
        exLenCum = new uint32 [nEx];
        for (uint32 iex=0; iex<nEx; iex++) {
            exinfo >> exSE[2*iex] >> exSE[2*iex+1] >> exLenCum[iex]; //reading all elements one after another
        };
        for (uint32 ii=0;ii<nTr;ii++) {
            uint32 iex1=trExI[ii]+trExN[ii]-1; //last exon of the transcript
            trLen[ii]=exLenCum[iex1]+exSE[2*iex1+1]-exSE[2*iex1]+1;
        };
        
        logOut() << "Loaded exon database, nEx="<<nEx<<endl;
        exinfo.close();
    };
    //load exon-gene structures
    if ( P.quant.geCount.yes ) {
        ifstream & exinfo = trOpen(trInfoDir+"/exonGeTrInfo.tab", ERROR_OUT, "SOLUTION: utilize --sjdbGTFfile /path/to/annotantions.gtf option at the genome generation step or mapping step");
        exinfo >> exG.nEx;
        exG.s=new uint64[exG.nEx];
        exG.e=new uint64[exG.nEx];
        exG.eMax=new uint64[exG.nEx];
        exG.str=new uint8[exG.nEx];
        exG.g=new uint32[exG.nEx];
        exG.t=new uint32[exG.nEx];
        for (uint ii=0;ii<exG.nEx;ii++) {
            int str1;
            exinfo >> exG.s[ii] >> exG.e[ii] >> str1 >> exG.g[ii] >> exG.t[ii];
            exG.str[ii] = (uint8) str1;
        };
        exinfo.close();
        //calculate eMax
        exG.eMax[0]=exG.e[0];
        for (uint iex=1;iex<exG.nEx;iex++) {
            exG.eMax[iex]=max(exG.eMax[iex-1],exG.e[iex]);
        };
    };

    if ( P.quant.geneFull.yes || P.quant.geneFull_ExonOverIntron.yes ) {
        ifstream & exinfo = trOpen(trInfoDir+"/exonGeTrInfo.tab", ERROR_OUT, "SOLUTION: utilize --sjdbGTFfile /path/to/annotantions.gtf option at the genome generation step or mapping step");
        exinfo >> exG.nEx;

        geneFull.s=new uint64[nGe];
        geneFull.e=new uint64[nGe];
        geneFull.eMax=new uint64[nGe];
        geneFull.g=new uint32[nGe];
        geneFull.str=new uint8[nGe];

        for (uint ig=0;ig<nGe;ig++) {
            geneFull.s[ig]=-1;//largest uint64
            geneFull.e[ig]=0;
        };

        for (uint ii=0;ii<exG.nEx;ii++) {
            uint64 s1,e1,str1,g1,t1;
            exinfo >> s1 >> e1 >> str1 >> g1 >> t1;
            geneFull.s[g1]=min(geneFull.s[g1],s1);
            geneFull.e[g1]=max(geneFull.e[g1],e1);
            geneFull.str[g1] = (uint8) str1;
        };
        exinfo.close();

        uint64 *gF=new uint64 [4*nGe];
        for (uint ii=0;ii<nGe;ii++) {
            gF[4*ii]   = geneFull.s[ii];
            gF[4*ii+1] = geneFull.e[ii];
            gF[4*ii+2] = geneFull.str[ii];
            gF[4*ii+3] = ii;
        };

        qsort((void*) gF, nGe, sizeof(uint64)*4, funCompareArrays<uint64,2>);

        for (uint ii=0;ii<nGe;ii++) {
            geneFull.s[ii]   = gF[4*ii];
            geneFull.e[ii]   = gF[4*ii+1];
            geneFull.str[ii] = gF[4*ii+2];
            geneFull.g[ii]   = gF[4*ii+3];
        };

        //calculate eMax
        geneFull.eMax[0]=geneFull.e[0];
        for (uint iex=1;iex<nGe;iex++) {
            geneFull.eMax[iex]=max(geneFull.eMax[iex-1],geneFull.e[iex]);
        };
    };

};

void Transcriptome::quantsAllocate() {
    if ( P.quant.geCount.yes ) {
        quants = new Quantifications (nGe);
    };
};

void Transcriptome::quantsOutput() {
    ofstream qOut(P.quant.geCount.outFile);
    qOut << "N_unmapped";
    for (int itype=0; itype<quants->geneCounts.nType; itype++) {
        qOut << "\t" <<g_statsAll.unmappedMismatch + g_statsAll.unmappedShort + g_statsAll.unmappedOther + g_statsAll.unmappedMulti;
    };
    qOut << "\n";

    qOut << "N_multimapping";
    for (int itype=0; itype<quants->geneCounts.nType; itype++){
        qOut << "\t" <<quants->geneCounts.cMulti;
    };
    qOut << "\n";

    qOut << "N_noFeature";
    for (int itype=0; itype<quants->geneCounts.nType; itype++){
        qOut << "\t" <<quants->geneCounts.cNone[itype];
    };
    qOut << "\n";

    qOut << "N_ambiguous";
    for (int itype=0; itype<quants->geneCounts.nType; itype++) {
        qOut << "\t" <<quants->geneCounts.cAmbig[itype];
    };
    qOut << "\n";

    for (uint32 ig=0; ig<nGe; ig++) {
        qOut << geID[ig];
        for (int itype=0; itype<quants->geneCounts.nType; itype++) {
            qOut << "\t" <<quants->geneCounts.gCount[itype][ig];
        };
        qOut << "\n";
    };
    qOut.close();
};

//below the constructor, so that its ERROR_OUT line numbers stay as they were
#include <thread>
#include <system_error>

class TranscriptomeEarly {//the Transcriptome constructed on a thread while genomeLoad reads the genome arrays
public:
    ostringstream logOut;//its Log.out lines, written by takeEarly
    string fatalError;//its fatal error, raised by takeEarly
    int fatalCode;
    Transcriptome *tr;
    std::thread thread;
    TranscriptomeEarly() : fatalCode(0), tr(NULL) {};
    ~TranscriptomeEarly() {
        if (thread.joinable())
            thread.join();
    };
};

struct TranscriptomeFatal {};//thrown on that thread once fatalError is stored

std::shared_ptr<TranscriptomeEarly> Transcriptome::loadEarly(Parameters &P)
{//reads only P and the files in trInfoDir, which nothing writes before the Transcriptome is constructed
    std::shared_ptr<TranscriptomeEarly> early=std::make_shared<TranscriptomeEarly>();
    TranscriptomeEarly *e=early.get();
    try {
        e->thread=std::thread([e, &P]() {
            try {
                e->tr=new Transcriptome(P, e);
            } catch (const TranscriptomeFatal &) {
            };
        });
    } catch (const std::system_error &) {//no thread: takeEarly is not called, the Transcriptome is constructed as usual
        return std::shared_ptr<TranscriptomeEarly>();
    };
    return early;
};

Transcriptome *Transcriptome::takeEarly(std::shared_ptr<TranscriptomeEarly> &early, Parameters &P)
{
    std::shared_ptr<TranscriptomeEarly> e=early;
    early.reset();
    e->thread.join();
    P.inOut->logMain << e->logOut.str() <<flush;
    if (e->fatalError!="")
        exitWithError(e->fatalError,std::cerr, P.inOut->logMain, e->fatalCode, P);
    e->tr->early=NULL;
    return e->tr;
};

std::ifstream & Transcriptome::trOpen(std::string fileName, std::string errorID, std::string solutionString)
{//ifstrOpen; on the early thread its error is stored instead of exiting
    if (early==NULL)
        return ifstrOpen(fileName, errorID, solutionString, P);

    std::ifstream & ifStream = *new std::ifstream(fileName.c_str());
    if (ifStream.fail()) {//ifstrOpen's message
        ostringstream errOut;
        errOut << errorID<<": exiting because of *INPUT FILE* error: could not open input file "<< fileName <<"\n";
        errOut << "Solution: check that the file exists and you have read permission for this file\n";
        if (solutionString.size()>0) {
            errOut << "          "<< solutionString <<"\n";
        };
        fatal(errOut.str(), EXIT_CODE_FILE_OPEN);
    };
    return ifStream;
};

void Transcriptome::fatal(const string &message, int errorInt)
{
    if (early==NULL)
        exitWithError(message,std::cerr, P.inOut->logMain, errorInt, P);
    early->fatalError=message;
    early->fatalCode=errorInt;
    throw TranscriptomeFatal();
};

ostream & Transcriptome::logOut()
{
    if (early==NULL)
        return P.inOut->logMain;
    return early->logOut;
};
