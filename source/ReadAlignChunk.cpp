#include "ReadAlignChunk.h"
#include <pthread.h>
#include <thread>
#include <system_error>
#include "ErrorWarning.h"
#include "GlobalVariables.h"
#ifdef __GLIBC__
#include <malloc.h>
#endif

void ReadAlignChunk::buildAll(Parameters &P, Genome &genomeIn, Transcriptome *TrIn, ReadAlignChunk **RAchunk)
{//new ReadAlignChunk for each of the P.runThreadN chunks; above 1 thread, one thread per chunk, and the chunks' Log.out lines written in chunk order after all are built
    int nChunks=P.runThreadN;
    if (nChunks<2) {
        for (int ii=0; ii<nChunks; ii++)
            RAchunk[ii]=new ReadAlignChunk(P, genomeIn, TrIn, ii);
        return;
    };

    vector<ostringstream> logBufs(nChunks);
    vector<std::thread> threads;
    threads.reserve(nChunks-1);
    vector<int> unstarted;//chunks whose thread could not be created: built on this thread
#ifdef M_TOP_PAD
    mallopt(M_TOP_PAD, 64<<20);//a new thread arena then starts with its whole heap writable; grown page by page, its mprotect calls serialize the other threads' page faults
#endif
    for (int ii=1; ii<nChunks; ii++) {
        try {
            threads.emplace_back([&P, &genomeIn, TrIn, RAchunk, &logBufs, ii]() {RAchunk[ii]=new ReadAlignChunk(P, genomeIn, TrIn, ii, &logBufs[ii]);});
        } catch (const std::system_error &) {
            unstarted.push_back(ii);
        };
    };
    RAchunk[0]=new ReadAlignChunk(P, genomeIn, TrIn, 0, &logBufs[0]);
    for (int ii : unstarted)
        RAchunk[ii]=new ReadAlignChunk(P, genomeIn, TrIn, ii, &logBufs[ii]);
    for (auto &th : threads)
        th.join();
#ifdef M_TOP_PAD
    mallopt(M_TOP_PAD, 128*1024);//glibc's default
#endif

    for (int ii=0; ii<nChunks; ii++) {
        P.inOut->logMain << logBufs[ii].str();
        RAchunk[ii]->logBuf=NULL;
    };
    P.inOut->logMain << flush;
};

ReadAlignChunk::ReadAlignChunk(Parameters& Pin, Genome &genomeIn, Transcriptome *TrIn, int iChunk, ostringstream *logBufIn) : P(Pin), mapGen(genomeIn) {//initialize chunk

    logBuf=logBufIn;
    iThread=iChunk;

    if ( P.quant.yes ) {//allocate transcriptome structures
        chunkTr=new Transcriptome(*TrIn);
        chunkTr->quantsAllocate();
    } else {
        chunkTr=NULL;
    };

    RA = new ReadAlign(P, mapGen, chunkTr, iChunk);//new local copy of RA for each chunk

    RA->iRead=0;

    chunkIn=new char* [P.readNends];
    readInStream=new istringstream* [P.readNends];
    
    for (uint ii=0;ii<P.readNends;ii++) {
       chunkIn[ii]=new char[P.chunkInSizeBytesArray];//reserve more space to finish loading one read
       memset(chunkIn[ii],'\n',P.chunkInSizeBytesArray);
       readInStream[ii] = new istringstream;
       readInStream[ii]->rdbuf()->pubsetbuf(chunkIn[ii],P.chunkInSizeBytesArray);
       RA->readInStream[ii]=readInStream[ii];
    };


    if (P.outSAMbool) {
        chunkOutBAM=new char [P.chunkOutBAMsizeBytes];
        RA->outBAMarray=chunkOutBAM;
        chunkOutBAMstream=new ostringstream;
        chunkOutBAMstream->rdbuf()->pubsetbuf(chunkOutBAM,P.chunkOutBAMsizeBytes);
        RA->outSAMstream=chunkOutBAMstream;
        RA->outSAMstream->seekp(0,ios::beg);
        chunkOutBAMtotal=0;
    };

    if (P.outBAMunsorted) {
        chunkOutBAMunsorted = new BAMoutput (P.inOut->outBAMfileUnsorted, P);
        RA->outBAMunsorted = chunkOutBAMunsorted;
    } else {
        chunkOutBAMunsorted=NULL;
        RA->outBAMunsorted=NULL;
    };

    if (P.outBAMcoord) {
        chunkOutBAMcoord = new BAMoutput (iChunk, P.outBAMsortTmpDir, P);
        RA->outBAMcoord = chunkOutBAMcoord;
    } else {
        chunkOutBAMcoord=NULL;
        RA->outBAMcoord=NULL;
    };

    if ( P.quant.trSAM.bamYes ) {
        chunkOutBAMquant = new BAMoutput (P.inOut->outQuantBAMfile,P);
        RA->outBAMquant = chunkOutBAMquant;
    } else {
        chunkOutBAMquant=NULL;
        RA->outBAMquant=NULL;
    };

    if (P.outSJ.yes) {
        chunkOutSJ  = new OutSJ (P.limitOutSJcollapsed, P, mapGen);
        RA->chunkOutSJ  = chunkOutSJ;
    } else {
        RA->chunkOutSJ  = NULL;
    };

    if (P.outFilterBySJoutStage == 1) {
        chunkOutSJ1 = new OutSJ (P.limitOutSJcollapsed, P, mapGen);
        RA->chunkOutSJ1 = chunkOutSJ1;
    } else {
        RA->chunkOutSJ1  = NULL;
    };

    
    

    if (P.pCh.segmentMin>0) {
       if (P.pCh.out.samOld) {
            chunkFstreamOpen(P.outFileTmp + "/Chimeric.out.sam.thread", iChunk, RA->chunkOutChimSAM);
       };
       if (P.pCh.out.junctions) {
            chunkFstreamOpen(P.outFileTmp + "/Chimeric.out.junction.thread", iChunk, *RA->chunkOutChimJunction);
       };
    };

    if (P.outReadsUnmapped=="Fastx" ) {
        for (uint32 imate=0; imate < P.readNends; imate++) 
            chunkFstreamOpen(P.outFileTmp + "/Unmapped.out.mate"+ to_string(imate) +".thread",iChunk, RA->chunkOutUnmappedReadsStream[imate]);
    };

    if (P.outFilterType=="BySJout") {
        chunkFstreamOpen(P.outFileTmp + "/FilterBySJoutFiles.mate1.thread",iChunk, RA->chunkOutFilterBySJoutFiles[0]);
        if (P.readNends==2) chunkFstreamOpen(P.outFileTmp + "/FilterBySJoutFiles.mate2.thread",iChunk, RA->chunkOutFilterBySJoutFiles[1]); //here we do not output barcode read
    };

    if (P.wasp.yes) {
        RA->waspRA= new ReadAlign(Pin,genomeIn,TrIn,iChunk);
    };
    if (P.peOverlap.yes) {
        RA->peMergeRA= new ReadAlign(Pin,genomeIn,TrIn,iChunk);
        delete RA->peMergeRA->chunkOutChimJunction;
        RA->peMergeRA->chunkOutChimJunction=RA->chunkOutChimJunction;//point to the same out-stream
        RA->peMergeRA->chimDet->ostreamChimJunction=RA->peMergeRA->chunkOutChimJunction;
        RA->peMergeRA->outBAMunsorted=RA->outBAMunsorted;
        RA->peMergeRA->outBAMcoord=RA->outBAMcoord;
    };
};

///////////////
void ReadAlignChunk::chunkFstreamOpen(string filePrefix, int iChunk, fstream &fstreamOut) {//open fstreams for chunks
    ostringstream fNameStream1;
    fNameStream1 << filePrefix << iChunk;
    string fName1=fNameStream1.str();
    ostream &logOut = logBuf!=NULL ? (ostream&) *logBuf : (ostream&) P.inOut->logMain;
    logOut << "Opening the file: " << fName1 << " ... " <<flush;

    remove(fName1.c_str()); //remove the file
    fstreamOut.open(fName1.c_str(),ios::out); //create empty file
    fstreamOut.close();
    fstreamOut.open(fName1.c_str(), ios::in | ios::out); //re-open the file in in/out mode

    if (fstreamOut.fail()) {
        logOut << "failed!\n";
        if (logBuf!=NULL) {//under the mutex exitWithError keeps until the exit, so never after another thread's exit began
            pthread_mutex_lock(&g_threadChunks.mutexError);
            P.inOut->logMain << logBuf->str() << flush;
            pthread_mutex_unlock(&g_threadChunks.mutexError);
        };
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not create output file "<< fName1 << "\n";
        errOut << "Solution: check that you have permission to write this file\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };
    logOut << "ok" <<endl;
};

void ReadAlignChunk::chunkFstreamCat (fstream &chunkOut, ofstream &allOut, bool mutexFlag, pthread_mutex_t &mutexVal){
    chunkOut.flush();
    chunkOut.seekg(0,ios::beg);
    if (mutexFlag) pthread_mutex_lock(&mutexVal);
    allOut << chunkOut.rdbuf();
    allOut.clear();
    allOut.flush();
    allOut.clear();
    if (mutexFlag) pthread_mutex_unlock(&mutexVal);
    chunkOut.clear();
    chunkOut.seekp(0,ios::beg); //set put pointer at the beginning
};


void ReadAlignChunk::chunkFilesCat(ostream *allOut, string filePrefix, uint &iC) {//concatenates a file into main output
            while (true) {
                ostringstream name1("");
                name1 << filePrefix <<iC;
                ifstream fileChunkIn(name1.str().c_str());
                if (fileChunkIn.good()) {
                    *allOut << fileChunkIn.rdbuf();
                    allOut->flush();
                    allOut->clear();
                    fileChunkIn.close();
                    fileChunkIn.clear();
                    remove(name1.str().c_str());
                    iC++;
                } else {
                    fileChunkIn.close();
                    break;
                };
            };
};

