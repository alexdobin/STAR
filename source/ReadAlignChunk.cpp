#include "ReadAlignChunk.h"
#include <pthread.h>
#include "ErrorWarning.h"

ReadAlignChunk::ReadAlignChunk(Parameters& Pin, Genome &genomeIn, Transcriptome *TrIn, int iChunk) : P(Pin), mapGen(genomeIn) {//initialize chunk

    iThread=iChunk;

    if ( P.quant.yes ) {//allocate transcriptome structures
        chunkTr=new Transcriptome(*TrIn);
        chunkTr->quantsAllocate();
    } else {
        chunkTr=NULL;
    };

    RA = new ReadAlign(P, mapGen, chunkTr, iChunk);//new local copy of RA for each chunk
    
    RA->iRead=0;

    chunkIn=new char* [P.readNmates];
    readInStream=new istringstream* [P.readNmates];
//     readInStream=new istringstream* [P.readNmates];
    for (uint ii=0;ii<P.readNmates;ii++) {
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

    if ( P.quant.trSAM.yes ) {
        chunkOutBAMquant = new BAMoutput (P.inOut->outQuantBAMfile,P);
        RA->outBAMquant = chunkOutBAMquant;
    } else {
        chunkOutBAMquant=NULL;
        RA->outBAMquant=NULL;
    };

    chunkOutSJ=new OutSJ (P.limitOutSJcollapsed, P, mapGen);
    chunkOutSJ1=new OutSJ (P.limitOutSJcollapsed, P, mapGen);

    RA->chunkOutSJ=chunkOutSJ;
    RA->chunkOutSJ1=chunkOutSJ1;

    if (P.pCh.segmentMin>0) {
       if (P.pCh.out.samOld) {        
            chunkFstreamOpen(P.outFileTmp + "/Chimeric.out.sam.thread", iChunk, RA->chunkOutChimSAM);
       };
       if (P.pCh.out.junctions) {       
            chunkFstreamOpen(P.outFileTmp + "/Chimeric.out.junction.thread", iChunk, *RA->chunkOutChimJunction);
       };
    };
    if (P.outReadsUnmapped=="Fastx" ) {
        chunkFstreamOpen(P.outFileTmp + "/Unmapped.out.mate1.thread",iChunk, RA->chunkOutUnmappedReadsStream[0]);
        if (P.readNmates==2) chunkFstreamOpen(P.outFileTmp + "/Unmapped.out.mate2.thread",iChunk, RA->chunkOutUnmappedReadsStream[1]);
    };
    if (P.outFilterType=="BySJout") {
        chunkFstreamOpen(P.outFileTmp + "/FilterBySJoutFiles.mate1.thread",iChunk, RA->chunkOutFilterBySJoutFiles[0]);
        if (P.readNmates==2) chunkFstreamOpen(P.outFileTmp + "/FilterBySJoutFiles.mate2.thread",iChunk, RA->chunkOutFilterBySJoutFiles[1]);
    };
    
    if (P.wasp.yes) {
        RA->waspRA= new ReadAlign(Pin,genomeIn,TrIn,iChunk);
    };
    if (P.peOverlap.yes) {
        RA->peMergeRA= new ReadAlign(Pin,genomeIn,TrIn,iChunk);
        delete RA->peMergeRA->chunkOutChimJunction;
        RA->peMergeRA->chunkOutChimJunction=RA->chunkOutChimJunction;//point to the same out-stream
        RA->peMergeRA->chimDet->ostreamChimJunction=RA->peMergeRA->chunkOutChimJunction;
    };    
};

///////////////
void ReadAlignChunk::chunkFstreamOpen(string filePrefix, int iChunk, fstream &fstreamOut) {//open fstreams for chunks
    ostringstream fNameStream1;
    fNameStream1 << filePrefix << iChunk;
    string fName1=fNameStream1.str();
    P.inOut->logMain << "Opening the file: " << fName1 << " ... " <<flush;

    remove(fName1.c_str()); //remove the file
    fstreamOut.open(fName1.c_str(),ios::out); //create empty file
    fstreamOut.close();
    fstreamOut.open(fName1.c_str(), ios::in | ios::out); //re-open the file in in/out mode

    if (fstreamOut.fail()) {
        P.inOut->logMain << "failed!\n";
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not create output file "<< fName1 << "\n";
        errOut << "Solution: check that you have permission to write this file\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };
    P.inOut->logMain << "ok" <<endl;
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

