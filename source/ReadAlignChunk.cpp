#include "ReadAlignChunk.h"
#include "sysRemoveDir.h"

#if !defined(_WIN32) && defined(USE_PTHREAD)
#include <pthread.h>
#else
#include <thread>
#endif
#include "ErrorWarning.h"
#ifdef _WIN32
#include "DirFunctions.h"
#endif

ReadAlignChunk::ReadAlignChunk(Parameters* Pin, Genome &genomeIn, Transcriptome *TrIn, int iChunk) : P(Pin) {//initialize chunk

	iThread = iChunk;

	if (P->quant.yes) {//allocate transcriptome structures
		chunkTr = new Transcriptome(*TrIn);
		chunkTr->quantsAllocate();
	}
	else {
		chunkTr = NULL;
	};

	RA = new ReadAlign(P, genomeIn, chunkTr, iChunk);//new local copy of RA for each chunk

	RA->iRead = 0;

	chunkIn = new char*[P->readNmates];
	readInStream = new istringstream*[P->readNmates];
	for (uint ii = 0; ii < P->readNmates; ii++) {
		chunkIn[ii] = new char[P->chunkInSizeBytesArray];//reserve more space to finish loading one read
		memset(chunkIn[ii], '\n', P->chunkInSizeBytesArray);

#ifndef _WIN32 // VS implementation of pubsetbuf does nothing, do this only for Non Windows
		readInStream[ii] = new istringstream;
		readInStream[ii]->rdbuf()->pubsetbuf(chunkIn[ii],P->chunkInSizeBytesArray);
		RA->readInStream[ii] = readInStream[ii];
#else
		// For VS, workaround using istringstream constructor added in ReadAlignChunk::mapChunk()
		// Here just Initialize to nullptr, we need to check readInStream pointer and delete in ReadAlignChunk::mapChunk()
		readInStream[ii] = nullptr;
#endif
	};

	if (P->outSAMbool) {
		chunkOutBAM = new char[P->chunkOutBAMsizeBytes];
		RA->outBAMarray = chunkOutBAM;

#ifndef _WIN32 // VS implementation of pubsetbuf does nothing, do this only for Non Windows.
		chunkOutBAMstream = new ostringstream;
		chunkOutBAMstream->rdbuf()->pubsetbuf(chunkOutBAM, P->chunkOutBAMsizeBytes);
#else // For VS, workaround using ostringstream constructor to set buffer.
		chunkOutBAMstream = new ostringstream(std::string(chunkOutBAM, P->chunkOutBAMsizeBytes), std::stringstream::in | std::stringstream::out);
#endif

		RA->outSAMstream = chunkOutBAMstream;
		RA->outSAMstream->seekp(0, ios::beg);
		chunkOutBAMtotal = 0;
	};

	if (P->outBAMunsorted) {
		chunkOutBAMunsorted = new BAMoutput(P->inOut->outBAMfileUnsorted, P);
		RA->outBAMunsorted = chunkOutBAMunsorted;
	}
	else {
		chunkOutBAMunsorted = NULL;
		RA->outBAMunsorted = NULL;
	};

	if (P->outBAMcoord) {
		chunkOutBAMcoord = new BAMoutput(iChunk, P->outBAMsortTmpDir, P);
		RA->outBAMcoord = chunkOutBAMcoord;
	}
	else {
		chunkOutBAMcoord = NULL;
		RA->outBAMcoord = NULL;
	};

	if (P->quant.trSAM.yes) {
		chunkOutBAMquant = new BAMoutput(P->inOut->outQuantBAMfile, P);
		RA->outBAMquant = chunkOutBAMquant;
	}
	else {
		chunkOutBAMquant = NULL;
		RA->outBAMquant = NULL;
	};

	chunkOutSJ = new OutSJ(P->limitOutSJcollapsed, P);
	chunkOutSJ1 = new OutSJ(P->limitOutSJcollapsed, P);

	RA->chunkOutSJ = chunkOutSJ;
	RA->chunkOutSJ1 = chunkOutSJ1;

	if (P->chimSegmentMin>0) {
		chunkFstreamOpen(P->outFileTmp + "/Chimeric.out.sam.thread", iChunk, RA->chunkOutChimSAM);
		chunkFstreamOpen(P->outFileTmp + "/Chimeric.out.junction.thread", iChunk, RA->chunkOutChimJunction);
	};
	if (P->outReadsUnmapped == "Fastx") {
		chunkFstreamOpen(P->outFileTmp + "/Unmapped.out.mate1.thread", iChunk, RA->chunkOutUnmappedReadsStream[0]);
		if (P->readNmates == 2) chunkFstreamOpen(P->outFileTmp + "/Unmapped.out.mate2.thread", iChunk, RA->chunkOutUnmappedReadsStream[1]);
	};
	if (P->outFilterType == "BySJout") {
		chunkFstreamOpen(P->outFileTmp + "/FilterBySJoutFiles.mate1.thread", iChunk, RA->chunkOutFilterBySJoutFiles[0]);
		if (P->readNmates == 2) chunkFstreamOpen(P->outFileTmp + "/FilterBySJoutFiles.mate2.thread", iChunk, RA->chunkOutFilterBySJoutFiles[1]);
	};
};


///////////////
void ReadAlignChunk::chunkFstreamOpen(string filePrefix, int iChunk, fstream &fstreamOut) {//open fstreams for chunks
    ostringstream fNameStream1;
    fNameStream1 << filePrefix << iChunk;
    string fName1=fNameStream1.str();
    P->inOut->logMain << "Opening the file: " << fName1 << " ... " <<flush;
    
    remove(fName1.c_str()); //remove the file
	fstreamOut.open(fName1.c_str(), ios::out | std::ios::binary); //create empty file
    fstreamOut.close();
	fstreamOut.open(fName1.c_str(), ios::in | ios::out | std::ios::binary); //re-open the file in in/out mode
    if (fstreamOut.fail()) {
        P->inOut->logMain << "failed!\n";
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not create output file "<< fName1 << "\n";
        errOut << "Solution: check that you have permission to write this file\n";        
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);        
    };
    P->inOut->logMain << "ok" <<endl;
};


#if !defined(_WIN32) && defined(USE_PTHREAD)

void ReadAlignChunk::chunkFstreamCat (fstream &chunkOut, ofstream &allOut, bool mutexFlag, pthread_mutex_t &mutexVal){
    chunkOut.flush();
    chunkOut.seekg(0,ios::beg);
    if (mutexFlag) pthread_mutex_lock(&mutexVal);
    allOut << chunkOut.rdbuf();
    allOut.clear();
    if (mutexFlag) pthread_mutex_unlock(&mutexVal);
    chunkOut.clear();
    chunkOut.seekp(0,ios::beg); //set put pointer at the beginning
};

#else

void ReadAlignChunk::chunkFstreamCat(fstream &chunkOut, ofstream &allOut, bool mutexFlag, std::mutex &mutexVal){
	chunkOut.flush();
	chunkOut.seekg(0, ios::beg);
	if (mutexFlag) mutexVal.lock();
	allOut << chunkOut.rdbuf();
	allOut.clear();
	if (mutexFlag) mutexVal.unlock();
	chunkOut.clear();
	chunkOut.seekp(0, ios::beg); //set put pointer at the beginning
};

#endif


void ReadAlignChunk::chunkFilesCat(ostream *allOut, string filePrefix, uint &iC) {//concatenates a file into main output
            while (true) {
                ostringstream name1("");
                name1 << filePrefix <<iC;
				ifstream fileChunkIn(name1.str().c_str(), ios_base::in | ios_base::binary);
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

void ReadAlignChunk::closeReadAlignFiles()
{
	if (!RA)
		return; 

	if (RA->chunkOutChimSAM.is_open())
		RA->chunkOutChimSAM.close(); 

	if (RA->chunkOutChimJunction.is_open())
		RA->chunkOutChimJunction.close();
	
	for (unsigned int i = 0; i < MAX_N_MATES; i++)
	{
		if (RA->chunkOutUnmappedReadsStream[i].is_open())
			RA->chunkOutUnmappedReadsStream[i].close();
	}

	for (unsigned int i = 0; i < MAX_N_MATES; i++)
	{
		if (RA->chunkOutFilterBySJoutFiles[i].is_open())
			RA->chunkOutFilterBySJoutFiles[i].close();
	}
}
