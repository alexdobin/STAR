#include "BAMoutput.h"
#include <sys/stat.h>
#include "GlobalVariables.h"
#include <pthread.h>
#include "serviceFuns.cpp"
#include "ThreadControl.h"

BAMoutput::BAMoutput (int iChunk, string tmpDir, Parameters &Pin) : P(Pin){//allocate bam array

    nBins=P.outBAMcoordNbins;
    binSize=P.chunkOutBAMsizeBytes/nBins;
    bamArraySize=binSize*nBins;
    bamArray = new char [bamArraySize];

    bamDir=tmpDir+to_string((uint) iChunk);//local directory for this thread (iChunk)

    mkdir(bamDir.c_str(),P.runDirPerm);
    binStart=new char* [nBins];
    binBytes=new uint64 [nBins];
    binStream=new ofstream* [nBins];
    binTotalN=new uint [nBins];
    binTotalBytes=new uint [nBins];
    for (uint ii=0;ii<nBins;ii++) {
        binStart[ii]=bamArray+bamArraySize/nBins*ii;
        binBytes[ii]=0;
        binStream[ii]=new ofstream((bamDir +"/"+to_string(ii)).c_str());    //open temporary files
        binTotalN[ii]=0;
        binTotalBytes[ii]=0;
    };

    binSize1=binStart[nBins-1]-binStart[0];
    nBins=1;//start with one bin to estimate genomic bin sizes
};

BAMoutput::BAMoutput (BGZF *bgzfBAMin, Parameters &Pin) : P(Pin){//allocate BAM array with one bin, streamed directly into bgzf file

    bamArraySize=P.chunkOutBAMsizeBytes;
    bamArray = new char [bamArraySize];
    binBytes1=0;
    bgzfBAM=bgzfBAMin;
    //not used
    binSize=0;
    binStream=NULL;
    binStart=NULL;
    binBytes=NULL;
    binTotalBytes=NULL;
    binTotalN=NULL;
    nBins=0;
};

void BAMoutput::unsortedOneAlign (char *bamIn, uint bamSize, uint bamSize2) {//record one alignment to the buffer, write buffer if needed

    if (bamSize==0) return; //no output, could happen if one of the mates is not mapped

    if (binBytes1+bamSize2 > bamArraySize) {//write out this buffer

        if (g_threadChunks.threadBool) pthread_mutex_lock(&g_threadChunks.mutexOutSAM);
        bgzf_write(bgzfBAM,bamArray,binBytes1);
        if (g_threadChunks.threadBool) pthread_mutex_unlock(&g_threadChunks.mutexOutSAM);

        binBytes1=0;//rewind the buffer
    };

    memcpy(bamArray+binBytes1, bamIn, bamSize);
    binBytes1 += bamSize;

};

void BAMoutput::unsortedFlush () {//flush all alignments
    if (g_threadChunks.threadBool) pthread_mutex_lock(&g_threadChunks.mutexOutSAM);
    bgzf_write(bgzfBAM,bamArray,binBytes1);
    if (g_threadChunks.threadBool) pthread_mutex_unlock(&g_threadChunks.mutexOutSAM);
    binBytes1=0;//rewind the buffer
};

void BAMoutput::coordOneAlign (char *bamIn, uint bamSize, uint iRead) {

    uint32 *bamIn32;
    uint alignG;
    uint32 iBin=0;

    if (bamSize==0) {
        return; //no output, could happen if one of the mates is not mapped
    } else {
        //determine which bin this alignment belongs to
        bamIn32=(uint32*) bamIn;
        alignG=( ((uint) bamIn32[1]) << 32 ) | ( (uint)bamIn32[2] );
        if (bamIn32[1] == ((uint32) -1) ) {//unmapped
            iBin=P.outBAMcoordNbins-1;
        } else if (nBins>1) {//bin starts have already been determined
            iBin=binarySearch1a <uint64> (alignG, P.outBAMsortingBinStart, (int32) (nBins-1));
        };
    };

//     if ( alignG == (uint32) -1 ) {//unmapped alignment, last bin
//         iBin=nBins-1;
//     } else {
//         iBin=(alignG + chrStart)/binGlen;
//     };

    //write buffer is filled
    if (binBytes[iBin]+bamSize+sizeof(uint) > ( (iBin>0 || nBins>1) ? binSize : binSize1) ) {//write out this buffer
        if ( nBins>1 || iBin==(P.outBAMcoordNbins-1) ) {//normal writing, bins have already been determined
            binStream[iBin]->write(binStart[iBin],binBytes[iBin]);
            binBytes[iBin]=0;//rewind the buffer
        } else {//the first chunk of reads was written in one bin, need to determine bin sizes, and re-distribute reads into bins
            coordBins();
            coordOneAlign (bamIn, bamSize, iRead);//record the current align into the new bins
            return;
        };
    };

    //record this alignment in its bin
    memcpy(binStart[iBin]+binBytes[iBin], bamIn, bamSize);
    binBytes[iBin] += bamSize;
    memcpy(binStart[iBin]+binBytes[iBin], &iRead, sizeof(uint));
    binBytes[iBin] += sizeof(uint);
    binTotalBytes[iBin] += bamSize+sizeof(uint);
    binTotalN[iBin] += 1;
    return;
};

void BAMoutput::coordBins() {//define genomic starts for bins
    nBins=P.outBAMcoordNbins;//this is the true number of bins

    //mutex here
    if (P.runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexBAMsortBins);
    if (P.outBAMsortingBinStart[0]!=0) {//it's set to 0 only after the bin sizes are determined
        //extract coordinates and sort
        uint *startPos = new uint [binTotalN[0]];//array of aligns start positions
        for (uint ib=0,ia=0;ia<binTotalN[0];ia++) {
            uint32 *bamIn32=(uint32*) (binStart[0]+ib);
            startPos[ia]  =( ((uint) bamIn32[1]) << 32) | ( (uint)bamIn32[2] );
            ib+=bamIn32[0]+sizeof(uint32)+sizeof(uint);//note that size of the BAM record does not include the size record itself
        };
        qsort((void*) startPos, binTotalN[0], sizeof(uint), funCompareUint1);

        //determine genomic starts of the bins
        P.inOut->logMain << "BAM sorting: "<<binTotalN[0]<< " mapped reads\n";
        P.inOut->logMain << "BAM sorting bins genomic start loci:\n";

        P.outBAMsortingBinStart[0]=0;
        for (uint32 ib=1; ib<(nBins-1); ib++) {
            P.outBAMsortingBinStart[ib]=startPos[binTotalN[0]/(nBins-1)*ib];
            P.inOut->logMain << ib <<"\t"<< (P.outBAMsortingBinStart[ib]>>32) << "\t" << ((P.outBAMsortingBinStart[ib]<<32)>>32) <<endl;
            //how to deal with equal boundaries???
        };
        delete [] startPos;
    };
    //mutex here
    if (P.runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexBAMsortBins);

    //re-allocate binStart
    uint binTotalNold=binTotalN[0];
    char *binStartOld=new char [binSize1];
    memcpy(binStartOld,binStart[0],binBytes[0]);

    binBytes[0]=0;
    binTotalN[0]=0;
    binTotalBytes[0]=0;

    //re-bin all aligns
    for (uint ib=0,ia=0;ia<binTotalNold;ia++) {
        uint32 *bamIn32=(uint32*) (binStartOld+ib);
        uint ib1=ib+bamIn32[0]+sizeof(uint32);//note that size of the BAM record does not include the size record itself
        coordOneAlign (binStartOld+ib, (uint) (bamIn32[0]+sizeof(uint32)), *((uint*) (binStartOld+ib1)) );
        ib=ib1+sizeof(uint);//iRead at the end of the BAM record
    };
    delete [] binStartOld;
    return;
};

void BAMoutput::coordFlush () {//flush all alignments
    if (nBins==1) {
        coordBins();
    };
    for (uint32 iBin=0; iBin<nBins; iBin++) {
        binStream[iBin]->write(binStart[iBin],binBytes[iBin]);
        binStream[iBin]->flush();
        binBytes[iBin]=0;//rewind the buffer
    };
};

void BAMoutput::coordUnmappedPrepareBySJout () {//flush all alignments
    uint iBin=P.outBAMcoordNbins-1;
    binStream[iBin]->write(binStart[iBin],binBytes[iBin]);
    binStream[iBin]->flush();
    binBytes[iBin]=0;//rewind the buffer
    binStream[iBin]->close();
    binStream[iBin]->open((bamDir +"/"+to_string(iBin)+".BySJout").c_str());
};
