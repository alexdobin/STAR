#include "BAMoutput.h"
#include <sys/stat.h>
#include "GlobalVariables.h"
#include <pthread.h>

BAMoutput::BAMoutput (uint64 bamArraySizeIn, uint32 nBinsIn, uint64 genomeSize, int iChunk, string tmpDir) {//allocate bam array

    nBins=nBinsIn;
    binGlen=genomeSize/(nBins-1)+1;
    binSize=bamArraySizeIn/nBins;
    bamArraySize=binSize*nBins;
    bamArray = new char [bamArraySize];

    mkdir((tmpDir+to_string((uint) iChunk)).c_str(),S_IRWXU);    
    binStart=new char* [nBins];
    binBytes=new uint64 [nBins];    
    binStream=new ofstream* [nBins];
    binTotalN=new uint [nBins];
    binTotalBytes=new uint [nBins];
    for (uint ii=0;ii<nBins;ii++) {
        binStart[ii]=bamArray+bamArraySize/nBins*ii;
        binBytes[ii]=0;
        binStream[ii]=new ofstream((tmpDir+to_string((uint) iChunk) +"/"+to_string(ii)).c_str());    //open temporary files
        binTotalN[ii]=0;
        binTotalBytes[ii]=0;
    };
};

BAMoutput::BAMoutput (uint64 bamArraySizeIn, BGZF *bgzfBAMin) {//allocate BAM array with one bin, streamed directly into bgzf file
    bamArraySize=bamArraySizeIn;
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

void BAMoutput::coordOneAlign (char *bamIn, uint bamSize, uint chrStart, uint iRead) {
    
    if (bamSize==0) return; //no output, could happen if one of the mates is not mapped
    
    uint iBin;
    uint32 alignG= *( (uint32*) (bamIn+2*sizeof(uint32)) );
    if ( alignG == (uint32) -1 ) {//unmapped alignment, last bin
        iBin=nBins-1;
    } else {
        iBin=(alignG + chrStart)/binGlen;
    };
        
    if (binBytes[iBin]+bamSize+sizeof(uint) > binSize) {//write out this buffer
        binStream[iBin]->write(binStart[iBin],binBytes[iBin]);
        binBytes[iBin]=0;//rewind the buffer
    };
    
    memcpy(binStart[iBin]+binBytes[iBin], bamIn, bamSize);
    binBytes[iBin] += bamSize;
    memcpy(binStart[iBin]+binBytes[iBin], &iRead, sizeof(uint));
    binBytes[iBin] += sizeof(uint);
    binTotalBytes[iBin] += bamSize+sizeof(uint);
    binTotalN[iBin] += 1;
    
};

void BAMoutput::coordFlush () {//flush all alignments
    for (uint32 iBin=0; iBin<nBins; iBin++) {
        binStream[iBin]->write(binStart[iBin],binBytes[iBin]);
        binStream[iBin]->flush();
        binBytes[iBin]=0;//rewind the buffer
    };
};
