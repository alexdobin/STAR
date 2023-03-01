#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ReadAlign.h"

ReadAlign::ReadAlign (Parameters& Pin, Genome &genomeIn, Transcriptome *TrIn, int iChunk)
                    : mapGen(genomeIn), genOut(*genomeIn.genomeOut.g), P(Pin), chunkTr(TrIn)
{
    readNmates=P.readNmates; //not readNends
    //RNGs
    rngMultOrder.seed(P.runRNGseed*(iChunk+1));
    rngUniformReal0to1=std::uniform_real_distribution<double> (0.0, 1.0);
    //transcriptome
    if ( P.quant.trSAM.yes ) {
        alignTrAll=new Transcript [P.alignTranscriptsPerReadNmax];
    };

    if (P.pGe.gType==101) {//SuperTranscriptome
        splGraph = new SpliceGraph(*mapGen.superTr, P, this);
    } else {//standard map algorithm:
        winBin = new uintWinBin* [2];
        winBin[0] = new uintWinBin [P.winBinN];
        winBin[1] = new uintWinBin [P.winBinN];
        memset(winBin[0],255,sizeof(winBin[0][0])*P.winBinN);
        memset(winBin[1],255,sizeof(winBin[0][0])*P.winBinN);
        //split
        splitR=new uint*[3];
        splitR[0]=new uint[P.maxNsplit]; splitR[1]=new uint[P.maxNsplit]; splitR[2]=new uint[P.maxNsplit];
        //alignments
        PC=new uiPC[P.seedPerReadNmax];
        WC=new uiWC[P.alignWindowsPerReadNmax];
        nWA=new uint[P.alignWindowsPerReadNmax];
        nWAP=new uint[P.alignWindowsPerReadNmax];
        WALrec=new uint[P.alignWindowsPerReadNmax];
        WlastAnchor=new uint[P.alignWindowsPerReadNmax];
    
        WA=new uiWA*[P.alignWindowsPerReadNmax];
        for (uint ii=0;ii<P.alignWindowsPerReadNmax;ii++)
            WA[ii]=new uiWA[P.seedPerWindowNmax];
        WAincl = new bool [P.seedPerWindowNmax];        

        #ifdef COMPILE_FOR_LONG_READS
        swWinCov = new uint[P.alignWindowsPerReadNmax];
        scoreSeedToSeed = new intScore [P.seedPerWindowNmax*(P.seedPerWindowNmax+1)/2];
        scoreSeedBest = new intScore [P.seedPerWindowNmax];
        scoreSeedBestInd = new uint [P.seedPerWindowNmax];
        scoreSeedBestMM = new uint [P.seedPerWindowNmax];
        seedChain = new uint [P.seedPerWindowNmax];
        #endif
    };

    //aligns a.k.a. transcripts
    trAll = new Transcript**[P.alignWindowsPerReadNmax+1];
    nWinTr = new uint[P.alignWindowsPerReadNmax];
    trArray = new Transcript[P.alignTranscriptsPerReadNmax];
    trArrayPointer =  new Transcript*[P.alignTranscriptsPerReadNmax];
    for (uint ii=0;ii<P.alignTranscriptsPerReadNmax;ii++)
        trArrayPointer[ii]= &(trArray[ii]);
    trInit = new Transcript;
    
    if (mapGen.genomeOut.convYes) {//allocate output transcripts
        alignsGenOut.alMult = new Transcript*[P.outFilterMultimapNmax];
        for (uint32 ii=0; ii<P.outFilterMultimapNmax; ii++) 
            alignsGenOut.alMult[ii]=new Transcript;
    };
    
    //read
    Read0 = new char*[P.readNends];
    Qual0 = new char*[P.readNends];
    readNameMates=new char* [P.readNends];
    for (uint32 ii=0; ii<P.readNends; ii++) {
        readNameMates[ii]= new char [DEF_readNameLengthMax];
        Read0[ii]        = new char [DEF_readSeqLengthMax+1];
        Qual0[ii]        = new char [DEF_readSeqLengthMax+1];        
    };
    readNameExtra.resize(P.readNends);
    readName = readNameMates[0];

    Read1 = new char*[3];
    Read1[0]=new char[DEF_readSeqLengthMax+1]; 
    Read1[1]=new char[DEF_readSeqLengthMax+1]; 
    Read1[2]=new char[DEF_readSeqLengthMax+1];
    
    for (auto &q: qualHist)
        q.fill(0);
    
    //outBAM
    outBAMoneAlignNbytes = new uint [P.readNmates+2]; //extra piece for chimeric reads //not readNends: this is alignment
    outBAMoneAlign = new char* [P.readNmates+2]; //extra piece for chimeric reads //not readNends: this is alignment
    for (uint ii=0; ii<P.readNmates+2; ii++) {//not readNends: this is alignment
        outBAMoneAlign[ii]=new char [BAMoutput_oneAlignMaxBytes];
    };
    resetN();
    
    //chim
    chunkOutChimJunction = new fstream;
    chimDet = new ChimericDetection(P, trAll, nWinTr, Read1, mapGen, chunkOutChimJunction, this);
    
    //solo
    soloRead = new SoloRead (P, iChunk);
    
    //clipping
    P.pClip.initializeClipMates(clipMates);

    //debug
    {
    #ifdef DEBUG_OutputLastRead
        lastReadStream.open((P.outFileTmp+"/lastRead_"+to_string(iChunk)).c_str());
    #endif
    };
};

void ReadAlign::resetN () {//reset resets the counters to 0 for a new read
    mapMarker=0;
    nA=0; nP=0; nW=0;
    nTr=0;
    nUM[0]=0; nUM[1]=0;
    storedLmin=0; uniqLmax=0; uniqLmaxInd=0; multLmax=0; multLmaxN=0; multNminL=0; multNmin=0; multNmax=0; multNmaxL=0;
    chimN=0;

    for (uint ii=0; ii<P.readNmates; ii++) {//not readNends: this is alignment
        maxScoreMate[ii]=0;
    };
};

