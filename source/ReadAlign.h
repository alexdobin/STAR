#ifndef CODE_READALIGN
#define CODE_READALIGN

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "Genome.h"
#include "Stats.h"
#include "OutSJ.h"
#include <time.h>
#include "Transcriptome.h"
#include "BAMoutput.h"
#include "Quantifications.h"
#include <random>

class ReadAlign : public Genome 
{
    public:
        Parameters* P; //pointer to the parameters, will be initialized on construction
          
        //mapping statistics
        Stats statsRA;
        
        //quantification
        Transcriptome *chunkTr;
        
        //mapping time
        time_t timeStart, timeFinish;
        
        //random number generators
        std::mt19937 rngMultOrder;//initialize in ReadAlign.cpp
        std::uniform_real_distribution<double> rngUniformReal0to1;//initialize in ReadAlign.cpp
        
        //input,output
        istream* readInStream[MAX_N_MATES];
        ostream* outSAMstream;
        OutSJ *chunkOutSJ, *chunkOutSJ1;
        fstream chunkOutChimSAM, chunkOutChimJunction, chunkOutUnmappedReadsStream[MAX_N_MATES], chunkOutFilterBySJoutFiles[MAX_N_MATES];
        uint outBAMbytes, outBAMbytes1; //number of bytes output to SAM/BAM with oneRead
        char *outBAMarray, *outBAMarray1;//pointer to the (last+1) position of the SAM/BAM output array
        BAMoutput *outBAMcoord, *outBAMunsorted, *outBAMquant;//sorted by coordinate, unsorted, transcriptomic BAM structure
//        char outBAMoneAlign[MAX_N_MATES][BAMoutput_oneAlignMaxBytes];//tmp array to store BAM alignmnent
//        uint outBAMoneAlignNbytes[MAX_N_MATES];//number of bytes in the tmp BAM array
        char** outBAMoneAlign;
        uint* outBAMoneAlignNbytes;
        
        ostringstream samStreamCIGAR, samStreamSJmotif, samStreamSJintron,samStreamSJannot;
        
        intScore maxScoreMate[MAX_N_MATES];
        intScore *scoreSeedToSeed, *scoreSeedBest;
        uint *scoreSeedBestInd, *seedChain, *scoreSeedBestMM;
        
        bool outFilterPassed; //true if alignment passed all filter and is output to SAM/BAM
        
//         StatsAll *statsRA;
        
        //transcript
        Transcript* trArray; //linear array of transcripts to store all of them from all windows
        Transcript** trArrayPointer; //linear array of transcripts to store all of them from all windows            
        
        //read
        uint iRead, iReadAll, iMate, readFilesIndex;
        char readFilter; //Illumina not passed Y/N
        bool revertStrand; //what to do with the strand, according to strandType and iMate
        uint Lread, readLength[MAX_N_MATES], readLengthOriginal[MAX_N_MATES], readLengthPair, readLengthPairOriginal;
        uint clip3pNtotal[MAX_N_MATES], clip5pNtotal[MAX_N_MATES], clip3pAdapterN[MAX_N_MATES]; //total number of trimmed bases from 5p,3p
        int readFileType; //file type: 1=fasta; 2=fastq
        uint outFilterMismatchNmaxTotal;
        
        char dummyChar[4096];
        char** Read0;
        char** Qual0;
        char** readNameMates;
        char* readName;
        char** Read1;
        char** Qual1; //modified QSs for scoring
        
        //split            
        uint** splitR;
        uint Nsplit;
        
//         uint fragLength[MAX_N_FRAG], fragStart[MAX_N_FRAG]; //fragment Lengths and Starts in read space
        
        //binned alignments
        uintWinBin **winBin; //binned genome: window ID (number) per bin
        
        //alignments
        uiPC *PC; //pieces coordinates
        uiWC *WC; //windows coordinates        
        uiWA **WA; //aligments per window
        
        int unmapType; //marker for why a read is unmapped
        
        uint mapMarker; //alignment marker (typically, if there is something wrong)
        uint nA, nP, nW, nWall, nUM[2]; //number of all alignments,  pieces, windows, U/M, 
        uint *nWA, *nWAP, *WALrec, *WlastAnchor; //number of alignments per window, per window per piece, min recordable length per window
        bool *WAincl; //alginment inclusion mask
        
        uint *swWinCov, *swWinGleft, *swWinGright, *swWinRleft, *swWinRright; //read coverage per window
        char *swT;
        
        uint storedLmin, uniqLmax, uniqLmaxInd, multLmax, multLmaxN, multNmin, multNminL, multNmax, multNmaxL;
        uint nTr, nTrMate; // number of transcripts called
        intScore maxScore, nextWinScore;//maximum alignment score, next best score
        
        uint chimN, chimRepeat, chimStr, chimMotif;
        Transcript trChim[MAX_N_CHIMERAS];
        
        Transcript trA, trA1, *trBest, *trNext, *trInit; //transcript, best tr, next best tr, initialized tr
        Transcript ***trAll; //all transcripts for all windows
        uint *nWinTr; //number of recorded transcripts per window
        
        Transcript *alignC, *extendC, *polyAtailC; //alignment rules/conditions
        
        Transcript* trMult[MAX_N_MULTMAP];//multimapping transcripts
        Transcript *alignTrAll;//alignments to transcriptome       
        
        ReadAlign (Parameters* Pin, const Genome &genomeIn, Transcriptome *TrIn, int iChunk);//allocate arrays
        void resetN();//resets the counters to 0
        void multMapSelect();
        int mapOneRead();
        uint maxMappableLength2strands(uint pieceStart, uint pieceLength, uint iDir, uint iSA1, uint iSA2, uint& maxL, uint iFrag);
        void storeAligns (uint iDir, uint Shift, uint Nrep, uint L, uint indStartEnd[2], uint iFrag);
        
        bool outputTranscript(Transcript *trOut, uint nTrOut, ofstream *outBED);
        uint outputTranscriptSAM(Transcript const &trOut, uint nTrOut, uint iTrOut, uint mateChr, uint mateStart, char mateStrand, int unmapType, bool *mateMapped, ostream *outStream);
        int alignBAM(Transcript const &trOut, uint nTrOut, uint iTrOut, uint trChrStart, uint mateChr, uint mateStart, char mateStrand, int unmapType, bool *mateMapped, vector<int> outSAMattrOrder, char** outBAMarray, uint* outBAMarrayN);
        void samAttrNM_MD (Transcript const &trOut, uint iEx1, uint iEx2, uint &tagNM, string &tagMD);
        
        void outputTranscriptSJ(Transcript const &trOut, uint nTrOut, OutSJ *outStream, uint sjReadStartN );
        string outputTranscriptCIGARp(Transcript const &trOut);
        void outTxtMain(ofstream*,Transcript&);
        int createExtendWindowsWithAlign(uint a1, uint aStr); //extends and windows with one alignment
        void assignAlignToWindow(uint a1, uint aLength, uint aStr, uint aNrep, uint aFrag, uint aRstart,bool aAnchor, uint sjA); //assigns one alignment to a window
        void stitchPieces(char **R, char **Q, char *G, PackedArray& SA, uint Lread);
        bool chimericDetection();
        void outputAlignments();
        void stitchWindowSeeds (uint iW, uint iWrec, char* R, char* Q, char* G);//stitches all seeds in one window: iW
        
        int oneRead();
        uint quantTranscriptome (Transcriptome *Tr, uint nAlignG, Transcript **alignG, Transcript *alignT);
};

#endif


