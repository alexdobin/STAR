#ifndef CODE_READALIGN
#define CODE_READALIGN

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "Genome.h"
#include "Stats.h"
#include "OutSJ.h"
#include "Transcriptome.h"
#include "BAMoutput.h"
#include "Quantifications.h"
#include "ChimericDetection.h"

#include <time.h>
#include <random>

class ReadAlign
{
    public:
         //methods
        ReadAlign (Parameters& Pin, Genome &genomeIn, Transcriptome *TrIn, int iChunk);//allocate arrays
        int oneRead();
        
        //vars

        Genome &mapGen; //mapped-to-genome structure

        uint iRead;
        char** Read1;
        
        Stats statsRA; //mapping statistics        
        
        istream* readInStream[MAX_N_MATES];
        BAMoutput *outBAMcoord, *outBAMunsorted, *outBAMquant;//sorted by coordinate, unsorted, transcriptomic BAM structure
        fstream chunkOutChimSAM, *chunkOutChimJunction, chunkOutUnmappedReadsStream[MAX_N_MATES], chunkOutFilterBySJoutFiles[MAX_N_MATES];
        OutSJ *chunkOutSJ, *chunkOutSJ1;

        ostream* outSAMstream;        
        uint outBAMbytes; //number of bytes output to SAM/BAM with oneRead
        char *outBAMarray;//pointer to the (last+1) position of the SAM/BAM output array
        
        uint outFilterMismatchNmaxTotal;
        uint Lread, readLength[MAX_N_MATES], readLengthOriginal[MAX_N_MATES], readLengthPair, readLengthPairOriginal;
        intScore maxScoreMate[MAX_N_MATES];

        uint readFilesIndex;
        
        ReadAlign *waspRA; //ReadAlign for alternative WASP alignment
        int waspType, waspType1; //alignment ASE-WASP type and

        ReadAlign *peMergeRA; //ReadAlign for merged PE mates
        
        ChimericDetection *chimDet;
        
    private:
        Parameters& P; //pointer to the parameters, will be initialized on construction

        //quantification
        Transcriptome *chunkTr;

        //mapping time
        time_t timeStart, timeFinish;

        //random number generators
        std::mt19937 rngMultOrder;//initialize in ReadAlign.cpp
        std::uniform_real_distribution<double> rngUniformReal0to1;//initialize in ReadAlign.cpp

        //input,output

        char** outBAMoneAlign;
        uint* outBAMoneAlignNbytes;

        ostringstream samStreamCIGAR, samStreamSJmotif, samStreamSJintron;
        vector <string> matesCIGAR;

        intScore *scoreSeedToSeed, *scoreSeedBest;
        uint *scoreSeedBestInd, *seedChain, *scoreSeedBestMM;

        bool outFilterPassed; //true if alignment passed all filter and is output to SAM/BAM

//         StatsAll *statsRA;

        //transcript
        Transcript* trArray; //linear array of transcripts to store all of them from all windows
        Transcript** trArrayPointer; //linear array of transcripts to store all of them from all windows

        //read
        uint iReadAll, iMate;
        char readFilter; //Illumina not passed Y/N
        bool revertStrand; //what to do with the strand, according to strandType and iMate
        uint clip3pNtotal[MAX_N_MATES], clip5pNtotal[MAX_N_MATES], clip3pAdapterN[MAX_N_MATES]; //total number of trimmed bases from 5p,3p
        int readFileType; //file type: 1=fasta; 2=fastq

        vector<string>readNameExtra;

        char dummyChar[4096];
        char** Read0;
        char** Qual0;
        char** readNameMates;
        char* readName;
        char** Qual1; //modified QSs for scoring

        uint readNmates;
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
        intScore maxScore;//maximum alignment score
        
        Transcript trA, trA1, *trBest, *trInit; //transcript, best tr, next best tr, initialized tr
        Transcript ***trAll; //all transcripts for all windows
        uint *nWinTr; //number of recorded transcripts per window

        //old chimeric detection
        uint chimN, chimRepeat, chimStr;
        int chimMotif;
        uint chimRepeat0, chimRepeat1, chimJ0, chimJ1;
        Transcript trChim[MAX_N_CHIMERAS];
        //new chimeric detection
        bool chimRecord; //true if chimeric aligment was detected
        
        Transcript *alignC, *extendC, *polyAtailC; //alignment rules/conditions

        Transcript* trMult[MAX_N_MULTMAP];//multimapping transcripts
        Transcript *alignTrAll;//alignments to transcriptome
        
        struct {
            bool yes;
            uint nOv;//number of overlapping bases 
            uint ovS;//first read base of the overlap
            uint mateStart[2];//mates starts in the merged read
        } peOv;//PE  mates overlap/merge/remap structure
        
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
        int createExtendWindowsWithAlign(uint a1, uint aStr); //extends and windows with one alignment
        void assignAlignToWindow(uint a1, uint aLength, uint aStr, uint aNrep, uint aFrag, uint aRstart,bool aAnchor, uint sjA); //assigns one alignment to a window
        
        void mappedFilter();
        void chimericDetection();
        bool chimericDetectionOld();
        void chimericDetectionOldOutput();
        bool chimericDetectionMult();        
        void chimericDetectionPEmerged(ReadAlign &seRa);
//         void chimericDetectionPEmergedTrim();
                
        void outputAlignments();
        void calcCIGAR(Transcript const &trOut, uint nMates, uint iExMate, uint leftMate);

        void stitchWindowSeeds (uint iW, uint iWrec, bool *WAexcl, char *R);//stitches all seeds in one window: iW
        void stitchPieces(char **R, uint Lread);

        uint quantTranscriptome (Transcriptome *Tr, uint nAlignG, Transcript **alignG, Transcript *alignT);
        
        void copyRead(ReadAlign&);
        void waspMap();
        void peOverlapMergeMap();
        void peMergeMates();
        void peOverlapSEtoPE(ReadAlign &seRA);

        
};

#endif


