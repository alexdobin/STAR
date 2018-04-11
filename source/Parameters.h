#ifndef PARAMETERS_DEF
#define PARAMETERS_DEF

#include "IncludeDefine.h"
#include "InOutStreams.h"
#include "ParameterInfo.h"
#include <map>
#include "TimeFunctions.h"
#include <unistd.h>
#include <signal.h>
#include "ParametersChimeric.h"
#include "ParametersGenome.h"

class Parameters {

    public:
        vector <ParameterInfoBase*> parArray, parArrayInitial;
        vector <string> parameterInputName;

        string commandLine, commandLineFull;

        //version
        uint versionSTAR;
        vector <uint> versionGenome;

        //system parameters
        string sysShell; //shell for executing system commands

        // run parameters
        string runMode;
        int runThreadN;
        mode_t runDirPerm;
        string runDirPermIn; //permission for directores created at run-time
        int runRNGseed; //random number generator seed

        //parameters
        vector <string> parametersFiles;

        //input
        string inputBAMfile;

        //genome
        char genomeNumToNT[6];
        ParametersGenome pGe;

        //binning,windows,anchors
        uint winBinChrNbits, winBinNbits, winAnchorDistNbins, winFlankNbins, winBinN;
        uint winAnchorMultimapNmax; //max number of alignments for anchors
        double winReadCoverageRelativeMin;
        uint winReadCoverageBasesMin;

        //read parameters
        vector <string> readFilesType;
        int readFilesTypeN;
        string readFilesPrefix;
        vector <string> readFilesIn, readFilesInTmp;
        uint32 readFilesN;
        vector <vector <string> > readFilesNames;
        vector <string> readFilesCommand;
        int readFilesIndex;
        pid_t readFilesCommandPID[MAX_N_MATES];

        uint readMapNumber;
        uint iReadAll;
        uint readNmates;
        string readMatesLengthsIn;

        vector <string> readNameSeparator;
        vector <char> readNameSeparatorChar;

        string outSAMreadID;

        vector <uint> clip5pNbases, clip3pNbases, clip3pAfterAdapterNbases;
        vector <double> clip3pAdapterMMp;
        vector <string> clip3pAdapterSeq;
        char *clip3pAdapterSeqNum[MAX_N_MATES];//adapter sequence - numerical
        bool readMatesEqualLengths; //whether or not the read mates have the same length, true if onyl one mate

        //align parameters
        uint alignSJoverhangMin,alignSJDBoverhangMin,alignSplicedMateMapLmin; //min SJ donor/acceptor length
        double alignSplicedMateMapLminOverLmate;
        uint alignWindowsPerReadNmax; //max number of alignment windows per read
        uint alignTranscriptsPerWindowNmax; //maximum number of transcripts recorded per window
        uint alignTranscriptsPerReadNmax;   //max number of alignments per read
        uint alignIntronMin;//min length to call a gap an intron
        uint alignIntronMax;//max length to call
        uint alignMatesGapMax;//max gap between the mates (if paired-end)
        string alignSoftClipAtReferenceEnds;
        vector <int32> alignSJstitchMismatchNmax;

        struct {
            string in;
            bool ext[2][2];
        } alignEndsType;

        struct {
            vector<string> in;
            int nBasesMax;
            bool concordantPair;
        } alignEndsProtrude;
        
        struct {
            string in;
            bool flushRight;
        } alignInsertionFlush;
        
        
        //seed parameters
        uint seedMultimapNmax; //max number of multiple alignments per piece
        uint seedSearchLmax; //max length of the seed
        uint seedPerReadNmax; //max number of pieces per Read
        uint seedPerWindowNmax; //max number of aligns per window
        uint seedNoneLociPerWindow; //max number of aligns from one piece per window
        uint seedSearchStartLmax;
        double seedSearchStartLmaxOverLread; //length of split start points
        uint seedSplitMin;

        //chunk parameters
        uint chunkInSizeBytes,chunkInSizeBytesArray,chunkOutBAMsizeBytes;

        //output
        string outFileNamePrefix, outStd;
        string outTmpDir, outTmpKeep;

        //SAM output
        string outBAMfileCoordName, outBAMfileUnsortedName, outQuantBAMfileName;
        string samHeader, samHeaderHD, samHeaderSortedCoord, samHeaderExtra;
        string outSAMmode, outSAMstrandField,  outSAMorder, outSAMprimaryFlag;
        vector<string> outSAMattributes, outSAMheaderHD, outSAMheaderPG;
        vector<string> outSAMattrRGline,outSAMattrRGlineSplit,outSAMattrRG;
        uint outSAMmultNmax,outSAMattrIHstart;
        string outSAMheaderCommentFile;
        int outSAMmapqUnique;

        int outSAMtlen;
        
        struct {bool NH,HI,AS,NM,MD,nM,jM,jI,RG,XS,rB,vG,vA,vW,ch,MC;} outSAMattrPresent, outSAMattrPresentQuant;

        vector <int> outSAMattrOrder, outSAMattrOrderQuant;
        int outBAMcompression;
        vector <string> outSAMtype;
        bool outBAMunsorted, outBAMcoord, outSAMbool;
        uint32 outBAMcoordNbins;
        uint32 outBAMsortingBinsN;//user-defined number of bins for sorting
        string outBAMsortTmpDir;
        
//         string bamRemoveDuplicatesType;
//         uint bamRemoveDuplicatesMate2basesN;
        struct
        {
            string mode;
            bool yes;
            bool markMulti;
            uint mate2basesN;
        } removeDuplicates;
        
        int outBAMsortingThreadN, outBAMsortingThreadNactual;
        uint64 *outBAMsortingBinStart; //genomic starts for bins for sorting BAM files
        uint16 outSAMflagOR, outSAMflagAND;

        struct
        {
            vector <string> mode;
            bool yes;
            bool within;//output unmapped reads within SAM/BAM files
            bool keepPairs;//keep mates together
        } outSAMunmapped;

        struct
        {
            vector <string> mode;
            bool yes;
            bool KeepOnlyAddedReferences;
            bool KeepAllAddedReferences;            
        } outSAMfilter;

        struct
        {
            string mode;
            bool random;
        } outMultimapperOrder;
        
        struct
        {
            bool yes;
            uint NbasesMin;
            double MMp;
        } peOverlap;

        string outReadsUnmapped;
        int outQSconversionAdd;
        string outFileTmp;

        //output filtering
        uint outFilterMismatchNmax;
        double outFilterMismatchNoverLmax, outFilterMismatchNoverReadLmax; //max proportion of all MM within all bases

        uint outFilterMatchNmin,outFilterMultimapNmax;//min number of matches
        double outFilterScoreMinOverLread, outFilterMatchNminOverLread;//normalzied to read length
        intScore outFilterScoreMin,outFilterMultimapScoreRange;//min score to output
        string outFilterIntronMotifs,outFilterIntronStrands;
        string outFilterType; //type of filtering
        int outFilterBySJoutStage; //indicates the stage of filtering by SJout

        //output filtering SJs
        string outSJfilterReads;
        vector <int32> outSJfilterCountUniqueMin, outSJfilterCountTotalMin;
        vector <int32> outSJfilterOverhangMin;
        vector <int32> outSJfilterDistToOtherSJmin; //min allowed distance to other SJ's donor/acceptor
        vector <int32> outSJfilterIntronMaxVsReadN;

        //wiggle output
        vector <string> outWigType, outWigStrand, outWigNorm;
        string outWigReferencesPrefix;
        struct {
            bool yes;
            bool strand;
            int type;
            int format;
            int norm;
        } outWigFlags;

        //2-pass
//         uint twoPass.pass1readsN, twoPass.sjLimit;
//         string twoPass.dir,twopassSJpass1file;
        struct {
            bool yes; //true in 2-pass mode
            bool pass2; //true if now running the 2nd pass
            uint pass1readsN;
            int pass1readsN_par;
            string dir;
            string pass1sjFile;
            string mode;
        } twoPass;

        //inserting junctions on the fly
        struct {
            bool yes; //insert?
            bool pass1;//insert on the 1st pass?
            bool pass2;//insert on the 2nd pass?
            string outDir;
        } sjdbInsert;

        //storage limits
        uint limitGenomeGenerateRAM;
        uint limitIObufferSize; //max size of the in/out buffer, bytes
        uint limitOutSAMoneReadBytes;
        uint limitOutSJoneRead, limitOutSJcollapsed;
        uint limitBAMsortRAM;
        uint limitSjdbInsertNsj;

        // penalties
        intScore scoreGap, scoreGapNoncan, scoreGapGCAG, scoreGapATAC, scoreDelBase, scoreDelOpen, scoreInsBase, scoreInsOpen;
        intScore scoreStitchSJshift;//Max negative score when
        double scoreGenomicLengthLog2scale;

        //old variables: CLEAN-up needed
        char outputBED[MAX_OUTPUT_FLAG]; //output flags

        //SW search
        uint swMode, swWinCoverageMinP;
        //SW penalties
        uint swPeoutFilterMatchNmin, swPenMismatch, swPenGapOpen, swPenGapExtend;
        uint swHsize;

        int annotScoreScale;//overall multiplication factor for the annotation
        string annotSignalFile;//binary file with annotation signal

        uint sjNovelN, *sjNovelStart, *sjNovelEnd; //novel junctions collapased and filtered        
        
        
        //quantification parameters
        //input

        struct
        {
          bool yes; //if any quantification is done
          vector <string> mode; //quantification mode input string

          struct
          {
              bool yes;
              bool indel;
              bool softClip;
              bool singleEnd;
              int bamCompression;
              string ban;
          } trSAM;

          struct
          {
              bool yes;
              string outFile;
          } geCount;

        } quant;
       
        //variation parameters
        struct
        {
            bool yes;
            string vcfFile;
        } var;
        
        struct
        {
            bool yes;
            bool SAMtag;
            string outputMode;
        } wasp;

        //chimeric
        ParametersChimeric pCh;


        //splitting
        char Qsplit;
        uint maxNsplit, minLsplit, minLmap;

        //limits


    ////////////////////// CLEAN-UP needed
    InOutStreams *inOut; //main input output streams

    uint Lread;

    Parameters();
    int readParsFromFile(ifstream*, ofstream*, int); //read parameters from one file
    int readPars(); // read parameters from all files
    int scanOneLine (string &lineIn, int inputLevel, int inputLevelRequested);
    void scanAllLines (istream &streamIn, int inputLevel, int inputLevelRequested);
    void inputParameters (int argInN, char* argIn[]); //input parameters: default, from files, from command line
    void openReadsFiles();
    void closeReadsFiles();
    void readSAMheader(const string readFilesCommandString, const vector<string> readFilesNames);

};
#endif  // Parameters.h
