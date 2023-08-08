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
#include "ParametersSolo.h"
#include "ParametersClip.h"
#include "ParametersGenome.h"
#include <vector>
#include <array>
#include <unordered_set>

class Parameters {

    public:
        vector <ParameterInfoBase*> parArray, parArrayInitial;
        vector <string> parameterInputName;

        string commandLine, commandLineFull;

        //version
        string versionGenome;

        //system parameters
        string sysShell; //shell for executing system commands

        // run parameters
        string runMode;
        vector<string> runModeIn;
        int runThreadN;
        mode_t runDirPerm;
        string runDirPermIn; //permission for directores created at run-time
        int runRNGseed; //random number generator seed

        struct {
            int32 type;//0 no restart, 1 no mapping - restart from _STARtmp files
        } runRestart; //restart options - in development
        
        //parameters
        vector <string> parametersFiles;

        //input
        string inputBAMfile;

        //genome
        char genomeNumToNT[6];
        ParametersGenome pGe, pGeOut;

        //binning,windows,anchors
        uint winBinChrNbits, winBinNbits, winAnchorDistNbins, winFlankNbins, winBinN;
        uint winAnchorMultimapNmax; //max number of alignments for anchors
        double winReadCoverageRelativeMin;
        uint winReadCoverageBasesMin;

        //read parameters
        vector <string> readFilesType;
        int readFilesTypeN;
        string readFilesPrefix, readFilesPrefixFinal;
        vector <string> readFilesIn, readFilesInTmp;
        uint32 readFilesN;
        vector <vector <string> > readFilesNames;
        vector <string> readFilesCommand;
        vector <string> readFilesManifest;
               
        string readFilesCommandString; //actual command string
        int readFilesIndex;
        pid_t readFilesCommandPID[MAX_N_MATES];

        uint readMapNumber;
        uint iReadAll;
        uint readNmates, readNends;
        string readMatesLengthsIn;
        uint32 readQualityScoreBase;

        vector <string> readNameSeparator;
        vector <char> readNameSeparatorChar;

        string outSAMreadID;
        bool outSAMreadIDnumber;
        
        //new: structure for readFiles parameters
        struct {
            vector<string> samAttrKeepIn; //input vector of SAM tags to keep, if readFilesType=SAMtag
            std::unordered_set<uint16_t> samAttrKeep;
            bool samAttrKeepAll, samAttrKeepNone;
        } readFiles;

        ParametersClip pClip;

        //align parameters
        uint alignSJoverhangMin,alignSJDBoverhangMin,alignSplicedMateMapLmin; //min SJ donor/acceptor length
        double alignSplicedMateMapLminOverLmate;
        uint alignWindowsPerReadNmax; //max number of alignment windows per read
        uint alignTranscriptsPerWindowNmax; //maximum number of transcripts recorded per window
        uint alignTranscriptsPerReadNmax;   //max number of alignments per read
        uint alignIntronMin;//min length to call a gap an intron
        uint alignIntronMax;//max length to call
        uint alignMatesGapMax;//max gap between the mates (if paired-end)
        vector <int32> alignSJstitchMismatchNmax;

        //         struct {
        //             string strandString;
        //             int32 strand;
        //         } pReads;

        struct {
            string in;
            bool yes;
        } alignSoftClipAtReferenceEnds;

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
        uint64 seedSplitMin, seedMapMin;

        //chunk parameters
        uint chunkInSizeBytes,chunkInSizeBytesArray,chunkOutBAMsizeBytes;

        //output
        string outFileNamePrefix, outStd;
        string outTmpDir, outTmpKeep;
        string outLogFileName;

        //SAM output
        string outBAMfileCoordName, outBAMfileUnsortedName, outQuantBAMfileName;
        string samHeader, samHeaderHD, samHeaderSortedCoord, samHeaderExtra;
        string outSAMmode,  outSAMorder, outSAMprimaryFlag;
        vector<string> outSAMattributes, outSAMheaderHD, outSAMheaderPG;
        vector<string> outSAMattrRGline,outSAMattrRGlineSplit,outSAMattrRG;
        uint outSAMmultNmax,outSAMattrIHstart;
        string outSAMheaderCommentFile;
        int outSAMmapqUnique;

        struct {
            string in;
            uint32 type;
        } outSAMstrandField;

        int outSAMtlen;

        struct {bool NH,HI,AS,NM,MD,nM,jM,jI,RG,XS,rB,vG,vA,vW,ha,ch,MC,CR,CY,UR,UY,CB,UB,GX,GN,gx,gn,sM,sS,sQ,cN,sF;} outSAMattrPresent, outSAMattrPresentQuant;

        vector <int> outSAMattrOrder, outSAMattrOrderQuant;
        int outBAMcompression;
        vector <string> outSAMtype;
        bool outBAMunsorted, outBAMcoord, outSAMbool;
        uint32 outBAMcoordNbins;
        uint32 outBAMsortingBinsN;//user-defined number of bins for sorting
        string outBAMsortTmpDir;

//         string bamRemoveDuplicatesType;
//         uint bamRemoveDuplicatesMate2basesN;
        struct {
            string mode;
            bool yes;
            bool markMulti;
            uint mate2basesN;
        } removeDuplicates;

        int outBAMsortingThreadN, outBAMsortingThreadNactual;
        uint64 *outBAMsortingBinStart; //genomic starts for bins for sorting BAM files
        uint16 outSAMflagOR, outSAMflagAND;

        struct {
            vector <string> mode;
            bool yes;
            bool within;//output unmapped reads within SAM/BAM files
            bool keepPairs;//keep mates together
        } outSAMunmapped;

        struct {
            vector <string> mode;
            bool yes;
            bool KeepOnlyAddedReferences;
            bool KeepAllAddedReferences;
        } outSAMfilter;

        struct {
            string mode;
            bool random;
        } outMultimapperOrder;

        struct {
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

        struct {
            vector<string> type;
            bool yes;
        } outSJ;
        
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
        uint64 limitGenomeGenerateRAM;
        vector<uint64> limitIObufferSize; //max size of the in/out buffer, bytes
        uint64 limitOutSAMoneReadBytes;
        uint64 limitOutSJoneRead, limitOutSJcollapsed;
        uint64 limitBAMsortRAM;
        uint64 limitSjdbInsertNsj;
        uint64 limitNreadsSoft;

        // penalties
        intScore scoreGap, scoreGapNoncan, scoreGapGCAG, scoreGapATAC, scoreDelBase, scoreDelOpen, scoreInsBase, scoreInsOpen;
        intScore scoreStitchSJshift;//Max negative score when
        double scoreGenomicLengthLog2scale;

        //quantification parameters
        //input

        struct {
          bool yes=false; //if any quantification is done
          vector <string> mode; //quantification mode input string

            struct {
                bool yes=false;
                bool bamYes;
                bool indel;
                bool softClip;
                bool singleEnd;
                int bamCompression;
                string ban;
            } trSAM;

            struct {
                bool yes=false;
                string outFile;
            } geCount;

            struct {
                bool yes=false;
            } geneFull;
          
            struct {
                bool yes=false;
            } geneFull_Ex50pAS;

            struct {
                bool yes=false;
            } geneFull_ExonOverIntron;

            struct {
                bool yes=false;
            } gene;
          
        } quant;

        //variation parameters
        struct {
            bool yes=false;
            bool heteroOnly=false;
            string vcfFile;
        } var;

        struct {
            bool yes=false;
            bool SAMtag;
            string outputMode;
        } wasp;

        //solo
        ParametersSolo pSolo;

        //chimeric
        ParametersChimeric pCh;

        //splitting
        uint maxNsplit;

        //not really parameters, but global variables:
        array<vector<uint64>,2> sjAll;
        uint64 sjNovelN, *sjNovelStart, *sjNovelEnd; //novel junctions collapased and filtered

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
    void readFilesInit();
    void closeReadsFiles();
    void readSAMheader(const string readFilesCommandString, const vector<string> readFilesNames);
    void samAttributes();
    void samAttrRequiresBAM(bool attrYes, string attrTag);
};
#endif  // Parameters.h
