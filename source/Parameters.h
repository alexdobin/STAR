#ifndef PARAMETERS_DEF
#define PARAMETERS_DEF

#include "IncludeDefine.h"
#include "InOutStreams.h"
#include "ParameterInfo.h"
#include <map>
#include "TimeFunctions.h"
#include <unistd.h>
#include <signal.h>

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
        
        //genome, SA, ...
        vector <uint> chrStart, chrLength;
        string genomeDir,genomeLoad;
        vector <string> genomeFastaFiles; 
        uint genomeSAsparseD;//sparsity=distance between indices
        uint genomeInsertL; //total length of the sequence to be inserted on the fly
        uint genomeInsertChrIndFirst; //index of the first inserted chromosome
        
        //binning,windows,anchors
        uint genomeChrBinNbits, genomeChrBinNbases, chrBinN, *chrBin;
        uint winBinChrNbits, winBinNbits, winAnchorDistNbins, winFlankNbins, winBinN;
        uint winAnchorMultimapNmax; //max number of alignments for anchors

        uint genomeSAindexNbases; //length of the SA pre-index strings
        uint *genomeSAindexStart;//starts of the L-mer indices in the SAindex, 1<=L<=genomeSAindexNbases
        
        char genomeNumToNT[6];
        //read parameters
        uint readMapNumber;
        uint iReadAll;
        int readFilesIndex;
        uint32 readFilesN;
        vector <string> readFilesIn, readFilesInTmp;
        vector <vector <string>> readFilesNames;
        uint readNmates;
        string readMatesLengthsIn;
        vector <string> readFilesCommand;
        pid_t readFilesCommandPID[MAX_N_MATES];
        
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
        
        struct 
        {
            string in;
            bool ext[2][2];
        } alignEndsType;
        
        
        //seed parameters
        uint seedMultimapNmax; //max number of multiple alignments per piece          
        uint seedSearchLmax; //max length of the seed
        uint seedPerReadNmax; //max number of pieces per Read
        uint seedPerWindowNmax; //max number of aligns per window
        uint seedNoneLociPerWindow; //max number of aligns from one piece per window
        uint seedSearchStartLmax;
        double seedSearchStartLmaxOverLread; //length of split start points

        //chunk parameters
        uint chunkInSizeBytes,chunkInSizeBytesArray,chunkOutBAMsizeBytes;
        
        //output
        string outFileNamePrefix, outStd, outTmpDir;
        
        //SAM output
        string outBAMfileCoordName, outBAMfileUnsortedName, outQuantBAMfileName;        
        string samHeader, samHeaderHD, samHeaderSortedCoord;
        string outSAMmode, outSAMstrandField,  outSAMunmapped, outSAMorder, outSAMprimaryFlag;
        vector<string> outSAMattributes, outSAMheaderHD, outSAMheaderPG;
        vector<string> outSAMattrRGline,outSAMattrRGlineSplit,outSAMattrRG;
        uint outSAMmultNmax,outSAMattrIHstart;
        string outSAMheaderCommentFile;
        int outSAMmapqUnique;
        struct {bool NH,HI,AS,NM,MD,nM,jM,jI,RG,XS;} outSAMattrPresent, outSAMattrPresentQuant;
        vector <int> outSAMattrOrder, outSAMattrOrderQuant;
        int outBAMcompression;
        vector <string> outSAMtype;
        bool outBAMunsorted, outBAMcoord, outSAMbool;
        uint32 outBAMcoordNbins;
        string outBAMsortTmpDir;
        string bamRemoveDuplicatesType;
        uint bamRemoveDuplicatesMate2basesN;
        int outBAMsortingThreadN, outBAMsortingThreadNactual;
        uint64 *outBAMsortingBinStart; //genomic starts for bins for sorting BAM files
        uint16 outSAMflagOR, outSAMflagAND;
        
        struct
        {
            vector <string> mode;
            bool yes;
            bool KeepOnlyAddedReferences;
        } outSAMfilter;
        
        struct
        {
            string mode;
            bool random;
        } outMultimapperOrder;
        
        string outReadsUnmapped;
        int outQSconversionAdd;
        string outFileTmp;
        
        //output filtering
        uint outFilterMismatchNmax;
        double outFilterMismatchNoverLmax, outFilterMismatchNoverReadLmax; //max proportion of all MM within all bases
        
        uint outFilterMatchNmin,outFilterMultimapNmax;//min number of matches
        double outFilterScoreMinOverLread, outFilterMatchNminOverLread;//normalzied to read length
        intScore outFilterScoreMin,outFilterMultimapScoreRange;//min score to output
        string outFilterIntronMotifs;
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
            string save;
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
    
        //SJ database parameters
        vector <string> sjdbFileChrStartEnd;
        string sjdbGTFfile, sjdbGTFchrPrefix, sjdbGTFfeatureExon, sjdbGTFtagExonParentTranscript, sjdbGTFtagExonParentGene;
        uint sjdbOverhang,sjdbLength; //length of the donor/acceptor, length of the sj "chromosome" =2*sjdbOverhang+1 including spacer
        int sjdbOverhang_par;
        int sjdbScore;        
        
        uint sjChrStart,sjdbN; //first sj-db chr
        uint sjGstart; //start of the sj-db genome sequence
        uint *sjDstart,*sjAstart,*sjStr, *sjdbStart, *sjdbEnd; //sjdb loci
        uint8 *sjdbMotif; //motifs of annotated junctions
        uint8 *sjdbShiftLeft, *sjdbShiftRight; //shifts of junctions
        uint8 *sjdbStrand; //junctions strand, not used yet
        
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
        
        //chimeric
        uint chimSegmentMin, chimJunctionOverhangMin; //min chimeric donor/acceptor length
        uint chimSegmentReadGapMax; //max read gap for stitching chimeric windows
        int chimScoreMin,chimScoreDropMax,chimScoreSeparation, chimScoreJunctionNonGTAG; //min chimeric score
        string chimOutType;
        vector <string> chimFilter;
        
        struct 
        {
            struct
            {
                bool genomicN;
            } filter;
        } chimPar;
        
        //splitting
        char Qsplit;
        uint maxNsplit, minLsplit, minLmap;
        
        //limits
    
        
    ////////////////////// CLEAN-UP needed
    InOutStreams *inOut; //main input output streams

    uint Lread;
            
    //Genome parameters
    uint nGenome, nSA, nSAbyte, nChrReal;//genome length, SA length, # of chromosomes, vector of chromosome start loci
    uint nGenome2, nSA2, nSAbyte2, nChrReal2; //same for the 2nd pass
    uint nSAi; //size of the SAindex
    vector <string> chrName;
    map <string,uint> chrNameIndex;
    unsigned char GstrandBit, SAiMarkNbit, SAiMarkAbsentBit; //SA index bit for strand information
    uint GstrandMask, SAiMarkAbsentMask, SAiMarkAbsentMaskC, SAiMarkNmask, SAiMarkNmaskC;//maske to remove strand bit from SA index, to remove mark from SAi index
    
   

    Parameters();
    void chrInfoLoad(); //find nChr and chrStart from genome
    void chrBinFill();//file chrBin array
    int readParsFromFile(ifstream*, ofstream*, int); //read parameters from one file
    int readPars(); // read parameters from all files
    int scanOneLine (string &lineIn, int inputLevel, int inputLevelRequested);
    void scanAllLines (istream &streamIn, int inputLevel, int inputLevelRequested);
    void inputParameters (int argInN, char* argIn[]); //input parameters: default, from files, from command line    
    void openReadsFiles();
    void closeReadsFiles();

};
#endif  // Parameters.h
