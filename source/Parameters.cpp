#include "IncludeDefine.h"
#include "Parameters.h"
#include "ErrorWarning.h"
#include "SequenceFuns.h"
#include "OutSJ.h"
#include "sysRemoveDir.h"
#include "stringSubstituteAll.h"
#include SAMTOOLS_BGZF_H
#include "GlobalVariables.h"
#include "signalFromBAM.h"
#include "bamRemoveDuplicates.h"

//for mkfifo
#include <sys/stat.h>

#define PAR_NAME_PRINT_WIDTH 30

Parameters::Parameters() {//initalize parameters info

    inOut = new InOutStreams;

    //versions
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "versionSTAR", &versionSTAR));
    parArray.push_back(new ParameterInfoVector <uint> (-1, -1, "versionGenome", &versionGenome));

    //parameters
    parArray.push_back(new ParameterInfoVector <string> (-1, 2, "parametersFiles", &parametersFiles));

    //system
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sysShell", &sysShell));

    //run
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "runMode", &runMode));
    parArray.push_back(new ParameterInfoScalar <int> (-1, -1, "runThreadN", &runThreadN));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "runDirPerm", &runDirPermIn));
    parArray.push_back(new ParameterInfoScalar <int> (-1, -1, "runRNGseed", &runRNGseed));

    //genome
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "genomeDir", &pGe.gDir));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "genomeLoad", &pGe.gLoad));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "genomeFastaFiles", &pGe.gFastaFiles));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "genomeChainFiles", &pGe.gChainFiles));
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "genomeSAindexNbases", &pGe.gSAindexNbases));
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "genomeChrBinNbits", &pGe.gChrBinNbits));
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "genomeSAsparseD", &pGe.gSAsparseD));
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "genomeSuffixLengthMax", &pGe.gSuffixLengthMax));
    parArray.push_back(new ParameterInfoVector <uint> (-1, -1, "genomeFileSizes", &pGe.gFileSizes));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "genomeConsensusFile", &pGe.gConsensusFile));

    //read
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "readFilesType", &readFilesType));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "readFilesIn", &readFilesIn));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "readFilesPrefix", &readFilesPrefix));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "readFilesCommand", &readFilesCommand));
    
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "readMatesLengthsIn", &readMatesLengthsIn));
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "readMapNumber", &readMapNumber));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "readNameSeparator", &readNameSeparator));


    //input from BAM
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "inputBAMfile", &inputBAMfile));

    //BAM processing
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "bamRemoveDuplicatesType", &removeDuplicates.mode));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "bamRemoveDuplicatesMate2basesN", &removeDuplicates.mate2basesN));

    //limits
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitGenomeGenerateRAM", &limitGenomeGenerateRAM));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitIObufferSize", &limitIObufferSize));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitOutSAMoneReadBytes", &limitOutSAMoneReadBytes));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitOutSJcollapsed", &limitOutSJcollapsed));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitOutSJoneRead", &limitOutSJoneRead));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitBAMsortRAM", &limitBAMsortRAM));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitSjdbInsertNsj", &limitSjdbInsertNsj));


    //output
    parArray.push_back(new ParameterInfoScalar <string>     (-1, 2, "outFileNamePrefix", &outFileNamePrefix));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, 2, "outTmpDir", &outTmpDir));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, 2, "outTmpKeep", &outTmpKeep));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, 2, "outStd", &outStd));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outReadsUnmapped", &outReadsUnmapped));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "outQSconversionAdd", &outQSconversionAdd));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outMultimapperOrder", &outMultimapperOrder.mode));

    //outSAM
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "outSAMtype", &outSAMtype));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMmode", &outSAMmode));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMstrandField", &outSAMstrandField));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "outSAMattributes", &outSAMattributes));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "outSAMunmapped", &outSAMunmapped.mode));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMorder", &outSAMorder));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMprimaryFlag", &outSAMprimaryFlag));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMreadID", &outSAMreadID));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "outSAMmapqUnique", &outSAMmapqUnique));
    parArray.push_back(new ParameterInfoScalar <uint16>        (-1, -1, "outSAMflagOR", &outSAMflagOR));
    parArray.push_back(new ParameterInfoScalar <uint16>        (-1, -1, "outSAMflagAND", &outSAMflagAND));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "outSAMattrRGline", &outSAMattrRGline));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "outSAMheaderHD", &outSAMheaderHD));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "outSAMheaderPG", &outSAMheaderPG));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMheaderCommentFile", &outSAMheaderCommentFile));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "outBAMcompression", &outBAMcompression));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "outBAMsortingThreadN", &outBAMsortingThreadN));
    parArray.push_back(new ParameterInfoScalar <uint32>        (-1, -1, "outBAMsortingBinsN", &outBAMsortingBinsN));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "outSAMfilter", &outSAMfilter.mode));
    parArray.push_back(new ParameterInfoScalar <uint>     (-1, -1, "outSAMmultNmax", &outSAMmultNmax));
    parArray.push_back(new ParameterInfoScalar <uint>     (-1, -1, "outSAMattrIHstart", &outSAMattrIHstart));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "outSAMtlen", &outSAMtlen));


   //output SJ filtering
    parArray.push_back(new ParameterInfoScalar <string>  (-1, -1, "outSJfilterReads", &outSJfilterReads));
    parArray.push_back(new ParameterInfoVector <int32>   (-1, -1, "outSJfilterCountUniqueMin", &outSJfilterCountUniqueMin));
    parArray.push_back(new ParameterInfoVector <int32>   (-1, -1, "outSJfilterCountTotalMin", &outSJfilterCountTotalMin));
    parArray.push_back(new ParameterInfoVector <int32>   (-1, -1, "outSJfilterOverhangMin", &outSJfilterOverhangMin));
    parArray.push_back(new ParameterInfoVector <int32>   (-1, -1, "outSJfilterDistToOtherSJmin", &outSJfilterDistToOtherSJmin));
    parArray.push_back(new ParameterInfoVector <int32>   (-1, -1, "outSJfilterIntronMaxVsReadN", &outSJfilterIntronMaxVsReadN));

    //output wiggle
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "outWigType", &outWigType));
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "outWigStrand", &outWigStrand));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "outWigReferencesPrefix", &outWigReferencesPrefix));
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "outWigNorm", &outWigNorm));

    //output filtering
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "outFilterType", &outFilterType) );

    parArray.push_back(new ParameterInfoScalar <uint>     (-1, -1, "outFilterMultimapNmax", &outFilterMultimapNmax));
    parArray.push_back(new ParameterInfoScalar <intScore> (-1, -1, "outFilterMultimapScoreRange", &outFilterMultimapScoreRange));

    parArray.push_back(new ParameterInfoScalar <intScore> (-1, -1, "outFilterScoreMin", &outFilterScoreMin));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "outFilterScoreMinOverLread", &outFilterScoreMinOverLread));

    parArray.push_back(new ParameterInfoScalar <uint>     (-1, -1, "outFilterMatchNmin", &outFilterMatchNmin));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "outFilterMatchNminOverLread", &outFilterMatchNminOverLread));

    parArray.push_back(new ParameterInfoScalar <uint>     (-1, -1, "outFilterMismatchNmax", &outFilterMismatchNmax));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "outFilterMismatchNoverLmax", &outFilterMismatchNoverLmax));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "outFilterMismatchNoverReadLmax", &outFilterMismatchNoverReadLmax));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "outFilterIntronMotifs", &outFilterIntronMotifs));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "outFilterIntronStrands", &outFilterIntronStrands));

    //clipping
    parArray.push_back(new ParameterInfoVector <uint>   (-1, -1, "clip5pNbases", &clip5pNbases));
    parArray.push_back(new ParameterInfoVector <uint>   (-1, -1, "clip3pNbases", &clip3pNbases));
    parArray.push_back(new ParameterInfoVector <uint>   (-1, -1, "clip3pAfterAdapterNbases", &clip3pAfterAdapterNbases));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "clip3pAdapterSeq", &clip3pAdapterSeq));
    parArray.push_back(new ParameterInfoVector <double> (-1, -1, "clip3pAdapterMMp", &clip3pAdapterMMp));

    //binning, anchors, windows
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "winBinNbits", &winBinNbits));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "winAnchorDistNbins", &winAnchorDistNbins));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "winFlankNbins", &winFlankNbins));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "winAnchorMultimapNmax", &winAnchorMultimapNmax));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "winReadCoverageRelativeMin", &winReadCoverageRelativeMin));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "winReadCoverageBasesMin", &winReadCoverageBasesMin));

    //scoring
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreGap", &scoreGap));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreGapNoncan", &scoreGapNoncan));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreGapGCAG", &scoreGapGCAG));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreGapATAC", &scoreGapATAC));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreStitchSJshift", &scoreStitchSJshift));
    parArray.push_back(new ParameterInfoScalar <double>     (-1, -1, "scoreGenomicLengthLog2scale", &scoreGenomicLengthLog2scale));

    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreDelBase", &scoreDelBase));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreDelOpen", &scoreDelOpen));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreInsOpen", &scoreInsOpen));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreInsBase", &scoreInsBase));

    //alignment
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedSearchLmax", &seedSearchLmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedSearchStartLmax", &seedSearchStartLmax));
    parArray.push_back(new ParameterInfoScalar <double>     (-1, -1, "seedSearchStartLmaxOverLread", &seedSearchStartLmaxOverLread));

    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedPerReadNmax", &seedPerReadNmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedPerWindowNmax", &seedPerWindowNmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedNoneLociPerWindow", &seedNoneLociPerWindow));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedMultimapNmax", &seedMultimapNmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedSplitMin", &seedSplitMin));

    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignIntronMin", &alignIntronMin));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignIntronMax", &alignIntronMax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignMatesGapMax", &alignMatesGapMax));

    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignTranscriptsPerReadNmax", &alignTranscriptsPerReadNmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignSJoverhangMin", &alignSJoverhangMin));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignSJDBoverhangMin", &alignSJDBoverhangMin));
    parArray.push_back(new ParameterInfoVector <int32>      (-1, -1, "alignSJstitchMismatchNmax", &alignSJstitchMismatchNmax));

    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignSplicedMateMapLmin", &alignSplicedMateMapLmin));
    parArray.push_back(new ParameterInfoScalar <double>     (-1, -1, "alignSplicedMateMapLminOverLmate", &alignSplicedMateMapLminOverLmate));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignWindowsPerReadNmax", &alignWindowsPerReadNmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignTranscriptsPerWindowNmax", &alignTranscriptsPerWindowNmax));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "alignEndsType", &alignEndsType.in));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "alignSoftClipAtReferenceEnds", &alignSoftClipAtReferenceEnds));

    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "alignEndsProtrude", &alignEndsProtrude.in));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "alignInsertionFlush", &alignInsertionFlush.in));

    //peOverlap
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "peOverlapNbasesMin", &peOverlap.NbasesMin));    
    parArray.push_back(new ParameterInfoScalar <double>     (-1, -1, "peOverlapMMp", &peOverlap.MMp));
    
    //chimeric
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "chimSegmentMin", &pCh.segmentMin));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "chimScoreMin", &pCh.scoreMin));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "chimScoreDropMax", &pCh.scoreDropMax));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "chimScoreSeparation", &pCh.scoreSeparation));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "chimScoreJunctionNonGTAG", &pCh.scoreJunctionNonGTAG));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "chimMainSegmentMultNmax", &pCh.mainSegmentMultNmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "chimJunctionOverhangMin", &pCh.junctionOverhangMin));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "chimOutType", &pCh.out.type));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "chimFilter", &pCh.filter.stringIn));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "chimSegmentReadGapMax", &pCh.segmentReadGapMax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "chimMultimapNmax", &pCh.multimapNmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "chimMultimapScoreRange", &pCh.multimapScoreRange));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "chimNonchimScoreDropMin", &pCh.nonchimScoreDropMin));

    //sjdb
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "sjdbFileChrStartEnd", &pGe.sjdbFileChrStartEnd));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbGTFfile", &pGe.sjdbGTFfile));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbGTFchrPrefix", &pGe.sjdbGTFchrPrefix));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbGTFfeatureExon", &pGe.sjdbGTFfeatureExon));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbGTFtagExonParentTranscript", &pGe.sjdbGTFtagExonParentTranscript));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbGTFtagExonParentGene", &pGe.sjdbGTFtagExonParentGene));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "sjdbOverhang", &pGe.sjdbOverhang));
    pGe.sjdbOverhang_par=parArray.size()-1;
    parArray.push_back(new ParameterInfoScalar <int>    (-1, -1, "sjdbScore", &pGe.sjdbScore));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbInsertSave", &pGe.sjdbInsertSave)); 
    
    //variation
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "varVCFfile", &var.vcfFile)); 

    //WASP
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "waspOutputMode", &wasp.outputMode)); 

    
    //quant
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "quantMode", &quant.mode));
    parArray.push_back(new ParameterInfoScalar <int>     (-1, -1, "quantTranscriptomeBAMcompression", &quant.trSAM.bamCompression));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "quantTranscriptomeBan", &quant.trSAM.ban));

    //2-pass
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "twopass1readsN", &twoPass.pass1readsN));
    twoPass.pass1readsN_par=parArray.size()-1;
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "twopassMode", &twoPass.mode));


//     //SW parameters
//     parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "swMode", &swMode));
//     parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "swWinCoverageMinP", &swWinCoverageMinP));

    parameterInputName.push_back("Default");
    parameterInputName.push_back("Command-Line-Initial");
    parameterInputName.push_back("Command-Line");
    parameterInputName.push_back("genomeParameters.txt");

};


void Parameters::inputParameters (int argInN, char* argIn[]) {//input parameters: default, from files, from command line

///////// Default parameters

    #include "parametersDefault.xxd"
    string parString( (const char*) parametersDefault,parametersDefault_len);
    stringstream parStream (parString);

    scanAllLines(parStream, 0, -1);
    for (uint ii=0; ii<parArray.size(); ii++) {
        if (parArray[ii]->inputLevel<0) {
            ostringstream errOut;
            errOut <<"BUG: DEFAULT parameter value not defined: "<<parArray[ii]->nameString;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

///////// Initial parameters from Command Line

    commandLine="";
    string commandLineFile="";

    if (argInN>1) {//scan parameters from command line
        commandLine += string(argIn[0]);
        for (int iarg=1; iarg<argInN; iarg++) {
            string oneArg=string(argIn[iarg]);

            if (oneArg=="--version") {//print version and exit
                std::cout << STAR_VERSION <<std::endl;
                exit(0);
            };

            size_t found = oneArg.find("=");
            if (found!=string::npos && oneArg.substr(0,2)=="--") {// --parameter=value
                string key = oneArg.substr(2, found - 2);
                string val = oneArg.substr(found + 1);
                if (val.find_first_of(" \t")!=std::string::npos) {//there is white space in the argument, put "" around
                    val ='\"' + val + '\"';
                };
                commandLineFile += '\n' + key + ' ' + val;
            }
            else if (oneArg.substr(0,2)=="--") {//parameter name, cut --
                commandLineFile +='\n' + oneArg.substr(2);
            } else {//parameter value
                if (oneArg.find_first_of(" \t")!=std::string::npos) {//there is white space in the argument, put "" around
                    oneArg ='\"'  + oneArg +'\"';
                };
                commandLineFile +=' ' + oneArg;
            };
            commandLine += ' ' + oneArg;
        };
        istringstream parStreamCommandLine(commandLineFile);
        scanAllLines(parStreamCommandLine, 1, 2); //read only initial Command Line parameters
    };

//     need to be careful since runMode and pGe.gDir are not Command-Line-Initial
//     if (runMode=="genomeGenerate" && outFileNamePrefix=="./") {// for genome generation, output into pGe.gDir
//         outFileNamePrefix=pGe.gDir;
//     };

    inOut->logMain.open((outFileNamePrefix + "Log.out").c_str());
    if (inOut->logMain.fail()) {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL ERROR: could not create output file: "<<outFileNamePrefix + "Log.out"<<"\n";
        errOut <<"Check if the path " << outFileNamePrefix << " exists and you have permissions to write there\n";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    inOut->logMain << "STAR version=" << STAR_VERSION << "\n"<<flush;
    inOut->logMain << "STAR compilation time,server,dir=" << COMPILATION_TIME_PLACE << "\n"<<flush;



    //define what goes to cout
    if (outStd=="Log") {
        inOut->logStdOut=& std::cout;
    } else if (outStd=="SAM" || outStd=="BAM_Unsorted" || outStd=="BAM_SortedByCoordinate" || outStd=="BAM_Quant") {
        inOut->logStdOutFile.open((outFileNamePrefix + "Log.std.out").c_str());
        inOut->logStdOut= & inOut->logStdOutFile;
    } else {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL PARAMETER error: outStd="<<outStd <<" is not a valid value of the parameter\n";
        errOut <<"SOLUTION: provide a valid value fot outStd: Log / SAM / BAM_Unsorted / BAM_SortedByCoordinate";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    inOut->logMain << "##### DEFAULT parameters:\n" <<flush;
    for (uint ii=0; ii<parArray.size(); ii++) {
        if (parArray[ii]->inputLevel==0) {
            inOut->logMain << setw(PAR_NAME_PRINT_WIDTH) << parArray[ii]->nameString <<"    "<< *(parArray[ii]) << endl;
        };
    };

    inOut->logMain <<"##### Command Line:\n"<<commandLine <<endl ;

    inOut->logMain << "##### Initial USER parameters from Command Line:\n";
    for (uint ii=0; ii<parArray.size(); ii++) {
        if (parArray[ii]->inputLevel==1) {
            inOut->logMain << setw(PAR_NAME_PRINT_WIDTH) << parArray[ii]->nameString <<"    "<< *(parArray[ii]) << endl;
        };
    };

///////// Parameters files

    if (parametersFiles.at(0) != "-") {//read parameters from a user-defined file
        for (uint ii=0; ii<parametersFiles.size(); ii++) {
            parameterInputName.push_back(parametersFiles.at(ii));
            ifstream parFile(parametersFiles.at(ii).c_str());
            if (parFile.good()) {
                inOut->logMain << "##### USER parameters from user-defined parameters file " <<parametersFiles.at(ii)<< ":\n" <<flush;
                scanAllLines(parFile, parameterInputName.size()-1, -1);
                parFile.close();
            } else {
                ostringstream errOut;
                errOut <<"EXITING because of fatal input ERROR: could not open user-defined parameters file " <<parametersFiles.at(ii)<< "\n" <<flush;
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            };
        };
    };

///////// Command Line Final

    if (argInN>1) {//scan all parameters from command line and override previuos values
        inOut->logMain << "###### All USER parameters from Command Line:\n" <<flush;
        istringstream parStreamCommandLine(commandLineFile);
        scanAllLines(parStreamCommandLine, 2, -1);
    };

    inOut->logMain << "##### Finished reading parameters from all sources\n\n" << flush;

    inOut->logMain << "##### Final user re-defined parameters-----------------:\n" << flush;

    ostringstream clFull;
    clFull << argIn[0];
    for (uint ii=0; ii<parArray.size(); ii++) {
        if (parArray[ii]->inputLevel>0) {
            inOut->logMain << setw(PAR_NAME_PRINT_WIDTH) << parArray[ii]->nameString <<"    "<< *(parArray[ii]) << endl;
            if (parArray[ii]->nameString != "parametersFiles" ) {
                clFull << "   --" << parArray[ii]->nameString << " " << *(parArray[ii]);
            };
        };
    };
    commandLineFull=clFull.str();
    inOut->logMain << "\n-------------------------------\n##### Final effective command line:\n" <<  clFull.str() << "\n";

//     parOut.close();
    inOut->logMain << "\n##### Final parameters after user input--------------------------------:\n" << flush;
//     parOut.open("Parameters.all.out");
    for (uint ii=0; ii<parArray.size(); ii++) {
        inOut->logMain << setw(PAR_NAME_PRINT_WIDTH) << parArray[ii]->nameString <<"    "<< *(parArray[ii]) << endl;
    };
//     parOut.close();

    inOut->logMain << "----------------------------------------\n\n" << flush;


    ///////////////////////////////////////// Old variables
    outputBED[0]=0; outputBED[1]=0; //controls the format of starBED output. Will be replaced with HDF5 output

    //annot scoring
    annotScoreScale=0;
    annotSignalFile="-";

    //splitting
    Qsplit=0;
    maxNsplit=10;
    minLsplit=seedSplitMin;
    minLmap=5;



////////////////////////////////////////////////////// Calculate and check parameters
    iReadAll=0;

    if (runDirPermIn=="User_RWX")
    {
        runDirPerm=S_IRWXU;
    } else if (runDirPermIn=="All_RWX")
    {
//         umask(0); //this should be defined by the user!
        runDirPerm= S_IRWXU | S_IRWXG | S_IRWXO;
    } else
    {
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT ERROR: unrecognized option in --runDirPerm=" << runDirPerm << "\n";
        errOut << "SOLUTION: use one of the allowed values of --runDirPerm : 'User_RWX' or 'All_RWX' \n";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (outTmpDir=="-") {
        outFileTmp=outFileNamePrefix +"_STARtmp/";
        sysRemoveDir (outFileTmp);
    } else {
        outFileTmp=outTmpDir;
    };

    if (mkdir (outFileTmp.c_str(),runDirPerm)!=0) {
        ostringstream errOut;
        errOut <<"EXITING because of fatal ERROR: could not make temporary directory: "<< outFileTmp<<"\n";
        errOut <<"SOLUTION: (i) please check the path and writing permissions \n (ii) if you specified --outTmpDir, and this directory exists - please remove it before running STAR\n"<<flush;
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    //threaded or not
    g_threadChunks.threadBool=(runThreadN>1);

    //wigOut parameters
    if (outWigType.at(0)=="None") {
        outWigFlags.yes=false;
    } else if (outWigType.at(0)=="bedGraph") {
        outWigFlags.yes=true;
        outWigFlags.format=0;
    } else if (outWigType.at(0)=="wiggle") {
        outWigFlags.yes=true;
        outWigFlags.format=1;
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT ERROR: unrecognized option in --outWigType=" << outWigType.at(0) << "\n";
        errOut << "SOLUTION: use one of the allowed values of --outWigType : 'None' or 'bedGraph' \n";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
    if (outWigStrand.at(0)=="Stranded") {
        outWigFlags.strand=true;
    } else if (outWigStrand.at(0)=="Unstranded") {
        outWigFlags.strand=false;
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT ERROR: unrecognized option in --outWigStrand=" << outWigStrand.at(0) << "\n";
        errOut << "SOLUTION: use one of the allowed values of --outWigStrand : 'Stranded' or 'Unstranded' \n";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (outWigType.size()==1) {//simple bedGraph
        outWigFlags.type=0;
    } else {
        if (outWigType.at(1)=="read1_5p") {
            outWigFlags.type=1;
        } else if (outWigType.at(1)=="read2") {
            outWigFlags.type=2;
        } else {
            ostringstream errOut;
            errOut << "EXITING because of FATAL INPUT ERROR: unrecognized second option in --outWigType=" << outWigType.at(1) << "\n";
            errOut << "SOLUTION: use one of the allowed values of --outWigType : 'read1_5p' \n";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

    //wigOut parameters
    if (outWigNorm.at(0)=="None") {
        outWigFlags.norm=0;
    } else if (outWigNorm.at(0)=="RPM") {
        outWigFlags.norm=1;
    } else {
        ostringstream errOut;
        errOut << "EXITING because of fatal parameter ERROR: unrecognized option in --outWigNorm=" << outWigNorm.at(0) << "\n";
        errOut << "SOLUTION: use one of the allowed values of --outWigNorm : 'None' or 'RPM' \n";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    
    //remove duplicates parameters
    if (removeDuplicates.mode=="UniqueIdentical")
    {
        removeDuplicates.yes=true;
        removeDuplicates.markMulti=true;
    } else if (removeDuplicates.mode=="UniqueIdenticalNotMulti")
    {
        removeDuplicates.yes=true;
        removeDuplicates.markMulti=false;
    } else if (removeDuplicates.mode!="-")
    {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in of --bamRemoveDuplicatesType="<<removeDuplicates.mode<<"\n";
            errOut << "SOLUTION: use allowed option: - or UniqueIdentical or UniqueIdenticalNotMulti";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);        
    };    
    
    if (runMode=="alignReads") {
        inOut->logProgress.open((outFileNamePrefix + "Log.progress.out").c_str());
    } else if (runMode=="inputAlignmentsFromBAM") {
        //at the moment, only wiggle output is implemented
        if (outWigFlags.yes) {
            *inOut->logStdOut << timeMonthDayTime() << " ..... reading from BAM, output wiggle\n" <<flush;
            inOut->logMain << timeMonthDayTime()    << " ..... reading from BAM, output wiggle\n" <<flush;
            string wigOutFileNamePrefix=outFileNamePrefix + "Signal";
            signalFromBAM(inputBAMfile, wigOutFileNamePrefix, *this);
            *inOut->logStdOut << timeMonthDayTime() << " ..... done\n" <<flush;
            inOut->logMain << timeMonthDayTime()    << " ..... done\n" <<flush;
        } else if (removeDuplicates.mode!="-") {
            *inOut->logStdOut << timeMonthDayTime() << " ..... reading from BAM, remove duplicates, output BAM\n" <<flush;
            inOut->logMain << timeMonthDayTime()    << " ..... reading from BAM, remove duplicates, output BAM\n" <<flush;
            bamRemoveDuplicates(inputBAMfile, (outFileNamePrefix+"Processed.out.bam").c_str(), *this);
            *inOut->logStdOut << timeMonthDayTime() << " ..... done\n" <<flush;
            inOut->logMain << timeMonthDayTime()    << " ..... done\n" <<flush;
        } else {
            ostringstream errOut;
            errOut <<"EXITING because of fatal INPUT ERROR: at the moment --runMode inputFromBAM only works with --outWigType bedGraph OR --bamRemoveDuplicatesType Identical"<<"\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
        sysRemoveDir (outFileTmp);
        exit(0);
    };

    outSAMbool=false;
    outBAMunsorted=false;
    outBAMcoord=false;
    if (runMode=="alignReads" && outSAMmode != "None") {//open SAM file and write header
        if (outSAMtype.at(0)=="BAM") {
            if (outSAMtype.size()<2) {
                ostringstream errOut;
                errOut <<"EXITING because of fatal PARAMETER error: missing BAM option\n";
                errOut <<"SOLUTION: re-run STAR with one of the allowed values of --outSAMtype BAM Unsorted OR SortedByCoordinate OR both\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            };
            for (uint32 ii=1; ii<outSAMtype.size(); ii++) {
                if (outSAMtype.at(ii)=="Unsorted") {
                    outBAMunsorted=true;
                } else if (outSAMtype.at(ii)=="SortedByCoordinate") {
                    outBAMcoord=true;
                } else {
                    ostringstream errOut;
                    errOut <<"EXITING because of fatal input ERROR: unknown value for the word " <<ii+1<<" of outSAMtype: "<< outSAMtype.at(ii) <<"\n";
                    errOut <<"SOLUTION: re-run STAR with one of the allowed values of --outSAMtype BAM Unsorted or SortedByCoordinate or both\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                };
            };
            //TODO check for conflicts
            if (outBAMunsorted) {
                if (outStd=="BAM_Unsorted") {
                    outBAMfileUnsortedName="-";
                } else {
                    outBAMfileUnsortedName=outFileNamePrefix + "Aligned.out.bam";
                };
                inOut->outBAMfileUnsorted = bgzf_open(outBAMfileUnsortedName.c_str(),("w"+to_string((long long) outBAMcompression)).c_str());
            };
            if (outBAMcoord) {
                if (outStd=="BAM_SortedByCoordinate") {
                    outBAMfileCoordName="-";
                } else {
                    outBAMfileCoordName=outFileNamePrefix + "Aligned.sortedByCoord.out.bam";
                };
                inOut->outBAMfileCoord = bgzf_open(outBAMfileCoordName.c_str(),("w"+to_string((long long) outBAMcompression)).c_str());
                if (outBAMsortingThreadN==0) {
                    outBAMsortingThreadNactual=min(6, runThreadN);
                } else {
                    outBAMsortingThreadNactual=outBAMsortingThreadN;
                };
                outBAMcoordNbins=max((uint32)outBAMsortingThreadNactual*3,outBAMsortingBinsN);
                outBAMsortingBinStart= new uint64 [outBAMcoordNbins];
                outBAMsortingBinStart[0]=1;//this initial value means that the bin sizes have not been determined yet

                outBAMsortTmpDir=outFileTmp+"/BAMsort/";
                mkdir(outBAMsortTmpDir.c_str(),runDirPerm);
            };
        } else if (outSAMtype.at(0)=="SAM") {
            if (outSAMtype.size()>1)
            {
                ostringstream errOut;
                errOut <<"EXITING because of fatal PARAMETER error: --outSAMtype SAM can cannot be combined with "<<outSAMtype.at(1)<<" or any other options\n";
                errOut <<"SOLUTION: re-run STAR with with '--outSAMtype SAM' only, or with --outSAMtype BAM Unsorted|SortedByCoordinate\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            };
            outSAMbool=true;
            if (outStd=="SAM") {
                inOut->outSAM = & std::cout;
            } else {
                inOut->outSAMfile.open((outFileNamePrefix + "Aligned.out.sam").c_str());
                inOut->outSAM = & inOut->outSAMfile;
            };
        } else if (outSAMtype.at(0)=="None") {
            //nothing to do, all flags are already false
        } else {
            ostringstream errOut;
            errOut <<"EXITING because of fatal input ERROR: unknown value for the first word of outSAMtype: "<< outSAMtype.at(0) <<"\n";
            errOut <<"SOLUTION: re-run STAR with one of the allowed values of outSAMtype: BAM or SAM \n"<<flush;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

    if (!outBAMcoord && outWigFlags.yes && runMode=="alignReads") {
        ostringstream errOut;
        errOut <<"EXITING because of fatal PARAMETER error: generating signal with --outWigType requires sorted BAM\n";
        errOut <<"SOLUTION: re-run STAR with with --outSAMtype BAM SortedByCoordinate, or, id you also need unsroted BAM, with --outSAMtype BAM SortedByCoordinate Unsorted\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    //versions
    for (uint ii=0;ii<2;ii++) {
        if (parArray[ii]->inputLevel>0) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal input ERROR: the version parameter "<< parArray[ii]->nameString << " cannot be re-defined by the user\n";
            errOut <<"SOLUTION: please remove this parameter from the command line or input files and re-start STAR\n"<<flush;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

    //run
    if (runThreadN<=0) {
        ostringstream errOut;
        errOut <<"EXITING: fatal input ERROR: runThreadN must be >0, user-defined runThreadN="<<runThreadN<<"\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    //
    if (outFilterType=="BySJout" && outSAMorder=="PairedKeepInputOrder") {
        ostringstream errOut;
        errOut <<"EXITING: fatal input ERROR: --outFilterType=BySJout is not presently compatible with --outSAMorder=PairedKeepInputOrder\n";
        errOut <<"SOLUTION: re-run STAR without setting one of those parameters.\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
    if (!outSAMbool && outSAMorder=="PairedKeepInputOrder") {
        ostringstream errOut;
        errOut <<"EXITING: fatal input ERROR: --outSAMorder=PairedKeepInputOrder is presently only compatible with SAM output, i.e. default --outSMAtype SAM\n";
        errOut <<"SOLUTION: re-run STAR without --outSAMorder=PairedKeepInputOrder, or with --outSAMorder=PairedKeepInputOrder --outSMAtype SAM .\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };    
    //SJ filtering
    for (int ii=0;ii<4;ii++) {
        if (outSJfilterOverhangMin.at(ii)<0) outSJfilterOverhangMin.at(ii)=numeric_limits<int32>::max();
        if (outSJfilterCountUniqueMin.at(ii)<0) outSJfilterCountUniqueMin.at(ii)=numeric_limits<int32>::max();
        if (outSJfilterCountTotalMin.at(ii)<0) outSJfilterCountTotalMin.at(ii)=numeric_limits<int32>::max();
        if (outSJfilterDistToOtherSJmin.at(ii)<0) outSJfilterDistToOtherSJmin.at(ii)=numeric_limits<int32>::max();

        if (alignSJstitchMismatchNmax.at(ii)<0) alignSJstitchMismatchNmax.at(ii)=numeric_limits<int32>::max();
    };

    if (limitGenomeGenerateRAM==0) {//must be >0
        inOut->logMain <<"EXITING because of FATAL PARAMETER ERROR: limitGenomeGenerateRAM=0\n";
        inOut->logMain <<"SOLUTION: please specify a >0 value for limitGenomeGenerateRAM\n"<<flush;
        exit(1);
    } else if (limitGenomeGenerateRAM>1000000000000) {//
        inOut->logMain <<"WARNING: specified limitGenomeGenerateRAM="<<limitGenomeGenerateRAM<<" bytes appears to be too large, if you do not have enough memory the code will crash!\n"<<flush;
    };


    {//read groups
        if (outSAMattrRGline.at(0)!="-") {
            string linefull;
            for (uint ii=0;ii<outSAMattrRGline.size(); ii++) {//concatenate into one line
                if (ii==0 || outSAMattrRGline.at(ii)==",") {//start new entry
                    if (ii>0) ++ii;//skip comma
                    outSAMattrRGlineSplit.push_back(outSAMattrRGline.at(ii)); //star new RG line with the first field which must be ID:xxx
                    if (outSAMattrRGlineSplit.back().substr(0,3)!="ID:") {
                        ostringstream errOut;
                        errOut <<"EXITING because of FATAL INPUT ERROR: the first word of a line from --outSAMattrRGline="<<outSAMattrRGlineSplit.back()<<" does not start with ID:xxx read group identifier\n";
                        errOut <<"SOLUTION: re-run STAR with all lines in --outSAMattrRGline starting with ID:xxx\n";
                        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                    };
                    outSAMattrRG.push_back(outSAMattrRGlineSplit.back().substr(3)); //this adds the ID field
                } else {//keep adding fields to this RG line, until the next comma
                    outSAMattrRGlineSplit.back()+="\t" + outSAMattrRGline.at(ii);
                };
            };
        };
    };

    outSAMfilter.KeepOnlyAddedReferences=false;
    outSAMfilter.KeepAllAddedReferences=false;
    outSAMfilter.yes=true;
    if (outSAMfilter.mode.at(0)=="KeepOnlyAddedReferences")
    {
        outSAMfilter.KeepOnlyAddedReferences=true;
    } else if (outSAMfilter.mode.at(0)=="KeepAllAddedReferences")  
    {
        outSAMfilter.KeepAllAddedReferences=true;
    } else if (outSAMfilter.mode.at(0)=="None")
    {
      outSAMfilter.yes=false;
    } else
    {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --outSAMfilter: "<<outSAMfilter.mode.at(0) <<"\n";
        errOut <<"SOLUTION: specify one of the allowed values: KeepOnlyAddedReferences or None\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (outMultimapperOrder.mode=="Old_2.4")
    {
        outMultimapperOrder.random=false;
    } else if (outMultimapperOrder.mode=="Random")
    {
        outMultimapperOrder.random=true;
    } else
    {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --outMultimapperOrder: "<<outMultimapperOrder.mode <<"\n";
        errOut <<"SOLUTION: specify one of the allowed values: Old_2.4 or SortedByCoordinate or Random\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    //read parameters
    if (readFilesType.at(0)=="Fastx") {
        readFilesTypeN=1;
    } else if (readFilesType.at(0)=="SAM"){
        readFilesTypeN=10;
    } else {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --readFilesType: "<<readFilesType.at(0) <<"\n";
        errOut <<"SOLUTION: specify one of the allowed values: Fastx or SAM\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
        
    if (readFilesTypeN==1) {
        readNmates=readFilesIn.size(); //for now the number of mates is defined by the number of input files
    } else if (readFilesTypeN==10) {//find the number of mates from the SAM file
        if (readFilesType.size()==2 && readFilesType.at(1)=="SE") {
            readNmates=1;
        } else if (readFilesType.size()==2 && readFilesType.at(1)=="PE") {
            readNmates=2;
        } else {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL INPUT ERROR: --readFilesType SAM requires specifying SE or PE reads"<<"\n";
            errOut <<"SOLUTION: specify --readFilesType SAM SE for single-end reads or --readFilesType SAM PE for paired-end reads\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

    if (runMode=="alignReads" && pGe.gLoad!="Remove" && pGe.gLoad!="LoadAndExit") {//open reads files to check if they are present
        openReadsFiles();

        //check sizes of the mate files, if not the same, assume mates are not the same length
        if (readNmates==1) {
            readMatesEqualLengths=true;
        } else if (readNmates > 2){
            ostringstream errOut;
            errOut <<"EXITING: because of fatal input ERROR: number of read mates files > 2: " <<readNmates << "\n";
            errOut <<"SOLUTION:specify only one or two files in the --readFilesIn option. If file names contain spaces, use quotes: \"file name\"\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        } else if (readMatesLengthsIn=="Equal") {
            readMatesEqualLengths=true;
        } else if (readMatesLengthsIn=="NotEqual") {
            readMatesEqualLengths=false;
        } else {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL input ERROR: the value of the parameter readMatesLengthsIn=" << readMatesLengthsIn <<" is not among the allowed values: Equal or NotEqual\n";
            errOut <<"SOLUTION: specify one of the allowed values: Equal or NotEqual\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };

        if ( runMode=="alignReads" && outReadsUnmapped=="Fastx" ) {//open unmapped reads file
            for (uint imate=0;imate<readNmates;imate++) {
                ostringstream ff;
                ff << outFileNamePrefix << "Unmapped.out.mate" << imate+1;
                inOut->outUnmappedReadsStream[imate].open(ff.str().c_str());
            };
        };


        if (outFilterType=="Normal") {
            outFilterBySJoutStage=0;
        } else if (outFilterType=="BySJout") {
            outFilterBySJoutStage=1;
        } else {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL input ERROR: unknown value of parameter outFilterType: " << outFilterType <<"\n";
            errOut <<"SOLUTION: specify one of the allowed values: Normal | BySJout\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

    if (outSAMmapqUnique<0 || outSAMmapqUnique>255) {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL input ERROR: out of range value for outSAMmapqUnique=" << outSAMmapqUnique <<"\n";
            errOut <<"SOLUTION: specify outSAMmapqUnique within the range of 0 to 255\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    // in/out buffers
    #define BUFFER_InSizeFraction 0.5
    if (limitIObufferSize<limitOutSJcollapsed*Junction::dataSize+1000000)
    {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: --limitIObufferSize="<<limitIObufferSize <<" is too small for ";
        errOut << "--limitOutSJcollapsed*"<<Junction::dataSize<<"="<< limitOutSJcollapsed<<"*"<<Junction::dataSize<<"="<<limitOutSJcollapsed*Junction::dataSize<<"\n";
        errOut <<"SOLUTION: re-run STAR with larger --limitIObufferSize or smaller --limitOutSJcollapsed\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
    chunkInSizeBytesArray=(uint) int((limitIObufferSize-limitOutSJcollapsed*Junction::dataSize)*BUFFER_InSizeFraction)/2;
    chunkOutBAMsizeBytes= (uint) int((1.0/BUFFER_InSizeFraction-1.0)*chunkInSizeBytesArray*2.0);
    chunkInSizeBytes=chunkInSizeBytesArray-2*(DEF_readSeqLengthMax+1)-2*DEF_readNameLengthMax;//to prevent overflow

    //basic trimming
    if (clip5pNbases.size()==1 && readNmates==2) clip5pNbases.push_back(clip5pNbases[0]);
    if (clip3pNbases.size()==1 && readNmates==2) clip3pNbases.push_back(clip3pNbases[0]);
    if (clip3pAfterAdapterNbases.size()==1 && readNmates==2) clip3pAfterAdapterNbases.push_back(clip3pAfterAdapterNbases[0]);

    //adapter clipping
    if (clip3pAdapterSeq.size()==1 && readNmates==2) clip3pAdapterSeq.push_back(clip3pAdapterSeq[0]);
    if (clip3pAdapterMMp.size()==1 && readNmates==2) clip3pAdapterMMp.push_back(clip3pAdapterMMp[0]);
    for (uint imate=0;imate<readNmates;imate++) {
        if (clip3pAdapterSeq.at(imate).at(0)=='-') {// no clipping
            clip3pAdapterSeq.at(imate).assign("");
        } else {//clipping
            clip3pAdapterSeqNum[imate]=new char [clip3pAdapterSeq.at(imate).length()];
            convertNucleotidesToNumbers(clip3pAdapterSeq.at(imate).data(),clip3pAdapterSeqNum[imate],clip3pAdapterSeq.at(imate).length());
            //inOut->fastaOutSeqs.open("Seqs.out.fasta");
        };
    };
    
    //variation
    var.yes=false;
    if (var.vcfFile!="-")
    {
        var.yes=true;
    };         
        
    //WASP
    wasp.yes=false;
    wasp.SAMtag=false;
    if (wasp.outputMode=="SAMtag") {
        wasp.yes=true;
        wasp.SAMtag=true;
    } else if (wasp.outputMode=="None") {
        //nothing to do
    } else {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: unknown/unimplemented --waspOutputMode option: "<<wasp.outputMode <<"\n";
        errOut <<"SOLUTION: re-run STAR with allowed --waspOutputMode options: None or SAMtag\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (wasp.yes && !var.yes) {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: --waspOutputMode option requires VCF file: "<<wasp.outputMode <<"\n";
        errOut <<"SOLUTION: re-run STAR with --waspOutputMode ... and --varVCFfile /path/to/file.vcf\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    
     if (wasp.yes && outSAMtype.at(0)!="BAM") {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: --waspOutputMode requires output to BAM file\n";
        errOut <<"SOLUTION: re-run STAR with --waspOutputMode ... and --outSAMtype BAM ... \n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
    
    //outSAMattributes
    outSAMattrPresent.NH=false;//TODO re-write as class with constructor?
    outSAMattrPresent.HI=false;
    outSAMattrPresent.AS=false;
    outSAMattrPresent.NM=false;
    outSAMattrPresent.MD=false;
    outSAMattrPresent.nM=false;
    outSAMattrPresent.jM=false;
    outSAMattrPresent.jI=false;
    outSAMattrPresent.RG=false;
    outSAMattrPresent.MC=false;    
    outSAMattrPresent.XS=false;
    outSAMattrPresent.vA=false;
    outSAMattrPresent.vG=false;
    outSAMattrPresent.vW=false;
    outSAMattrPresent.ch=false;
    outSAMattrPresent.rB=false;
            
    //for quant SAM output only NH and HI flags
    outSAMattrPresentQuant=outSAMattrPresent;
    outSAMattrPresentQuant.NH=true;
    outSAMattrPresentQuant.HI=true;
    outSAMattrOrderQuant.push_back(ATTR_NH);
    outSAMattrOrderQuant.push_back(ATTR_HI);

    vector<string> vAttr1;
    if (outSAMattributes.at(0)=="None") {
    } else if (outSAMattributes.at(0)=="All"){
        vAttr1={"NH","HI","AS","nM","NM","MD","jM","jI","rB","MC","ch"};
    } else if (outSAMattributes.at(0)=="Standard"){
        vAttr1={"NH","HI","AS","nM"};
    } else {
        vAttr1=outSAMattributes;
    };

    for (uint ii=0;ii<vAttr1.size();ii++) {
        if        (vAttr1.at(ii)== "NH") {
            outSAMattrOrder.push_back(ATTR_NH);
            outSAMattrPresent.NH=true;
        } else if (vAttr1.at(ii)== "HI") {
            outSAMattrOrder.push_back(ATTR_HI);
            outSAMattrPresent.HI=true;
        } else if (vAttr1.at(ii)== "AS") {
            outSAMattrOrder.push_back(ATTR_AS);
            outSAMattrPresent.AS=true;
        } else if (vAttr1.at(ii)== "NM") {
            outSAMattrOrder.push_back(ATTR_NM);
            outSAMattrPresent.NM=true;
        } else if (vAttr1.at(ii)== "MD") {
            outSAMattrOrder.push_back(ATTR_MD);
            outSAMattrPresent.MD=true;
        } else if (vAttr1.at(ii)== "nM") {
            outSAMattrOrder.push_back(ATTR_nM);
            outSAMattrPresent.nM=true;
        } else if (vAttr1.at(ii)== "jM") {
            outSAMattrOrder.push_back(ATTR_jM);
            outSAMattrPresent.jM=true;
        } else if (vAttr1.at(ii)== "jI") {
            outSAMattrOrder.push_back(ATTR_jI);
            outSAMattrPresent.jI=true;
        } else if (vAttr1.at(ii)== "vA") {
            outSAMattrOrder.push_back(ATTR_vA);
            outSAMattrPresent.vA=true;
        } else if (vAttr1.at(ii)== "vG") {
            outSAMattrOrder.push_back(ATTR_vG);
            outSAMattrPresent.vG=true; 
        } else if (vAttr1.at(ii)== "vW") {
            outSAMattrOrder.push_back(ATTR_vW);
            outSAMattrPresent.vW=true;              
        } else if (vAttr1.at(ii)== "RG") {
            outSAMattrOrder.push_back(ATTR_RG);
            outSAMattrOrderQuant.push_back(ATTR_RG);
            outSAMattrPresent.RG=true;
        } else if (vAttr1.at(ii)== "rB") {
            outSAMattrOrder.push_back(ATTR_rB);
            outSAMattrOrderQuant.push_back(ATTR_rB);
            outSAMattrPresent.rB=true;            
        } else if (vAttr1.at(ii)== "ch") {
            outSAMattrOrder.push_back(ATTR_ch);
            outSAMattrOrderQuant.push_back(ATTR_ch);
            outSAMattrPresent.ch=true;      
        } else if (vAttr1.at(ii)== "MC") {
            outSAMattrOrder.push_back(ATTR_MC);
            outSAMattrOrderQuant.push_back(ATTR_MC);
            outSAMattrPresent.MC=true;                        
        } else if (vAttr1.at(ii)== "XS") {
            outSAMattrOrder.push_back(ATTR_XS);
            outSAMattrPresent.XS=true;
            if (outSAMstrandField!="intronMotif") {
                inOut->logMain << "WARNING --outSAMattributes contains XS, therefore STAR will use --outSAMstrandField intronMotif" <<endl;
                outSAMstrandField="intronMotif";
            };
        } else {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL INPUT ERROR: unknown/unimplemented SAM atrribute (tag): "<<vAttr1.at(ii) <<"\n";
            errOut <<"SOLUTION: re-run STAR with --outSAMattributes that contains only implemented attributes\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };
    
    if  (!var.yes && (outSAMattrPresent.vA | outSAMattrPresent.vG)) {
        ostringstream errOut;
        errOut <<"EXITING because of fatal PARAMETER error: --outSAMattributes contains vA and/or vG tag(s), but --varVCFfile is not set\n";
        errOut <<"SOLUTION: re-run STAR with a --varVCFfile option, or without vA/vG tags in --outSAMattributes\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);        
    };
    
    if (!wasp.yes && outSAMattrPresent.vW) {
        ostringstream errOut;
        errOut <<"EXITING because of fatal PARAMETER error: --outSAMattributes contains vW tag, but --waspOutputMode is not set\n";
        errOut <<"SOLUTION: re-run STAR with a --waspOutputMode option, or without vW tags in --outSAMattributes\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);        
    };
    
    if (outSAMattrRG.size()>0 && !outSAMattrPresent.RG) {
        outSAMattrOrder.push_back(ATTR_RG);
        outSAMattrOrderQuant.push_back(ATTR_RG);
        outSAMattrPresent.RG=true;
        inOut->logMain << "WARNING --outSAMattrRG defines a read group, therefore STAR will output RG attribute" <<endl;
    } else if (outSAMattrRG.size()==0 && outSAMattrPresent.RG) {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL INPUT ERROR: --outSAMattributes contains RG tag, but --outSAMattrRGline is not set\n";
            errOut <<"SOLUTION: re-run STAR with a valid read group parameter --outSAMattrRGline\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (outSAMstrandField=="intronMotif" && !outSAMattrPresent.XS) {
        outSAMattrOrder.push_back(ATTR_XS);
        inOut->logMain << "WARNING --outSAMstrandField=intronMotif, therefore STAR will output XS attribute" <<endl;
    };

    if (wasp.yes && !outSAMattrPresent.vW) {
        outSAMattrOrder.push_back(ATTR_vW);
        outSAMattrOrderQuant.push_back(ATTR_vW);
        outSAMattrPresent.vW=true;
        inOut->logMain << "WARNING --waspOutputMode is set, therefore STAR will output vW attribute" <<endl;
    };    
    
    //chimeric
    pCh.out.bam=false;
    pCh.out.junctions=false;
    pCh.out.samOld=false;
    pCh.out.bamHardClip=true;//default
    if (pCh.segmentMin>0) {//only if chimeric detection is activated
        for (const auto& type1 : pCh.out.type) {
            if (type1=="WithinBAM") {
                pCh.out.bam=true;
            } else if (type1=="SeparateSAMold") {
                pCh.out.samOld=true;
            } else if (type1=="Junctions") {
                pCh.out.junctions=true;
            } else if (type1=="HardClip") {
                pCh.out.bamHardClip=true;
            } else if (type1=="SoftClip") {
                pCh.out.bamHardClip=false;     
            } else {
                ostringstream errOut;
                errOut <<"EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --chimOutType: "<<type1 <<"\n";
                errOut <<"SOLUTION: re-run STAR with --chimOutType Junctions , SeparateSAMold  , WithinBAM , HardClip \n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            };
        };

        if (pCh.out.bam && !outBAMunsorted && !outBAMcoord) {
                ostringstream errOut;
                errOut <<"EXITING because of fatal PARAMETERS error: --chimOutType WithinBAM requires BAM output\n";
                errOut <<"SOLUTION: re-run with --outSAMtype BAM Unsorted/SortedByCoordinate\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };

        if (pCh.multimapNmax>0 && (pCh.out.bam || pCh.out.samOld)) {       
                ostringstream errOut;
                errOut <<"EXITING because of fatal PARAMETERS error: --chimMultimapNmax > 0 (new chimeric detection) presently only works with --chimOutType Junctions\n";
                errOut <<"SOLUTION: re-run with --chimOutType Junctions\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };

        if (peOverlap.NbasesMin > 0) {
            if (pCh.multimapNmax == 0 && (pCh.out.junctions || pCh.out.samOld)) {
                    ostringstream errOut;
                    errOut <<"EXITING because of fatal PARAMETERS error: --chimMultimapNmax 0 (default old chimeric detection) and --peOverlapNbasesMin > 0 (merging ovelrapping mates) presently only works with --chimOutType WithinBAM\n";
                    errOut <<"SOLUTION: re-run with --chimOutType WithinBAM\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            };
        };



        if (pCh.out.bam && !outSAMattrPresent.NM) {
           outSAMattrOrder.push_back(ATTR_NM);
           inOut->logMain << "WARNING --chimOutType=WithinBAM, therefore STAR will output NM attribute" <<endl;
        };

        pCh.filter.genomicN=false;
        for (uint ii=0; ii<pCh.filter.stringIn.size(); ii++)
        {
            if (pCh.filter.stringIn.at(ii)=="banGenomicN")
            {
                pCh.filter.genomicN=true;
            }
            else if (pCh.filter.stringIn.at(ii)=="None")
            {//nothing to do
            }
            else
            {
                ostringstream errOut;
                errOut << "EXITING because of fatal PARAMETERS error: unrecognized value of --chimFilter="<<pCh.filter.stringIn.at(ii)<<"\n";
                errOut << "SOLUTION: use allowed values: banGenomicN || None";
                exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            };
        };
    };

    alignEndsType.ext[0][0]=false;
    alignEndsType.ext[0][1]=false;
    alignEndsType.ext[1][0]=false;
    alignEndsType.ext[1][1]=false;

    if (alignEndsType.in=="EndToEnd")
    {
        alignEndsType.ext[0][0]=true;
        alignEndsType.ext[0][1]=true;
        alignEndsType.ext[1][0]=true;
        alignEndsType.ext[1][1]=true;
    } else if (alignEndsType.in=="Extend5pOfRead1" )
    {
        alignEndsType.ext[0][0]=true;
    } else if (alignEndsType.in=="Extend5pOfReads12" )
    {
        alignEndsType.ext[0][0]=true;        
        alignEndsType.ext[1][0]=true;  
    } else if (alignEndsType.in=="Extend3pOfRead1" )
    {
        alignEndsType.ext[0][1]=true;
    } else if (alignEndsType.in=="Local")
    {
        //nothing to do for now
    } else
    {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --alignEndsType: "<<alignEndsType.in <<"\n";
        errOut <<"SOLUTION: re-run STAR with --alignEndsType Local OR EndToEnd OR Extend5pOfRead1 OR Extend3pOfRead1\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

//     #ifdef COMPILE_NO_SHM
//         if (pGe.gLoad!="NoSharedMemory") {
//             ostringstream errOut;
//             errOut <<"EXITING because of FATAL INPUT ERROR: The code was compiled with NO SHARED MEMORY support, but pGe.gLoad="<<pGe.gLoad<<"\n";
//             errOut <<"SOLUTION: run STAR with    --genomeLoad NoSharedMemory    option\n";
//             exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
//         };
//     #endif

    //open compilation-dependent streams
    #ifdef OUTPUT_localChains
            inOut->outLocalChains.open((outFileNamePrefix + "LocalChains.out.tab").c_str());
    #endif

//     genomeNumToNT={'A','C','G','T','N'};
    strcpy(genomeNumToNT,"ACGTN");

    if (pGe.gLoad!="LoadAndKeep" && pGe.gLoad!="LoadAndRemove" && pGe.gLoad!="Remove" && pGe.gLoad!="LoadAndExit" && pGe.gLoad!="NoSharedMemory") {// find shared memory fragment
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT ERROR: --genomeLoad=" << pGe.gLoad << "\n" <<flush;
        errOut << "SOLUTION: use one of the allowed values of --genomeLoad : NoSharedMemory,LoadAndKeep,LoadAndRemove,LoadAndExit,Remove.\n" <<flush;
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    //quantification parameters
    quant.yes=false;
    quant.geCount.yes=false;
    quant.trSAM.yes=false;
    if (quant.mode.at(0) != "-") {
        quant.yes=true;
        for (uint32 ii=0; ii<quant.mode.size(); ii++) {
            if (quant.mode.at(ii)=="TranscriptomeSAM") {
                quant.trSAM.yes=true;
                if (outStd=="BAM_Quant") {
                    outFileNamePrefix="-";
                } else {
                    outQuantBAMfileName=outFileNamePrefix + "Aligned.toTranscriptome.out.bam";
                };
                inOut->outQuantBAMfile=bgzf_open(outQuantBAMfileName.c_str(),("w"+to_string((long long) quant.trSAM.bamCompression)).c_str());
                if (quant.trSAM.ban=="IndelSoftclipSingleend")
                {
                    quant.trSAM.indel=false;
                    quant.trSAM.softClip=false;
                    quant.trSAM.singleEnd=false;
                } else if (quant.trSAM.ban=="Singleend")
                {
                    quant.trSAM.indel=true;
                    quant.trSAM.softClip=true;
                    quant.trSAM.singleEnd=false;
                };
            } else if  (quant.mode.at(ii)=="GeneCounts") {
                quant.geCount.yes=true;
                quant.geCount.outFile=outFileNamePrefix + "ReadsPerGene.out.tab";
            } else {
                ostringstream errOut;
                errOut << "EXITING because of fatal INPUT error: unrecognized option in --quant.mode=" << quant.mode.at(ii) << "\n";
                errOut << "SOLUTION: use one of the allowed values of --quant.mode : TranscriptomeSAM or - .\n";
                exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            };
        };
    };


    //two-pass
    if (parArray.at(twoPass.pass1readsN_par)->inputLevel>0  && twoPass.mode=="None")
    {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: --twopass1readsN is defined, but --twoPassMode is not defined\n";
        errOut << "SOLUTION: to activate the 2-pass mode, use --twopassMode Basic";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };



    twoPass.yes=false;
    if (twoPass.mode!="None") {//2-pass parameters
        if (runMode!="alignReads")
        {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: 2-pass mapping option  can only be used with --runMode alignReads\n";
            errOut << "SOLUTION: remove --twopassMode option";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };

        if (twoPass.mode!="Basic")
        {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: unrecognized value of --twopassMode="<<twoPass.mode<<"\n";
            errOut << "SOLUTION: for the 2-pass mode, use allowed values --twopassMode: Basic";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };

        if (twoPass.pass1readsN==0)
        {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: --twopass1readsN = 0 in the 2-pass mode\n";
            errOut << "SOLUTION: for the 2-pass mode, specify --twopass1readsN > 0. Use a very large number or -1 to map all reads in the 1st pass.\n";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };

        if (pGe.gLoad!="NoSharedMemory") {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: 2-pass method is not compatible with --genomeLoad "<<pGe.gLoad<<"\n";
            errOut << "SOLUTION: re-run STAR with --genomeLoad NoSharedMemory ; this is the only option compatible with --twopassMode Basic .\n";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
        twoPass.yes=true;
        twoPass.dir=outFileNamePrefix+"_STARpass1/";
        sysRemoveDir (twoPass.dir);
        if (mkdir (twoPass.dir.c_str(),runDirPerm)!=0) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal ERROR: could not make pass1 directory: "<< twoPass.dir<<"\n";
            errOut <<"SOLUTION: please check the path and writing permissions \n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

   //sjdb insert on the fly

    sjdbInsert.pass1=false;
    sjdbInsert.pass2=false;
    sjdbInsert.yes=false;
    if (pGe.sjdbFileChrStartEnd.at(0)!="-" || pGe.sjdbGTFfile!="-")
    {//will insert annotated sjdb on the fly
       sjdbInsert.pass1=true;
       sjdbInsert.yes=true;
    };
    if (twoPass.yes)
    {
       sjdbInsert.pass2=true;
       sjdbInsert.yes=true;
    };

    if (pGe.gLoad!="NoSharedMemory" && sjdbInsert.yes )
    {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: on the fly junction insertion and 2-pass mappng cannot be used with shared memory genome \n" ;
        errOut << "SOLUTION: run STAR with --genomeLoad NoSharedMemory to avoid using shared memory\n" <<flush;
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (runMode=="alignReads" && sjdbInsert.yes )
    {//run-time genome directory, this is needed for genome files generated on the fly
        if (pGe.sjdbOverhang<=0) {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: pGe.sjdbOverhang <=0 while junctions are inserted on the fly with --sjdbFileChrStartEnd or/and --sjdbGTFfile\n";
            errOut << "SOLUTION: specify pGe.sjdbOverhang>0, ideally readmateLength-1";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
        sjdbInsert.outDir=outFileNamePrefix+"_STARgenome/";
        sysRemoveDir (sjdbInsert.outDir);
        if (mkdir (sjdbInsert.outDir.c_str(),runDirPerm)!=0) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal ERROR: could not make run-time genome directory directory: "<< sjdbInsert.outDir<<"\n";
            errOut <<"SOLUTION: please check the path and writing permissions \n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

    if (outBAMcoord && limitBAMsortRAM==0) {//check limitBAMsortRAM
        if (pGe.gLoad!="NoSharedMemory") {
            ostringstream errOut;
            errOut <<"EXITING because of fatal PARAMETERS error: limitBAMsortRAM=0 (default) cannot be used with --genomeLoad="<<pGe.gLoad <<", or any other shared memory options\n";
            errOut <<"SOLUTION: please use default --genomeLoad NoSharedMemory, \n        OR specify --limitBAMsortRAM the amount of RAM (bytes) that can be allocated for BAM sorting in addition to shared memory allocated for the genome.\n        --limitBAMsortRAM typically has to be > 10000000000 (i.e 10GB).\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
        inOut->logMain<<"WARNING: --limitBAMsortRAM=0, will use genome size as RAM limit for BAM sorting\n";
    };

    for (uint ii=0; ii<readNameSeparator.size(); ii++)
    {
        if (readNameSeparator.at(ii)=="space")
        {
            readNameSeparatorChar.push_back(' ');
        }
        else if (readNameSeparator.at(ii)=="none")
        {
            //nothing to do
        }
        else if (readNameSeparator.at(ii).size()==1)
        {
            readNameSeparatorChar.push_back(readNameSeparator.at(ii).at(0));
        }
        else
        {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: unrecognized value of --readNameSeparator="<<readNameSeparator.at(ii)<<"\n";
            errOut << "SOLUTION: use allowed values: space OR single characters";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };


    //outSAMunmapped
    if (outSAMunmapped.mode.at(0)=="None" && outSAMunmapped.mode.size()==1) {
        outSAMunmapped.yes=false;
        outSAMunmapped.within=false;
        outSAMunmapped.keepPairs=false;
    } else if (outSAMunmapped.mode.at(0)=="Within" && outSAMunmapped.mode.size()==1) {
        outSAMunmapped.yes=true;
        outSAMunmapped.within=true;
        outSAMunmapped.keepPairs=false;   
    } else if (outSAMunmapped.mode.at(0)=="Within" && outSAMunmapped.mode.at(1)=="KeepPairs") {
        outSAMunmapped.yes=true;
        outSAMunmapped.within=true;
        outSAMunmapped.keepPairs=true;   
    } else {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option for --outSAMunmapped=";
        for (uint ii=0; ii<outSAMunmapped.mode.size(); ii++) errOut <<" "<< outSAMunmapped.mode.at(ii);
        errOut << "\n";
        errOut << "SOLUTION: use allowed options: None OR Within OR Within KeepPairs";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    alignEndsProtrude.nBasesMax=stoi(alignEndsProtrude.in.at(0),nullptr);
    if (alignEndsProtrude.nBasesMax>0) {//allow ends protrusion
        if (alignEndsProtrude.in.at(1)=="ConcordantPair") {
            alignEndsProtrude.concordantPair=true;
        } else if (alignEndsProtrude.in.at(1)=="DiscordantPair") {
            alignEndsProtrude.concordantPair=false;
        } else  {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in of --alignEndsProtrude="<<alignEndsProtrude.in.at(1)<<"\n";
            errOut << "SOLUTION: use allowed option: ConcordantPair or DiscordantPair";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };
    
    if (alignInsertionFlush.in=="None") {
        alignInsertionFlush.flushRight=false;
    } else if (alignInsertionFlush.in=="Right") {
        alignInsertionFlush.flushRight=true;
    } else  {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in of --alignInsertionFlush="<<alignInsertionFlush.in<<"\n";
        errOut << "SOLUTION: use allowed option: None or Right";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
        
    
    //peOverlap
    if (peOverlap.NbasesMin>0) {
        peOverlap.yes=true;
    } else {
        peOverlap.yes=false;
    };
    
    ////////////////////////////////////////////////
    inOut->logMain << "Finished loading and checking parameters\n" <<flush;
};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Parameters::scanAllLines (istream &streamIn, int inputLevel,  int inputLevelRequested) {//scan
//     istringstream stringInStream (stringIn);
    string lineIn;
    while (getline(streamIn,lineIn)) {
        scanOneLine(lineIn, inputLevel, inputLevelRequested);
    };
};

int Parameters::scanOneLine (string &lineIn, int inputLevel, int inputLevelRequested) {//scan one line and load the parameters,
                                                             //0 if comment, 1 if OK
    if (lineIn=="") return 0; //empty line

    istringstream lineInStream (lineIn);

    if (inputLevel==0 && ( lineIn.substr(0,1)==" " || lineIn.substr(0,1)=="\t" ) ) return 0;//for Default input spaces also mark comments, for nice formatting

    string parIn("");
    lineInStream >> parIn;
    if (parIn=="" || parIn.substr(0,2)=="//" || parIn.substr(0,1)=="#") return 0; //this is a comment

    uint iPar;
    for (iPar=0; iPar<parArray.size(); iPar++) {
        if (parIn==parArray[iPar]->nameString) {//
            if (inputLevelRequested < 0 || inputLevelRequested == parArray[iPar]->inputLevelAllowed) {
                break;//will read this parameter values
            } else {
                return 1; //do not read inputs not requested at this level
            };
        };
    };

    string parV("");
    lineInStream >> parV;
    if (parV=="") {//parameter value cannot be empty
        ostringstream errOut;
        errOut << "EXITING: FATAL INPUT ERROR: empty value for parameter \""<< parIn << "\" in input \"" << parameterInputName.at(inputLevel) <<"\"\n";
        errOut << "SOLUTION: use non-empty value for this parameter\n"<<flush;
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    lineInStream.str(lineIn); lineInStream.clear(); lineInStream >> parIn; //get the correct state of stream, past reading parIn

    if (iPar==parArray.size()) {//string is not identified
        ostringstream errOut;
        errOut << "EXITING: FATAL INPUT ERROR: unrecoginzed parameter name \""<< parIn << "\" in input \"" << parameterInputName.at(inputLevel) <<"\"\n";
        errOut << "SOLUTION: use correct parameter name (check the manual)\n"<<flush;
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    } else {//found the corresponding parameter
        if (inputLevel==0 && parArray[iPar]->inputLevel>0) {//this is one of the initial parameters, it was read from Command Line and should not be re-defined
            getline(lineInStream,parV);
            inOut->logMain << setiosflags(ios::left) << setw(PAR_NAME_PRINT_WIDTH) << parArray[iPar]->nameString <<parV<<" ... is RE-DEFINED on Command Line as: " << *(parArray[iPar]) <<"\n";
        } else if (parArray[iPar]->inputLevelAllowed>0 && parArray[iPar]->inputLevelAllowed < inputLevel) {//this is initial parameter and cannot be redefined
            ostringstream errOut;
            errOut << "EXITING: FATAL INPUT ERROR: parameter \""<< parIn << "\" cannot be defined at the input level \"" << parameterInputName.at(inputLevel) << "\"\n";
            errOut << "SOLUTION: define parameter \""<< parIn << "\" in \"" << parameterInputName.at(parArray[iPar]->inputLevelAllowed) <<"\"\n" <<flush;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        } else if (parArray[iPar]->inputLevel==inputLevel) {//this parameter was already defined at this input level
            ostringstream errOut;
            errOut << "EXITING: FATAL INPUT ERROR: duplicate parameter \""<< parIn << "\" in input \"" << parameterInputName.at(inputLevel) << "\"\n";
            errOut << "SOLUTION: keep only one definition of input parameters in each input source\n"<<flush;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        } else {//read values
            parArray[iPar]->inputValues(lineInStream);
            parArray[iPar]->inputLevel=inputLevel;
            if ( inOut->logMain.good() ) {
                inOut->logMain << setiosflags(ios::left) << setw(PAR_NAME_PRINT_WIDTH) << parArray[iPar]->nameString << *(parArray[iPar]);
                if ( parArray[iPar]->inputLevel > 0 ) inOut->logMain <<"     ~RE-DEFINED";
                inOut->logMain << endl;
            };
        };
    };
    return 0;
};



