#include "ParametersSolo.h"
#include "Parameters.h"
#include "ErrorWarning.h"
#include "streamFuns.h"
#include "SequenceFuns.h"
#include "serviceFuns.cpp"

#include <stdlib.h>

void ParametersSolo::initialize(Parameters *pPin)
{
    pP=pPin;

    cellFiltering();
    if (pP->runMode=="soloCellFiltering") {//only filtering happens, do not need any other parameters
        yes=true;
        umiDedup.typesIn = {"NoDedup"}; //this does not affect the results - the dedup had been done when raw matrix was generated
        umiDedup.initialize(this);
        return;
    };
    
    //constants - may turn into parameters in the future
    redistrReadsNfiles = 3*pP->runThreadN;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////--soloType
    barcodeStart=barcodeEnd=0; //this means that the entire barcodeRead is considered barcode. Will change it for simple barcodes.    
    yes = true;
    if (typeStr=="None" || typeStr=="SmartSeq") {//solo SAM attributes not allowed
            
        if (   pP->outSAMattrPresent.CR  || pP->outSAMattrPresent.CY  || pP->outSAMattrPresent.UR 
            || pP->outSAMattrPresent.UY  || pP->outSAMattrPresent.CB  || pP->outSAMattrPresent.UB
            || pP->outSAMattrPresent.sS  || pP->outSAMattrPresent.sQ  || pP->outSAMattrPresent.sM
            || pP->outSAMattrPresent.sF
           ) {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL INPUT ERROR: --outSAMattributes contains CR/CY/UR/UY/CB/UB tags which are not allowed for --soloType " << typeStr <<'\n';
            errOut <<"SOLUTION: re-run STAR without these attribures\n";
            exitWithError(errOut.str(), std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };        
    };
    
    if (typeStr=="None") {
        type = SoloTypes::None;
        yes = false;
        samAttrYes = false;     
        return;        
    } else if (typeStr=="CB_UMI_Simple" || typeStr=="Droplet") {
        type=SoloTypes::CB_UMI_Simple;        
        if (umiL > 16) {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: UMI length is too long: --soloUMIlen="<<umiL<<"\n";
            errOut << "SOLUTION: UMI length cannot be longer than 16";
            exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };
        if (cbL > 31) {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: CB length is too long: --soloCBlen="<<cbL<<"\n";
            errOut << "SOLUTION: CB length cannot be longer than 31";
            exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };
        
        cbumiL=cbL+umiL;
        if (bL==1)
            bL=cbumiL;
        
        barcodeStart=min(cbS,umiS)-1;
        barcodeEnd=max(cbS+cbL,umiS+umiL)-2;
        
    } else if (typeStr=="CB_UMI_Complex") {
        type=SoloTypes::CB_UMI_Complex;
        bL=0;//no fixed length for barcode sequence
        cbumiL=0;

    } else if (typeStr=="CB_samTagOut") {
        type=SoloTypes::CB_samTagOut;
        cbumiL=cbL+umiL;
        if (bL==1)
            bL=cbumiL;

        barcodeStart=min(cbS,umiS)-1;
        barcodeEnd=max(cbS+cbL,umiS+umiL)-2;
        
    } else if (typeStr=="SmartSeq") {
        type=SoloTypes::SmartSeq;     
       
    } else  {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in --soloType="<<typeStr<<"\n";
        errOut << "SOLUTION: use allowed option: None OR CB_UMI_Simple OR CB_UMI_Complex\n";
        errOut << "Obsolete option Droplet should be replaced with CB_UMI_Simple";
        exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };

    /////////////////////////// CB type
    if (CBtype.typeString=="Sequence") {
        CBtype.type = 1;
    } else if (CBtype.typeString=="String") {
        if (soloCBwhitelist[0]!="None") {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: --soloCBtype String cannot be used with passlist\n";
            errOut << "SOLUTION: use --soloCBwhitelist None\n";
            exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);

        };
        CBtype.type = 2;
        CBtype.strMtx = new std::mutex;
        CBtype.strMap.reserve(100000);        
    } else  {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in --soloCBtype="<< CBtype.typeString <<"\n";
        errOut << "SOLUTION: use allowed option: Sequence OR String\n";
        exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////input read files
    barcodeRead=-1; 
    barcodeReadSeparate=false;
    
    if (pP->readFilesTypeN != 10) {//input from FASTQ
        if (type==SoloTypes::SmartSeq) {//no barcode read
            //TODO: a lot of parameters should not be defined for SmartSeq option - check it here
        } else {//all other types require barcode read
            if (barcodeReadIn == 0) {//barcode read is separate - needs to be the last read in --readFilesIn
                if (pP->readNends < 2) {
                    exitWithError("EXITING because of fatal PARAMETERS error: --soloType (except SmartSeq) with --soloBarcodeMate 0 (default) require 2 reads or 3 reads, where the last read is the barcode read.\n"
                                  "SOLUTION: if barcode is in a separate mate, specify it as the last file in --readFilesIn. If barcode sequence is a part of one of the mates, specify that mate with --soloBarcodeMate 1 (or 2 or 3)"
                                  ,std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
                };
                pP->readNmates = pP->readNends-1; //true mates, excluding barcode read
                barcodeRead = pP->readNends-1;//the barcode read is always the last one 
            } else if (barcodeReadIn > pP->readNends) {
                exitWithError("EXITING because of fatal PARAMETERS error: --soloBarcodeMate " +to_string(barcodeReadIn)+ "is larger than number of mates " + to_string(pP->readNends) +
                                "\nSOLUTION: specify --soloBarcodeMate <= than the number of mates.",std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
            } else {//barcode sequence is in one of the mates
                if (type != SoloTypes::CB_UMI_Simple) {
                    exitWithError("EXITING because of fatal PARAMETERS error: --soloBarcodeMate " +to_string(barcodeReadIn)+ ">0 for is not allowed for --soloType " + typeStr +
                                    "\nSOLUTION: specify --soloBarcodeMate 0   or   --soloType CB_UMI_Simple",std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
                };
                
                barcodeRead = barcodeReadIn-1;
                barcodeReadSeparate = true;
                bL = 0; //>0 only for separate barcode read
                    
                if ( pP->pClip.in[0].N[barcodeRead] == 0 && pP->pClip.in[1].N[barcodeRead] == 0 ){//clipping not specified for the mate with barcodes
                    exitWithError("EXITING because of fatal PARAMETERS error: --soloBarcodeMate " +to_string(barcodeReadIn)+ " specifies that barcode sequence is a part of the mate " +to_string(barcodeReadIn)+
                                    ", which requires clipping the barcode off this mate."
                                    "\nSOLUTION: clip the barcode sequence from 5' or/and 3' with --clip5pNbases   or/and --clip3pNbases . The values for mate1 and mate2 have to specified, specify 0 for no clipping.",std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
                };
            };
        };
    } else if (pP->readFilesTypeN == 10) {//input from SAM
        if (typeStr=="SmartSeq") {
            exitWithError("EXITING because of fatal PARAMETERS error: --readFilesType SAM SE/PE cannot be used with --soloType SmartSeq\n"
                          "SOLUTION: for Smart-seq input from BAM files, use --soloType CB_UMI_Simple , create whitelist of SmartSeq file names, and specify the SAM tag that records these file names in --soloInputSAMattrBarcodeSeq"
                          , std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };            
            
        
        if (samAtrrBarcodeSeq.at(0) == "-") {
            exitWithError("EXITING because of fatal PARAMETERS error: --readFilesType SAM SE/PE requires --soloInputSAMattrBarcodeSeq.\n"
                          "SOLUTION: specify input SAM attributes for barcode sequence in --soloInputSAMattrBarcodeSeq, and (optionally) quality with --soloInputSAMattrBarcodeQual"
                          , std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };
        
        if (samAtrrBarcodeQual.at(0) == "-") {
            warningMessage(" since --readFilesType SAM SE/PE --soloInputSAMattrBarcodeQual - : qualities for barcode read will be replaced with 'H'", pP->inOut->logMain,std::cerr, *pP);
            samAtrrBarcodeQual.clear();
        };
        
        for (auto &tag: samAtrrBarcodeSeq) {
            if (tag.size()!=2) {
                exitWithError("EXITING because of fatal PARAMETERS error: --soloInputSAMattrBarcodeSeq attributes have to be two-letter strings.\n"
                              "SOLUTION: specify correct two-letter strings in --soloInputSAMattrBarcodeSeq", 
                              std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
            };
            tag.insert(0,1,'\t');
        };
        
        for (auto &tag: samAtrrBarcodeQual) {
            if (tag.size()!=2) {
                exitWithError("EXITING because of fatal PARAMETERS error: --soloInputSAMattrBarcodeQual attributes have to be two-letter strings.\n"
                              "SOLUTION: specify correct two-letter strings in --soloInputSAMattrBarcodeQual", 
                              std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
            };
            tag.insert(0,1,'\t');
        };        
    };
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////// soloStrand
    if (strandStr=="Unstranded") {
        strand=-1;
    } else if (strandStr=="Forward") {
        strand=0;
    } else if (strandStr=="Reverse") {
        strand=1;
    } else  {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in --soloStrand="<<strandStr<<"\n";
        errOut << "SOLUTION: use allowed option: Unstranded OR Forward OR Reverse";
        exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////// features
    featureInd.fill(-1);
    featureYes.fill(false);
    featureFirst = -1;
    for (auto &fin : featureIn) {
        for (uint32 ii=0; ii<SoloFeatureTypes::Names.size(); ii++) {
            if (fin==SoloFeatureTypes::Names[ii]) {
                featureYes[ii]=true;
                features.push_back(ii);
                if (featureFirst == -1)
                    featureFirst = ii; //recored the first feature on the list
                break;
            };
        };
    };
    nFeatures=features.size();
    std::sort(features.begin(), features.end());
    for (uint32 ii=0; ii<nFeatures; ii++)
        featureInd[features[ii]]=ii;

    if (nFeatures != featureIn.size()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option(s) in --soloFeatures ";
        for (auto &fname : featureIn)
                errOut << fname <<" ";
        errOut << "\nSOLUTION: use allowed option: ";
        errOut <<SoloFeatureTypes::Names[0];
        for (auto &fname : SoloFeatureTypes::Names)
            errOut <<"   OR   "<< fname;
        exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };

    if (featureYes[SoloFeatureTypes::Velocyto] && type==SoloTypes::SmartSeq) {
        string errOut = "EXITING because of fatal PARAMETERS error: --soloFeatures Velocyto is presently not compatible with --soloType SmartSeq .\n";
        errOut       += "SOLUTION: re-run without --soloFeatures Velocyto .";
        exitWithError(errOut, std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
    
    if (featureYes[SoloFeatureTypes::Gene]) {
        pP->quant.gene.yes = true;
        pP->quant.yes = true;
    };
    if (featureYes[SoloFeatureTypes::GeneFull]) {
        pP->quant.geneFull.yes = true;
        pP->quant.yes = true;
        if (!featureYes[SoloFeatureTypes::Gene]) 
            pP->quant.gene.yes=false; //if GeneFull is requested, but Gene is not, turn it off - it could have been turned on because of GX/GN attributes
    };
    if (featureYes[SoloFeatureTypes::GeneFull_Ex50pAS]) {
        pP->quant.geneFull_Ex50pAS.yes = true;
        pP->quant.yes = true;
        if (!featureYes[SoloFeatureTypes::Gene]) 
            pP->quant.gene.yes=false; //if GeneFull is requested, but Gene is not, turn it off - it could have been turned on because of GX/GN attributes
    };
    if (featureYes[SoloFeatureTypes::GeneFull_ExonOverIntron]) {
        pP->quant.geneFull_ExonOverIntron.yes = true;
        pP->quant.gene.yes = true; //needed to prioritize exons over introns
        pP->quant.yes = true;
    };          
    
    //initialize CB match to WL types
    init_CBmatchWL();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////umiDedup
    umiDedup.initialize(this);
    
    ///////////// finished parameters input
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //make output directory if needed
    if ( outFileNames[0].find_last_of("/") < outFileNames[0].size() ) {//need to create dir
        string dir1=pP->outFileNamePrefix+outFileNames[0].substr(0,outFileNames[0].find_last_of("/"));
        if (mkdir(dir1.c_str(),pP->runDirPerm)!=0 && errno!=EEXIST) {
            ostringstream errOut;
            errOut << "EXITING because of fatal OUTPUT FILE error: could not create Solo output directory"<<dir1<<"\n";
            errOut << "SOLUTION: check the path and permisssions";
            exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };
    };

    QSbase=33;//TODO make these user-definable
    QSmax=33;
    cbMinP=0.975;

    umiMaskLow=(uint32) ( (((uint64)1)<<umiL) - 1);
    umiMaskHigh=~umiMaskLow;

    //////////////////////////////////////////////////////CB whitelist
    if (type==SoloTypes::CB_UMI_Simple || type==SoloTypes::CB_samTagOut) {//simple whitelist
        if (soloCBwhitelist.size()>1) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in INPUT parameters: --soloCBwhitelist contains more than one file which is not allowed with --soloType CB_UMI_Simple \n";
            errOut << "SOLUTION: in --soloCBwhitelist specify only one whitelist file \n";
            exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_INPUT_FILES, *pP);
        } else if (soloCBwhitelist[0]=="-") {
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR in INPUT parameters: --soloCBwhitelist is not defined\n";
                errOut << "SOLUTION: in --soloCBwhitelist specify path to and name of the whitelist file, or None for CB demultiplexing without whitelist \n";
                exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_INPUT_FILES, *pP);
        } else if (soloCBwhitelist[0]=="None") {
            cbWLyes=false;
        } else {
            cbWLyes=true;
            ifstream & cbWlStream = ifstrOpen(soloCBwhitelist[0], ERROR_OUT, "SOLUTION: check the path and permissions of the CB whitelist file: " + soloCBwhitelist[0], *pP);
            string seq1;
            while (cbWlStream >> seq1) {
                if (seq1.size() != cbL) {
                    ostringstream errOut;
                    errOut << "EXITING because of FATAL ERROR in input CB whitelist file: "<< soloCBwhitelist[0] <<" the total length of barcode sequence is "  << seq1.size() << " not equal to expected " <<bL <<"\n"  ;
                    errOut << "SOLUTION: make sure that the barcode read is the second in --readFilesIn and check that is has the correct formatting\n";
                    exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_INPUT_FILES, *pP);
                };
                uint64 cb1;
                if (convertNuclStrToInt64(seq1,cb1)) {//convert to 2-bit format
                    cbWL.push_back(cb1);
                } else {
                    pP->inOut->logMain << "WARNING: CB whitelist sequence contains non-ACGT base and is ignored: " << seq1 <<endl;
                };
            };
            if (cbWL.size()==0) {//empty whitelist
                exitWithError("EXITING because of FATAL ERROR: CB whitelist file " + soloCBwhitelist[0] + \
                               " is empty. \nSOLUTION: provide non-empty whitelist.\n" , \
                               std::cerr, pP->inOut->logMain, EXIT_CODE_INPUT_FILES, *pP);
            };
        };

        std::sort(cbWL.begin(),cbWL.end());//sort
        auto un1=std::unique(cbWL.begin(),cbWL.end());//collapse identical
        cbWL.resize(std::distance(cbWL.begin(),un1));        
        cbWLsize=cbWL.size();
        pP->inOut->logMain << "Number of CBs in the whitelist = " << cbWLsize <<endl;
        
        cbWLstr.resize(cbWLsize);
        for (uint64 ii=0; ii<cbWLsize; ii++)
             cbWLstr[ii] = convertNuclInt64toString(cbWL[ii],cbL);        
        
    //////////////////////////////////////////////////////////////////////////////////
    } else if (type==SoloTypes::SmartSeq) {
        cbWLstr=pP->outSAMattrRG;
        cbWLsize=cbWLstr.size();
        cbWLyes=true; 
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    } else if (type==SoloTypes::CB_UMI_Complex) {//complex barcodes: multiple whitelist (one for each CB), varying CB length
        cbWLyes=true; //for complex barcodes, no-whitelist option is not allowed for now
        
        adapterYes=false;
        if (adapterSeq!="-")
            adapterYes=true;
        
        if (cbPositionStr.size() != soloCBwhitelist.size()) {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETER error: number of barcodes in --soloCBposition : "<< cbPositionStr.size() <<" is not equal to the number of WhiteLists in --soloCBwhitelist : " << soloCBwhitelist.size() <<"\n"  ;
            errOut << "SOLUTION: make sure that the number of CB whitelists and CB positions are the same\n";
            exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_INPUT_FILES, *pP);
        };
        cbV.resize(cbPositionStr.size());
        for (uint32 ii=0; ii<cbPositionStr.size(); ii++) {
            cbV[ii].extractPositionsFromString(cbPositionStr[ii]);
        };
        
        umiV.extractPositionsFromString(umiPositionStr);
        umiL = 0; //this will be defined when the first barcode is processed
              
        umiV.adapterLength=adapterSeq.size();//one adapter for all
        cbWLsize=1;
        for (uint32 icb=0; icb<cbV.size(); icb++) {//cycle over WL files
            cbV[icb].adapterLength=adapterSeq.size();//one adapter for all
            
            ifstream & cbWlStream = ifstrOpen(soloCBwhitelist[icb], ERROR_OUT, "SOLUTION: check the path and permissions of the CB whitelist file: " + soloCBwhitelist[icb], *pP);
            
            string seq1;
            while (cbWlStream >> seq1) {//cycle over one WL file
                uint64 cb1;
                if (!convertNuclStrToInt64(seq1,cb1)) {//convert to 2-bit format
                    pP->inOut->logMain << "WARNING: CB whitelist sequence contains non-ACGT base and is ignored: " << seq1 <<endl;
                    continue;
                };
                
                uint32 len1=seq1.size();
                if (len1>=cbV[icb].wl.size())
                    cbV[icb].wl.resize(len1+1);//add new possible lengths to this CB
                cbV[icb].wl.at(len1).push_back(cb1);
            };
            
            cbV[icb].sortWhiteList(this);
            cbV[icb].wlFactor=cbWLsize;
            cbWLsize *= cbV[icb].totalSize;
        };
        
        complexWLstrings();
    };

    time_t rawTime;
    time(&rawTime);
    pP->inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished reading, sorting and deduplicating CB whitelist sequences." <<endl;

    //////////////////////////////////////////////////////////////SAM attributes
    samAttrYes=false;
    if ( (pP->outSAMattrPresent.CB || pP->outSAMattrPresent.UB) && type!=SoloTypes::CB_samTagOut) {
        samAttrYes=true;
        if (!pP->outBAMcoord) {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: CB and/or UB attributes in --outSAMattributes can only be output in the sorted BAM file.\n";
            errOut << "SOLUTION: re-run STAR with --outSAMtype BAM SortedByCoordinate ...\n";
            exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };
    } else if ( pP->outSAMattrPresent.UB && type==SoloTypes::CB_samTagOut) {
        exitWithError("EXITING because of fatal PARAMETERS error: UB attribute (corrected UMI) in --outSAMattributes cannot be used with --soloType CB_samTagOut \n" \
                      "SOLUTION: instead, use UR (uncorrected UMI) in --outSAMattributes\n",   std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
    
    ////////////////////////////////////////////////////////////////readInfoYes: which feature is used to fill readInfo. Only one feature is allowed
    readInfoYes.fill(false);
    if (featureYes[SoloFeatureTypes::VelocytoSimple] || featureYes[SoloFeatureTypes::Velocyto]) {//turn readInfo on for Gene needed by VelocytoSimple
        readInfoYes[SoloFeatureTypes::Gene]=true;
    };
    samAttrFeature = featureFirst;
    if (samAttrYes){//pSolo.samAttrFeature=0 by default, so need to check samAttrYes
        if (   featureFirst == SoloFeatureTypes::Gene || featureFirst == SoloFeatureTypes::GeneFull ||
               featureFirst == SoloFeatureTypes::GeneFull_Ex50pAS || featureFirst == SoloFeatureTypes::GeneFull_ExonOverIntron ) {
            //all good
        } else {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: CB and/or UB attributes in --outSAMattributes require --soloFeatures Gene OR/AND GeneFull OR/AND GeneFull_Ex50pAS.\n";
            errOut << "SOLUTION: re-run STAR adding Gene AND/OR GeneFull OR/AND GeneFull_Ex50pAS OR/AND GeneFull_ExonOverIntron to --soloFeatures\n";
            exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };
        readInfoYes[samAttrFeature]=true;
    };
    readIndexYes = readInfoYes;

    /////////////////////////////////////////////////////////////////// readFlag output
    readStats.yes = false;
    readStatsYes.fill(false);
    if (readStats.type == "Standard") {
        readStats.yes = true;
        readStatsYes.fill(true);
        readStatsYes[SoloFeatureTypes::VelocytoSimple] = false; //this could be allowed, but it will have the same info as Gene
        readStatsYes[SoloFeatureTypes::Velocyto] = false;
        readStatsYes[SoloFeatureTypes::SJ] = false; //output for SJ requires careful consideration for SoloReadFeature_record.cpp, for output of reads that do not have splice junction
        for (uint32 ff=0; ff<readIndexYes.size(); ff++) {//merge with previous values
            readIndexYes[ff] |= readStatsYes[ff];
        };
    } else if (readStats.type != "None" ) {
        exitWithError("EXITING because of fatal PARAMETERS error: unrecognized option in --soloCellReadStats" + readStats.type
                      + "\nSOLUTION: use allowed options: None OR Standard \n",
                      std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
       
    ///////////////////////////////////////////////////////////////////umi filtering
    if (umiFiltering.type[0]=="MultiGeneUMI") {
        umiFiltering.MultiGeneUMI = true;
        umiFiltering.yes = true;
    } else if (umiFiltering.type[0]=="MultiGeneUMI_All") {
        umiFiltering.MultiGeneUMI_All = true;
        umiFiltering.yes = true;
    } else if (umiFiltering.type[0]=="MultiGeneUMI_CR") {
        umiFiltering.MultiGeneUMI_CR = true;
        if (umiDedup.typesIn.size()>1 || umiDedup.typesIn.at(0) != "1MM_CR")
            exitWithError("EXITING because of fatal PARAMETERS error: --soloUMIfiltering MultiGeneUMI_CR only works with --soloUMIdedup 1MM_CR"
                            "\nSOLUTION: rerun with --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \n",
                            std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    } else if (umiFiltering.type[0]=="-") {
        //nothing to do
    } else {
        exitWithError("EXITING because of fatal PARAMETERS error: unrecognized option in --soloUMIfiltering=" + umiFiltering.type[0] 
                      + "\nSOLUTION: use allowed options: - or MultiGeneUMI or MultiGeneUMI_CR \n",
                      std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
    

    
    //////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////// MultiMappers
    multiMap.initialize(this);
    if (multiMap.yes.multi) {
        if (type==SoloTypes::CB_samTagOut || type==SoloTypes::SmartSeq) {
            exitWithError("EXITING because of fatal PARAMETERS error: multimapping options do not work for --soloType " +typeStr+ 
                          "\nSOLUTION: use default option --soloMultiMappers Unique\n",
                            std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };
        readIndexYes[SoloFeatureTypes::Gene]=true;
        readIndexYes[SoloFeatureTypes::GeneFull]=true;
        readIndexYes[SoloFeatureTypes::GeneFull_Ex50pAS]=true;
        readIndexYes[SoloFeatureTypes::GeneFull_ExonOverIntron]=true;
    };
};

/////////////////////////////////
void ParametersSolo::umiSwapHalves(uint32 &umi) {
    uint32 high=umi>>(umiL);
    umi &= umiMaskLow; //remove high
    umi <<= (umiL); //move low to high
    umi |= high; //add high
};

void ParametersSolo::complexWLstrings() {
    
    cbWLstr.resize(cbWLsize);
    
    for (auto &cb : cbV) {//initialize
        cb.iCB=0;
        cb.iLen=cb.minLen;
    };

    for (uint32 ii=0; ii<cbWLsize; ii++) {//cycle over full WL
        for (uint64 ii=0; ii<cbV.size(); ii++) {//check for overflow and re-calculate all indexes
            SoloBarcode &cb=cbV[ii];
            if (cb.iCB == cb.wl[cb.iLen].size()) {//advance length
                cb.iLen++;
                cb.iCB=0;//reset iCB
            };
            if (cb.iLen == cb.wl.size()){//advance barcode
                cbV[ii+1].iCB++;
                cb.iLen=cb.minLen;//reset length
            };
        };
        
        for (auto &cb : cbV) 
            cbWLstr[ii] += convertNuclInt64toString(cb.wl[cb.iLen][cb.iCB], cb.iLen) + "_";
        cbWLstr[ii].pop_back();
        
        cbV[0].iCB++;//shift by one for the next CB
    };
};

void ParametersSolo::cellFiltering()
{//cell filtering
    string pars1;
    for (auto s=cellFilter.type.begin()+1; s!=cellFilter.type.end(); s++)
        pars1 += ' ' + *s; //concatenate parameters into one string - easier to process that way
    
    if (cellFilter.type[0]=="CellRanger2.2") {           
        if (cellFilter.type.size()==1) {
            pP->inOut->logMain << "ParametersSolo: using hardcoded filtering parameters for --soloCellFilterType CellRanger2.2" <<endl;
            pars1="3000 0.99 10";
        } else if (cellFilter.type.size()<4) {            
            string errOut="EXITING because of fatal PARAMETERS error: --soloCellFilterType CellRanger2.2 requires exactly 3 numerical parameters";
            errOut +=     "\nSOLUTION: re-run with --soloCellFilterType CellRanger2.2 <nExpectedCells> <maxPercentile> <maxMinRatio>\n";
            exitWithError(errOut, std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };
        
        pP->inOut->logMain << "ParametersSolo: --soloCellFilterType CellRanger2.2 filtering parameters: " << pars1 <<endl;
        istringstream parsStream(pars1);
        parsStream >> cellFilter.knee.nExpectedCells >> cellFilter.knee.maxPercentile >> cellFilter.knee.maxMinRatio;

    } else if (cellFilter.type[0]=="EmptyDrops_CR") {
        if (cellFilter.type.size()==1) {
            pP->inOut->logMain << "ParametersSolo: using hardcoded filtering parameters for --soloCellFilterType EmptyDrops_CR\n";
            pars1="3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000";
        } else if (cellFilter.type.size()<11) {            
            string errOut="EXITING because of fatal PARAMETERS error: --soloCellFilterType EmptyDrops_CR requires exactly 10 numerical parameters";
            errOut +=     "\nSOLUTION: re-run with --soloCellFilterType EmptyDrops_CR ";
            errOut +=     "<nExpectedCells> <maxPercentile> <maxMinRatio> <indMin> <indMax> <umiMin> <umiMinFracMedian> <candMaxN> <FDR> <simN>\n";
            exitWithError(errOut, std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };

        pP->inOut->logMain << "ParametersSolo: --soloCellFilterType EmptyDrops_CR filtering parameters: " << pars1 <<endl;
        istringstream parsStream(pars1);
        parsStream >> cellFilter.knee.nExpectedCells >> cellFilter.knee.maxPercentile >> cellFilter.knee.maxMinRatio;
        parsStream >> cellFilter.eDcr.indMin >> cellFilter.eDcr.indMax >> cellFilter.eDcr.umiMin >> cellFilter.eDcr.umiMinFracMedian;
        parsStream >> cellFilter.eDcr.candMaxN >> cellFilter.eDcr.FDR >> cellFilter.eDcr.simN;
        
    } else if (cellFilter.type[0]=="TopCells") {
        if (cellFilter.type.size()<2) {
            string errOut="EXITING because of fatal PARAMETERS error: number of cells not specified for --soloCellFilterType TopCells";
            errOut +=     "\nSOLUTION: --soloCellFilterType TopCells <NumberOfCells>\n";
            exitWithError(errOut, std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };
        cellFilter.topCells=stoi(cellFilter.type[1]);
    } else if (cellFilter.type[0]=="None") {
        //nothing to do
    } else {
        exitWithError("EXITING because of fatal PARAMETERS error: unrecognized option in --soloCellFilterType=" + cellFilter.type[0] + "\nSOLUTION: use allowed options: CellRanger2.2 or None\n",
                    std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
};

void UMIdedup::initialize(ParametersSolo *pS)
{
    yes.N = 0;
    countInd.I.fill((uint32_t) -1); //marks types not used
    yes.B.fill(false);

    for (uint32_t iin=0; iin<typesIn.size(); iin++) {
        uint32_t itype;
        for (itype=0; itype<tN; itype++) {
            if (typesIn[iin] == typeNames[itype])
                break; //found match
        };
        
        if (itype==tN) {//no match
            std::string tall;
            for (auto &t: typeNames)
                tall +=" " + t; // concatenate allowed values
            
            exitWithError("EXITING because of fatal PARAMETERS error: unrecognzied option --soloUMIdedup = " + typesIn[iin] + '\n'
                          + "SOLUTION: use allowed values: " + tall + '\n'
                          ,std::cerr, pS->pP->inOut->logMain, EXIT_CODE_PARAMETER, *pS->pP);
        };
        
        types.push_back(itype);
        yes.B[itype] = true;
        yes.N++;
        countInd.I[itype] = iin + 1; //for each type, which column itype's recorded in
            
        if (pS->type == pS->SoloTypes::SmartSeq && (yes.All || yes.Directional || yes.CR) )
            exitWithError("EXITING because of fatal PARAMETERS error: --soloUMIdedup = " + typesIn[iin] + " is not allowed for --soloType SmartSeq\n"
                    + "SOLUTION: use allowed options: Exact and/or NoDedup\n"
                    ,std::cerr, pS->pP->inOut->logMain, EXIT_CODE_PARAMETER, *pS->pP);            
    };
       
    //hard-coded for now
    typeMain = types[0]; //main is the 0th entry in typesIn
    countInd.main = 1;   //hard-coded - 1 column is always main
};

void MultiMappers::initialize(ParametersSolo* pS)
{
    yes.N = 0;
    countInd.I.fill((uint32) -1); //marks types not used
    yes.B.fill(false);

    for (uint32 iin=0; iin<typesIn.size(); iin++) {
        uint32 itype;
        for (itype=0; itype<tN; itype++) {
            if (typesIn[iin] == typeNames[itype])
                break; //found match
        };
        
        if (itype==tN) {//no match
            std::string tall;
            for (auto &t: typeNames)
                tall +=" " + t; // concatenate allowed values
            
            exitWithError("EXITING because of fatal PARAMETERS error: unrecognzied option --soloMultiMappers = " + typesIn[iin] + '\n'
                          + "SOLUTION: use allowed values: " + tall + '\n'
                          ,std::cerr, pS->pP->inOut->logMain, EXIT_CODE_PARAMETER, *pS->pP);
        };
        
        if (itype == typeI::Unique)
            continue; //Unique type does not have to be recorded TODO use it to perform filtering and stats

        
        types.push_back(itype);
        yes.B[itype] = true;
        yes.N++;
    };
    
    if (yes.N==0) {//only Unique, no multimappers
        yes.multi=false;
        return;
    };
    
    uint32 ind1=1; //start
    for (const auto &itype : types) {
        countInd.I[itype] = ind1;
        ind1 += pS->umiDedup.yes.N;
    };
    
    //hard-coded for now
    typeMain = types[0]; //main is the 0th entry in typesIn
    countInd.main = 1;   //hard-coded - 1 column is always main
    
    yes.multi = yes.Uniform | yes.Rescue | yes.PropUnique | yes.EM;
};

////////////////////////////////////////////////////
void ParametersSolo::init_CBmatchWL()
{//CBmatchWL
    bool incomp1 =        typeStr=="CB_UMI_Complex" && (CBmatchWL.type!="Exact" && CBmatchWL.type!="1MM" && CBmatchWL.type!="EditDist_2");
    incomp1 = incomp1 || (typeStr=="CB_samTagOut"   && (CBmatchWL.type!="Exact" && CBmatchWL.type!="1MM"));
    incomp1 = incomp1 || (typeStr!="CB_UMI_Complex" &&  CBmatchWL.type=="EditDist_2");
    
    if ( incomp1 ) {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: --soloCBmatchWLtype "<< CBmatchWL.type << " does not work with --soloType " << typeStr << "\n";
        errOut << "SOLUTION: use allowed option: use --soloCBmatchWLtype Exact (exact matches only) OR 1MM (one match with 1 mismatched base)\n";
        exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
    
    
    CBmatchWL.mm1 = false;
    CBmatchWL.mm1_multi = false;
    CBmatchWL.mm1_multi_pc = false;        
    CBmatchWL.mm1_multi_Nbase = false;
    CBmatchWL.oneExact = false; //if true, for a CB matching with 1 mismatch to a WL-CB, requires at least one other read to match this WL-CB. 
                                //this is true for all options except pseudocount
    CBmatchWL.EditDist_2 = false;
        
    if (CBmatchWL.type=="Exact") {    
        CBmatchWL.oneExact=true;
    } else if (CBmatchWL.type=="1MM") {
        CBmatchWL.mm1=true;
        CBmatchWL.oneExact=true;        
    } else if (CBmatchWL.type=="1MM_multi") {
        CBmatchWL.mm1=true;
        CBmatchWL.mm1_multi=true;
        CBmatchWL.oneExact=true;    
    } else if (CBmatchWL.type=="1MM_multi_pseudocounts") {
        CBmatchWL.mm1=true;
        CBmatchWL.mm1_multi=true;
        CBmatchWL.mm1_multi_pc=true;
    } else if (CBmatchWL.type=="1MM_multi_Nbase_pseudocounts") {
        CBmatchWL.mm1=true;
        CBmatchWL.mm1_multi=true;
        CBmatchWL.mm1_multi_pc=true;
        CBmatchWL.mm1_multi_Nbase = true;
    } else if (CBmatchWL.type=="EditDist_2") {
        CBmatchWL.EditDist_2=true;
    } else {
        exitWithError("EXITING because of fatal PARAMETERS error: unrecognized option in --soloCBmatchWLtype " +CBmatchWL.type + "\nSOLUTION: use allowed options: Exact or 1MM or 1MM_multi or 1MM_multi_pseudocounts 1MM_multi_Nbase_pseudocounts\n",
                        std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
};
