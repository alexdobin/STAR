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

    //constants - may turn into parameters in the future
    redistrReadsNfiles = 3*pP->runThreadN;
    
    barcodeReadYes=false; //will select true later if needed
    barcodeStart=barcodeEnd=0; //this means that the entire barcodeRead is considered barcode. Will change it for simple barcodes.

    yes = true;
    if (typeStr=="None") {
        type = SoloTypes::None;
        yes = false;
        samAttrYes = false;
        //solo SAM attributes not allowed
        if ( pP->outSAMattrPresent.CR || pP->outSAMattrPresent.CY || pP->outSAMattrPresent.UR 
          || pP->outSAMattrPresent.UY  || pP->outSAMattrPresent.CB  || pP->outSAMattrPresent.UB
          || pP->outSAMattrPresent.sS  || pP->outSAMattrPresent.sQ  || pP->outSAMattrPresent.sM
           ) {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL INPUT ERROR: --outSAMattributes contains CR/CY/UR/UY/CB/UB tags, but --soloType is not set\n";
            errOut <<"SOLUTION: re-run STAR without these attribures, or with --soloType set\n";
            exitWithError(errOut.str(), std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };        
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
        bL=0;
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

    if (type==SoloTypes::SmartSeq) {//no barcode read
        pP->readNmates=pP->readNmatesIn; 
        barcodeRead=-1; 
        barcodeReadYes=false;
        //TODO: a lot of parameters should not be defined for SmartSeq option - check it here
        //umiDedup
    } else {//all other types require barcode read
        if (pP->readNmatesIn<2) {
            exitWithError("solo* options require 2 reads or 3 reads, where the last read is the barcode read.\n"
                          "SOLUTION: specify 2 or 3 files in --readFilesIn",std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };
        pP->readNmates=pP->readNmatesIn-1; //output mates    
        barcodeRead=pP->readNmates;//TODO make this flexible
        barcodeReadYes=true;
    };
    
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

    //////////////////////////////////// features
    featureInd.fill(-1);
    featureYes.fill(false);
    for (uint32 ii=0; ii<SoloFeatureTypes::Names.size(); ii++) {
        for (auto &fin : featureIn) {
            if (fin==SoloFeatureTypes::Names[ii]) {
                featureYes[ii]=true;
                features.push_back(ii);
                featureInd[ii]=features.size()-1;
                break;
            };
        };
    };
    nFeatures=features.size();

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

    if (featureYes[SoloFeatureTypes::Gene]) {
        pP->quant.gene.yes=true;
        pP->quant.yes = true;
    };
    if (featureYes[SoloFeatureTypes::GeneFull]) {
        pP->quant.geneFull.yes=true;
        pP->quant.yes = true;
    };
    
    ////////////////////////////////////////////umiDedup
    umiDedupColumns.resize(umiDedup.size());
    if (type==SoloTypes::SmartSeq) {
        for (uint32 ii=0; ii<umiDedup.size(); ii++) {
            if (umiDedup[ii]=="Exact") {
                umiDedupColumns[ii]=0;
            } else if (umiDedup[ii]=="NoDedup") {
                if (type!=SoloTypes::SmartSeq) {
                    ostringstream errOut;
                    errOut << "EXITING because of fatal PARAMETERS error: --soloUMIdedup NoDedup is allowed only for --soloType SmartSeq\n";
                    errOut << "SOLUTION: used allowed options: Exact or 1MM_All or 1MM_Directional";
                    exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
                };
                umiDedupColumns[ii]=1;                
            } else {
                ostringstream errOut;
                errOut << "EXITING because of fatal PARAMETERS error: --soloUMIdedup="<<umiDedup[ii] <<" is not allowed for --soloType SmartSeq\n";
                errOut << "SOLUTION: used allowed options: Exact or NoDedup";
                exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
            };
        };        
    } else {
        for (uint32 ii=0; ii<umiDedup.size(); ii++) {
            if (umiDedup[ii]=="Exact") {
                umiDedupColumns[ii]=0;
            } else if (umiDedup[ii]=="1MM_All") {
                umiDedupColumns[ii]=1;
            } else if (umiDedup[ii]=="1MM_Directional") {
                umiDedupColumns[ii]=2;                
            } else {
                ostringstream errOut;
                errOut << "EXITING because of fatal PARAMETERS error: --soloUMIdedup="<<umiDedup[ii] <<" is not allowed for --soloType " << typeStr <<"\n";
                errOut << "SOLUTION: used allowed options: Exact or 1MM_All or 1MM_Directional";
                exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
            };
        };
    };
    ///////////// finished parameters input

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
        
        /* debug
        cbV[0].length=16;
        cbV[0].anchorType=0;
        cbV[0].anchorDist=0;
        cbV[0].pos=0;
        
        umiV.length=10;
        umiV.anchorType=0;
        umiV.anchorDist=0;
        umiV.pos=16;

        //hard-coded inDrop
        cbV.resize(2);
        cbV[0].anchorType={0,2};
        cbV[0].anchorDist={0,-1};//CB ends one base before the start of the anchor
        cbV[1].anchorType={3,3};
        cbV[1].anchorDist={1,8};//CB starts 1 base after the end of the anchor  
        umiV.anchorType={3,3};
        umiV.anchorDist={9,14};
        ///////////////////
        */
        
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
            
            cbV[icb].sortWhiteList();
            cbV[icb].wlFactor=cbWLsize;
            cbWLsize *= cbV[icb].totalSize;
        };
        
        complexWLstrings();
    };

    time_t rawTime;
    time(&rawTime);
    pP->inOut->logMain << timeMonthDayTime(rawTime) << "Finished reading, sorting and deduplicating CB whitelist sequences." <<endl;

    //SAM attributes
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
    
    readInfoYes.fill(false);
    if (samAttrYes) {//pSolo.samAttrFeature=0 by default, so need to check samAttrYes
    	if (featureYes[SoloFeatureTypes::Gene]) {
    		samAttrFeature=SoloFeatureTypes::Gene;
    	} else if (featureYes[SoloFeatureTypes::GeneFull]) {//if Gene is not defined
    		samAttrFeature=SoloFeatureTypes::GeneFull;
    	} else {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: CB and/or UB attributes in --outSAMattributes require --soloFeatures Gene OR/AND GeneFull.\n";
            errOut << "SOLUTION: re-run STAR adding Gene AND/OR GeneFull to --soloFeatures\n";
            exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    	};

        readInfoYes[samAttrFeature]=true;
    };
    if (featureYes[SoloFeatureTypes::VelocytoSimple] || featureYes[SoloFeatureTypes::Velocyto]) //turn readInfo on for Gene needed by VelocytoSimple
        readInfoYes[SoloFeatureTypes::Gene]=true;    
    
    //cell filtering
    if (cellFilter.type[0]=="CellRanger2.2") {
        cellFilter.cr2expectedCells=stod(cellFilter.type[1]);
        cellFilter.cr2maxPercentile=stod(cellFilter.type[2]);
        cellFilter.cr2maxMinRatio=stod(cellFilter.type[3]);
        cellFilter.cr2maxCellInd=(uint64) (cellFilter.cr2expectedCells*(1.0-cellFilter.cr2maxPercentile)+0.5);//cell number for robust max
    } else if (cellFilter.type[0]=="TopCells") {
        cellFilter.topCells=stoi(cellFilter.type[1]);
    } else if (cellFilter.type[0]=="None") {
    } else {
        exitWithError("EXITING because of fatal PARAMETERS error: unrecognized option in --soloCellFilterType=" + cellFilter.type[0] + "\nSOLUTION: use allowed options: CellRanger2.2 or None\n",
                       std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
    
    //umi filtering
    umiFiltering.MultiGeneUMI=false;
    if (umiFiltering.type[0]=="MultiGeneUMI") {
        umiFiltering.MultiGeneUMI=true;
    } else if (umiFiltering.type[0]=="-") {
        //nothing to do
    } else {
        exitWithError("EXITING because of fatal PARAMETERS error: unrecognized option in --soloUMIfiltering=" + umiFiltering.type[0] + "\nSOLUTION: use allowed options: - OR MultiGeneUMI\n",
                      std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
    
    //CBmatchWL
    if ( (typeStr=="CB_UMI_Complex" || typeStr=="CB_samTagOut") && (CBmatchWL.type!="Exact" && CBmatchWL.type!="1MM") ) {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: --soloCBmatchWLtype "<< CBmatchWL.type << " does not work with --soloType " << typeStr << "\n";
        errOut << "SOLUTION: use allowed option: use --soloCBmatchWLtype Exact (exact matches only) OR 1MM (one match with 1 mismatched base)\n";
        exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
    
    if (CBmatchWL.type=="Exact") {
        CBmatchWL.mm1=false;
        CBmatchWL.mm1_multi=false;
        CBmatchWL.mm1_multi_pc=false;        
        CBmatchWL.oneExact=true;
    } else if (CBmatchWL.type=="1MM") {
        CBmatchWL.mm1=true;
        CBmatchWL.mm1_multi=false;
        CBmatchWL.mm1_multi_pc=false;
        CBmatchWL.oneExact=true;        
    } else if (CBmatchWL.type=="1MM_multi") {
        CBmatchWL.mm1=true;
        CBmatchWL.mm1_multi=true;
        CBmatchWL.mm1_multi_pc=false;
        CBmatchWL.oneExact=true;    
    } else if (CBmatchWL.type=="1MM_multi_pseudocounts") {
        CBmatchWL.mm1=true;
        CBmatchWL.mm1_multi=true;
        CBmatchWL.mm1_multi_pc=true;
        CBmatchWL.oneExact=false;    
    } else {
        exitWithError("EXITING because of fatal PARAMETERS error: unrecognized option in --soloCBmatchWL.type " +CBmatchWL.type + "\nSOLUTION: use allowed options: Exact -or- 1MM -or- 1MM_multi -or- 1MM_multi_pseudocounts\n",
                      std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
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
        cb.icb=0;
        cb.ilen=cb.minLen;
    };

    for (uint32 ii=0; ii<cbWLsize; ii++) {//cycle over full WL
        for (uint64 ii=0; ii<cbV.size(); ii++) {//check for overflow and re-calculate all indexes
            SoloBarcode &cb=cbV[ii];
            if (cb.icb == cb.wl[cb.ilen].size()) {//advance length
                cb.ilen++;
                cb.icb=0;//reset icb
            };
            if (cb.ilen == cb.wl.size()){//advance barcode
                cbV[ii+1].icb++;
                cb.ilen=cb.minLen;//reset length
            };
        };
        
        for (auto &cb : cbV) 
            cbWLstr[ii] += convertNuclInt64toString(cb.wl[cb.ilen][cb.icb], cb.ilen) + "_";
        cbWLstr[ii].pop_back();
        
        cbV[0].icb++;//shift by one for the next CB
    };
};
