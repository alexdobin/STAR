#include "ParametersSolo.h"
#include "Parameters.h"
#include "ErrorWarning.h"
#include "streamFuns.h"
#include "SequenceFuns.h"
#include "serviceFuns.cpp"

#include <stdlib.h>

const vector<string> ParametersSolo::featureNames={"Gene","SJ","GeneFull","Transcript3p"};

void ParametersSolo::initialize(Parameters *pPin)
{
    pP=pPin;

    if (typeStr=="None") {
        type=0;
        //solo SAM attributes not allowed
        if (pP->outSAMattrPresent.CR || pP->outSAMattrPresent.CY || pP->outSAMattrPresent.UR || pP->outSAMattrPresent.UY  || pP->outSAMattrPresent.CB  || pP->outSAMattrPresent.UB) {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL INPUT ERROR: --outSAMattributes contains CR/CY/UR/UY tags, but --soloType is not set\n";
            errOut <<"SOLUTION: re-run STAR without these attribures, or with --soloType set\n";
            exitWithError(errOut.str(), std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };        
        return;
    } else if (typeStr=="CB_UMI_Simple") {
        type=1;        
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
        if (bL==1)
            bL=cbL+umiL;
        pP->readNmates=1; //output mates TODO: check that readNmatesIn==2       
    } else if (typeStr=="CB_UMI_Complex") {
        type=2;
    } else  {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in --soloType="<<typeStr<<"\n";
        errOut << "SOLUTION: use allowed option: None OR CB_UMI_Simple OR CB_UMI_Complex\n";
        errOut << "Obsolete option Droplet should be replaced with CB_UMI_Simple";
        exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };

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

    featureYes=new bool [featureNames.size()];
    featureInd.resize(featureNames.size(),-1);
    for (auto &fin : featureIn) {
        bool finGood=false;
        for (uint32 ii=0; ii<featureNames.size(); ii++) {
            if (fin==featureNames[ii]) {
                finGood=true;
                featureYes[ii]=true;
                features.push_back(ii);
                featureInd[ii]=features.size()-1;
                break;
            };
        };
        if (!finGood) {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in --soloFeatures="<<fin<<"\n";
            errOut << "SOLUTION: use allowed option: ";
            errOut <<featureNames[0]<< "   OR   ";
            for (auto &fname : featureNames)
                errOut << fname <<" ";
            exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };
    };

    nFeatures=features.size();

    umiDedupYes.resize(3,false);
    umiDedupColumns.resize(umiDedup.size());
    for (uint32 ii=0; ii<umiDedup.size(); ii++) {
        if (umiDedup[ii]=="1MM_NotCollapsed") {
            umiDedupYes[0]=true;
            umiDedupColumns[ii]=0;
        } else if (umiDedup[ii]=="1MM_All") {
            umiDedupYes[1]=true;
            umiDedupColumns[ii]=1;
        } else if (umiDedup[ii]=="1MM_Directional") {
            umiDedupYes[2]=true;
            umiDedupColumns[ii]=2;
        } else {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in --soloUMIdedup="<<umiDedup[ii]<<"\n";
            exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
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

    //load the CB whitelist
    if (type==1) {//simple whitelist
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
        };

        std::sort(cbWL.begin(),cbWL.end());//sort
        auto un1=std::unique(cbWL.begin(),cbWL.end());//collapse identical
        cbWL.resize(std::distance(cbWL.begin(),un1));        
        
        pP->inOut->logMain << "Number of CBs in the whitelist = " << cbWL.size() <<endl;
        
    } else if (type==2) {//complex barcodes: multiple whitelist (one for each CB), varying CB length
        cbWLyes=true; //for complex barcodes, no-whitelist option is not allowed for now
        cbWLv.resize(soloCBwhitelist.size());
        for (uint32 icb=0; icb<soloCBwhitelist.size(); icb++) {//cycle over WL files
            ifstream & cbWlStream = ifstrOpen(soloCBwhitelist[icb], ERROR_OUT, "SOLUTION: check the path and permissions of the CB whitelist file: " + soloCBwhitelist[icb], *pP);
            
            string seq1;
            while (cbWlStream >> seq1) {//cycle over one WL file
                uint64 cb1;
                if (!convertNuclStrToInt64(seq1,cb1)) {//convert to 2-bit format
                    pP->inOut->logMain << "WARNING: CB whitelist sequence contains non-ACGT base and is ignored: " << seq1 <<endl;
                    
                };
                
                uint32 len1=seq1.size();
                if (len1>cbWLv.back().size())
                    cbWLv[icb].resize(len1);//add new possible lengths to this CB
                cbWLv[icb].at(len1).push_back(cb1);
            };
            
            uint32 lenMin=(uint32)-1;
            for (uint32 ilen=1; ilen<cbWLv[icb].size(); ilen++) {//scan through different lengths for this CB
                if (cbWLv[icb][ilen].size()>0) {
                    if (ilen<lenMin)
                        lenMin=ilen;
                    std::sort(cbWLv[icb][ilen].begin(),cbWLv[icb][ilen].end());//sort
                    auto un1=std::unique(cbWLv[icb][ilen].begin(),cbWLv[icb][ilen].end());//collapse identical
                    cbWLv[icb][ilen].resize(std::distance(cbWLv[icb][ilen].begin(),un1));        
                };
            };
        };
    };

    time_t rawTime;
    time(&rawTime);
    pP->inOut->logMain << timeMonthDayTime(rawTime) << "Finished reading, sorting and deduplicating CB whitelist sequences." <<endl;

    
    if (!pP->quant.trSAM.yes) {
        pP->quant.yes = true;
        pP->quant.trSAM.yes = true;
        pP->quant.trSAM.bamYes = false;
        pP->quant.trSAM.bamCompression = -2;
        pP->quant.trSAM.indel = true;
        pP->quant.trSAM.softClip = true;
        pP->inOut->logMain << "Turning on Genomic->Transcriptomic coordinate conversion for STARsolo\n";
    };
    
    if (featureYes[2])
        pP->quant.geneFull.yes=true;

    samAttrYes=false;
    samAttrFeature=0;//hard-coded for now to follow Cell Ranger - error correction only done for feature=Gene
    if (pP->outSAMattrPresent.CB || pP->outSAMattrPresent.UB) {
        samAttrYes=true;
        if (!pP->outBAMcoord) {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: CB and/or UB attributes in --outSAMattributes can only be output in the sorted BAM file.\n";
            errOut << "SOLUTION: re-run STAR with --outSAMtype BAM SortedByCoordinate ...\n";
            exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };
    };  
};

void ParametersSolo::umiSwapHalves(uint32 &umi) {
    uint32 high=umi>>(umiL);
    umi &= umiMaskLow; //remove high
    umi <<= (umiL); //move low to high
    umi |= high; //add high
};