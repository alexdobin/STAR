#include "ParametersChimeric.h"
#include "Parameters.h"
#include "ErrorWarning.h"
#include "GlobalVariables.h"

void ParametersChimeric::initialize(Parameters *pPin)
{
    out.bam=false;
    out.junctions=false;
    out.samOld=false;
    out.bamHardClip=true;//default
    
    if (segmentMin==0)
        return;
    
    pP=pPin;

    pthread_mutex_init(&g_threadChunks.mutexOutChimSAM, NULL);
    pthread_mutex_init(&g_threadChunks.mutexOutChimJunction, NULL);

    for (const auto& type1 : out.type) {
        if (type1=="WithinBAM") {
            out.bam=true;
        } else if (type1=="SeparateSAMold") {
            out.samOld=true;
        } else if (type1=="Junctions") {
            out.junctions=true;
        } else if (type1=="HardClip") {
            out.bamHardClip=true;
        } else if (type1=="SoftClip") {
            out.bamHardClip=false;
        } else {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --chimOutType: "<<type1 <<"\n";
            errOut <<"SOLUTION: re-run STAR with --chimOutType Junctions , SeparateSAMold  , WithinBAM , HardClip \n";
            exitWithError(errOut.str(), std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };
    };
    
    if (out.samOld) {
        pP->inOut->outChimSAM.open((pP->outFileNamePrefix + "Chimeric.out.sam").c_str());
        pP->inOut->outChimSAM << pP->samHeader;
    };
    
    if (out.junctions) {
        pP->inOut->outChimJunction.open((pP->outFileNamePrefix + "Chimeric.out.junction").c_str());
        
        if (multimapNmax>0)
            // column headers for Chimeric.out.junction file
            pP->inOut->outChimJunction <<
                                        "chr_donorA" <<"\t"<<
                                        "brkpt_donorA" <<"\t"<<
                                        "strand_donorA" <<"\t"<<
                                        "chr_acceptorB" <<"\t"<<
                                        "brkpt_acceptorB" <<"\t"<<
                                        "strand_acceptorB" <<"\t"<<
                                        "junction_type" <<"\t"<<
                                        "repeat_left_lenA" <<"\t"<<
                                        "repeat_right_lenB" <<"\t"<<
                                        "read_name" <<"\t"<<
                                        "start_alnA" <<"\t"<<
                                        "cigar_alnA" <<"\t"<<
                                        "start_alnB" <<"\t"<<
                                        "cigar_alnB" <<"\t"<<
                                        "num_chim_aln" <<"\t"<<
                                        "max_poss_aln_score" <<"\t"<<
                                        "non_chim_aln_score" <<"\t"<<
                                        "this_chim_aln_score" <<"\t"<<
                                        "bestall_chim_aln_score" <<"\t"<<
                                        "PEmerged_bool" <<"\t"<<
                                        "readgrp" <<"\n";
    };
    


    if (out.bam && !pP->outBAMunsorted && !pP->outBAMcoord) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal PARAMETERS error: --chimOutType WithinBAM requires BAM output\n";
            errOut <<"SOLUTION: re-run with --outSAMtype BAM Unsorted/SortedByCoordinate\n";
            exitWithError(errOut.str(), std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };

    if (multimapNmax>0 && out.samOld) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal PARAMETERS error: --chimMultimapNmax > 0 (new chimeric detection) presently only works with --chimOutType Junctions/WithinBAM\n";
            errOut <<"SOLUTION: re-run with --chimOutType Junctions/WithinBAM\n";
            exitWithError(errOut.str(), std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };

    if (pP->peOverlap.NbasesMin > 0) {
        if (multimapNmax == 0 && (out.junctions || out.samOld)) {
                ostringstream errOut;
                errOut <<"EXITING because of fatal PARAMETERS error: --chimMultimapNmax 0 (default old chimeric detection) and --peOverlapNbasesMin > 0 (merging ovelrapping mates) presently only works with --chimOutType WithinBAM\n";
                errOut <<"SOLUTION: re-run with --chimOutType WithinBAM\n";
                exitWithError(errOut.str(), std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };
    };

    if (out.bam && !pP->outSAMattrPresent.NM) {
       pP->outSAMattrOrder.push_back(ATTR_NM);
       pP->inOut->logMain << "WARNING --chimOutType=WithinBAM, therefore STAR will output NM attribute" <<endl;
    };

    filter.genomicN=false;
    for (uint ii=0; ii<filter.stringIn.size(); ii++) {
        if (filter.stringIn.at(ii)=="banGenomicN") {
            filter.genomicN=true;
        } else if (filter.stringIn.at(ii)=="None") {
            //nothing to do
        }
        else{
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: unrecognized value of --chimFilter="<<filter.stringIn.at(ii)<<"\n";
            errOut << "SOLUTION: use allowed values: banGenomicN || None";
            exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        };
    };
};
