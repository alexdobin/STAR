#include "ReadAlign.h"
#include "SequenceFuns.h"
#include "ErrorWarning.h"

uint64 ReadAlign::outputSpliceGraphSAM(Transcript const &trOut, uint nTrOut, uint iTrOut, ostream *outStream) 
{

    uint64 outStreamPos0=(uint64)outStream->tellp();

    uint16 samFLAG=0;
    if (readFilter=='Y') 
        samFLAG |= 0x200; //not passing quality control
        
    if (unmapType>=0) {//unmapped reads: SAM
        samFLAG |= 0x4;
        
        *outStream << readName+1 <<"\t"<< samFLAG <<"\t"<< '*' <<"\t"<< '0' <<"\t"<< '0' <<"\t"<< '*';
        //mate
        *outStream <<"\t"<< '*' <<"\t"<< '0' <<"\t"<< '0' ;
        
        *outStream <<"\t"<< Read0[0] <<"\t"<< (readFileType==2 ? Qual0[0]:"*") \
                <<"\tNH:i:0" <<"\tHI:i:0" <<"\tAS:i:"<<trOut.maxScore <<"\tnM:i:"<<trOut.nMM<<"\tuT:A:" <<unmapType;
        if (!P.outSAMattrRG.empty()) *outStream<< "\tRG:Z:" <<P.outSAMattrRG.at(readFilesIndex);

        if (P.readFilesTypeN==10 && !readNameExtra[0].empty()) {//SAM files as input - output extra attributes
            *outStream << "\t" <<readNameExtra[0];
        };
        *outStream <<"\n";
        
        return (uint)outStream->tellp()-outStreamPos0;
    };
    
    //mapped
    if (!trOut.primaryFlag)
        samFLAG |= 0x100;
    if (trOut.Str==1)
        samFLAG |= 0x10;

    string cigar;
    const static char* cigarChars="MIDNS";
    for (auto &cc : trOut.cigar) {
        cigar += to_string(cc[1]) + cigarChars[cc[0]];
    };

    char seqRev[DEF_readSeqLengthMax+1], qualRev[DEF_readSeqLengthMax+1];
    char *seqOut=NULL, *qualOut=NULL;

    if ( trOut.Str==0 )  {//seq strand is correct
        seqOut=Read0[0];
        qualOut=Qual0[0];
    } else {
        revComplementNucleotides(Read0[0], seqRev, Lread);
        for (uint ii=0;ii<Lread; ii++) 
            qualRev[ii]=Qual0[0][Lread-1-ii];
        seqOut=seqRev;
        qualOut=qualRev;
    };
    seqOut[Lread]=0;//to ensure string termination
    qualOut[Lread]=0;
    
    int MAPQ=P.outSAMmapqUnique;
    if (nTrOut>=5) {
        MAPQ=0;
    } else if (nTrOut>=3) {
        MAPQ=1;
    } else if (nTrOut==2) {
        MAPQ=3;
    };

    *outStream << readName+1 <<"\t"<< ((samFLAG & P.outSAMflagAND) | P.outSAMflagOR) <<"\t"<< genOut.chrName[trOut.Chr] <<"\t"<< trOut.gStart + 1 - genOut.chrStart[trOut.Chr]
                <<"\t"<< MAPQ <<"\t"<< cigar;

    *outStream <<"\t"<< "*" <<"\t"<< 0 <<"\t"<< 0;

    *outStream <<"\t"<< seqOut;

    if (readFileType==2 && P.outSAMmode != "NoQS") {//fastq
        *outStream <<"\t"<< qualOut ;
    } else {
        *outStream <<"\t"<< "*";
    };
    
    uint32 tagNM = trOut.nMM + trOut.lIns + trOut.lDel;

        for (uint ii=0;ii<P.outSAMattrOrder.size();ii++) {
            switch (P.outSAMattrOrder[ii]) {
                case ATTR_NH:
                    *outStream <<"\tNH:i:" << nTrOut;
                    break;
                case ATTR_HI:
                    *outStream <<"\tHI:i:"<<iTrOut+P.outSAMattrIHstart;
                    break;
                case ATTR_AS:
                    *outStream<<"\tAS:i:"<<trOut.maxScore;
                    break;
                case ATTR_nM:
                    *outStream<<"\tnM:i:"<<trOut.nMM;
                    break;
                case ATTR_NM:
                    *outStream<< "\tNM:i:" <<tagNM;
                    break;
                    
//                 case ATTR_jM:
//                     *outStream<<"\tjM:B:c"<< SJmotif;
//                     break;
//                 case ATTR_jI:
//                     *outStream<<"\tjI:B:i"<< SJintron;
//                     break;
//                 case ATTR_XS:
//                     if (trOut.sjMotifStrand==1) {
//                         *outStream<<"\tXS:A:+";
//                     } else if (trOut.sjMotifStrand==2) {
//                         *outStream<<"\tXS:A:-";
//                     };
//                     break;
//                 case ATTR_MD:
//                     *outStream<< "\tMD:Z:" <<tagMD;
//                     break;
                case ATTR_RG:
                    *outStream<< "\tRG:Z:" <<P.outSAMattrRG.at(readFilesIndex);
                    break;
//                 case ATTR_MC:
//                     if (nMates>1) {
//                         *outStream<< "\tMC:Z:" <<matesCIGAR[1-imate];
//                     };
//                     break;
                    
                //do nothing - this attributes only work for BAM output
                case ATTR_jM:
                case ATTR_jI:
                case ATTR_XS:
                case ATTR_MD:
                    
                case ATTR_ch:
                case ATTR_CR:
                case ATTR_CY:
                case ATTR_UR:
                case ATTR_UY:
                case ATTR_CB:
                case ATTR_UB:
                case ATTR_sM:
                case ATTR_sS:
                case ATTR_sQ:
                case ATTR_rB:
                case ATTR_vG:
                case ATTR_vA:
                case ATTR_vW:
                    break;
                default:
                    ostringstream errOut;
                    errOut <<"EXITING because of FATAL error: unknown/unimplemented SAM atrribute (tag): "<<P.outSAMattrOrder[ii] <<"\n";
                    errOut <<"SOLUTION: contact Alex Dobin at dobin@cshl.edu\n";
                    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
            };
        };

    if (P.readFilesTypeN==10 && !readNameExtra[0].empty()) {//SAM files as input - output extra attributes
         *outStream << "\t" << readNameExtra.at(0);
    };

    *outStream << "\n"; //done with one SAM line

    return (uint)outStream->tellp()-outStreamPos0;
};
