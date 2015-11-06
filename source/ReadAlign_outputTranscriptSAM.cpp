#include "ReadAlign.h"
#include "SequenceFuns.h"
#include "ErrorWarning.h"

uint ReadAlign::outputTranscriptSAM(Transcript const &trOut, uint nTrOut, uint iTrOut, uint mateChr, uint mateStart, char mateStrand, int unmapType, bool *mateMapped, ostream *outStream) {    
    
    if (P->outSAMmode=="None") return 0; //no SAM output  
    
    uint outStreamPos0=(uint)outStream->tellp();
    
    if (unmapType>=0) {//unmapped reads: SAM
        for (uint imate=0;imate<P->readNmates;imate++) {//cycle over mates
            if (!mateMapped[imate]) {
                uint16 samFLAG=0x4;
                if (P->readNmates==2) {//paired read
                    samFLAG+=0x1 + (imate==0 ? 0x40 : 0x80);
                    if (mateMapped[1-imate]) {//mate mapped
                        if (trBest->Str!=1-imate) samFLAG+=0x20;//mate strand reverted
                    } else {//mate unmapped
                        samFLAG+=0x8;
                    };
                };

                if (readFilter=='Y') samFLAG+=0x200; //not passing quality control
                
                *outStream << readName+1 <<"\t"<< samFLAG \
                        <<"\t"<< '*' <<"\t"<< '0' <<"\t"<< '0' <<"\t"<< '*';

                if (mateMapped[1-imate]) {//mate is mapped
                    *outStream <<"\t"<< P->chrName[trBest->Chr] <<"\t"<< trBest->exons[0][EX_G] + 1 - P->chrStart[trBest->Chr];
                } else {
                    *outStream <<"\t"<< '*' <<"\t"<< '0';
                };

                *outStream <<"\t"<< '0' <<"\t"<< Read0[imate] <<"\t"<< (readFileType==2 ? Qual0[imate]:"*") \
                        <<"\tNH:i:0" <<"\tHI:i:0" <<"\tAS:i:"<<trBest->maxScore <<"\tnM:i:"<<trBest->nMM<<"\tuT:A:" <<unmapType;
                if (!P->outSAMattrRG.empty()) *outStream<< "\tRG:Z:" <<P->outSAMattrRG.at(readFilesIndex); 
                *outStream <<"\n";
                
            };
        };
        return (uint)outStream->tellp()-outStreamPos0;
    };//if (unmapType>=0 && outStream != NULL) //unmapped reads: SAM
    
    
    bool flagPaired = P->readNmates==2;
    string CIGAR;

    //for SAM output need to split mates
    uint iExMate; //last exon of the first mate
    uint nMates=1;
    for (iExMate=0;iExMate<trOut.nExons-1;iExMate++) {
        if (trOut.canonSJ[iExMate]==-3){
            nMates=2;
            break;
        };
    };

    uint samFLAG=0;
    uint leftMate=0; //the mate (0 or 1) which is on the left
    for (uint imate=0;imate<nMates;imate++) {
        if (flagPaired) {//paired reads
            samFLAG=0x0001;
            if (iExMate==trOut.nExons-1) {//single mate
                if (mateChr>P->nChrReal) samFLAG+=0x0008; //not mapped as pair
            } else {//properly paired
                samFLAG+=0x0002; //mapped as pair
            };                                 
        } else {//single end
            samFLAG=0;
        };

        if (readFilter=='Y') samFLAG+=0x200; //not passing quality control

        uint iEx1 = (imate==0 ? 0 : iExMate+1);
        uint iEx2 = (imate==0 ? iExMate : trOut.nExons-1);
        uint Mate=trOut.exons[iEx1][EX_iFrag];        
        uint Str= trOut.Str;//note that Strand = the mate on the left
        
        if (Mate==0) {
            samFLAG += Str*0x10;
            if (nMates==2) samFLAG += (1-Str)*0x20;
        } else {//second mate strand need to be reverted
            samFLAG += (1-Str)*0x10;
            if (nMates==2) samFLAG += Str*0x20;            
        };        
        
        if (flagPaired) {
            leftMate=Str;
            samFLAG += (Mate==0 ? 0x0040 : 0x0080);
            if (flagPaired && nMates==1 && mateStrand==1) samFLAG +=0x20;//revert strand using inout value of mateStrand (e.g. for chimeric aligns)
        };        

        //not primary align?
        if (!trOut.primaryFlag) samFLAG +=0x100;

        //empty streams
        samStreamCIGAR.str(std::string());
        samStreamSJmotif.str(std::string());
        samStreamSJintron.str(std::string());
//         samStreamSJannot.str(std::string());

        uint trimL;
        if (Str==0 && Mate==0) {
            trimL=clip5pNtotal[Mate];
        } else if (Str==0 && Mate==1) {
            trimL=clip3pNtotal[Mate];
        } else if (Str==1 && Mate==0) {
            trimL=clip3pNtotal[Mate];
        } else {
            trimL=clip5pNtotal[Mate];
        };
 
        uint trimL1 = trimL + trOut.exons[iEx1][EX_R] - (trOut.exons[iEx1][EX_R]<readLength[leftMate] ? 0 : readLength[leftMate]+1);
        if (trimL1>0) {
            samStreamCIGAR << trimL1 << "S"; //initial trimming
        };                      
        
        for (uint ii=iEx1;ii<=iEx2;ii++) {
            if (ii>iEx1) {//record gaps
                uint gapG=trOut.exons[ii][EX_G]-(trOut.exons[ii-1][EX_G]+trOut.exons[ii-1][EX_L]);
                uint gapR=trOut.exons[ii][EX_R]-trOut.exons[ii-1][EX_R]-trOut.exons[ii-1][EX_L];
                //it's possible to have a D or N and I at the same time
                if (gapR>0){
                    samStreamCIGAR << gapR;
                    samStreamCIGAR << "I";                            
                };                
                if (trOut.canonSJ[ii-1]>=0 || trOut.sjAnnot[ii-1]==1) {//junction: N
                    samStreamCIGAR << gapG;                        
                    samStreamCIGAR << "N";
                    samStreamSJmotif <<','<< trOut.canonSJ[ii-1] + (trOut.sjAnnot[ii-1]==0 ? 0 : SJ_SAM_AnnotatedMotifShift); //record junction type
//                     samStreamSJannot <<','<< (int) trOut.sjAnnot[ii-1]; //record annotation type                    
                    samStreamSJintron <<','<< trOut.exons[ii-1][EX_G] + trOut.exons[ii-1][EX_L] + 1 - P->chrStart[trOut.Chr] <<','\
                                   << trOut.exons[ii][EX_G] - P->chrStart[trOut.Chr]; //record intron loci
                } else if (gapG>0) {//deletion: N
                    samStreamCIGAR << gapG;
                    samStreamCIGAR << "D";
                };
            };                
            samStreamCIGAR << trOut.exons[ii][EX_L] << "M";
        };

        string SJmotif = samStreamSJmotif.str();
        string SJintron = samStreamSJintron.str();
//         string SJannot = samStreamSJannot.str();
        
        if (SJmotif.length()==0) {//no junctions recorded, mark with -1
            SJmotif=",-1"; 
            SJintron=",-1";
//             SJannot=",-1";
        };
        
        uint trimR1=(trOut.exons[iEx1][EX_R]<readLength[leftMate] ? \
            readLengthOriginal[leftMate] : readLength[leftMate]+1+readLengthOriginal[Mate]) \
            - trOut.exons[iEx2][EX_R]-trOut.exons[iEx2][EX_L] - trimL;
        if ( trimR1 > 0 ) {
            samStreamCIGAR << trimR1 << "S"; //final trimming
        };
        CIGAR=samStreamCIGAR.str();

           
        char seqMate[DEF_readSeqLengthMax+1], qualMate[DEF_readSeqLengthMax+1];
        char *seqOut=NULL, *qualOut=NULL;

        if ( Mate==Str )  {//seq strand is correct
            seqOut=Read0[Mate];
            qualOut=Qual0[Mate];
        } else {
            revComplementNucleotides(Read0[Mate], seqMate, readLengthOriginal[Mate]);
            seqMate[readLengthOriginal[Mate]]=0;
            for (uint ii=0;ii<readLengthOriginal[Mate]; ii++) qualMate[ii]=Qual0[Mate][readLengthOriginal[Mate]-1-ii];
            qualMate[readLengthOriginal[Mate]]=0;
            seqOut=&seqMate[0];
            qualOut=&qualMate[0];            
        };
        
//         return;

        int MAPQ=P->outSAMmapqUnique;
        if (nTrOut>=5) {
            MAPQ=0;
        } else if (nTrOut>=3) {
            MAPQ=1;
        } else if (nTrOut==2) {
            MAPQ=3;
        };
        
        *outStream << readName+1 <<"\t"<< ((samFLAG & P->outSAMflagAND) | P->outSAMflagOR) <<"\t"<< P->chrName[trOut.Chr] <<"\t"<< trOut.exons[iEx1][EX_G] + 1 - P->chrStart[trOut.Chr]
                <<"\t"<< MAPQ <<"\t"<< CIGAR;

        if (nMates>1) {
            *outStream <<"\t"<< "=" <<"\t"<< trOut.exons[(imate==0 ? iExMate+1 : 0)][EX_G]+  1 - P->chrStart[trOut.Chr]
                     <<"\t"<< (imate==0? "":"-") << trOut.exons[trOut.nExons-1][EX_G]+trOut.exons[trOut.nExons-1][EX_L]-trOut.exons[0][EX_G];
        } else if (mateChr<P->nChrReal){//mateChr is given in the function parameters
            *outStream <<"\t"<< P->chrName[mateChr] <<"\t"<< mateStart+1-P->chrStart[mateChr] <<"\t"<< 0;
        } else {
            *outStream <<"\t"<< "*" <<"\t"<< 0 <<"\t"<< 0;
        };

        
        *outStream <<"\t"<< seqOut;
        
        if (readFileType==2 && P->outSAMmode != "NoQS") {//fastq
            *outStream <<"\t"<< qualOut ;
        } else {
            *outStream <<"\t"<< "*";
        };   
        
//         vector<string> customAttr(outSAMattrN,"");
   
        uint tagNM=0;
        string tagMD("");
        if (P->outSAMattrPresent.NM || P->outSAMattrPresent.MD) {
            char* R=Read1[trOut.roStr==0 ? 0:2];
            uint matchN=0;
            for (uint iex=iEx1;iex<=iEx2;iex++) {
                for (uint ii=0;ii<trOut.exons[iex][EX_L];ii++) {
                    char r1=R[ii+trOut.exons[iex][EX_R]];
                    char g1=G[ii+trOut.exons[iex][EX_G]];
                    if ( r1!=g1 || r1==4 || g1==4) {
                        ++tagNM;
//                         if (matchN>0 || (ii==0 && iex>0 && trOut.canonSJ[iex]==-1) ) {
                        tagMD+=to_string(matchN);
//                         };
                        tagMD+=P->genomeNumToNT[(uint8) g1];
                        matchN=0;
                    } else {
                        matchN++;
                    };
                };
                if (iex<iEx2) {
                    if (trOut.canonSJ[iex]==-1) {//deletion
                        tagNM+=trOut.exons[iex+1][EX_G]-(trOut.exons[iex][EX_G]+trOut.exons[iex][EX_L]);
                        tagMD+=to_string(matchN) + "^";
                        for (uint ii=trOut.exons[iex][EX_G]+trOut.exons[iex][EX_L];ii<trOut.exons[iex+1][EX_G];ii++) {
                            tagMD+=P->genomeNumToNT[(uint8) G[ii]];
                        };
                        matchN=0;
                    } else if (trOut.canonSJ[iex]==-2) {//insertion
                        tagNM+=trOut.exons[iex+1][EX_R]-trOut.exons[iex][EX_R]-trOut.exons[iex][EX_L];
                    };
                };
            };      
            tagMD+=to_string(matchN);
        };
        for (int ii=0;ii<P->outSAMattrOrder.size();ii++) {
            switch (P->outSAMattrOrder[ii]) {
                case ATTR_NH:
                    *outStream <<"\tNH:i:" << nTrOut;
                    break;
                case ATTR_HI:
                    *outStream <<"\tHI:i:"<<iTrOut+P->outSAMattrIHstart;
                    break;
                case ATTR_AS:
                    *outStream<<"\tAS:i:"<<trOut.maxScore;
                    break;                    
                case ATTR_nM:
                    *outStream<<"\tnM:i:"<<trOut.nMM;   
                    break;                            
                case ATTR_jM:
                    *outStream<<"\tjM:B:c"<< SJmotif;
                    break;                            
                case ATTR_jI:
                    *outStream<<"\tjI:B:i"<< SJintron;
                    break;
                case ATTR_XS:
                    if (trOut.sjMotifStrand==1) {
                        *outStream<<"\tXS:A:+";
                    } else if (trOut.sjMotifStrand==2) {
                        *outStream<<"\tXS:A:-";
                    };                    
                    break;
                case ATTR_NM:
                    *outStream<< "\tNM:i:" <<tagNM;
                    break;
                case ATTR_MD:
                    *outStream<< "\tMD:Z:" <<tagMD;
                    break;
                case ATTR_RG:
                    *outStream<< "\tRG:Z:" <<P->outSAMattrRG.at(readFilesIndex);                    
                    break;
                default:
                    ostringstream errOut;
                    errOut <<"EXITING because of FATAL BUG: unknown/unimplemented SAM atrribute (tag): "<<P->outSAMattrOrder[ii] <<"\n";
                    errOut <<"SOLUTION: contact Alex Dobin at dobin@cshl.edu\n";
                    exitWithError(errOut.str(), std::cerr, P->inOut->logMain, EXIT_CODE_PARAMETER, *P);                         
            };
        };
     
        
//         for (uint ii=0;ii<customAttr.size();ii++) *outStream <<"\t"<< customAttr.at(ii); //output all attributes in the right order
        
        *outStream << "\n"; //done with one SAM line
    };//for (uint imate=0;imate<nMates;imate++)
    
    return (uint)outStream->tellp()-outStreamPos0;
};
