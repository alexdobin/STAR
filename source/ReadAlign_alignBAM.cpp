#include "ReadAlign.h"
#include "SequenceFuns.h"
#include "ErrorWarning.h"
#include "IncludeDefine.h"
# include "BAMfunctions.h"

#include <type_traits>

void ReadAlign::samAttrNM_MD (Transcript const &trOut, uint iEx1, uint iEx2, uint &tagNM, string &tagMD) {
    tagNM=0;
    tagMD="";
    char* R=Read1[trOut.roStr==0 ? 0:2];
    uint matchN=0;
    uint32 nMM=0, nI=0, nD=0;
    for (uint iex=iEx1;iex<=iEx2;iex++) {
        for (uint ii=0;ii<trOut.exons[iex][EX_L];ii++) {
            char r1=R[ii+trOut.exons[iex][EX_R]];
            char g1=genOut.G[ii+trOut.exons[iex][EX_G]];
            if ( r1!=g1 || r1==4 || g1==4) {
                ++nMM;
                tagMD+=to_string(matchN);
                tagMD+=P.genomeNumToNT[(uint8) g1];
                matchN=0;
            } else {
                matchN++;
            };
        };
        if (iex<iEx2) {
            if (trOut.canonSJ[iex] < 0) {//indels, junction are not included in edit distance
                nD+=trOut.exons[iex+1][EX_G]-(trOut.exons[iex][EX_G]+trOut.exons[iex][EX_L]);//deletion
            };
            nI+=trOut.exons[iex+1][EX_R]-trOut.exons[iex][EX_R]-trOut.exons[iex][EX_L];//insertion
            if (trOut.canonSJ[iex]==-1) {//deletion. TODO: This does not work if there is both deletion and insertion!!!
                tagMD+=to_string(matchN) + "^";
                for (uint ii=trOut.exons[iex][EX_G]+trOut.exons[iex][EX_L];ii<trOut.exons[iex+1][EX_G];ii++) {
                    tagMD+=P.genomeNumToNT[(uint8) genOut.G[ii]];
                };
                matchN=0;
            };
        };
    };
    tagMD+=to_string(matchN);
    //cout << nMM<<" "<<nD<<" "<<nI<<endl;
    tagNM=nMM+nI+nD;
};

int ReadAlign::alignBAM(Transcript const &trOut, uint nTrOut, uint iTrOut, uint trChrStart, uint mateChr, uint mateStart, char mateStrand, int alignType, bool *mateMap, vector<int> outSAMattrOrder, char** outBAMarray, uint* outBAMarrayN) {
    //return: number of lines (mates)

    //alignType>=0: unmapped reads
    //          -1: normal mapped reads
    //          -10: chimeric alignment, not supplemental (like -11,-12,-13)
    //          -11: chimeric alignment, supplemental, hard-clipping, chimeric junction on the left
    //          -12: chimeric alignment, supplemental, hard-clipping, chimeric junction on the right
    //          -13: chimeric alignment, supplemental, soft-clipping


    if (P.outSAMmode=="None") return 0; //no SAM/BAM output

    uint32 recSize=0; //record size - total for both mates
    outBAMarrayN[0]=0;
    outBAMarrayN[1]=0;

    //for SAM output need to split mates
    uint iExMate=0; //last exon of the first mate

    uint16 samFLAG=0;


    bool flagPaired = P.readNmates==2; //not readNends: this is about alignments

    uint nMates=1;
    if (alignType<0) {//mapped reads: SAM
        for (iExMate=0;iExMate<trOut.nExons-1;iExMate++) {
            if (trOut.canonSJ[iExMate]==-3){
                nMates=2;
                break;
            };
        };
    } else {
        nMates=0;
    };

    uint tLen=0,leftMostMate=0;
    if (nMates>1 && P.outSAMtlen==2) {
        tLen=max(trOut.exons[trOut.nExons-1][EX_G]+trOut.exons[trOut.nExons-1][EX_L],trOut.exons[iExMate][EX_G]+trOut.exons[iExMate][EX_L])-min(trOut.exons[0][EX_G],trOut.exons[iExMate+1][EX_G]);
        leftMostMate=(trOut.exons[0][EX_G]<=trOut.exons[iExMate+1][EX_G] ? 0 : 1);
    };

    uint leftMate=0; //the mate (0 or 1) which is on the left
    if (flagPaired) {
        leftMate=trOut.Str;
    };

    if (P.outSAMattrPresent.MC) {
        calcCIGAR(trOut, nMates, iExMate, leftMate);
    };

    int32 gxgnGene = -1; //gene to output in GX/GN
    if (P.outSAMattrPresent.GX || P.outSAMattrPresent.GN) {
        auto annFeat = readAnnot.annotFeatures[P.pSolo.samAttrFeature];
        if (annFeat.fSet.size()==1 && iTrOut < annFeat.fAlign.size() && annFeat.fAlign[iTrOut].size()==1) {
            //                        2nd condition needed when this function is called to output transcriptomic alignments, where GX/GN are not allowed, and iTrOut can be large
            gxgnGene = *annFeat.fAlign[iTrOut].begin();
        };
    };

    for (uint imate=0;imate < (alignType<0 ? nMates:P.readNmates);imate++) { //not readNends: this is about alignments

        uint iEx1=0;
        uint iEx2=0;
        uint Mate=0;
        uint Str=0;
        uint32 packedCIGAR[BAM_CIGAR_MaxSize];
        uint32 nCIGAR=0; //number of CIGAR operations
        int MAPQ=0;
        uint32 attrN=0;
        char attrOutArray[BAM_ATTR_MaxSize];
        uint trimL1=0, trimR1=0;

        if (alignType>=0) {//this mate is unmapped
            if (mateMap!=NULL && mateMap[imate]) continue; //this mate was mapped, do not record it as unmapped
            samFLAG=0x4;
            if (P.readNmates==2) {//paired read, //not readNends: this is about alignments
                samFLAG|=0x1 + (imate==0 ? 0x40 : 0x80);
                if (mateMap[1-imate]) {//mate mapped
                    if (trOut.Str!=(1-imate))
                    {//mate strand reverted
                       samFLAG|=0x20;
                    };
                    mateChr=trOut.Chr;
                    trChrStart=genOut.chrStart[mateChr];
                    mateStart=trOut.exons[0][EX_G] - trChrStart;
                    mateStrand= trOut.Str == (1-imate) ? 0 : 1;

                    if (!trOut.primaryFlag && P.outSAMunmapped.keepPairs)
                    {//mapped mate is not primary, keep unmapped mate for each pair, hence need to mark some as not primary
                        samFLAG|=0x100;
                    };

                } else {//mate unmapped
                    samFLAG|=0x8;
                };
            };

            if (readFilter=='Y') samFLAG|=0x200; //not passing quality control

            if (mateMap[1-imate])
            {//mate is mapped, fill the infromation from trOut

            };

            Mate=imate;
            Str=Mate;

            attrN=0;
            attrN+=bamAttrArrayWriteInt(0,"NH",attrOutArray+attrN,P);
            attrN+=bamAttrArrayWriteInt(0,"HI",attrOutArray+attrN,P);
            attrN+=bamAttrArrayWriteInt(trOut.maxScore,"AS",attrOutArray+attrN,P);
            attrN+=bamAttrArrayWriteInt(trOut.nMM,"nM",attrOutArray+attrN,P);
            attrN+=bamAttrArrayWrite((to_string((uint) alignType)).at(0), "uT",attrOutArray+attrN); //cast to uint is only necessary for old compilers

            if (P.outSAMattrPresent.RG)
                attrN+=bamAttrArrayWrite(P.outSAMattrRG.at(readFilesIndex),"RG",attrOutArray+attrN);
            if (P.outSAMattrPresent.CR)
                attrN+=bamAttrArrayWrite(soloRead->readBar->cbSeq,"CR",attrOutArray+attrN);
            if (P.outSAMattrPresent.CY)
                attrN+=bamAttrArrayWrite(soloRead->readBar->cbQual,"CY",attrOutArray+attrN);
            if (P.outSAMattrPresent.CB && soloRead->readBar->cbSeqCorrected!="")
                attrN+=bamAttrArrayWrite(soloRead->readBar->cbSeqCorrected,"CB",attrOutArray+attrN);            
            if (P.outSAMattrPresent.UR)
                attrN+=bamAttrArrayWrite(soloRead->readBar->umiSeq,"UR",attrOutArray+attrN);
            if (P.outSAMattrPresent.UY)
                attrN+=bamAttrArrayWrite(soloRead->readBar->umiQual,"UY",attrOutArray+attrN);
            if (P.outSAMattrPresent.sM)
                attrN+=bamAttrArrayWrite(soloRead->readBar->cbMatch,"sM",attrOutArray+attrN);
            if (P.outSAMattrPresent.sS)
                attrN+=bamAttrArrayWrite(soloRead->readBar->bSeq,"sS",attrOutArray+attrN);
            if (P.outSAMattrPresent.sQ)
                attrN+=bamAttrArrayWrite(soloRead->readBar->bQual,"sQ",attrOutArray+attrN);
            if (P.outSAMattrPresent.cN) {
                vector <int32> v1={(int32)clipMates[imate][0].clippedN, (int32)clipMates[imate][1].clippedN};
                attrN+=bamAttrArrayWrite(v1,"cN",attrOutArray+attrN);
            };
            
        } else {//this mate is mapped
            if (flagPaired) {//paired reads
                samFLAG=0x0001;
                if (iExMate==trOut.nExons-1) {//single mate
                    if (mateChr>genOut.nChrReal) samFLAG|=0x0008; //not mapped as pair
                } else {//properly paired
                    samFLAG|=0x0002; //mapped as pair
                };
            } else {//single end
                samFLAG=0;
            };

            if (readFilter=='Y') samFLAG|=0x200; //not passing quality control

            if (alignType==-11 || alignType==-12 || alignType==-13) samFLAG|=0x800; //mark chimeric alignments
            if (!trOut.primaryFlag) samFLAG|=0x100; //mark not primary align

            iEx1 = (imate==0 ? 0 : iExMate+1);
            iEx2 = (imate==0 ? iExMate : trOut.nExons-1);
            Mate=trOut.exons[iEx1][EX_iFrag];
            Str= trOut.Str;//note that Strand = the mate on the left

            if (Mate==0) {
                samFLAG|= Str*0x10;
                if (nMates==2) samFLAG|= (1-Str)*0x20;
            } else {//second mate strand need to be reverted
                samFLAG|= (1-Str)*0x10;
                if (nMates==2) samFLAG|= Str*0x20;
            };

            if (flagPaired) {
                samFLAG|= (Mate==0 ? 0x0040 : 0x0080);
                if (flagPaired && nMates==1 && mateStrand==1) samFLAG|=0x20;//revert strand using inout value of mateStrand (e.g. for chimeric aligns)
            };


            uint trimL;
            if (Str==0 && Mate==0) {
                trimL=clipMates[Mate][0].clippedN;
            } else if (Str==0 && Mate==1) {
                trimL=clipMates[Mate][1].clippedN;
            } else if (Str==1 && Mate==0) {
                trimL=clipMates[Mate][1].clippedN;
            } else {
                trimL=clipMates[Mate][0].clippedN;
            };

            nCIGAR=0; //number of CIGAR operations

            trimL1 = trimL + trOut.exons[iEx1][EX_R] - (trOut.exons[iEx1][EX_R]<readLength[leftMate] ? 0 : readLength[leftMate]+1);
            if (trimL1>0) {
                packedCIGAR[nCIGAR++]=trimL1<<BAM_CIGAR_OperationShift | (alignType==-11 ? BAM_CIGAR_H : BAM_CIGAR_S);
            };

            vector<int32> SJintron;
            vector<char> SJmotif;

            for (uint ii=iEx1;ii<=iEx2;ii++) {
                if (ii>iEx1) {//record gaps
                    uint gapG=trOut.exons[ii][EX_G]-(trOut.exons[ii-1][EX_G]+trOut.exons[ii-1][EX_L]);
                    uint gapR=trOut.exons[ii][EX_R]-trOut.exons[ii-1][EX_R]-trOut.exons[ii-1][EX_L];
                    //it's possible to have a D or N and I at the same time
                    if (gapR>0){
                        packedCIGAR[nCIGAR++]=gapR<<BAM_CIGAR_OperationShift | BAM_CIGAR_I;
                    };
                    if (trOut.canonSJ[ii-1]>=0 || trOut.sjAnnot[ii-1]==1) {//junction: N
                        packedCIGAR[nCIGAR++]=gapG<<BAM_CIGAR_OperationShift | BAM_CIGAR_N;
                        SJmotif.push_back(trOut.canonSJ[ii-1] + (trOut.sjAnnot[ii-1]==0 ? 0 : SJ_SAM_AnnotatedMotifShift)); //record junction type
                        SJintron.push_back((int32) (trOut.exons[ii-1][EX_G] + trOut.exons[ii-1][EX_L] + 1 - trChrStart) );//record intron start
                        SJintron.push_back((int32) (trOut.exons[ii][EX_G] - trChrStart)); //record intron end
                    } else if (gapG>0) {//deletion: D
                        packedCIGAR[nCIGAR++]=gapG<<BAM_CIGAR_OperationShift | BAM_CIGAR_D;
                    };
                };
                if (trOut.exons[ii][EX_L]>0) //do not record 0-length blocks (which could appear near junctions followed by deletions)
                    packedCIGAR[nCIGAR++]=trOut.exons[ii][EX_L]<<BAM_CIGAR_OperationShift | BAM_CIGAR_M;
            };

            if (SJmotif.size()==0) {//no junctions recorded, mark with -1
                SJmotif.push_back(-1);
                SJintron.push_back(-1);
            };

            trimR1=(trOut.exons[iEx1][EX_R]<readLength[leftMate] ? \
                readLengthOriginal[leftMate] : readLength[leftMate]+1+readLengthOriginal[Mate]) \
                - trOut.exons[iEx2][EX_R]-trOut.exons[iEx2][EX_L] - trimL;
            if ( trimR1 > 0 ) {
                packedCIGAR[nCIGAR++]=trimR1<<BAM_CIGAR_OperationShift | (alignType==-12 ? BAM_CIGAR_H : BAM_CIGAR_S);
            };

            MAPQ=P.outSAMmapqUnique;
            if (nTrOut>=5) {
                MAPQ=0;
            } else if (nTrOut>=3) {
                MAPQ=1;
            } else if (nTrOut==2) {
                MAPQ=3;
            };

            //attribute string
            uint tagNM=(uint) -1;
            string tagMD("");

            attrN=0;
            for (uint ii=0;ii<outSAMattrOrder.size();ii++) {
                switch (outSAMattrOrder[ii]) {
                    case ATTR_NH:
                        attrN+=bamAttrArrayWriteInt(nTrOut,"NH",attrOutArray+attrN,P);
                        break;
                    case ATTR_HI:
                        attrN+=bamAttrArrayWriteInt(iTrOut+P.outSAMattrIHstart,"HI",attrOutArray+attrN,P);
                        break;
                    case ATTR_AS:
                        attrN+=bamAttrArrayWriteInt(trOut.maxScore,"AS",attrOutArray+attrN,P);
                        break;
                    case ATTR_nM:
                        attrN+=bamAttrArrayWriteInt(trOut.nMM,"nM",attrOutArray+attrN,P);
                        break;
                    case ATTR_jM:
                        attrN+=bamAttrArrayWrite(SJmotif,"jM",attrOutArray+attrN);
                        break;
                    case ATTR_jI:
                        attrN+=bamAttrArrayWrite(SJintron,"jI",attrOutArray+attrN);
                        break;
                    case ATTR_XS:
                        if (trOut.sjMotifStrand==1) {
                            attrN+=bamAttrArrayWrite('+',"XS",attrOutArray+attrN);
                        } else if (trOut.sjMotifStrand==2) {
                            attrN+=bamAttrArrayWrite('-',"XS",attrOutArray+attrN);
                        };
                        break;
                        
                    case ATTR_cN:
                        {
                            vector <int32> v1={(int32)clipMates[imate][0].clippedN, (int32)clipMates[imate][1].clippedN};
                            attrN+=bamAttrArrayWrite(v1,"cN",attrOutArray+attrN);
                            break;
                        };
                    case ATTR_NM:
                        if ( tagNM == (uint) -1 ) samAttrNM_MD (trOut, iEx1, iEx2, tagNM, tagMD);
                        attrN+=bamAttrArrayWriteInt(tagNM,"NM",attrOutArray+attrN,P);
                        break;
                    case ATTR_MD:
                        if ( tagMD.size()==0 ) samAttrNM_MD (trOut, iEx1, iEx2, tagNM, tagMD);
                        attrN+=bamAttrArrayWrite(tagMD,"MD",attrOutArray+attrN);
                        break;
                    case ATTR_RG:
                        attrN+=bamAttrArrayWrite(P.outSAMattrRG.at(readFilesIndex),"RG",attrOutArray+attrN);
                        break;
                    case ATTR_rB:
                        {
                            vector <int32> rb;
                            for (uint ii=iEx1;ii<=iEx2;ii++) {
                                rb.push_back( (int32) trOut.exons[ii][EX_R]+1 );
                                rb.push_back( (int32) trOut.exons[ii][EX_R]+trOut.exons[ii][EX_L]);
                                rb.push_back( (int32) (trOut.exons[ii][EX_G]-genOut.chrStart[trOut.Chr]+1) );
                                rb.push_back( (int32) (trOut.exons[ii][EX_G]-genOut.chrStart[trOut.Chr]+trOut.exons[ii][EX_L]) );
                            };
                            attrN+=bamAttrArrayWrite(rb,"rB",attrOutArray+attrN);
                        };
                        break;
                    case ATTR_vG:
                    {
                        const vector <int32> &v1=trOut.varGenCoord;
                        if (v1.size()>0)
                            attrN+=bamAttrArrayWrite(v1,"vG",attrOutArray+attrN);
                        break;
                    };
                    case ATTR_vA:
                    {
                        const vector <char> &v1=trOut.varAllele;
                        if (v1.size()>0)
                            attrN+=bamAttrArrayWrite(v1,"vA",attrOutArray+attrN);
                        break;
                    };
                    case ATTR_vW:
                    {
                        if (waspType!=-1)
                            attrN+=bamAttrArrayWrite( (int32) waspType, "vW", attrOutArray+attrN );
                        break;
                    };

                    case ATTR_ha:
                    {
                        if (mapGen.pGe.transform.type==2) {
                            attrN+=bamAttrArrayWrite( (int32) trOut.haploType, "ha", attrOutArray+attrN );
                        };
                        break; 
                    };

                    case ATTR_ch:
                        if (alignType<=-10)
                        {//chimeric alignment
                            attrN+=bamAttrArrayWrite('1',"ch",attrOutArray+attrN);
                        };
                        break;
                    case ATTR_MC:
                        if (nMates>1)
                        {//chimeric alignment
                            attrN+=bamAttrArrayWrite(matesCIGAR[1-imate],"MC",attrOutArray+attrN);
                        };
                        break;
                        
                    //solo
                    case ATTR_CR:
                        attrN+=bamAttrArrayWrite(soloRead->readBar->cbSeq,"CR",attrOutArray+attrN);
                        break;
                    case ATTR_CY:
                        attrN+=bamAttrArrayWrite(soloRead->readBar->cbQual,"CY",attrOutArray+attrN);
                        break;
                    case ATTR_UR:
                        attrN+=bamAttrArrayWrite(soloRead->readBar->umiSeq,"UR",attrOutArray+attrN);
                        break;
                    case ATTR_UY:
                        attrN+=bamAttrArrayWrite(soloRead->readBar->umiQual,"UY",attrOutArray+attrN);
                        break;

                    case ATTR_GX:
                        {
                            string algenes = gxgnGene<0 ? "-" : chunkTr->geID[gxgnGene];
                            attrN+=bamAttrArrayWrite(algenes, "GX",attrOutArray+attrN);
                        };
                        break;
                            
                    case ATTR_GN:
                        {
                            string algenes = gxgnGene<0 ? "-" : chunkTr->geName[gxgnGene];
                            attrN+=bamAttrArrayWrite(algenes, "GN",attrOutArray+attrN);
                        };
                        break;       

                    case ATTR_gx:
                        {
                            string algenes;
                            for (auto &gg: readAnnot.annotFeatures[P.pSolo.samAttrFeature].fAlign[iTrOut])
                                algenes += chunkTr->geID[gg] + ";";
                            if (algenes.size()>0) {
                                algenes.pop_back();
                            } else {
                                algenes = "-";
                            };
                            attrN+=bamAttrArrayWrite(algenes,"gx",attrOutArray+attrN);                            
                        };
                        break;

                    case ATTR_gn:
                        {
                            string algenes;
                            for (auto &gg: readAnnot.annotFeatures[P.pSolo.samAttrFeature].fAlign[iTrOut])
                                algenes += chunkTr->geName[gg] + ";";
                            if (algenes.size()>0) {
                                algenes.pop_back();
                            } else {
                                algenes = "-";
                            };
                            attrN+=bamAttrArrayWrite(algenes,"gn",attrOutArray+attrN);                            
                        };
                        break;                                                       
                        
                    case ATTR_sF:
                        {
                            vector<int32_t> sf(2);
                            sf[0] = readAnnot.annotFeatures[P.pSolo.samAttrFeature].ovType;
                            sf[1] = readAnnot.annotFeatures[P.pSolo.samAttrFeature].fSet.size();
                            if ( (sf[0]==1 || sf[0]==3 || sf[0]==5) && readAnnot.annotFeatures[P.pSolo.samAttrFeature].fAlign[iTrOut].size()==0) {
                                sf={-1,-1};
                            };
                            attrN += bamAttrArrayWrite(sf,"sF",attrOutArray+attrN);

                        };
                        break;

                    case ATTR_sM:
                        attrN+=bamAttrArrayWrite(soloRead->readBar->cbMatch,"sM",attrOutArray+attrN);
                        break;
                    case ATTR_sS:
                        attrN+=bamAttrArrayWrite(soloRead->readBar->bSeq,"sS",attrOutArray+attrN);
                        break;
                    case ATTR_sQ:
                        attrN+=bamAttrArrayWrite(soloRead->readBar->bQual,"sQ",attrOutArray+attrN);
                        break;
                        
                    case ATTR_CB:
                        if (soloRead->readBar->cbSeqCorrected!="") //only if corrected CB is defined (e.g. CB_samTagOut), otherwise it will be added at the BAM sorting stage
                        	attrN+=bamAttrArrayWrite(soloRead->readBar->cbSeqCorrected,"CB",attrOutArray+attrN);
                        break;
                    case ATTR_UB: //will be added at the BAM sorting stage
                        break;
                        
                    default:
                        ostringstream errOut;
                        errOut <<"EXITING because of FATAL BUG: unknown/unimplemented SAM/BAM atrribute (tag): "<<outSAMattrOrder[ii] <<"\n";
                        errOut <<"SOLUTION: please make sure that SAM attributes in --outSAMattributes are compatible with the --outSAMtype option\n";
                        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
                };
            };
        };

        if (P.readFilesTypeN==10 && !P.readFiles.samAttrKeepNone) {
            attrN+=bamAttrArrayWriteSAMtags(readNameExtra[Mate], attrOutArray+attrN, P);
        };
////////////////////////////// prepare sequence and qualities
        char seqMate[DEF_readSeqLengthMax+1], qualMate[DEF_readSeqLengthMax+1];
        char *seqOut=NULL, *qualOut=NULL;

        if ( Mate==Str)  {//seq strand is correct, or mate is unmapped
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

        uint seqMateLength=readLengthOriginal[Mate];
        if (alignType==-11) {//hard-clip on the left
            seqMateLength-=trimL1;
            seqOut+=trimL1;
            qualOut+=trimL1;
        } else if (alignType==-12) {
            seqMateLength-=trimR1;
        } else {//no-chimeric alignment
        };

        //pack sequence
        nuclPackBAM(seqOut,seqMate,seqMateLength);

/////////////////////////////////// write BAM
        uint32 *pBAM=(uint32*) (outBAMarray[imate]);
        recSize=0;

        //1: refID: Reference sequence ID, -1 <= refID <= n ref; -1 for a read without a mapping position.
        if (alignType<0) {
            pBAM[1]=trOut.Chr;
        } else {
            pBAM[1]=(uint32) -1;
        };

        //2: pos: 0-based leftmost coordinate (= POS - 1): int32_t
        if (alignType<0) {
            pBAM[2]=trOut.exons[iEx1][EX_G] - trChrStart;
        } else {
            pBAM[2]=(uint32) -1;
        };

        //3: bin mq nl bin<<16|MAPQ<<8|l read name; bin is computed by the > reg2bin() function in Section 4.3; l read name is the length> of read name below (= length(QNAME) + 1).> uint32 t
        if (alignType<0) {
            pBAM[3]=( ( reg2bin(trOut.exons[iEx1][EX_G] - trChrStart,trOut.exons[iEx2][EX_G] + trOut.exons[iEx2][EX_L] - trChrStart) << 16 ) \
                   |( MAPQ<<8 ) | ( strlen(readName) ) ); //note:read length includes 0-char
        } else {
            pBAM[3]=( reg2bin(-1,0) << 16 |  strlen(readName) );//4680=reg2bin(-1,0)
        };

        //4: FLAG<<16|n cigar op; n cigar op is the number of operations in CIGAR.
        pBAM[4]=( ( ((samFLAG & P.outSAMflagAND) | P.outSAMflagOR) << 16 ) | (nCIGAR) );

        //5: l seq Length of SEQ
        pBAM[5]=seqMateLength;

        //6: next refID Ref-ID of the next segment (􀀀1  mate refID < n ref)
        if (nMates>1) {
            pBAM[6]=trOut.Chr;
        } else if (mateChr<genOut.nChrReal){
            pBAM[6]=mateChr;
        } else {
            pBAM[6]=-1;
        };

        //7: next pos 0-based leftmost pos of the next segment (= PNEXT 􀀀 1)
        if (nMates>1) {
            pBAM[7]=trOut.exons[(imate==0 ? iExMate+1 : 0)][EX_G] - trChrStart;
        } else if (mateChr<genOut.nChrReal){
            pBAM[7]=mateStart;
        } else {
            pBAM[7]=-1;
        };

        //8: tlen Template length (= TLEN)
        if (nMates>1) {
            if (P.outSAMtlen==1) {
                int32 tlen=trOut.exons[trOut.nExons-1][EX_G]+trOut.exons[trOut.nExons-1][EX_L]-trOut.exons[0][EX_G];
                pBAM[8]=(imate==0 ? tlen : -tlen);
            } else if (P.outSAMtlen==2) {
                int32 tlen=(int32)tLen;
                pBAM[8]=(imate==leftMostMate ? tlen : -tlen);
            };
        } else {
            pBAM[8]=0;
        };

        recSize+=9*sizeof(int32); //core record size

        //Read name1, NULL terminated (QNAME plus a tailing `\0')
        strcpy(outBAMarray[imate]+recSize,readName+1);
        recSize+=strlen(readName);

        //CIGAR: op len<<4|op. `MIDNSHP=X'!`012345678'
        memcpy(outBAMarray[imate]+recSize,packedCIGAR, nCIGAR*sizeof(int32));
        recSize+=nCIGAR*sizeof(int32);

        //4-bit encoded read: `=ACMGRSVTWYHKDBN'! [0; 15]; other characters mapped to `N'; high nybble rst (1st base in the highest 4-bit of the 1st byte)
        memcpy(outBAMarray[imate]+recSize,seqMate,(seqMateLength+1)/2);
        recSize+=(seqMateLength+1)/2;

        //Phred base quality (a sequence of 0xFF if absent)
        if (readFileType==2 && P.outSAMmode != "NoQS") {//output qualtiy
            for (uint32 ii=0; ii<seqMateLength; ii++) {
                outBAMarray[imate][recSize+ii]=qualOut[ii]-33;
            };
//             memcpy(outBAM+recSize,qualOut,readLengthOriginal[Mate]);
        } else {
            memset(outBAMarray[imate]+recSize,0xFF,seqMateLength);
        };
        recSize+=seqMateLength;

        //atributes
        memcpy(outBAMarray[imate]+recSize,attrOutArray,attrN);
        recSize+=attrN;

        //total size of the record
        pBAM[0]=recSize-sizeof(uint32);//record size excluding the size entry itself
        outBAMarrayN[imate]=recSize;
    };//for (uint imate=0;imate<nMates;imate++)

    return ( outBAMarrayN[1]==0 ? 1 : 2);
};
