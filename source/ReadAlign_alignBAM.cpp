#include "ReadAlign.h"
#include "SequenceFuns.h"
#include "ErrorWarning.h"
#include "IncludeDefine.h"
#include <type_traits> // C++0x

void ReadAlign::samAttrNM_MD (Transcript const &trOut, uint iEx1, uint iEx2, uint &tagNM, string &tagMD) {
    tagNM=0;
    tagMD="";
    char* R=Read1[trOut.roStr==0 ? 0:2];
    uint matchN=0;
    for (uint iex=iEx1;iex<=iEx2;iex++) {
        for (uint ii=0;ii<trOut.exons[iex][EX_L];ii++) {
            char r1=R[ii+trOut.exons[iex][EX_R]];
            char g1=mapGen.G[ii+trOut.exons[iex][EX_G]];
            if ( r1!=g1 || r1==4 || g1==4) {
                ++tagNM;
                tagMD+=to_string(matchN);
                tagMD+=P.genomeNumToNT[(uint8) g1];
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
                    tagMD+=P.genomeNumToNT[(uint8) mapGen.G[ii]];
                };
                matchN=0;
            } else if (trOut.canonSJ[iex]==-2) {//insertion
                tagNM+=trOut.exons[iex+1][EX_R]-trOut.exons[iex][EX_R]-trOut.exons[iex][EX_L];
            };
        };
    };
    tagMD+=to_string(matchN);
};
// calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open)
int reg2bin(int beg, int end)
{
    --end;
    if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
    return 0;
};

int bamAttrArrayWrite(int32 attr, const char* tagName, char* attrArray ) {
    attrArray[0]=tagName[0];attrArray[1]=tagName[1];
    attrArray[2]='i';
    *( (int32*) (attrArray+3))=attr;
    return 3+sizeof(int32);
};
int bamAttrArrayWrite(float attr, const char* tagName, char* attrArray ) {
    attrArray[0]=tagName[0];attrArray[1]=tagName[1];
    attrArray[2]='f';
    *( (float*) (attrArray+3))=attr;
    return 3+sizeof(int32);
};
int bamAttrArrayWrite(char attr, const char* tagName, char* attrArray ) {
    attrArray[0]=tagName[0];attrArray[1]=tagName[1];
    attrArray[2]='A';
    attrArray[3]=attr;
    return 3+sizeof(char);
};
int bamAttrArrayWrite(string &attr, const char* tagName, char* attrArray ) {
    attrArray[0]=tagName[0];attrArray[1]=tagName[1];
    attrArray[2]='Z';
    memcpy(attrArray+3,attr.c_str(),attr.size()+1);//copy string data including \0
    return 3+attr.size()+1;
};
int bamAttrArrayWrite(const vector<char> &attr, const char* tagName, char* attrArray ) {
    attrArray[0]=tagName[0];attrArray[1]=tagName[1];
    attrArray[2]='B';
    attrArray[3]='c';
    *( (int32*) (attrArray+4))=attr.size();
    memcpy(attrArray+4+sizeof(int32),attr.data(),attr.size());//copy array data
    return 4+sizeof(int32)+attr.size();
};
int bamAttrArrayWrite(const vector<int32> &attr, const char* tagName, char* attrArray ) {
    attrArray[0]=tagName[0];attrArray[1]=tagName[1];
    attrArray[2]='B';
    attrArray[3]='i';
    *( (int32*) (attrArray+4))=attr.size();
    memcpy(attrArray+4+sizeof(int32),attr.data(),sizeof(int32)*attr.size());//copy array data
    return 4+sizeof(int32)+sizeof(int32)*attr.size();
};

int bamAttrArrayWriteSAMtags(string &attrStr, char *attrArray) {//write bam record into attrArray for string attribute attString
    size_t pos1=0, pos2=0;
    int nattr=0;
    do {//cycle over multiple tags separated by tab
        pos2 = attrStr.find('\t',pos1);
        string attr1 = attrStr.substr(pos1, pos2-pos1);
        pos1=pos2+1;
        
        if (attr1.empty())
            continue; //extra tab at the beginning, or consecutive tabs
        
        switch (attr1.at(3)) {
            case 'i':
            {
                int32 a1=stol(attr1.substr(5));
                nattr += bamAttrArrayWrite(a1,attr1.c_str(),attrArray+nattr);
                break;
            };
            case 'A':
            {
                char a1=attr1.at(5);
                nattr += bamAttrArrayWrite(a1,attr1.c_str(),attrArray+nattr);                
                break;
            };                
                break;
            case 'Z':
            {
                string a1=attr1.substr(5);
                nattr += bamAttrArrayWrite(a1,attr1.c_str(),attrArray+nattr);                
                break;
            };
            case 'f':
            {
                float a1=stof(attr1.substr(5));
                nattr += bamAttrArrayWrite(a1,attr1.c_str(),attrArray+nattr);                
                break;
            };       
        };
    } while (pos2!= string::npos);
    
    return nattr;
};

template <typename intType>
int bamAttrArrayWriteInt(intType xIn, const char* tagName, char* attrArray, Parameters &P) {//adapted from samtools
    attrArray[0]=tagName[0];attrArray[1]=tagName[1];
    #define ATTR_RECORD_INT(_intChar,_intType,_intValue) attrArray[2] = _intChar; *(_intType*)(attrArray+3) = (_intType) _intValue; return 3+sizeof(_intType)
    int64 x = (int64) xIn;
    if (x < 0) {
        if (x >= -127) {
            ATTR_RECORD_INT('c',int8_t,x);
        } else if (x >= -32767) {
            ATTR_RECORD_INT('s',int16_t,x);
        } else {
            ATTR_RECORD_INT('i',int32_t,x);
            if (!(x>=-2147483647)) {
                ostringstream errOut;
                errOut <<"EXITING because of FATAL BUG: integer out of range for BAM conversion: "<< x <<"\n";
                errOut <<"SOLUTION: contact Alex Dobin at dobin@cshl.edu\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);
            };
        };
    } else {
        if (x <= 255) {
            ATTR_RECORD_INT('C',uint8_t,x);
        } else if (x <= 65535) {
            ATTR_RECORD_INT('S',uint16_t,x);
        } else {
            ATTR_RECORD_INT('I',uint32_t,x);
            if (!(x<=4294967295)) {
                ostringstream errOut;
                errOut <<"EXITING because of FATAL BUG: integer out of range for BAM conversion: "<< x <<"\n";
                errOut <<"SOLUTION: contact Alex Dobin at dobin@cshl.edu\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);
            };
        };
    };
};




int ReadAlign::alignBAM(Transcript const &trOut, uint nTrOut, uint iTrOut, uint trChrStart, uint mateChr, uint mateStart, char mateStrand, int alignType, bool *mateMapped, vector<int> outSAMattrOrder, char** outBAMarray, uint* outBAMarrayN) {
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
    

    bool flagPaired = P.readNmates==2;

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

    for (uint imate=0;imate < (alignType<0 ? nMates:P.readNmates);imate++) {

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
            if (mateMapped!=NULL && mateMapped[imate]) continue; //this mate was mapped, do not record it as unmapped
            samFLAG=0x4;
            if (P.readNmates==2) {//paired read
                samFLAG|=0x1 + (imate==0 ? 0x40 : 0x80);
                if (mateMapped[1-imate]) {//mate mapped
                    if (trOut.Str!=(1-imate))
                    {//mate strand reverted
                       samFLAG|=0x20;
                    };
                    mateChr=trOut.Chr;
                    trChrStart=mapGen.chrStart[mateChr];
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

            if (mateMapped[1-imate])
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

            if (!P.outSAMattrRG.empty()) attrN+=bamAttrArrayWrite(P.outSAMattrRG.at(readFilesIndex),"RG",attrOutArray+attrN);
            
        } else {//this mate is mapped
            if (flagPaired) {//paired reads
                samFLAG=0x0001;
                if (iExMate==trOut.nExons-1) {//single mate
                    if (mateChr>mapGen.nChrReal) samFLAG|=0x0008; //not mapped as pair
                } else {//properly paired
                    samFLAG|=0x0002; //mapped as pair
                };
            } else {//single end
                samFLAG=0;
            };

            if (readFilter=='Y') samFLAG|=0x200; //not passing quality control

            if (alignType==-11 || alignType==-12 || alignType==-13) {
                samFLAG|=0x800; //mark chimeric alignments
            } else {//only non-chimeric alignments will be marked as non-primary, since chimeric are already marked with 0x800
                if (!trOut.primaryFlag) samFLAG|=0x100;//mark not primary align
            };

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
                trimL=clip5pNtotal[Mate];
            } else if (Str==0 && Mate==1) {
                trimL=clip3pNtotal[Mate];
            } else if (Str==1 && Mate==0) {
                trimL=clip3pNtotal[Mate];
            } else {
                trimL=clip5pNtotal[Mate];
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
                    } else if (gapG>0) {//deletion: N
                        packedCIGAR[nCIGAR++]=gapG<<BAM_CIGAR_OperationShift | BAM_CIGAR_D;
                    };
                };
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
                                rb.push_back( (int32) (trOut.exons[ii][EX_G]-mapGen.chrStart[trOut.Chr]+1) );
                                rb.push_back( (int32) (trOut.exons[ii][EX_G]-mapGen.chrStart[trOut.Chr]+trOut.exons[ii][EX_L]) );                                
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
                    default:
                        ostringstream errOut;
                        errOut <<"EXITING because of FATAL BUG: unknown/unimplemented SAM/BAM atrribute (tag): "<<outSAMattrOrder[ii] <<"\n";
                        errOut <<"SOLUTION: contact Alex Dobin at dobin@cshl.edu\n";
                        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
                };
            };
        };
            
        if (P.readFilesTypeN==10) {
//             if (readNameExtra[Mate].size()<1)
//                 cout << iReadAll <<" " <<readName <<endl; 
            attrN+=bamAttrArrayWriteSAMtags(readNameExtra[Mate],attrOutArray+attrN);
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
        } else if (mateChr<mapGen.nChrReal){
            pBAM[6]=mateChr;
        } else {
            pBAM[6]=-1;
        };

        //7: next pos 0-based leftmost pos of the next segment (= PNEXT 􀀀 1)
        if (nMates>1) {
            pBAM[7]=trOut.exons[(imate==0 ? iExMate+1 : 0)][EX_G] - trChrStart;
        } else if (mateChr<mapGen.nChrReal){
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
