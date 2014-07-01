#include "signalFromBAM.h"
#include <iostream>
void signalFromBAM(const string bamFileName, const string sigFileName, const bool flagStranded, const int mapQunique, const Stats &mapStat, const uint nMates) {
    BGZF *bamIn=bgzf_open(bamFileName.c_str(),"r");
    bam_hdr_t *bamHeader=bam_hdr_read(bamIn);
    
    int sigN=flagStranded ? 4 : 2;
    double *normFactor=new double[sigN];
    
    ofstream **sigOutAll=new ofstream* [sigN];
    
    sigOutAll[0]=new ofstream ( (sigFileName+".Unique.str1.out.bg").c_str() );
    normFactor[0]=1.0e6 / mapStat.mappedReadsU / nMates;
    sigOutAll[1]=new ofstream ( (sigFileName+".UniqueMultiple.str1.out.bg").c_str() );
    normFactor[1]=1.0e6 / (mapStat.mappedReadsU+mapStat.mappedReadsM) / nMates;    
    if (flagStranded) {
        sigOutAll[2]=new ofstream( (sigFileName+".Unique.str2.out.bg").c_str() );
        normFactor[2]=normFactor[0];
        sigOutAll[3]=new ofstream( (sigFileName+".UniqueMultiple.str2.out.bg").c_str() );
        normFactor[3]=normFactor[1];
    };

    for (uint32_t is=0;is<sigN;is++) {//formatting double output
        *sigOutAll[is]<<setiosflags(ios::fixed) << setprecision(5);
    };
    bam1_t *bamA;
    bamA=bam_init1();

    
    int iChr=-1;
    double *sigAll=NULL;
    while ( true ) {//until the end of file
        int bamBytes1=bam_read1(bamIn, bamA);
        uint32_t chrLen;
        if (bamA->core.tid!=iChr || bamBytes1<0) {
            //output to file
            if (iChr!=-1) {//otherwise nothing was processed yet
                for (uint32_t is=0;is<sigN;is++) {
                    double prevSig=0;
                    for (uint32_t ig=0;ig<chrLen;ig++) { 
                        double newSig=sigAll[sigN*ig+is];
                        if (newSig!=prevSig) {
                            if (prevSig!=0) {//finish previous record
                                *sigOutAll[is] <<ig<<"\t"<<prevSig*normFactor[is] <<"\n"; //1-based end
                            };
                            if (newSig!=0) {
                                *sigOutAll[is] << bamHeader->target_name[iChr] <<"\t"<< ig <<"\t"; //0-based beginning
                            };
                            prevSig=newSig;
                        };
                    };
                };
            };
            if (bamBytes1<0) {//no more reads
                break;
            };
            
            iChr=bamA->core.tid;
            chrLen=bamHeader->target_len[iChr]+1;//one extra base at the end which sohuld always be 0            
            delete [] sigAll;        
            sigAll= new double[sigN*chrLen];
            memset(sigAll, 0, sizeof(*sigAll)*sigN*chrLen);
        };
        
//         uint32_t nCigar =(bamA->core.flag<<16)>>16;
//         uint32_t mapFlag=bamA->core.flag>>16;
//         uint32_t mapQ=(bamA->core.flag<<16)>>24;
        
        #define BAM_CIGAR_OperationShift 4
        #define BAM_CIGAR_LengthBits 28
        #define BAM_CIGAR_M 0
        #define BAM_CIGAR_I 1
        #define BAM_CIGAR_D 2
        #define BAM_CIGAR_N 3
        #define BAM_CIGAR_S 4
        #define BAM_CIGAR_H 5
        #define BAM_CIGAR_P 6
        #define BAM_CIGAR_EQ 7
        #define BAM_CIGAR_X 8
        
        //NH attribute
        uint32_t aNH=bam_aux2i(bam_aux_get(bamA,"NH")); //write a safer function allowing for lacking NH tag

        uint32_t* cigar=(uint32_t*) (bamA->data+bamA->core.l_qname);
        uint32_t aG=bamA->core.pos;
        for (uint32_t ic=0; ic<bamA->core.n_cigar; ic++) {
            uint32_t cigOp=(cigar[ic]<<BAM_CIGAR_LengthBits)>>BAM_CIGAR_LengthBits;
            uint32_t cigL=cigar[ic]>>BAM_CIGAR_OperationShift;
            switch (cigOp) {
                case(BAM_CIGAR_D):
                case(BAM_CIGAR_N):
                    aG+=cigL;
                    break;
                case(BAM_CIGAR_M):
                    for (uint32_t ig=0;ig<cigL;ig++) {
                        if (aG>=chrLen) {
                            cerr << "BUG: alignment extends past chromosome in signalFromBAM.cpp\n";
                            exit(-1);
                        };
                        uint32_t iStrand=0;
                        if (flagStranded) {//strand for stranded data from SAM flag
                            iStrand=2*( ( (bamA->core.flag & 0x10) > 0 ) == ( (bamA->core.flag & 0x80) == 0 ) );
                        };                       
                        
//                         if (bamA->core.qual>=mapQunique) {//unique mappers
                        if (aNH==1) {//unique mappers
                            sigAll[aG*sigN+0+iStrand]++;
                        };
                        sigAll[aG*sigN+1+iStrand]+=1.0/aNH;//U+M, need to normalize by the number of loci
                        
                        aG++;
                    };
            };
        };
    };
    delete [] sigAll;        

    for (int is=0; is<sigN; is++) {// flush/close all signal files
        sigOutAll[is]->flush();
        sigOutAll[is]->close();
    };
};
