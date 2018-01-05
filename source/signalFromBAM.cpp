#include "signalFromBAM.h"
#include <iostream>
#include <regex>

void signalFromBAM(const string bamFileName, const string sigFileName, Parameters P) {

    bam1_t *bamA;
    bamA=bam_init1();

    double nMult=0, nUniq=0;

    if (P.outWigFlags.norm==1) {//count reads in the BAM file
        BGZF *bamIn=bgzf_open(bamFileName.c_str(),"r");
        bam_hdr_t *bamHeader=bam_hdr_read(bamIn);
        while ( true ) {//until the end of file
            int bamBytes1=bam_read1(bamIn, bamA);
            if (bamBytes1<0) break; //end of file
            if (bamA->core.tid<0) continue; //unmapped read
//             if ( !std::regex_match(chrName.at(bamA->core.tid),std::regex(P.outWigReferencesPrefix))) continue; //reference does not mathc required references
            if ( P.outWigReferencesPrefix!="-" && (P.outWigReferencesPrefix.compare(0,P.outWigReferencesPrefix.size(),bamHeader->target_name[bamA->core.tid],P.outWigReferencesPrefix.size())!=0) ) continue; //reference does not match required references

            uint8_t* aNHp=bam_aux_get(bamA,"NH");
            if (aNHp!=NULL) {
                uint32_t aNH=bam_aux2i(aNHp);
                if (aNH==1) {//unique mappers
                    ++nUniq;
                } else if (aNH>1) {
                    nMult+=1.0/aNH;
                };
            };
        };
        bgzf_close(bamIn);
    };

    BGZF *bamIn=bgzf_open(bamFileName.c_str(),"r");
    bam_hdr_t *bamHeader=bam_hdr_read(bamIn);

    int sigN=P.outWigFlags.strand ? 4 : 2;

    double *normFactor=new double[sigN];

    ofstream **sigOutAll=new ofstream* [sigN];

    string* sigOutFileName=new string[sigN];
    sigOutFileName[0]=sigFileName+".Unique.str1.out";
    sigOutFileName[1]=sigFileName+".UniqueMultiple.str1.out";
    if (P.outWigFlags.strand) {
        sigOutFileName[2]=sigFileName+".Unique.str2.out";
        sigOutFileName[3]=sigFileName+".UniqueMultiple.str2.out";
    };

    for (int ii=0; ii<sigN; ii++) {
        sigOutFileName[ii]+= (P.outWigFlags.format==0 ? ".bg" : ".wig");
        sigOutAll[ii]=new ofstream ( sigOutFileName[ii].c_str() );
    };

    if (P.outWigFlags.norm==0) {//raw counts
        normFactor[0]=1;
        normFactor[1]=1;
    } else if (P.outWigFlags.norm==1) {//normlaized
        normFactor[0]=1.0e6 / nUniq;
        normFactor[1]=1.0e6 / (nUniq+nMult);
        for (int is=0;is<sigN;is++) {//formatting double output
            *sigOutAll[is]<<setiosflags(ios::fixed) << setprecision(5);
        };
    };
    if (P.outWigFlags.strand) {
        normFactor[2]=normFactor[0];
        normFactor[3]=normFactor[1];
    };


    int iChr=-999;
    double *sigAll=NULL;
    uint32_t chrLen=0;
    while ( true ) {//until the end of file
        int bamBytes1=bam_read1(bamIn, bamA);
        if (bamA->core.tid!=iChr || bamBytes1<0) {
            //output to file
            if (iChr!=-999) {//iChr=-999 marks chromosomes that are not output, including unmapped reads
                for (int is=0;is<sigN;is++) {
                    if (P.outWigFlags.format==1) {
                        *sigOutAll[is] <<"variableStep chrom="<<bamHeader->target_name[iChr] <<"\n";
                    };
                    double prevSig=0;
                    for (uint32_t ig=0;ig<chrLen;ig++) {
                        double newSig=sigAll[sigN*ig+is];
                        if (P.outWigFlags.format==0) {//bedGraph
                            if (newSig!=prevSig) {
                                if (prevSig!=0) {//finish previous record
                                    *sigOutAll[is] <<ig<<"\t"<<prevSig*normFactor[is] <<"\n"; //1-based end
                                };
                                if (newSig!=0) {
                                    *sigOutAll[is] << bamHeader->target_name[iChr] <<"\t"<< ig <<"\t"; //0-based beginning
                                };
                                prevSig=newSig;
                            };
                        } else if (P.outWigFlags.format==1){//wiggle
                            if (newSig!=0) {
                                *sigOutAll[is] <<ig+1<<"\t"<<newSig*normFactor[is] <<"\n";
                            };
                        };
                    };
                };
            };
            if (bamBytes1<0) {//no more reads
                break;
            };

            iChr=bamA->core.tid;
            if ( iChr==-1 || (P.outWigReferencesPrefix!="-" && (P.outWigReferencesPrefix.compare(0,P.outWigReferencesPrefix.size(),bamHeader->target_name[bamA->core.tid],P.outWigReferencesPrefix.size())!=0) ) ) {
                iChr=-999;
                continue; //reference does not match required references
            };

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

        //by default, alignments marked as duplicate are not processed
        if ( (bamA->core.flag & 0x400) > 0 ) continue;

        //NH attribute
        uint8_t* aNHp=bam_aux_get(bamA,"NH");
        uint32_t aNH; 
        if (aNHp==NULL) {
            aNH=1; //no NH tag: assume NH=1
            //continue; //do not process lines without NH field
        } else {
            aNH=bam_aux2i(bam_aux_get(bamA,"NH")); //write a safer function allowing for lacking NH tag
        };
        if (aNH==0) continue; //do not process lines without NH=0
        uint32_t aG=bamA->core.pos;
        uint32_t iStrand=0;
        if (P.outWigFlags.strand) {//strand for stranded data from SAM flag
            iStrand= ( (bamA->core.flag & 0x10) > 0 ) == ( (bamA->core.flag & 0x80) == 0 );//0/1 for +/-
        };
        if (P.outWigFlags.type==1) {//5' of the1st read signal only, RAMPAGE/CAGE
            if ( (bamA->core.flag & 0x80)>0) continue; //skip if this the second mate
            if (iStrand==0) {
                if (aNH==1) {//unique mappers
                    sigAll[aG*sigN+0+2*iStrand]++;
                };
                sigAll[aG*sigN+1+2*iStrand]+=1.0/aNH;//U+M, normalized by the number of multi-mapping loci
                continue; //record only the first position
            };
        };

        uint32_t* cigar=(uint32_t*) (bamA->data+bamA->core.l_qname);

        for (uint32_t ic=0; ic<bamA->core.n_cigar; ic++) {
            uint32_t cigOp=(cigar[ic]<<BAM_CIGAR_LengthBits)>>BAM_CIGAR_LengthBits;
            uint32_t cigL=cigar[ic]>>BAM_CIGAR_OperationShift;
            switch (cigOp) {
                case(BAM_CIGAR_D):
                case(BAM_CIGAR_N):
                    aG+=cigL;
                    break;
                case(BAM_CIGAR_M):
                    if (P.outWigFlags.type==0 || (P.outWigFlags.type==2 && (bamA->core.flag & 0x80)>0 )) {//full signal, or second mate onyl signal
                        for (uint32_t ig=0;ig<cigL;ig++) {
                            if (aG>=chrLen) {
                                cerr << "BUG: alignment extends past chromosome in signalFromBAM.cpp\n";
                                exit(-1);
                            };
                            if (aNH==1) {//unique mappers
                                sigAll[aG*sigN+0+2*iStrand]++;
                            };
                            sigAll[aG*sigN+1+2*iStrand]+=1.0/aNH;//U+M, normalized by the number of multi-mapping loci
                            aG++;
                        };
                    } else {
                        aG+=cigL;
                    };
            };
        };
        if (P.outWigFlags.type==1) {//full signal
            --aG;
            if (aNH==1) {//unique mappers
                sigAll[aG*sigN+0+2*iStrand]++;
            };
            sigAll[aG*sigN+1+2*iStrand]+=1.0/aNH;//U+M, normalized by the number of multi-mapping loci
        };
    };
    delete [] sigAll;

    for (int is=0; is<sigN; is++) {// flush/close all signal files
        sigOutAll[is]->flush();
        sigOutAll[is]->close();
    };
};
