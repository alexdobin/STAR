#include "sjdbPrepare.h"
#include "ErrorWarning.h"
#include "serviceFuns.cpp"

void sjdbPrepare (SjdbClass &sjdbLoci, Parameters *P, uint nGenomeReal, string outDir, char *G, char *Gsj) {
    uint *sjdbS=new uint [sjdbLoci.chr.size()];
    uint *sjdbE=new uint [sjdbLoci.chr.size()];

    uint8 *sjdbMotif=new uint8 [sjdbLoci.chr.size()];
    uint8 *sjdbShiftLeft=new uint8 [sjdbLoci.chr.size()];
    uint8 *sjdbShiftRight=new uint8 [sjdbLoci.chr.size()];


    string chrOld="";
    uint iChr=0;
    for (uint ii=0;ii<sjdbLoci.chr.size();ii++) {
        if (chrOld!=sjdbLoci.chr.at(ii)) {//find numeric value of the chr
            for (iChr=0;iChr<P->nChrReal;iChr++) {
                if (sjdbLoci.chr.at(ii)==P->chrName[iChr]) break;
            };
            if (iChr>=P->nChrReal) {
                ostringstream errOut;
                errOut << "EXITING because of FATAL error, the sjdb chromosome " << sjdbLoci.chr.at(ii) << " is not found among the genomic chromosomes\n";
                errOut << "SOLUTION: fix your file(s) --sjdbFileChrStartEnd, offending junciton:" <<sjdbLoci.chr.at(ii)<<"\t"<<sjdbLoci.start.at(ii)<<"\t"<<sjdbLoci.end.at(ii)<<"\n";
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
            };
            chrOld=sjdbLoci.chr.at(ii);
        };

        sjdbS[ii] = sjdbLoci.start.at(ii) + P->chrStart[iChr] - 1;//sj names contain 1-based intron loci
        sjdbE[ii] = sjdbLoci.end.at(ii)   + P->chrStart[iChr] - 1;

        //motifs
        if ( G[sjdbS[ii]]==2 && G[sjdbS[ii]+1]==3 && G[sjdbE[ii]-1]==0 && G[sjdbE[ii]]==2 ) {//GTAG
            sjdbMotif[ii]=1;
        } else if ( G[sjdbS[ii]]==1 && G[sjdbS[ii]+1]==3 && G[sjdbE[ii]-1]==0 && G[sjdbE[ii]]==1 ) {//CTAC
            sjdbMotif[ii]=2;
        } else if ( G[sjdbS[ii]]==2 && G[sjdbS[ii]+1]==1 && G[sjdbE[ii]-1]==0 && G[sjdbE[ii]]==2 ) {//GCAG
            sjdbMotif[ii]=3;
        } else if ( G[sjdbS[ii]]==1 && G[sjdbS[ii]+1]==3 && G[sjdbE[ii]-1]==2 && G[sjdbE[ii]]==1 ) {//CTGC
            sjdbMotif[ii]=4;
        } else if ( G[sjdbS[ii]]==0 && G[sjdbS[ii]+1]==3 && G[sjdbE[ii]-1]==0 && G[sjdbE[ii]]==1 ) {//ATAC
            sjdbMotif[ii]=5;
        } else if ( G[sjdbS[ii]]==2 && G[sjdbS[ii]+1]==3 && G[sjdbE[ii]-1]==0 && G[sjdbE[ii]]==3 ) {//GTAT
            sjdbMotif[ii]=6;
        } else {
            sjdbMotif[ii]=0;
        };
        //repeat length: go back and forth around jR to find repeat length
        uint jjL=0,jjR=0;
        while ( jjL <= sjdbS[ii]-1 && G[sjdbS[ii]-1-jjL]==G[sjdbE[ii]-jjL] && G[sjdbS[ii]-1-jjL]<4 && jjL<255) {//go back
            jjL++;
        };
        sjdbShiftLeft[ii]=jjL;

        while ( sjdbS[ii]+jjR < nGenomeReal && G[sjdbS[ii]+jjR]==G[sjdbE[ii]+1+jjR] && G[sjdbS[ii]+jjR]<4 && jjR<255) {//go forward
            jjR++;
        };
        sjdbShiftRight[ii]=jjR;


        if (jjR==255 || jjL==255) {
            P->inOut->logMain << "WARNING: long repeat for junction # " << ii+1 <<" : " \
                    << sjdbLoci.chr.at(ii) <<" "<<sjdbS[ii] - P->chrStart[iChr] + 1 <<" "<< sjdbE[ii] - P->chrStart[iChr] + 1 \
                    << "; left shift = "<< (int) sjdbShiftLeft[ii] <<"; right shift = "<< (int) sjdbShiftRight[ii] <<"\n";
        };

        sjdbS[ii]-=sjdbShiftLeft[ii];
        sjdbE[ii]-=sjdbShiftLeft[ii];
    };

    //sort sjdb
    uint *sjdbSort=new uint [sjdbLoci.chr.size()*3];
    for (uint ii=0;ii<sjdbLoci.chr.size();ii++) {
        uint shift1=0;
        switch (sjdbLoci.str.at(ii)) {
            case '+':
                shift1=0;
                break;
            case '-':
                shift1=nGenomeReal;
                break;
            default:
                shift1=2*nGenomeReal;
        };
        sjdbSort[ii*3]=sjdbS[ii]+shift1; //separate sorting of +/- strand
        sjdbSort[ii*3+1]=sjdbE[ii]+shift1;
        sjdbSort[ii*3+2]=ii;
    };

    qsort((void *) sjdbSort, sjdbLoci.chr.size(), sizeof(uint)*3, funCompareUint2);

    uint *I=new uint [sjdbLoci.chr.size()];
    uint nsj=0;
    for (uint ii=0;ii<sjdbLoci.chr.size();ii++) {
        uint isj=sjdbSort[ii*3+2];//index of the next sorted junction

        uint isj0;
        if (nsj>0)
        {
            isj0=I[nsj-1]; //index of the last recorded junctions
        };

        if (nsj==0 || sjdbS[isj]!=sjdbS[isj0] || sjdbE[isj]!=sjdbE[isj0])
        {//different intron coordinates
            I[nsj++]=isj;// add new junction
        } else if (sjdbLoci.priority.at(isj)<sjdbLoci.priority.at(isj0))
        {//new junction has lower priority
            //do nothing, i.e. keep the old junction
        } else if (sjdbLoci.priority.at(isj)>sjdbLoci.priority.at(isj0))
        {//new junction has higher priority
            I[nsj-1]=isj;//replace the old junction
        } else if ( (sjdbMotif[isj]>0 && sjdbMotif[isj0]==0) \
              || ( ((sjdbMotif[isj]>0) == (sjdbMotif[isj0]>0)) && sjdbShiftLeft[isj]<sjdbShiftLeft[isj0]) ) {
            //new and old junctions have the same priority
            //new junction is canonical or left-most, so it wins
            I[nsj-1]=isj;//replace the old junction
            // else - keep the old junction, do not add the new one
        };
    };

    //sort again, after returning canonical junctions back to original loci:
    for (uint ii=0;ii<nsj;ii++) {
        sjdbSort[ii*3]  =sjdbS[I[ii]] + (sjdbMotif[I[ii]]==0 ? 0 : sjdbShiftLeft[I[ii]]);
        sjdbSort[ii*3+1]=sjdbE[I[ii]] + (sjdbMotif[I[ii]]==0 ? 0 : sjdbShiftLeft[I[ii]]);
        sjdbSort[ii*3+2]=I[ii];
    };

    qsort((void *) sjdbSort, nsj, sizeof(uint)*3, funCompareUint2);

    P->sjdbStart=new uint [nsj];
    P->sjdbEnd=new uint [nsj];
    P->sjdbMotif=new uint8 [nsj];
    P->sjdbShiftLeft=new uint8 [nsj];
    P->sjdbShiftRight=new uint8 [nsj];
    P->sjdbStrand=new uint8 [nsj];

    uint nsj1=0;
    for (uint ii=0;ii<nsj;ii++) {
        uint isj=sjdbSort[ii*3+2];

        if ( nsj1>0 && P->sjdbStart[nsj1-1]==sjdbSort[ii*3] && P->sjdbEnd[nsj1-1]==sjdbSort[ii*3+1] ) {//same loci on opposite strands
            uint isj0=sjdbSort[(ii-1)*3+2];

            if (sjdbLoci.priority.at(isj)<sjdbLoci.priority.at(isj0))
            {//new junction has lower priority
                continue;//keep old junction, do not add new
            } else if (sjdbLoci.priority.at(isj)>sjdbLoci.priority.at(isj0))
            {//new junction has higher priority
                nsj1--;//replace the old junction with the new one
            } else if (P->sjdbStrand[nsj1-1]>0 && sjdbLoci.str.at(isj)=='.')
            {//new junction strand is not defined
                continue;
            } else if (P->sjdbStrand[nsj1-1]==0 && sjdbLoci.str.at(isj)!='.')
            {//old junction strand is not defined
                nsj1--; //replace old with new
            } else if (P->sjdbMotif[nsj1-1]==0 && sjdbMotif[isj]==0)
            {//both are non-canonical (on opposite strand)
                P->sjdbStrand[nsj1-1]=0;//do not record new junction, keep old with undefined strand
                continue;
            } else if ( (P->sjdbMotif[nsj1-1]>0 && sjdbMotif[isj]==0) ||(P->sjdbMotif[nsj1-1]%2 == (2-P->sjdbStrand[nsj1-1])) ){//both strands defined, both junctions canonical
                //old junction is canonical, new is not, OR old junction is on correct strand
                continue;
            } else {
                //new junction is on correct strand, replace the old one
                nsj1--;
            };
        };

        //record junction
        P->sjdbStart[nsj1]=sjdbSort[ii*3];
        P->sjdbEnd[nsj1]=sjdbSort[ii*3+1];
        P->sjdbMotif[nsj1]=sjdbMotif[isj];
        P->sjdbShiftLeft[nsj1]=sjdbShiftLeft[isj];
        P->sjdbShiftRight[nsj1]=sjdbShiftRight[isj];
        if (sjdbLoci.str.at(isj)=='+') {
            P->sjdbStrand[nsj1]=1;
        } else if (sjdbLoci.str.at(isj)=='-') {
            P->sjdbStrand[nsj1]=2;
        } else {
            if (P->sjdbMotif[nsj1]==0) {//strand un-defined
                P->sjdbStrand[nsj1]=0;
            } else {
                P->sjdbStrand[nsj1]=2-P->sjdbMotif[nsj1]%2;
            };
        };
        nsj1++;
    };
    P->sjdbN=nsj1;
    P->sjDstart = new uint [P->sjdbN];
    P->sjAstart = new uint [P->sjdbN];

    ofstream sjdbInfo((outDir+"/sjdbInfo.txt").c_str());
    ofstream sjdbList ((outDir+"/sjdbList.out.tab").c_str());
    char strandChar[3]={'.','+','-'};
    //first line is some general useful information
    sjdbInfo << P->sjdbN <<"\t"<< P->sjdbOverhang <<"\n";
    uint sjGstart=0;

    for (uint ii=0;ii<P->sjdbN;ii++)
    {//add sjdb sequence to genome
        P->sjDstart[ii]   = P->sjdbStart[ii]  - P->sjdbOverhang;
        P->sjAstart[ii]   = P->sjdbEnd[ii] + 1;
        if (P->sjdbMotif[ii]==0) {//shift non-canonical junctions back to their true coordinates
            P->sjDstart[ii] += P->sjdbShiftLeft[ii];
            P->sjAstart[ii] += P->sjdbShiftLeft[ii];
        };
        memcpy(Gsj+sjGstart,G+P->sjDstart[ii],P->sjdbOverhang);//sjdbStart contains 1-based intron loci
        memcpy(Gsj+sjGstart+P->sjdbOverhang,G+P->sjAstart[ii],P->sjdbOverhang);//sjdbStart contains 1-based intron loci
        sjGstart += P->sjdbLength;
        Gsj[sjGstart-1]=GENOME_spacingChar;//spacing char between the sjdb seqs
        sjdbInfo << P->sjdbStart[ii] <<"\t"<< P->sjdbEnd[ii] <<"\t"<<(int) P->sjdbMotif[ii] <<"\t"<<(int) P->sjdbShiftLeft[ii] <<"\t"<<(int) P->sjdbShiftRight[ii]<<"\t"<<(int) P->sjdbStrand[ii] <<"\n";
        uint chr1=P->chrBin[ P->sjdbStart[ii] >> P->genomeChrBinNbits];
        sjdbList << P->chrName[chr1]<< "\t" << P->sjdbStart[ii]-P->chrStart[chr1] + 1 + (P->sjdbMotif[ii]>0 ? 0:P->sjdbShiftLeft[ii]) \
                                    << "\t"<<  P->sjdbEnd[ii]-P->chrStart[chr1] + 1 + (P->sjdbMotif[ii]>0 ? 0:P->sjdbShiftLeft[ii]) \
                                    << "\t"<< strandChar[P->sjdbStrand[ii]]<<"\n";
    };
    sjdbInfo.close();
    sjdbList.close();

};

