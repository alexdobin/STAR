#include "GTF.h"
#include "serviceFuns.cpp"
#include "streamFuns.h"

//#include <ctime>
#include <map>

#define GTF_extrLoci_size 6
#define GTF_extrTrStart(ii) ((ii)*GTF_extrLoci_size)
#define GTF_extrTrEnd(ii) ((ii)*GTF_extrLoci_size+1)
#define GTF_extrTrID(ii) ((ii)*GTF_extrLoci_size+2)
#define GTF_extrExStart(ii) ((ii)*GTF_extrLoci_size+3)
#define GTF_extrExEnd(ii) ((ii)*GTF_extrLoci_size+4)
#define GTF_extrGeID(ii) ((ii)*GTF_extrLoci_size+5)

#define GTF_exgeLoci_size 5
#define GTF_exgeExStart(ii) ((ii)*GTF_exgeLoci_size+0)
#define GTF_exgeExEnd(ii) ((ii)*GTF_exgeLoci_size+1)
#define GTF_exgeExStrand(ii) ((ii)*GTF_exgeLoci_size+2)
#define GTF_exgeGeID(ii) ((ii)*GTF_exgeLoci_size+3)
#define GTF_exgeTrID(ii) ((ii)*GTF_exgeLoci_size+4)

uint64 GTF::transcriptGeneSJ(const string &dirOut) 
{//sort exonLoci by transcript ID and exon coordinates
 //fills sjdbLoci from GTF junctions
    
    if (!gtfYes)
        return 0;
    
    exonN=exonLoci.size();
    qsort((void*) exonLoci.data(), exonN, sizeof(uint64)*exL, funCompareUint2);

    {//exon-gene data structures: exon start/end/strand/gene/transcript
        //re-sort exons by exons loci
        uint64* exgeLoci=new uint64 [exonN*GTF_exgeLoci_size]; //this also contains transcripts start and end

        for (uint64 iex=0; iex<exonN; iex++) {
            exgeLoci[GTF_exgeExStart(iex)]=exonLoci[iex][exS];
            exgeLoci[GTF_exgeExEnd(iex)]=exonLoci[iex][exE];
            exgeLoci[GTF_exgeExStrand(iex)]=transcriptStrand[exonLoci[iex][exT]];
            exgeLoci[GTF_exgeGeID(iex)]=exonLoci[iex][exG];
            exgeLoci[GTF_exgeTrID(iex)]=exonLoci[iex][exT];
        };

        qsort((void*) exgeLoci, exonN, sizeof(uint64)*GTF_exgeLoci_size, funCompareArrays<uint64,5>);

        ofstream & exgeOut = ofstrOpen(dirOut+"/exonGeTrInfo.tab",ERROR_OUT,P);
        exgeOut<<exonN<<"\n";
        for (uint64 iex=0; iex<exonN; iex++) {
             exgeOut<<exgeLoci[GTF_exgeExStart(iex)] <<"\t"<<  exgeLoci[GTF_exgeExEnd(iex)] <<"\t"<< exgeLoci[GTF_exgeExStrand(iex)] \
              <<"\t"<< exgeLoci[GTF_exgeGeID(iex)] <<"\t"<< exgeLoci[GTF_exgeTrID(iex)] <<"\n"; //the last value, transript-number, is worng here since tranascripts are re-sorted later
        };
        exgeOut.close();

        ofstream & geOut = ofstrOpen(dirOut+"/geneInfo.tab",ERROR_OUT,P);
        geOut << geneID.size() << "\n";
        for (uint64 ig=0; ig<geneID.size(); ig++) {//just geneID for now
            geOut << geneID[ig] <<"\t"<< geneAttr[ig][0] <<"\t"<< geneAttr[ig][1] <<"\n";
        };
        geOut.close();

    };

    {//exon-transcript data structures
        //re-sort transcripts by transcript start/end
        uint64* extrLoci=new uint64 [exonN*GTF_extrLoci_size]; //this also contains transcripts start and end

        uint64 trex1=0;
        for (uint64 iex=0; iex<=exonN; iex++) {
            if (iex==exonN || exonLoci[iex][exT] != exonLoci[trex1][exT]) {
                for (uint64 iex1=trex1; iex1<iex; iex1++) {//go back and fill the trend
                    extrLoci[GTF_extrTrEnd(iex1)]=exonLoci[iex-1][exE];
                };
                if (iex==exonN) break;
                trex1=iex;
            };
            extrLoci[GTF_extrTrStart(iex)]=exonLoci[trex1][exS];
            extrLoci[GTF_extrTrID(iex)]=exonLoci[iex][exT];
            extrLoci[GTF_extrExStart(iex)]=exonLoci[iex][exS];
            extrLoci[GTF_extrExEnd(iex)]=exonLoci[iex][exE];
            extrLoci[GTF_extrGeID(iex)]=exonLoci[iex][exG];
        };

        qsort((void*) extrLoci, exonN, sizeof(uint64)*GTF_extrLoci_size, funCompareArrays<uint64,5>);

        ofstream trOut ((dirOut+"/transcriptInfo.tab").c_str());
        trOut<<transcriptID.size() << "\n";
        ofstream exOut ((dirOut+"/exonInfo.tab").c_str());
        exOut<<exonN<<"\n";

        uint64 trid=extrLoci[GTF_extrTrID(0)];
        uint64 trex=0;
        uint64 trstart=extrLoci[GTF_extrTrStart(0)];
        uint64 trend=extrLoci[GTF_extrTrEnd(0)];
        uint64 exlen=0;
        for (uint64 iex=0;iex<=exonN; iex++) {
            if (iex==exonN || extrLoci[GTF_extrTrID(iex)] != trid) {//start of the new transcript
                //write out previous transcript
                trOut << transcriptID.at(trid) <<"\t"<< extrLoci[GTF_extrTrStart(iex-1)]<<"\t"<< extrLoci[GTF_extrTrEnd(iex-1)] \
                       <<"\t"<< trend << "\t"<< (uint64) transcriptStrand[trid]  <<"\t"<< iex-trex <<"\t"<<trex<<"\t"<<extrLoci[GTF_extrGeID(iex-1)]<<"\n";
                if (iex==exonN) break;
                trid=extrLoci[GTF_extrTrID(iex)];
                trstart=extrLoci[GTF_extrTrStart(iex)];
                trex=iex;
                trend=max(trend,extrLoci[GTF_extrTrEnd(iex-1)]);
                exlen=0;
            };
            exOut << extrLoci[GTF_extrExStart(iex)]-trstart <<"\t"<< extrLoci[GTF_extrExEnd(iex)]-trstart <<"\t"<< exlen <<"\n";
            exlen+=extrLoci[GTF_extrExEnd(iex)]-extrLoci[GTF_extrExStart(iex)]+1;
        };
        trOut.close();
        exOut.close();
    };

    //make junctions
    const uint64 sjStride=4;
    uint64* sjLoci = new uint64 [exonN*sjStride];
    uint64 trIDn=exonLoci[0][exT];
    uint64 sjN=0;
    for (uint64 iex=1; iex<exonN; iex++) {
        if (trIDn==exonLoci[iex][exT]) {
                uint64 chr1=genome.chrBin[exonLoci[iex][exS] >> genome.pGe.gChrBinNbits];
            if ( exonLoci[iex][exS]<=exonLoci[iex-1][exE]+1 ) {//touching - nothing to do
            } else if ( exonLoci[iex][exS]<=exonLoci[iex-1][exE] ) {//overlapping
                P.inOut->logMain << "WARNING: while processing pGe.sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<": overlapping exons:\n";
                P.inOut->logMain << genome.chrName[chr1] <<"\t"<< exonLoci[iex-1][exS]+1-genome.chrStart[chr1] << "\t"<< exonLoci[iex-1][exE]+1-genome.chrStart[chr1]  <<"\n";
                P.inOut->logMain << genome.chrName[chr1] <<"\t"<< exonLoci[iex][exS]  +1-genome.chrStart[chr1] << "\t"<< exonLoci[iex][exE]  +1-genome.chrStart[chr1]  <<"\n";
            } else {
                sjLoci[sjN*sjStride]=exonLoci[iex-1][exE]+1;
                sjLoci[sjN*sjStride+1]=exonLoci[iex][exS]-1;
                sjLoci[sjN*sjStride+2]=(uint64) transcriptStrand[trIDn];
                sjLoci[sjN*sjStride+3]=exonLoci[iex][exG]+1;//genes are numbered from 1
                sjN++;
            };
        } else {
            trIDn=exonLoci[iex][exT];
        };
    };

    qsort((void*) sjLoci, sjN, sizeof(uint64)*sjStride, funCompareUint2);

    char strandChar[3]={'.','+','-'};
    uint64 sjdbN1=sjdbLoci.chr.size();
    sjdbLoci.gene.resize(sjdbN1); //need to resize in case sjdbLoci was loaded from files without gene attribute. TODO make sure gene is always present
    for (uint64 ii=0;ii<sjN;ii++) {
        if ( ii==0 || (sjLoci[ii*sjStride]!=sjLoci[(ii-1)*sjStride])
                   || (sjLoci[ii*sjStride+1]!=sjLoci[(ii-1)*sjStride+1])
                   || (sjLoci[ii*sjStride+2]!=sjLoci[(ii-1)*sjStride+2]) ) {
            uint64 chr1=genome.chrBin[sjLoci[ii*sjStride] >> genome.pGe.gChrBinNbits];
            sjdbLoci.chr.push_back(genome.chrName[chr1]);
            sjdbLoci.start.push_back(sjLoci[ii*sjStride]+1-genome.chrStart[chr1]);
            sjdbLoci.end.push_back(sjLoci[ii*sjStride+1]+1-genome.chrStart[chr1]);
            sjdbLoci.str.push_back(strandChar[sjLoci[ii*sjStride+2]]);
            sjdbLoci.gene.push_back({sjLoci[ii*sjStride+3]});
        } else {
            sjdbLoci.gene.back().insert(sjLoci[ii*sjStride+3]);
        };
    };

    ofstream sjdbList ((dirOut+"/sjdbList.fromGTF.out.tab").c_str());
    for (uint64 ii=sjdbN1;ii<sjdbLoci.chr.size(); ii++) {
        sjdbList << sjdbLoci.chr.at(ii)<<"\t"<< sjdbLoci.start.at(ii) << "\t"<< sjdbLoci.end.at(ii)  <<"\t"<< sjdbLoci.str.at(ii);

        auto gg=sjdbLoci.gene[ii].cbegin();//iterator for genes
        sjdbList <<"\t"<< *gg;
        ++gg;
        for (; gg!=sjdbLoci.gene[ii].cend(); gg++)
            sjdbList <<","<< *gg;
        sjdbList<<"\n";
    };
    sjdbList.close();

    sjdbLoci.priority.resize(sjdbLoci.chr.size(),20);

    P.inOut->logMain << "Processing pGe.sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<", found:\n";
    P.inOut->logMain << "\t\t"  << transcriptID.size() <<" transcripts\n" << "\t\t"  << exonN << " exons (non-collapsed)\n" << "\t\t"  << sjdbLoci.chr.size()-sjdbN1 << " collapsed junctions\n";
    P.inOut->logMain << "Total junctions: " <<sjdbLoci.chr.size()<<"\n";
    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ..... finished GTF processing\n\n" <<flush;

    return sjdbLoci.chr.size()-sjdbN1;
};
