#include "loadGTF.h"

#include "IncludeDefine.h"
#include "ErrorWarning.h"
#include "serviceFuns.cpp"
#include "streamFuns.h"
#include "TimeFunctions.h"

#include <ctime>
#include <map>


#define GTF_exonLoci_size 4
#define GTF_exonTrID(ii) ((ii)*GTF_exonLoci_size)
#define GTF_exonStart(ii) ((ii)*GTF_exonLoci_size+1)
#define GTF_exonEnd(ii) ((ii)*GTF_exonLoci_size+2)
#define GTF_exonGeID(ii) ((ii)*GTF_exonLoci_size+3)

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



uint64 loadGTF(SjdbClass &sjdbLoci, Parameters &P, string dirOut, Genome &mapGen) {//load gtf file, add junctions to P.sjdb
    //returns number of added junctions
    if (mapGen.sjdbOverhang==0 || mapGen.pGe.sjdbGTFfile=="-") //no GTF
        return 0;
        
    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ..... processing annotations GTF\n" <<flush;
    *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ..... processing annotations GTF\n" <<flush;

    vector<array<string,2>> geneAttr;
    
    ifstream sjdbStreamIn ( mapGen.pGe.sjdbGTFfile.c_str() );
    if (sjdbStreamIn.fail()) {
        ostringstream errOut;
        errOut << "FATAL error, could not open file pGe.sjdbGTFfile=" << mapGen.pGe.sjdbGTFfile <<"\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    if (mapGen.chrNameIndex.size()==0)
    {
        for (uint64 ii=0;ii<mapGen.nChrReal;ii++) {
            mapGen.chrNameIndex[mapGen.chrName[ii]]=ii;
        };
    };

    std::map <string,uint64> transcriptIDnumber, geneIDnumber;

    uint64 exonN=0;
    while (sjdbStreamIn.good()) {//count the number of exons
        string chr1,ddd2,featureType;
        sjdbStreamIn >> chr1 >> ddd2 >> featureType;
        if (chr1.substr(0,1)!="#" && featureType==mapGen.pGe.sjdbGTFfeatureExon) {
            exonN++;
        };
        sjdbStreamIn.ignore(1000000000,'\n'); //ignore the rest of the line
    };

    if (exonN==0)         {
        ostringstream errOut;
        errOut << "Fatal INPUT FILE error, no ""exon"" lines in the GTF file: " << mapGen.pGe.sjdbGTFfile <<"\n";
        errOut << "Solution: check the formatting of the GTF file, it must contain some lines with ""exon"" in the 3rd column.\n";
        errOut << "          Make sure the GTF file is unzipped.\n";
        errOut << "          If exons are marked with a different word, use --sjdbGTFfeatureExon .\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    uint64* exonLoci=new uint64 [exonN*GTF_exonLoci_size];
    char* transcriptStrand = new char [exonN];
    vector <string> transcriptID, geneID;

    exonN=0;//re-calculate
    sjdbStreamIn.clear();
    sjdbStreamIn.seekg(0,ios::beg);
    while (sjdbStreamIn.good()) {

        string oneLine,chr1,ddd2,featureType;
        getline(sjdbStreamIn,oneLine);
        istringstream oneLineStream (oneLine);

        oneLineStream >> chr1 >> ddd2 >> featureType;
        if (chr1.substr(0,1)!="#" && featureType==mapGen.pGe.sjdbGTFfeatureExon) {//exonic line, process

            if (mapGen.pGe.sjdbGTFchrPrefix!="-") chr1=mapGen.pGe.sjdbGTFchrPrefix + chr1;

            if (mapGen.chrNameIndex.count(chr1)==0) {//chr not in Genome
                P.inOut->logMain << "WARNING: while processing pGe.sjdbGTFfile=" << mapGen.pGe.sjdbGTFfile <<": chromosome '"<<chr1<<"' not found in Genome fasta files for line:\n";
                P.inOut->logMain << oneLine <<"\n"<<flush;
                continue; //do not process exons/transcripts on missing chromosomes
            };

            uint64 ex1,ex2;
            char str1;
            oneLineStream >> ex1 >> ex2 >> ddd2 >> str1 >> ddd2; //read all fields except the last

            string oneLine1;
            getline(oneLineStream, oneLine1);//get the last field
            replace(oneLine1.begin(),oneLine1.end(),';',' ');//to separate attributes
            replace(oneLine1.begin(),oneLine1.end(),'=',' ');//for GFF3 processing
            oneLineStream.str(oneLine1);
            oneLineStream.clear();

            //string trID(""), gID(""), attr1(""),gName(""),gBiotype("");
            array<string,4> exAttrNames({mapGen.pGe.sjdbGTFtagExonParentTranscript,mapGen.pGe.sjdbGTFtagExonParentGene,"gene_name", "gene_biotype"}); //trID, gID, gName, gBiotype
            array<string,4> exAttr; //trID, gID, gName, gBiotype
            while (oneLineStream.good()) {
                string attrName1("");
                oneLineStream >> attrName1;
                for (uint32 ii=0; ii<exAttrNames.size(); ii++) {
                    if (attrName1==exAttrNames[ii]) {
                        string attr1;
                        oneLineStream >> attr1;
                        attr1.erase(remove(attr1.begin(),attr1.end(),'"'),attr1.end());
                        attr1.erase(remove(attr1.begin(),attr1.end(),';'),attr1.end());
                        exAttr[ii]=attr1;
                    };
                };
            };

            if (exAttr[0]=="") {//no transcript ID
                P.inOut->logMain << "WARNING: while processing pGe.sjdbGTFfile=" << mapGen.pGe.sjdbGTFfile <<": no transcript_id for line:\n";
                P.inOut->logMain << oneLine <<"\n"<<flush;
            } else {
                transcriptIDnumber.insert(std::pair <string,uint64> (exAttr[0],(uint64) transcriptIDnumber.size()));//insert new element if necessary with a new numeric value
                if (transcriptID.size() < transcriptIDnumber.size()) transcriptID.push_back(exAttr[0]);
                if (str1=='+') {
                   transcriptStrand[transcriptIDnumber[exAttr[0]]]=1;
                } else if (str1=='-') {
                   transcriptStrand[transcriptIDnumber[exAttr[0]]]=2;
                } else {
                   transcriptStrand[transcriptIDnumber[exAttr[0]]]=0;
                };
            };

            if (exAttr[1]=="") {//no gene ID
                P.inOut->logMain << "WARNING: while processing pGe.sjdbGTFfile=" << mapGen.pGe.sjdbGTFfile <<": no gene_id for line:\n";
                P.inOut->logMain << oneLine <<"\n"<<flush;
            } else {//add gene ID if necessary
                geneIDnumber.insert(std::pair <string,uint64> (exAttr[1],(uint64) geneIDnumber.size()));//insert new element if necessary with a $
                if (geneID.size() < geneIDnumber.size()) {//new gene is added
                    geneID.push_back(exAttr[1]);
                    geneAttr.push_back({exAttr[2],exAttr[3]});
                };
            };

            exonLoci[GTF_exonTrID(exonN)]=transcriptIDnumber[exAttr[0]];
            exonLoci[GTF_exonStart(exonN)]=ex1+mapGen.chrStart[mapGen.chrNameIndex[chr1]]-1;
            exonLoci[GTF_exonEnd(exonN)]=ex2+mapGen.chrStart[mapGen.chrNameIndex[chr1]]-1;
            exonLoci[GTF_exonGeID(exonN)]=geneIDnumber[exAttr[1]];
            ++exonN;
        };//if (chr1.substr(0,1)!="#" && featureType=="exon")
    };//

    if (exonN==0) {
        ostringstream errOut;
        errOut << "Fatal INPUT FILE error, no valid ""exon"" lines in the GTF file: " << mapGen.pGe.sjdbGTFfile <<"\n";
        errOut << "Solution: check the formatting of the GTF file. Most likely cause is the difference in chromosome naming between GTF and FASTA file.\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };
    //sort exonLoci by transcript ID and exon coordinates
    qsort((void*) exonLoci, exonN, sizeof(uint64)*GTF_exonLoci_size, funCompareUint2);

    {//exon-gene data structures: exon start/end/strand/gene/transcript
        //re-sort exons by exons loci
        uint64* exgeLoci=new uint64 [exonN*GTF_exgeLoci_size]; //this also contains transcripts start and end

        for (uint64 iex=0; iex<exonN; iex++) {
            exgeLoci[GTF_exgeExStart(iex)]=exonLoci[GTF_exonStart(iex)];
            exgeLoci[GTF_exgeExEnd(iex)]=exonLoci[GTF_exonEnd(iex)];
            exgeLoci[GTF_exgeExStrand(iex)]=transcriptStrand[exonLoci[GTF_exonTrID(iex)]];
            exgeLoci[GTF_exgeGeID(iex)]=exonLoci[GTF_exonGeID(iex)];
            exgeLoci[GTF_exgeTrID(iex)]=exonLoci[GTF_exonTrID(iex)];
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
            if (iex==exonN || exonLoci[GTF_exonTrID(iex)] != exonLoci[GTF_exonTrID(trex1)]) {
                for (uint64 iex1=trex1; iex1<iex; iex1++) {//go back and fill the trend
                    extrLoci[GTF_extrTrEnd(iex1)]=exonLoci[GTF_exonEnd(iex-1)];
                };
                if (iex==exonN) break;
                trex1=iex;
            };
            extrLoci[GTF_extrTrStart(iex)]=exonLoci[GTF_exonStart(trex1)];
            extrLoci[GTF_extrTrID(iex)]=exonLoci[GTF_exonTrID(iex)];
            extrLoci[GTF_extrExStart(iex)]=exonLoci[GTF_exonStart(iex)];
            extrLoci[GTF_extrExEnd(iex)]=exonLoci[GTF_exonEnd(iex)];
            extrLoci[GTF_extrGeID(iex)]=exonLoci[GTF_exonGeID(iex)];
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
    uint64 trIDn=exonLoci[0];
    uint64 sjN=0;
    for (uint64 exI=1; exI<exonN; exI++) {
        if (trIDn==exonLoci[GTF_exonTrID(exI)]) {
            uint64 chr1=mapGen.chrBin[exonLoci[GTF_exonStart(exI)] >> mapGen.pGe.gChrBinNbits];
            if ( exonLoci[GTF_exonStart(exI)]<=exonLoci[GTF_exonEnd(exI-1)]+1 ) {
                P.inOut->logMain << "WARNING: while processing pGe.sjdbGTFfile=" << mapGen.pGe.sjdbGTFfile <<": overlapping or touching exons:\n";
                P.inOut->logMain << mapGen.chrName[chr1] <<"\t"<< exonLoci[GTF_exonStart(exI-1)]+1-mapGen.chrStart[chr1] << "\t"<< exonLoci[GTF_exonEnd(exI-1)]+1-mapGen.chrStart[chr1]  <<"\n";
                P.inOut->logMain << mapGen.chrName[chr1] <<"\t"<< exonLoci[GTF_exonStart(exI)]+1-mapGen.chrStart[chr1] << "\t"<< exonLoci[GTF_exonEnd(exI)]+1-mapGen.chrStart[chr1]  <<"\n";
            } else {
                sjLoci[sjN*sjStride]=exonLoci[GTF_exonEnd(exI-1)]+1;
                sjLoci[sjN*sjStride+1]=exonLoci[GTF_exonStart(exI)]-1;
                sjLoci[sjN*sjStride+2]=(uint64) transcriptStrand[trIDn];
                sjLoci[sjN*sjStride+3]=exonLoci[GTF_exonGeID(exI)]+1;//genes are numbered from 1
                sjN++;
            };
        } else {
            trIDn=exonLoci[GTF_exonTrID(exI)];
        };
    };

    qsort((void*) sjLoci, sjN, sizeof(uint64)*sjStride, funCompareUint2);

    char strandChar[3]={'.','+','-'};
    uint64 sjdbN1=sjdbLoci.chr.size();
    for (uint64 ii=0;ii<sjN;ii++) {
        if ( ii==0 || (sjLoci[ii*sjStride]!=sjLoci[(ii-1)*sjStride]) 
                   || (sjLoci[ii*sjStride+1]!=sjLoci[(ii-1)*sjStride+1]) 
                   || (sjLoci[ii*sjStride+2]!=sjLoci[(ii-1)*sjStride+2]) ) {
            uint64 chr1=mapGen.chrBin[sjLoci[ii*sjStride] >> mapGen.pGe.gChrBinNbits];
            sjdbLoci.chr.push_back(mapGen.chrName[chr1]);
            sjdbLoci.start.push_back(sjLoci[ii*sjStride]+1-mapGen.chrStart[chr1]);
            sjdbLoci.end.push_back(sjLoci[ii*sjStride+1]+1-mapGen.chrStart[chr1]);
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

    P.inOut->logMain << "Processing pGe.sjdbGTFfile=" << mapGen.pGe.sjdbGTFfile <<", found:\n";
    P.inOut->logMain << "\t\t"  << transcriptIDnumber.size() <<" transcripts\n" << "\t\t"  << exonN << " exons (non-collapsed)\n" << "\t\t"  << sjdbLoci.chr.size()-sjdbN1 << " collapsed junctions\n";
    time(&rawTime);
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ..... finished GTF processing\n" <<flush;

    return sjdbLoci.chr.size()-sjdbN1;
};
