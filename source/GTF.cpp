#include "GTF.h"

#include "ErrorWarning.h"
#include "streamFuns.h"
#include "TimeFunctions.h"

GTF::GTF(Genome &genome, Parameters &P, const string &dirOut, SjdbClass &sjdbLoci) 
                 : genome(genome), P(P), dirOut(dirOut), sjdbLoci(sjdbLoci), superTrome(P)
 {//initialize; load gtf file; returns number of added junctions
                     
    if (genome.sjdbOverhang==0 || genome.pGe.sjdbGTFfile=="-") {//no GTF
        gtfYes=false;
        return;
    };
    gtfYes=true;
    mkdir(dirOut.c_str(), P.runDirPerm);
    
    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ..... processing annotations GTF\n" <<flush;
    *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ..... processing annotations GTF\n" <<flush;

    std::map <string,uint64> transcriptIDnumber, geneIDnumber;    
    
    ifstream sjdbStreamIn ( genome.pGe.sjdbGTFfile.c_str() );
    if (sjdbStreamIn.fail()) {
        ostringstream errOut;
        errOut << "FATAL error, could not open file pGe.sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<"\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    if (genome.chrNameIndex.size()==0) {
        for (uint64 ii=0;ii<genome.nChrReal;ii++) {
            genome.chrNameIndex[genome.chrName[ii]]=ii;
        };
    };

    exonN=0;
    while (sjdbStreamIn.good()) {//count the number of exons
        
        string oneLine,chr1,ddd2,featureType;
        getline(sjdbStreamIn,oneLine);
        istringstream oneLineStream (oneLine);

        oneLineStream >> chr1 >> ddd2 >> featureType;
        if (chr1.substr(0,1)!="#" && featureType==genome.pGe.sjdbGTFfeatureExon) {
            exonN++;
        };
    };

    if (exonN==0)         {
        ostringstream errOut;
        errOut << "Fatal INPUT FILE error, no ""exon"" lines in the GTF file: " << genome.pGe.sjdbGTFfile <<"\n";
        errOut << "Solution: check the formatting of the GTF file, it must contain some lines with ""exon"" in the 3rd column.\n";
        errOut << "          Make sure the GTF file is unzipped.\n";
        errOut << "          If exons are marked with a different word, use --sjdbGTFfeatureExon .\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    exonLoci.resize(exonN);
    
    exonN=0;//will re-calculate
    sjdbStreamIn.clear();
    sjdbStreamIn.seekg(0,ios::beg);
    while (sjdbStreamIn.good()) {

        string oneLine,chr1,ddd2,featureType;
        getline(sjdbStreamIn,oneLine);
        istringstream oneLineStream (oneLine);

        oneLineStream >> chr1 >> ddd2 >> featureType;
        if (chr1.substr(0,1)!="#" && featureType==genome.pGe.sjdbGTFfeatureExon) {//exonic line, process

            if (genome.pGe.sjdbGTFchrPrefix!="-") chr1=genome.pGe.sjdbGTFchrPrefix + chr1;

            if (genome.chrNameIndex.count(chr1)==0) {//chr not in Genome
                P.inOut->logMain << "WARNING: while processing sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<": chromosome '"<<chr1<<"' not found in Genome fasta files for line:\n";
                P.inOut->logMain << oneLine <<"\n"<<flush;
                continue; //do not process exons/transcripts on missing chromosomes
            };

            uint64 ex1,ex2;
            char str1;
            oneLineStream >> ex1 >> ex2 >> ddd2 >> str1 >> ddd2; //read all fields except the last
            if ( ex2 > genome.chrLength[genome.chrNameIndex[chr1]] ) {
            	warningMessage("while processing sjdbGTFfile=" + genome.pGe.sjdbGTFfile + ", line:\n" + oneLine + "\n exon end = " + to_string(ex2) +  \
            			       " is larger than the chromosome " + chr1 + " length = " + to_string(genome.chrLength[genome.chrNameIndex[chr1]] ) + " , will skip this exon\n", \
							   std::cerr, P.inOut->logMain, P);
            	continue;
            };
            
            string oneLine1;
            getline(oneLineStream, oneLine1);//get the last field     
            replace(oneLine1.begin(),oneLine1.end(),';',' ');//to separate attributes
            replace(oneLine1.begin(),oneLine1.end(),'=',' ');//for GFF3 processing
            replace(oneLine1.begin(),oneLine1.end(),'\t',' ');//replace tabs
            replace(oneLine1.begin(),oneLine1.end(),'\"',' ');//now the only separator is space
            
            //string trID(""), gID(""), attr1(""),gName(""),gBiotype("");
            vector<vector<string>> exAttrNames({ {genome.pGe.sjdbGTFtagExonParentTranscript}, {genome.pGe.sjdbGTFtagExonParentGene}, genome.pGe.sjdbGTFtagExonParentGeneName, genome.pGe.sjdbGTFtagExonParentGeneType }); //trID, gID, gName, gBiotype
            vector<string> exAttr; //trID, gID, gName, gBiotype
            exAttr.resize(exAttrNames.size());
            
            for (uint32 ii=0; ii<exAttrNames.size(); ii++) {
                for (auto &attr1 : exAttrNames[ii]) {//scan through possible names
                    size_t pos1=oneLine1.find(" " + attr1 + " "); //attribute name is separated by spaces
                    if (pos1!=string::npos)
                        pos1=oneLine1.find_first_not_of(" ", pos1+attr1.size()+1);
                    if (pos1!=string::npos) {
                        exAttr[ii]=oneLine1.substr(pos1, oneLine1.find_first_of(" ",pos1)-pos1);
                    };
                };
            };

            if (exAttr[0]=="") {//no transcript ID
                P.inOut->logMain << "WARNING: while processing pGe.sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<": no transcript_id for line:\n";
                P.inOut->logMain << oneLine <<"\n"<<flush;
                exAttr[0]="tr_" + chr1 +"_"+ to_string(ex1) +"_"+ to_string(ex2) +"_"+ to_string(exonN); //unique name for the transcript
            };
            
            if (exAttr[1]=="") {//no gene ID
                P.inOut->logMain << "WARNING: while processing pGe.sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<": no gene_id for line:\n";
                P.inOut->logMain << oneLine <<"\n"<<flush;
                exAttr[1]="MissingGeneID";
            };
            
            if (exAttr[2]=="") {//no gene name
                exAttr[2]=exAttr[1];
            };
            
            if (exAttr[3]=="") {//no gene name
                exAttr[3]="MissingGeneType";
            };
            
            transcriptIDnumber.insert(std::pair <string,uint64> (exAttr[0],(uint64) transcriptIDnumber.size()));//insert new element if necessary with a new numeric value
            if (transcriptID.size() < transcriptIDnumber.size()) {//new transcript
                transcriptID.push_back(exAttr[0]);
                if (str1=='+') {
                   transcriptStrand.push_back(1);
                } else if (str1=='-') {
                   transcriptStrand.push_back(2);
                } else {
                   transcriptStrand.push_back(0);
                };
            };

            geneIDnumber.insert(std::pair <string,uint64> (exAttr[1],(uint64) geneIDnumber.size()));//insert new element if necessary with a $
            if (geneID.size() < geneIDnumber.size()) {//new gene is added
                geneID.push_back(exAttr[1]);
                geneAttr.push_back({exAttr[2],exAttr[3]});
            };

            exonLoci[exonN][exT]=transcriptIDnumber[exAttr[0]];
            exonLoci[exonN][exS]=ex1+genome.chrStart[genome.chrNameIndex[chr1]]-1;
            exonLoci[exonN][exE]=ex2+genome.chrStart[genome.chrNameIndex[chr1]]-1;
            exonLoci[exonN][exG]=geneIDnumber[exAttr[1]];
            ++exonN;
        };//if (chr1.substr(0,1)!="#" && featureType=="exon")
    };//

    if (exonN==0) {
        ostringstream errOut;
        errOut << "Fatal INPUT FILE error, no valid ""exon"" lines in the GTF file: " << genome.pGe.sjdbGTFfile <<"\n";
        errOut << "Solution: check the formatting of the GTF file. One likely cause is the difference in chromosome naming between GTF and FASTA file.\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };
    
    exonLoci.resize(exonN); //previous exonLoci.size() was an estimate, need to reszie to the actual value
    
    return;
};
