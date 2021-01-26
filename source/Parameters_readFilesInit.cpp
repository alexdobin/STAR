#include "Parameters.h"
#include "ErrorWarning.h"
#include "streamFuns.h"
#include <fstream>
#include <sys/stat.h>
#include "serviceFuns.cpp"

void Parameters::readFilesInit() 
{//initialize read files - but do not open yet

    if (readFilesType.at(0) == "Fastx") {
        readFilesTypeN=1;
    } else if (readFilesType.at(0) == "SAM"){
        readFilesTypeN=10;
        readFiles.samAttrKeepAll = false;
        readFiles.samAttrKeepNone = false;
        if (readFiles.samAttrKeepIn.at(0) == "All") {
            readFiles.samAttrKeepAll = true;
        } else if (readFiles.samAttrKeepIn.at(0) == "None") {
            readFiles.samAttrKeepNone = true;
        } else {
            for (auto &tag: readFiles.samAttrKeepIn) {
                if (tag.size()!=2) {
                    exitWithError("EXITING because of FATAL PARAMETER ERROR: each SAM tags in --readFilesSAMtagsKeep should contain two letters\n\
                                  SOLUTION: specify only two-letter tags in --readFilesSAMtagsKeep.",
                                  std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                };
                //array<char,2> taga = {tag[0], tag[1]};
                uint16_t tagn = * ( (uint16_t*) tag.c_str() );
                readFiles.samAttrKeep.insert(tagn);
            };
        };
    } else {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --readFilesType: "<<readFilesType.at(0) <<"\n";
        errOut <<"SOLUTION: specify one of the allowed values: Fastx or SAM\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    readFilesPrefixFinal=(readFilesPrefix=="-" ? "" : readFilesPrefix);
    
    if (readFilesManifest[0]=="-") {//no manifest, file names in readFilesIn
        readFilesNames.resize(readFilesIn.size());
        
        for (uint32 imate=0; imate<readFilesNames.size(); imate++) {
            splitString(readFilesIn[imate], ',', readFilesNames[imate]);
            if (readFilesNames[imate].back().empty()) {//extra comma at the end
                readFilesNames[imate].pop_back();
            };
        
            if (imate>0 && readFilesNames[imate].size() != readFilesNames[imate-1].size() ) {
                ostringstream errOut;
                errOut <<"EXITING: because of fatal INPUT ERROR: number of input files for mate" << imate+1 <<"="<< readFilesNames[imate].size()  <<" is not equal to that for mate"<< imate-1 <<"="<< readFilesNames[imate-1].size() <<"\n";
                errOut <<"Make sure that the number of files in --readFilesIn is the same for both mates\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            };
            
            for ( auto &fn : readFilesNames[imate] )
                fn = readFilesPrefixFinal + fn; //add prefix
        };

        readFilesN = readFilesNames[0].size();

        //read groups
        if (outSAMattrRGline.at(0)!="-") {
            string linefull;
            for (uint ii=0;ii<outSAMattrRGline.size(); ii++) {//concatenate into one line
                if (ii==0 || outSAMattrRGline.at(ii)==",") {//start new entry
                    if (ii>0) ++ii;//skip comma
                    outSAMattrRGlineSplit.push_back(outSAMattrRGline.at(ii)); //start new RG line with the first field which must be ID:xxx
                    if (outSAMattrRGlineSplit.back().substr(0,3)!="ID:") {
                        ostringstream errOut;
                        errOut <<"EXITING because of FATAL INPUT ERROR: the first word of a line from --outSAMattrRGline="<<outSAMattrRGlineSplit.back()<<" does not start with ID:xxx read group identifier\n";
                        errOut <<"SOLUTION: re-run STAR with all lines in --outSAMattrRGline starting with ID:xxx\n";
                        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                    };
                    outSAMattrRG.push_back(outSAMattrRGlineSplit.back().substr(3)); //this adds the ID field
                } else {//keep adding fields to this RG line, until the next comma
                    outSAMattrRGlineSplit.back()+="\t" + outSAMattrRGline.at(ii);
                };
            };
        };
        
        if (outSAMattrRG.size()>1 && outSAMattrRG.size()!=readFilesN) {
            ostringstream errOut;
            errOut <<"EXITING: because of fatal INPUT ERROR: number of input read files: "<< readFilesN << " does not agree with number of read group RG entries: "<< outSAMattrRG.size() <<"\n";
            errOut <<"Make sure that the number of RG lines in --outSAMattrRGline is equal to either 1, or the number of input read files in --readFilesIn\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        } else if (outSAMattrRG.size()==1) {//use the same read group for all files
            for (uint32 ifile=1; ifile<readFilesN; ifile++) {
                outSAMattrRG.push_back(outSAMattrRG.at(0));
            };
        };           
        
    } else {//read file names from manifest
        //TODO check that outSAMattrRGline and readFilesIn are not set, throw an error
        
        ifstream & rfM = ifstrOpen(readFilesManifest[0], ERROR_OUT, "SOLUTION: check the path and permissions for readFilesManifest = " + readFilesManifest[0], *this);
        inOut->logMain << "Reading input file names and read groups from readFileManifest " << readFilesManifest[0] << endl;

        readFilesNames.resize(2);
        string rfMline;
        while (getline(rfM, rfMline)) {
        	if (rfMline.find_first_not_of(" \t")>=rfMline.size())
        		continue; //skip blank lines

            uint32 itab1=0, itab2=0;
            for (uint32 imate=0; imate<2; imate++) {//SE manifest 2nd column contains "-"
                itab2=rfMline.find('\t',itab1);
                if (itab2>=rfMline.size()) {
                    ostringstream errOut;
                    errOut <<"EXITING because of FATAL INPUT FILE error: readFileManifest file " << readFilesManifest[0] <<  " has to contain at least 3 tab separated columns\n";
                    errOut <<"SOLUTION: fix the formatting of the readFileManifest file: Read1 <tab> Read2 <tab> ReadGroup. For single-end reads, use - in the 2nd column.\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_INPUT_FILES, *this);
                };
                readFilesNames[imate].push_back( readFilesPrefixFinal + rfMline.substr(itab1,itab2-itab1) );
                itab1=itab2+1;
                
                inOut->logMain << readFilesNames[imate].back() <<'\t';
            };
            
            outSAMattrRGlineSplit.push_back(rfMline.substr(itab2+1));
            
            if (outSAMattrRGlineSplit.back().substr(0,3)!="ID:")
                outSAMattrRGlineSplit.back().insert(0,"ID:");
            
            itab2=outSAMattrRGlineSplit.back().find('\t');
            outSAMattrRG.push_back(outSAMattrRGlineSplit.back().substr(3,itab2-3));
            
            inOut->logMain <<  outSAMattrRGlineSplit.back() <<'\n';
            
        };
        rfM.close();
        
        readNends = ( readFilesNames[1][0].back()=='-' ? 1 : 2);
        readFilesNames.resize(readNends);//resize if readFilesN=1
        readFilesN = readFilesNames[0].size();
    };

    inOut->logMain << "Number of fastq files for each mate = " << readFilesN << endl;
    
    readFilesCommandString="";
    if (readFilesCommand.at(0)=="-") {
        if (readFilesN>1)
            readFilesCommandString="cat   ";//concatenate multiple files
    } else {
        for (uint ii=0; ii<readFilesCommand.size(); ii++) 
            readFilesCommandString+=readFilesCommand.at(ii)+"   "; //concatenate into one string
    };    
    
    if (readFilesTypeN==1) {
        readNends=readFilesNames.size(); //for now the number of mates is defined by the number of input files
    } else if (readFilesTypeN==10) {//find the number of mates from the SAM file
        if (readFilesType.size()==2 && readFilesType.at(1)=="SE") {
            readNends=1;
        } else if (readFilesType.size()==2 && readFilesType.at(1)=="PE") {
            readNends=2;
        } else {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL INPUT ERROR: --readFilesType SAM requires specifying SE or PE reads"<<"\n";
            errOut <<"SOLUTION: specify --readFilesType SAM SE for single-end reads or --readFilesType SAM PE for paired-end reads\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };
    
    readNmates=readNends; //this may be changed later if one of the reads is barcode rea
};
