#include "genomeScanFastaFiles.h"
#include "ErrorWarning.h"
#include "SequenceFuns.h"

uint genomeScanFastaFiles (Parameters &P, char* G, bool flagRun, Genome &mapGen) {//scans fasta files. flagRun=false: check and find full size, flaRun=true: collect all the data

    uint N=0;//total number of bases in the genome, including chr "spacers"
    if (!flagRun && mapGen.chrLength.size()>0)
    {//previous chr records exist
       mapGen.chrStart.pop_back();//remove last record, it will be recorded again
       N =  mapGen.chrStart.back()+mapGen.chrLength.back();
       mapGen.chrLength.pop_back();//remove last record, it will be recorded again
    };

    ifstream fileIn;
    for (uint ii=0;ii<mapGen.pGe.gFastaFiles.size();ii++) {//all the input files
        fileIn.open(mapGen.pGe.gFastaFiles.at(ii).c_str());
        if ( !fileIn.good() )
        {//
            ostringstream errOut;
            errOut << "EXITING because of INPUT ERROR: could not open genomeFastaFile: " <<mapGen.pGe.gFastaFiles.at(ii) <<"\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };
        char cc=fileIn.peek();
        if ( !fileIn.good() )
        {//
            ostringstream errOut;
            errOut << "EXITING because of INPUT ERROR: could not read from genomeFastaFile: " <<mapGen.pGe.gFastaFiles.at(ii) <<"\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };
        if (cc!='>')
        {
            ostringstream errOut;
            errOut << "EXITING because of INPUT ERROR: the file format of the genomeFastaFile: " <<mapGen.pGe.gFastaFiles.at(ii) << " is not fasta:";
            errOut << " the first character is '" <<cc<<"' ("<< (cc+0) << "), not '>'.\n";
            errOut << " Solution: check formatting of the fasta file. Make sure the file is uncompressed (unzipped).\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };         while(!fileIn.eof()) {//read each file until eof
            string lineIn (4096,'.');
            getline(fileIn,lineIn);
            if (lineIn[0]=='>') {//new chromosome
                if (!flagRun) {
                    istringstream lineInStream (lineIn);
                    lineInStream.ignore(1,' ');
                    string chrName1;
                    lineInStream >> chrName1;
                    mapGen.chrName.push_back(chrName1);
                };

                if (!flagRun && mapGen.chrStart.size()>0) mapGen.chrLength.push_back(N-mapGen.chrStart.at(mapGen.chrStart.size()-1)); //true length of the chr

                if (N>0) {//pad the chromosomes to bins boudnaries
                    N = ( (N+1)/mapGen.genomeChrBinNbases+1 )*mapGen.genomeChrBinNbases;
                };

                if (!flagRun) {
                    mapGen.chrStart.push_back(N);
                    P.inOut->logMain << mapGen.pGe.gFastaFiles.at(ii)<<" : chr # " << mapGen.chrStart.size()-1 << "  \""<<mapGen.chrName.at(mapGen.chrStart.size()-1)<<"\" chrStart: "<<N<<"\n"<<flush;
                };
            } else {//char lines
                if (flagRun) {
                    N += convertNucleotidesToNumbersRemoveControls(lineIn.c_str(),G+N,lineIn.size());
                } else {
                    for (uint ii=0; ii<lineIn.size();ii++) {
                        if (int(lineIn.at(ii))>=32)
                            ++N;
                    };
                };
                
            };
        };
        fileIn.close();
    };


    if (!flagRun) 
        mapGen.chrLength.push_back(N-mapGen.chrStart.at(mapGen.chrStart.size()-1)); //true length of the last chr

    N = ( (N+1)/mapGen.genomeChrBinNbases+1)*mapGen.genomeChrBinNbases;

    if (!flagRun) {
        mapGen.nChrReal=mapGen.chrStart.size();
        mapGen.chrStart.push_back(N); //last chromosome end+1
        for (uint ii=0;ii<mapGen.nChrReal;ii++) {
            mapGen.chrNameIndex[mapGen.chrName[ii]]=ii;
        };
    };

    return N;
};

