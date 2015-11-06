#include "genomeScanFastaFiles.h"
#include "ErrorWarning.h"


uint genomeScanFastaFiles (Parameters *P, char* G, bool flagRun) {//scans fasta files. flagRun=false: check and find full size, flaRun=true: collect all the data

    uint N=0;//total number of bases in the genome, including chr "spacers"
    if (!flagRun && P->chrLength.size()>0) 
    {//previous chr records exist
       P->chrStart.pop_back();//remove last record
       N =  P->chrStart.back()+P->chrLength.back(); 
    };

    ifstream fileIn;
    for (uint ii=0;ii<P->genomeFastaFiles.size();ii++) {//all the input files
        fileIn.open(P->genomeFastaFiles.at(ii).c_str());
        if ( !fileIn.good() ) 
        {//
            ostringstream errOut;
            errOut << "EXITING because of INPUT ERROR: could not open genomeFastaFile: " <<P->genomeFastaFiles.at(ii) <<"\n";
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);            
        };
        char cc=fileIn.peek();
        if ( !fileIn.good() ) 
        {//
            ostringstream errOut;
            errOut << "EXITING because of INPUT ERROR: could not read from genomeFastaFile: " <<P->genomeFastaFiles.at(ii) <<"\n";
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);            
        };      
        if (cc!='>')
        {
            ostringstream errOut;
            errOut << "EXITING because of INPUT ERROR: the file format of the genomeFastaFile: " <<P->genomeFastaFiles.at(ii) << "is not fasta:";
            errOut << " the first character is " <<cc<<" , not > .\n";
            errOut << " Solution: check formatting of the fasta file. Make sure the file is uncompressed (unzipped).\n";
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);            
        };         while(!fileIn.eof()) {//read each file until eof
            string lineIn (4096,'.');
            getline(fileIn,lineIn);
            if (lineIn[0]=='>') {//new chromosome
                if (!flagRun) {
                    istringstream lineInStream (lineIn);
                    lineInStream.ignore(1,' ');
                    string chrName1;
                    lineInStream >> chrName1;
                    P->chrName.push_back(chrName1);
                };
                
                if (!flagRun && P->chrStart.size()>0) P->chrLength.push_back(N-P->chrStart.at(P->chrStart.size()-1)); //true length of the chr  
                
                if (N>0) {//pad the chromosomes to bins boudnaries
                    N = ( (N+1)/P->genomeChrBinNbases+1 )*P->genomeChrBinNbases;
                };

                if (!flagRun) {
                    P->chrStart.push_back(N);    
                    P->inOut->logMain << P->genomeFastaFiles.at(ii)<<" : chr # " << P->chrStart.size()-1 << "  \""<<P->chrName.at(P->chrStart.size()-1)<<"\" chrStart: "<<N<<"\n"<<flush;
                };
            } else {//char lines
                if (flagRun) lineIn.copy(G+N,lineIn.size(),0);
                N += lineIn.size();
            };
        };
        fileIn.close();        
    };
    
   
    if (!flagRun) P->chrLength.push_back(N-P->chrStart.at(P->chrStart.size()-1)); //true length of the chr  

    N = ( (N+1)/P->genomeChrBinNbases+1)*P->genomeChrBinNbases;
        
    if (!flagRun) 
    { 
        P->nChrReal=P->chrStart.size();
        P->chrStart.push_back(N); //last chromosome end+1
        for (uint ii=0;ii<P->nChrReal;ii++) {
            P->chrNameIndex[P->chrName[ii]]=ii;
        };
    } else    
    {//convert the genome to 0,1,2,3,4
        for (uint jj=0;jj<N;jj++) {
            switch (int(G[jj])){
                case(65): case(97):  G[jj]=char(0);break;//A
                case(67): case(99):  G[jj]=char(1);break;//C           
                case(71): case(103): G[jj]=char(2);break;//G                       
                case(84): case(116): G[jj]=char(3);break;//T                                
                case(78): case(110): G[jj]=char(4);break;//N
                case(48):            G[jj]=GENOME_spacingChar;break;//chromosomal breaks within the sequences
                default:              //anything else
                    if (G[jj]!=GENOME_spacingChar) {
                         //P->inOut->logMain << "Unexpected character: char="<< G[jj] << "   int="<<int(G[jj])<<"   at " << jj << " , replacing with N\n";
                         G[jj]=char(4);                                 
                    };
            };
        }; 
    };
    
    return N;
};        

