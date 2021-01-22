#include "Genome.h"
#include "SuffixArrayFuns.h"
#include "PackedArray.h"
#include "ErrorWarning.h"
#include "streamFuns.h"


void Genome::genomeOutLoad(){//allocate and load *output* Genome
    
    Parameters P1;

    ifstream parFile((pGe.gDir+("/genomeParameters.txt")).c_str());
    if (parFile.good()) {
        P.inOut->logMain << "Reading output genome generation parameters:\n";
        P1.inOut = P.inOut;
        P1.scanAllLines(parFile,3,-1);
        parFile.close();
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not open genome file "<< pGe.gDir+("/genomeParameters.txt") << endl;
        errOut << "SOLUTION: check that the path to genome files, specified in --genomeDir is correct and the files are present, and have user read permsissions\n" <<flush;
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };

    //find chr starts from files
    chrInfoLoad();
    
    
    ifstream GenomeIn;
    nGenome = OpenStream("Genome", GenomeIn, 0);
    G=new char[nGenome];
    //uint64 genomeReadBytesN = 
    fstreamReadBig(GenomeIn,G,nGenome);
    GenomeIn.close();
    
    Genome::loadSJDB(pGe.gDir);
    
    //record required genome parameters in P
    pGe.gSAindexNbases=P1.pGe.gSAindexNbases;
    pGe.gChrBinNbits=P1.pGe.gChrBinNbits;
    genomeChrBinNbases=1LLU<<pGe.gChrBinNbits;
    pGe.gSAsparseD=P1.pGe.gSAsparseD;
   
    pGe.gType=1;

    chrBinFill();

    ifstream &convStream= ifstrOpen(genomeOut.convFile, ERROR_OUT, "SOLUTION: regenerate genome.", P);
    uint32 nconv;
    convStream >> nconv >> genomeOut.nMinusStrandOffset;
    genomeOut.convBlocks.resize(nconv+1);
    for (uint32 ii=0; ii<nconv; ii++)
        convStream >> genomeOut.convBlocks[ii][0] >> genomeOut.convBlocks[ii][1] >> genomeOut.convBlocks[ii][2];
    
    //genomeOut.convBlocks[nconv][0]=genomeOut.convBlocks[nconv-1][0]+genomeOut.convBlocks[nconv-1][1];
    genomeOut.convBlocks[nconv-1][1] +=1;//increase the length of the last block so that we never reach the last base
    genomeOut.convBlocks[nconv][0]=(uint64)-1;//start of the block after last is infinity
};
