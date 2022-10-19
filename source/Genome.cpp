#include "Genome.h"
#include "SuffixArrayFuns.h"
#include "PackedArray.h"
#include "ErrorWarning.h"
#include "streamFuns.h"
#include "SharedMemory.h"
#include "genomeScanFastaFiles.h"
#include "systemFunctions.h"

#include <time.h>
#include <cmath>
#include <unistd.h>
#include <sys/stat.h>

Genome::Genome (Parameters &P, ParametersGenome &pGe): shmStart(NULL), P(P), pGe(pGe), sharedMemory(NULL)
{
    struct stat stbuf;
    stat(pGe.gDir.c_str(), &stbuf);
    shmKey=stbuf.st_ino;
    genomeOut.g=this;//will change if genomeOut is different from genomeMain
    genomeOut.convYes=false;
    sjdbOverhang = pGe.sjdbOverhang; //will be re-defined later if another value was used for the generated genome
    sjdbLength = pGe.sjdbOverhang==0 ? 0 : pGe.sjdbOverhang*2+1;
};

// Genome::~Genome()
// {
//     if (sharedMemory != NULL)
//         delete sharedMemory;
//     sharedMemory = NULL;
// }

void Genome::freeMemory(){//free big chunks of memory used by genome and suffix array

    if (pGe.gLoad=="NoSharedMemory") {//can deallocate only for non-shared memory
        delete[] G1;
        G1=NULL;
        SA.deallocateArray();
        SApass2.deallocateArray();
        SAi.deallocateArray();

        P.inOut->logMain << "RAM after freeing genome index memory:\n"
                         <<  linuxProcMemory() << flush;

    };
};

uint Genome::OpenStream(string name, ifstream & stream, uint size)
{
    stream.open((pGe.gDir+ "/" +name).c_str(), ios::binary);
    if (!stream.good()) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not open genome file: "<< pGe.gDir << "/" << name <<"\n";
        errOut << "SOLUTION: check that the path to genome files, specified in --genomeDir is correct and the files are present, and have user read permissions\n" <<flush;
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };


    if (size>0) {
        P.inOut->logMain << name << ": size given as a parameter = " << size <<"\n";
    } else {
        P.inOut->logMain << "Checking " << name << " size";
        stream.seekg (0, ios::end);
        int64 size1 = stream.tellg();
        if (size1<=0) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR: failed reading from genome file: "<< pGe.gDir << "/" << name <<"\n";
            errOut << "SOLUTION: re-generate the genome index\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, 1, P);
        };
        size=(uint) size1;
        stream.clear();
        stream.seekg (0, ios::beg);
        P.inOut->logMain << "file size: "<< size <<" bytes; state: good=" <<stream.good()\
                <<" eof="<<stream.eof()<<" fail="<<stream.fail()<<" bad="<<stream.bad()<<"\n"<<flush;
    };

    return size;
};



void Genome::HandleSharedMemoryException(const SharedMemoryException & exc, uint64 shmSize)
{
    ostringstream errOut;
    errOut << "Shared memory error: " << exc.GetErrorCode() << ", errno: " << strerror(exc.GetErrorDetail()) << "(" << errno << ")" << endl;

    int exitCode = EXIT_CODE_SHM;
    switch (exc.GetErrorCode())
    {
        case EOPENFAILED:
            errOut << "EXITING because of FATAL ERROR: problems with shared memory: error from shmget() or shm_open()." << endl << flush;
            errOut << "SOLUTION: check shared memory settings as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory" << endl << flush;
            break;
        case EEXISTS:
            errOut << "EXITING: fatal error from shmget() trying to allocate shared memory piece." << endl;
            errOut << "Possible cause 1: not enough RAM. Check if you have enough RAM of at least " << shmSize+2000000000 << " bytes" << endl;
            errOut << "Possible cause 2: not enough virtual memory allowed with ulimit. SOLUTION: run ulimit -v " <<  shmSize+2000000000 << endl;
            errOut << "Possible cause 3: allowed shared memory size is not large enough. SOLUTIONS: (i) consult STAR manual on how to increase shared memory allocation; " \
            "(ii) ask your system administrator to increase shared memory allocation; (iii) run STAR with --genomeLoad NoSharedMemory" << endl<<flush;
            break;
        case EFTRUNCATE:
            errOut << "EXITING: fatal error from ftruncate() error shared memory."  << endl;
            errOut << "Possible cause 1: not enough RAM. Check if you have enough RAM of at least " << shmSize+2000000000 << " bytes" << endl << flush;
            exitCode = EXIT_CODE_MEMORY_ALLOCATION;
            break;
        case EMAPFAILED:
            errOut << "EXITING because of FATAL ERROR: problems with shared memory: error from shmat() while trying to get address of the shared memory piece." << endl << flush;
            errOut << "SOLUTION: check shared memory settings as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory" << endl << flush;
            break;
        case ECLOSE:
            errOut << "EXITING because of FATAL ERROR: could not close the shared memory object." << endl << flush;
            break;
        case EUNLINK:
            #ifdef POSIX_SHARED_MEM
            errOut << "EXITING because of FATAL ERROR:  could not delete the shared object." << endl << flush;
            #else
            errOut << "EXITING because of FATAL ERROR: problems with shared memory: error from shmctl() while trying to remove shared memory piece." << endl << flush;
            errOut << "SOLUTION: check shared memory settings as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory" << endl << flush;
            #endif
            break;
        default:
            errOut << "EXITING because of FATAL ERROR: There was an issue with the shared memory allocation. Try running STAR with --genomeLoad NoSharedMemory to avoid using shared memory.";
            break;
    }

    try
    {
        if (sharedMemory != NULL)
            sharedMemory->Clean();
    }
    catch(...)
    {}

    exitWithError(errOut.str(),std::cerr, P.inOut->logMain, exitCode, P);
};

//////////////////////////////////////////////////////////////////////////////////////////
void Genome::chrInfoLoad() {//find chrStart,Length,nChr from Genome G

    //load chr names
    ifstream chrStreamIn ( (pGe.gDir+"/chrName.txt").c_str() );
    if (chrStreamIn.fail()) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL error, could not open file " << (pGe.gDir+"/chrName.txt") <<"\n";
        errOut << "SOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    char chrInChar[1000];

    while (chrStreamIn.good()) {
        string chrIn;
        chrStreamIn.getline(chrInChar,1000);
        chrIn=chrInChar;
        if (chrIn=="") break;
        chrName.push_back(chrIn);
    };
    chrStreamIn.close();
    nChrReal=chrName.size();

    P.inOut->logMain << "Number of real (reference) chromosomes= " << nChrReal <<"\n"<<flush;
    chrStart.resize(nChrReal+1);
    chrLength.resize(nChrReal);

    //load chr lengths
    chrStreamIn.open( (pGe.gDir+"/chrLength.txt").c_str() );
    if (chrStreamIn.fail()) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL error, could not open file " << (pGe.gDir+"/chrLength.txt") <<"\n";
        errOut << "SOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    for  (uint ii=0;ii<nChrReal;ii++) {
        chrStreamIn >> chrLength[ii];
    };
    chrStreamIn.close();

    //load chr starts
    chrStreamIn.open( (pGe.gDir+"/chrStart.txt").c_str() );
    if (chrStreamIn.fail()) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL error, could not open file " << (pGe.gDir+"/chrStart.txt") <<"\n";
        errOut << "SOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    for  (uint ii=0;ii<=nChrReal;ii++) {
        chrStreamIn >> chrStart[ii];
    };
    chrStreamIn.close();

    //log
    for (uint ii=0; ii<nChrReal;ii++) {
        P.inOut->logMain << ii+1 <<"\t"<< chrName[ii] <<"\t"<<chrLength[ii]<<"\t"<<chrStart[ii]<<"\n"<<flush;
        chrNameIndex[chrName[ii]]=ii;
    };

    //chr sets
    for (auto &cm: pGe.chrSet.mitoStrings) {
        uint64 ind1 = std::find(chrName.begin(), chrName.end(), cm) - chrName.begin();
        pGe.chrSet.mito.insert(ind1);
    };

};

//////////////////////////////////////////////////////////
void Genome::chrBinFill() {
    chrBinN = chrStart[nChrReal]/genomeChrBinNbases+1;
    chrBin = new uint [chrBinN];
    for (uint ii=0, ichr=1; ii<chrBinN; ++ii) {
        if (ii*genomeChrBinNbases>=chrStart[ichr]) ichr++;
        chrBin[ii]=ichr-1;
    };
};

//////////////////////////////////////////////////////////
void Genome::genomeSequenceAllocate(uint64 nGenomeIn, uint64 &nG1allocOut, char*& Gout, char*& G1out)
{
    nG1allocOut=(nGenomeIn + 100)*2; //extra 100 bytes at the beginning, just in case
    
    if (P.limitGenomeGenerateRAM < (nG1allocOut+nG1allocOut/3)) {//allocate nG1alloc/3 for SA generation
        ostringstream errOut;
        errOut <<"EXITING because of FATAL PARAMETER ERROR: limitGenomeGenerateRAM="<< (P.limitGenomeGenerateRAM) <<"is too small for your genome\n";
        errOut <<"SOLUTION: please specify --limitGenomeGenerateRAM not less than "<< nG1allocOut+nG1allocOut/3 <<" and make that much RAM available \n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };    
    
    G1out=new char[nG1allocOut];
    Gout=G1out+100;

    memset(G1out,GENOME_spacingChar,nG1allocOut);//initialize to K-1 all bytes
};
