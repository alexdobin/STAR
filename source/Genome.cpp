#include "Genome.h"
#include "Parameters.h"
#include "SuffixArrayFuns.h"
#include "PackedArray.h"
#include "ErrorWarning.h"
#include "streamFuns.h"
#include "SharedMemory.h"

#include <time.h>
#include <cmath>
#include <unistd.h>

//addresses with respect to shmStart of several genome values
#define SHM_sizeG 0
#define SHM_sizeSA 8
#define SHM_startG 16
// #define SHM_startSA 24
// 
// //first available byt of the shm
// #define SHM_startSHM 32


//arbitrary number for ftok function
#define SHM_projectID 23

Genome::Genome (Parameters* Pin ): P(Pin), shmStart(NULL) {
            shmKey=ftok(P->genomeDir.c_str(),SHM_projectID);
        };

void Genome::freeMemory(){//free big chunks of memory used by genome and suffix array

    if (P->genomeLoad=="NoSharedMemory") {//can deallocate only for non-shared memory
        delete[] G1;
        G1=NULL;
        SA.deallocateArray();
        SA2.deallocateArray();
        SAi.deallocateArray();
    };
};

uint Genome::OpenStream(string name, ifstream & stream)
{
    stream.open((P->genomeDir+ "/" +name).c_str(), ios::binary); 
    if (!stream.good()) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not open genome file "<< P->genomeDir << "/" << name <<"\n" << endl;
        errOut << "SOLUTION: check that the path to genome files, specified in --genomDir is correct and the files are present, and have user read permsissions\n" <<flush;     
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_GENOME_FILES, *P);
    };

    uint size = 0;

    P->inOut->logMain << "checking " << name << " size";
    stream.seekg (0, ios::end);
    size=(uint) stream.tellg();
    stream.clear();        
    stream.seekg (0, ios::beg);
    P->inOut->logMain << "file size: "<< size <<" bytes; state: good=" <<stream.good()\
            <<" eof="<<stream.eof()<<" fail="<<stream.fail()<<" bad="<<stream.bad()<<"\n"<<flush;

    return size;
}


void Genome::genomeLoad(){//allocate and load Genome
    uint *shmNG=NULL, *shmNSA=NULL;   //pointers to shm stored values , *shmSG, *shmSSA
    uint64 shmSize=0;//, shmStartG=0; shmStartSA=0;
    
    uint L=200,K=6;    
    
    time_t rawtime;
    time ( &rawtime );
    
    
    //protect some parameters
    uint sjdbOverhang1=P->sjdbOverhang;

    vector <uint> versionGenomeMin=P->versionGenome;
    P->versionGenome[0]=0;
    
    ifstream parFile((P->genomeDir+("/genomeParameters.txt")).c_str());
    if (parFile.good()) {
        P->inOut->logMain << "Reading genome generation parameters:\n";
        P->scanAllLines(parFile,3,-1);
        parFile.close();
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not open genome file "<< P->genomeDir+("/genomeParameters.txt") << endl;
        errOut << "SOLUTION: check that the path to genome files, specified in --genomDir is correct and the files are present, and have user read permsissions\n" <<flush;     
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_GENOME_FILES, *P);
    };            
    
    //check genome version
    if (P->versionGenome[0]==0) {//
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: read no value for the versionGenome parameter from genomeParameters.txt file\n";
        errOut << "SOLUTION: please re-generate genome from scratch with the latest version of STAR\n";
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_GENOME_FILES, *P);
    } else if (P->sjdbFileChrStartEnd.at(0)=="-" && P->versionGenome.at(0) >= versionGenomeMin.at(0)) {//
        P->inOut->logMain << "Genome version is compatible with current STAR version\n";
    } else if (P->sjdbFileChrStartEnd.at(0)!="-" && P->versionGenome.at(0) >= versionGenomeMin.at(1)) {//
        P->inOut->logMain << "Genome version is compatible with current STAR version\n";        
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: Genome version is INCOMPATIBLE with current STAR version\n";
        errOut << "SOLUTION: please re-generate genome from scratch with the latest version of STAR\n";
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_GENOME_FILES, *P);
    };


    if (P->sjdbOverhang==0) P->sjdbOverhang=sjdbOverhang1; //for twopass sjdbOverhang may be defined at the mapping stage
    //TODO make a more careful check, if the values are different, break!

    P->inOut->logMain << "Started loading the genome: " << asctime (localtime ( &rawtime ))<<"\n"<<flush;    
  
    ifstream GenomeIn, SAin, SAiIn;

    P->nGenome = OpenStream("Genome",GenomeIn);
    P->nSAbyte = OpenStream("SA",SAin);
    OpenStream("/SAindex",SAiIn);

    uint SAiInBytes=0;
    SAiInBytes += fstreamReadBig(SAiIn,(char*) &P->genomeSAindexNbases, sizeof(P->genomeSAindexNbases));
    P->genomeSAindexStart = new uint[P->genomeSAindexNbases+1];
    SAiInBytes += fstreamReadBig(SAiIn,(char*) P->genomeSAindexStart, sizeof(P->genomeSAindexStart[0])*(P->genomeSAindexNbases+1));  
    P->nSAi=P->genomeSAindexStart[P->genomeSAindexNbases];
    P->inOut->logMain << "Read from SAindex: genomeSAindexNbases=" << P->genomeSAindexNbases <<"  nSAi="<< P->nSAi <<endl <<flush;


    /////////////////////////////////// at this point all array sizes should be known: calculate packed array lengths
    P->GstrandBit = (uint) floor(log(P->nGenome)/log(2))+1;
    if (P->GstrandBit<32) P->GstrandBit=32; //TODO: use simple access function for SA

    P->GstrandMask = ~(1LLU<<P->GstrandBit);
    P->nSA=(P->nSAbyte*8)/(P->GstrandBit+1);
    SA.defineBits(P->GstrandBit+1,P->nSA);  
    
    P->SAiMarkNbit=P->GstrandBit+1;
    P->SAiMarkAbsentBit=P->GstrandBit+2;
    
    P->SAiMarkNmaskC=1LLU << P->SAiMarkNbit;
    P->SAiMarkNmask=~P->SAiMarkNmaskC;
    P->SAiMarkAbsentMaskC=1LLU << P->SAiMarkAbsentBit;
    P->SAiMarkAbsentMask=~P->SAiMarkAbsentMaskC;
    
    SAi.defineBits(P->GstrandBit+3,P->nSAi); 

    P->inOut->logMain << "nGenome=" << P->nGenome << ";  nSAbyte=" << P->nSAbyte <<endl<< flush;       
    P->inOut->logMain <<"GstrandBit="<<int(P->GstrandBit)<<"   SA number of indices="<<P->nSA<<endl<<flush;      
    
    shmSize=SA.lengthByte + P->nGenome+L+L+SHM_startG+8;
    shmSize+= SAi.lengthByte;                
    if (P->annotScoreScale>0) shmSize+=P->nGenome;


    if ((P->genomeLoad=="LoadAndKeep" ||
         P->genomeLoad=="LoadAndRemove" ||
         P->genomeLoad=="LoadAndExit" ||
         P->genomeLoad=="Remove") && sharedMemory == NULL)
    {
        bool unloadLast = P->genomeLoad=="LoadAndRemove";
        try
        {
            sharedMemory = new SharedMemory(shmKey, unloadLast);

            if (!sharedMemory->NeedsAllocation())
                P->inOut->logMain <<"Found genome in shared memory\n"<<flush;

            if (P->genomeLoad=="Remove") {//kill the genome and exit
                if (sharedMemory->NeedsAllocation()) {//did not find genome in shared memory, nothing to kill
                    ostringstream errOut;
                    errOut << "EXITING: Did not find the genome in memory, did not remove any genomes from shared memory\n";
                    exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_GENOME_FILES, *P);
                } else {
                    sharedMemory->Clean();
                    P->inOut->logMain <<"DONE: removed the genome from shared memory\n"<<flush;            
                    return;
                };
            }

            if (sharedMemory->NeedsAllocation())
                sharedMemory->Allocate(shmSize);
        }
        catch (const SharedMemoryException & exc)
        {
            ostringstream errOut;
            errOut << "SharedMemory Exception: " << exc.GetErrorCode() << ", errno: " << exc.GetErrorDetail() << endl;
            int err = exc.GetErrorDetail();

            int exitCode = EXIT_CODE_SHM;
            switch (exc.GetErrorCode())
            {
                case EOPENFAILED:            
                    errOut << "EXITING because of FATAL ERROR: problems with shared memory: error from shmget() or shm_open(): " << strerror(err) << "\n" << flush;
                    errOut << "SOLUTION: check shared memory settings as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory\n" <<flush;     
                case EEXISTS:
                    errOut << "EXITING: fatal error from shmget() trying to allocate shared memory piece: error type: " << strerror(err) <<"\n";
                    errOut << "Possible cause 1: not enough RAM. Check if you have enough RAM of at least " << shmSize+2000000000 << " bytes\n";
                    errOut << "Possible cause 2: not enough virtual memory allowed with ulimit. SOLUTION: run ulimit -v " <<  shmSize+2000000000 <<"\n";
                    errOut << "Possible cause 3: allowed shared memory size is not large enough. SOLUTIONS: (i) consult STAR manual on how to increase shared memory allocation; " \
                    "(ii) ask your system administrator to increase shared memory allocation; (iii) run STAR with --genomeLoad NoSharedMemory\n"<<flush;
                case EFTRUNCATE:
                    errOut << "EXITING: fatal error from ftruncate() error shared memory: error type: " << strerror(err) << endl;
                    errOut << "Possible cause 1: not enough RAM. Check if you have enough RAM of at least " << shmSize+2000000000 << " bytes\n";
                    exitCode = EXIT_CODE_MEMORY_ALLOCATION;
                case EMAPFAILED:
                    errOut << "EXITING because of FATAL ERROR: problems with shared memory: error from shmat() while trying to get address of the shared memory piece: " << strerror(err) << "\n" <<flush;
                    errOut << "SOLUTION: check shared memory settings as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory\n" <<flush;     
                case ECLOSE:
                    errOut << "EXITING because of FATAL ERROR: could not close the shared memory object: " << strerror(err) << "\n" <<flush;     
                case EUNLINK:
                    #ifdef POSIX_SHARED_MEM
                    errOut << "EXITING because of FATAL ERROR:  could not delete the shared object: " << strerror(err) <<flush;
                    #else
                    errOut << "EXITING because of FATAL ERROR: problems with shared memory: error from shmctl() while trying to remove shared memory piece: " << strerror(err) << "\n" <<flush;
                    errOut << "SOLUTION: check shared memory settings as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory\n" <<flush;     
                    #endif
                default:
                    errOut << "EXITING because of FATAL ERROR: There was an issue with the shared memory allocation.";
            }      

            try 
            {
                sharedMemory->Clean();
            }
            catch(...)
            {}

            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, exitCode, *P);
        }

        shmStart = (char*) sharedMemory->GetMapped();
        shmNG= (uint*) (shmStart+SHM_sizeG);
        shmNSA= (uint*) (shmStart+SHM_sizeSA);                          
        
        if (!sharedMemory->IsAllocator())
        {
            // genome is in shared memory or being loaded
            // wait for the process that will populate it
            // and record the sizes
        
            uint iwait=0;
            while (*shmNG != P->nGenome) {
                iwait++;
                P->inOut->logMain <<"Another job is still loading the genome, sleeping for 1 min\n" <<flush;
                sleep(60);                    
                if (iwait==100) {
                    ostringstream errOut;
                    errOut << "EXITING because of FATAL ERROR: waited too long for the other job to finish loading the genome" << strerror(errno) << "\n" <<flush;
                    errOut << "SOLUTION: remove the shared memory chunk by running STAR with --genomeLoad Remove, and restart STAR" <<flush;     
                    exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_GENOME_LOADING_WAITED_TOO_LONG, *P);                
                };
            };

            if (P->nSAbyte!=*shmNSA)
            {
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: the SA file size did not match what we found in shared memory" << "\n" << flush;
                errOut << "SOLUTION: remove the shared memory chunk by running STAR with --genomeLoad Remove, and restart STAR" << flush;     
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, *P);  
            }
            
            P->inOut->logMain << "Using shared memory for genome. key=0x" <<hex<<shmKey<<dec<< ";   shmid="<< sharedMemory->GetId() <<endl<<flush;
        }

        G1=shmStart+SHM_startG;
        SA.pointArray(G1+P->nGenome+L+L);
        char* shmNext=SA.charArray+P->nSAbyte;
        
        SAi.pointArray(shmNext);
        shmNext += SAi.lengthByte;

        if (P->annotScoreScale>0) {//optional allocation
            sigG = shmNext;
            shmNext += P->nGenome;
        }    
    }
    else if (P->genomeLoad=="NoSharedMemory")// simply allocate memory, do not use shared memory
    {
        try {
            if (P->twopass1readsN==0) {//1-pass, no extra memory
                G1=new char[P->nGenome+L+L];        
                SA.allocateArray();
            } else {//2-pass: reserve extra memory
                P->nGenome2=P->nGenome+P->twopassSJlimit*P->sjdbLength;
                SA2.defineBits(P->GstrandBit+1,P->nSA+2*P->twopassSJlimit*P->sjdbLength);
                G1=new char[P->nGenome2+L+L];        
                SA2.allocateArray();
                SA.pointArray(SA2.charArray+SA2.lengthByte-SA.lengthByte);
            };            
            SAi.allocateArray();
            P->inOut->logMain <<"Shared memory is not used for genomes. Allocated a private copy of the genome.\n"<<flush;                
        } catch (exception & exc) {
            ostringstream errOut;           
            errOut <<"EXITING: fatal error trying to allocate genome arrays, exception thrown: "<<exc.what()<<endl;
            errOut <<"Possible cause 1: not enough RAM. Check if you have enough RAM " << P->nGenome+L+L+SA.lengthByte+SAi.lengthByte+2000000000 << " bytes\n";
            errOut <<"Possible cause 2: not enough virtual memory allowed with ulimit. SOLUTION: run ulimit -v " <<  P->nGenome+L+L+SA.lengthByte+SAi.lengthByte+2000000000<<endl <<flush;
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_MEMORY_ALLOCATION, *P);            
        };

    }


//     if (twopass1readsN==0) {//not 2-pass
//         shmStartG=SHM_startSHM;
//         shmStartSA=0;
//     } else {//2-pass
//         ostringstream errOut;
//         errOut << "EXITING because of FATAL ERROR: 2-pass procedure cannot be used with genome already loaded im memory'  "\n" ;
//         errOut << "SOLUTION: check shared memory settings as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory\n" <<flush;     
//         exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_SHM, *P);
//     };
    

    G=G1+L;

    bool isAllocatorProcess = sharedMemory != NULL && sharedMemory->IsAllocator();

    if (P->genomeLoad=="NoSharedMemory" || isAllocatorProcess) {//load genome and SAs from files
        //load genome
        P->inOut->logMain <<"Genome file size: "<<P->nGenome <<" bytes; state: good=" <<GenomeIn.good()\
                <<" eof="<<GenomeIn.eof()<<" fail="<<GenomeIn.fail()<<" bad="<<GenomeIn.bad()<<"\n"<<flush;        
        P->inOut->logMain <<"Loading Genome ... " << flush;        
        uint genomeReadBytesN=fstreamReadBig(GenomeIn,G,P->nGenome);    
        P->inOut->logMain <<"done! state: good=" <<GenomeIn.good()\
                <<" eof="<<GenomeIn.eof()<<" fail="<<GenomeIn.fail()<<" bad="<<GenomeIn.bad()<<"; loaded "<<genomeReadBytesN<<" bytes\n" << flush;            
        GenomeIn.close();
        
        for (uint ii=0;ii<L;ii++) {// attach a tail with the largest symbol
            G1[ii]=K-1;
            G[P->nGenome+ii]=K-1;        
        };    
      
        //load SAs
        P->inOut->logMain <<"SA file size: "<<SA.lengthByte <<" bytes; state: good=" <<SAin.good()\
                <<" eof="<<SAin.eof()<<" fail="<<SAin.fail()<<" bad="<<SAin.bad()<<"\n"<<flush;        
        P->inOut->logMain <<"Loading SA ... " << flush;               
        genomeReadBytesN=fstreamReadBig(SAin,SA.charArray, SA.lengthByte);
        P->inOut->logMain <<"done! state: good=" <<SAin.good()\
                <<" eof="<<SAin.eof()<<" fail="<<SAin.fail()<<" bad="<<SAin.bad()<<"; loaded "<<genomeReadBytesN<<" bytes\n" << flush;            
        SAin.close();
        
        P->inOut->logMain <<"Loading SAindex ... " << flush;             
        SAiInBytes +=fstreamReadBig(SAiIn,SAi.charArray, SAi.lengthByte);
        P->inOut->logMain <<"done: "<<SAiInBytes<<" bytes\n" << flush;       
    };
    
    SAiIn.close();            

    if ((P->genomeLoad=="LoadAndKeep" || 
         P->genomeLoad=="LoadAndRemove" || 
         P->genomeLoad=="LoadAndExit") && isAllocatorProcess ) 
    {
        //record sizes. This marks the end of genome loading
        *shmNG=P->nGenome;
        *shmNSA=P->nSAbyte;
    };
    
    time ( &rawtime );
    P->inOut->logMain << "Finished loading the genome: " << asctime (localtime ( &rawtime )) <<"\n"<<flush;    
      
    #ifdef COMPILE_FOR_MAC
    {
        uint sum1=0;
        for (uint ii=0;ii<P->nGenome; ii++) sum1 +=  (uint) (unsigned char) G[ii];
        P->inOut->logMain << "Sum of all Genome bytes: " <<sum1 <<"\n"<<flush;  
        sum1=0;        
        for (uint ii=0;ii<SA.lengthByte; ii++) sum1 +=  (uint) (unsigned char) SA.charArray[ii];
        P->inOut->logMain << "Sum of all SA bytes: " <<sum1 <<"\n"<<flush;
        sum1=0;        
        for (uint ii=0;ii<SAi.lengthByte; ii++) sum1 +=  (uint) (unsigned char) SAi.charArray[ii];
        P->inOut->logMain << "Sum of all SAi bytes: " <<sum1 <<"\n"<<flush;
    };
    #endif
    
    if (P->genomeLoad=="LoadAndExit") {
    	uint shmSum=0;
    	for (uint ii=0;ii<shmSize;ii++) shmSum+=shmStart[ii];
        P->inOut->logMain << "genomeLoad=LoadAndExit: completed, the genome is loaded and kept in RAM, EXITING now.\n"<<flush;
//         system("echo `date` ..... Finished genome loading >> Log.timing.out");
//        exit(0);
    return;
    };
    
    //find chr starts from files
    P->chrInfoLoad();

    P->chrBinFill();
 
    //splice junctions database
    if (P->nGenome==P->chrStart[P->nChrReal]) {//no sjdb
        P->sjdbN=0;
        P->sjGstart=P->chrStart[P->nChrReal]+1; //not sure why I need that
    } else {//there are sjdb chromosomes
        ifstream sjdbInfo((P->genomeDir+"/sjdbInfo.txt").c_str());
        if (sjdbInfo.fail()) {
            ostringstream errOut;                            
            errOut << "EXITING because of FATAL error, could not open file " << (P->genomeDir+"/sjdbInfo.txt") <<"\n";
            errOut << "SOLUTION: check that the path to genome files, specified in --genomDir is correct and the files are present, and have user read permsissions\n" <<flush;     
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
        };
    
        
        sjdbInfo >> P->sjdbN >> P->sjdbOverhang;
        P->sjdbLength=P->sjdbOverhang*2+1;
        P->inOut->logMain << "Processing splice junctions database sjdbN=" <<P->sjdbN<<",   sjdbOverhang=" <<P->sjdbOverhang <<" \n";    
        
        P->sjChrStart=P->nChrReal;
        P->sjGstart=P->chrStart[P->sjChrStart];

        //fill the sj-db to genome translation array
        P->sjDstart=new uint [P->sjdbN];
        P->sjAstart=new uint [P->sjdbN];
        P->sjdbStart=new uint [P->sjdbN];
        P->sjdbEnd=new uint [P->sjdbN];
        
        P->sjdbMotif=new uint8 [P->sjdbN];
        P->sjdbShiftLeft=new uint8 [P->sjdbN];
        P->sjdbShiftRight=new uint8 [P->sjdbN];
        P->sjdbStrand=new uint8 [P->sjdbN];

        for (uint ii=0;ii<P->sjdbN;ii++) {//get the info about junctions from sjdbInfo.txt       
            {
                uint16 d1,d2,d3,d4;
                sjdbInfo >> P->sjdbStart[ii] >> P->sjdbEnd[ii] >> d1 >> d2 >> d3 >> d4;
                P->sjdbMotif[ii]      = (uint8) d1;
                P->sjdbShiftLeft[ii]  = (uint8) d2;
                P->sjdbShiftRight[ii] = (uint8) d3;
                P->sjdbStrand[ii] = (uint8) d4;
            };
            P->sjDstart[ii]   = P->sjdbStart[ii]  - P->sjdbOverhang; 
            P->sjAstart[ii]   = P->sjdbEnd[ii] + 1;     
            if (P->sjdbMotif[ii]==0) {//shinon-canonical junctions back to their true coordinates
                P->sjDstart[ii] += P->sjdbShiftLeft[ii];
                P->sjAstart[ii] += P->sjdbShiftLeft[ii];
            };
        };
    };     
    
    //check and redefine some parameters
    //max intron size
    if (P->alignIntronMax==0 && P->alignMatesGapMax==0) {
        P->inOut->logMain << "alignIntronMax=alignMatesGapMax=0, the max intron size will be approximately determined by (2^winBinNbits)*winAnchorDistNbins=" \
                << (1LLU<<P->winBinNbits)*P->winAnchorDistNbins <<endl;
    } else {
        //redefine winBinNbits
        P->winBinNbits=max( (uint) floor(log2(P->nGenome/40000)+0.5), (uint) floor(log2(max(max(4LLU,P->alignIntronMax),P->alignMatesGapMax)/4)+0.5) );
        P->inOut->logMain << "To accomodate alignIntronMax="<<P->alignIntronMax<<" redefined winBinNbits="<< P->winBinNbits <<endl;
    };
    
    if (P->winBinNbits > P->genomeChrBinNbits) {
       P->inOut->logMain << "winBinNbits=" <<P->winBinNbits <<" > " << "genomeChrBinNbits=" << P->genomeChrBinNbits << "   redefining:\n";
       P->winBinNbits=P->genomeChrBinNbits;
       P->inOut->logMain << "winBinNbits=" <<P->winBinNbits <<endl;
    };    
    
    
    if (P->alignIntronMax==0 && P->alignMatesGapMax==0) {
    } else {
        //redefine winFlankNbins,winAnchorDistNbins
        P->winFlankNbins=max(P->alignIntronMax,P->alignMatesGapMax)/(1LLU<<P->winBinNbits)+1;
        P->winAnchorDistNbins=2*P->winFlankNbins;
        P->inOut->logMain << "To accomodate alignIntronMax="<<P->alignIntronMax<<" and alignMatesGapMax="<<P->alignMatesGapMax<<\
                ", redefined winFlankNbins="<<P->winFlankNbins<<" and winAnchorDistNbins="<<P->winAnchorDistNbins<<endl;
    };
    
    P->winBinChrNbits=P->genomeChrBinNbits-P->winBinNbits;
    
};

