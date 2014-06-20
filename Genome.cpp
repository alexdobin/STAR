#include "Genome.h"
#include "Parameters.h"
#include "SuffixArraysFuns.h"
#include "PackedArray.h"
#include "ErrorWarning.h"
#include "streamFuns.h"
#include <time.h>
#include <cmath>
#include <unistd.h>

#define SHM_sizeG 0
#define SHM_sizeSA 8
#define SHM_startG 16
#define SHM_projectID 23

// Genome::Genome(Parameters* Pin) {
//     P=Pin;
// };

Genome::~Genome() {
    P->inOut->logMain << "--genomeLoad=" << P->genomeLoad <<" ."<<endl;
    if (P->genomeLoad=="LoadAndRemove") {//mark genome for removal after the jobs complete, if there are no other jobs attached to it
        struct shmid_ds shmStat;
        shmctl(shmID,IPC_STAT,&shmStat);
        if (shmStat.shm_nattch>1) {
            P->inOut->logMain << shmStat.shm_nattch-1 << " other job(s) are attached to the shared memory segment, will not remove it." <<endl;
        } else {
            shmctl(shmID,IPC_RMID,&shmStat);
            P->inOut->logMain <<"No other jobs are attached to the shared memory segement, removing it."<<endl;
        };
    }; 
};

void Genome::genomeLoad(){//allocate and load Genome
    shmID=0;
    bool shmLoad=false;   
    key_t shmKey=ftok(P->genomeDir.c_str(),SHM_projectID);;    
    char *shmStart=NULL;
    uint *shmNG=NULL, *shmNSA=NULL;   
    uint64 shmSize=0;

    uint L=200,K=6;    
    
    time_t rawtime;
    time ( &rawtime );
    
    
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
    } else if (P->sjdbFileChrStartEnd=="-" && P->versionGenome.at(0) >= versionGenomeMin.at(0)) {//
        P->inOut->logMain << "Genome version is compatible with current STAR version\n";
    } else if (P->sjdbFileChrStartEnd!="-" && P->versionGenome.at(0) >= versionGenomeMin.at(1)) {//
        P->inOut->logMain << "Genome version is compatible with current STAR version\n";        
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: Genome version is INCOMPATIBLE with current STAR version\n";
        errOut << "SOLUTION: please re-generate genome from scratch with the latest version of STAR\n";
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_GENOME_FILES, *P);
    };

    P->inOut->logMain << "Started loading the genome: " << asctime (localtime ( &rawtime ))<<"\n"<<flush;    
  
    
    ifstream GenomeIn((P->genomeDir+"/Genome").c_str(), ios::binary); 
    if (!GenomeIn.good()) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not open genome file "<< P->genomeDir<<"/Genome" <<"\n" << endl;
        errOut << "SOLUTION: check that the path to genome files, specified in --genomDir is correct and the files are present, and have user read permsissions\n" <<flush;     
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_GENOME_FILES, *P);
    };
    
    ifstream SAin((P->genomeDir + "/SA").c_str(), ios::binary);
    if (!SAin.good()) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not open genome file "<< P->genomeDir<<"/SA" <<"\n" << endl;
        errOut << "SOLUTION: check that the path to genome files, specified in --genomDir is correct and the files are present, and have user read permsissions\n" <<flush;     
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_GENOME_FILES, *P);
    };
    
    ifstream SAiIn((P->genomeDir+"/SAindex").c_str(),ios::binary);
    if (!SAiIn.good()) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not open genome file "<< P->genomeDir<<"/SAindex" <<"\n" << endl;
        errOut << "SOLUTION: check that the path to genome files, specified in --genomDir is correct and the files are present, and have user read permsissions\n" <<flush;     
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_GENOME_FILES, *P);
    };

    uint SAiInBytes=0;
    SAiInBytes += fstreamReadBig(SAiIn,(char*) &P->genomeSAindexNbases, sizeof(P->genomeSAindexNbases));
    P->genomeSAindexStart = new uint[P->genomeSAindexNbases+1];
    SAiInBytes += fstreamReadBig(SAiIn,(char*) P->genomeSAindexStart, sizeof(P->genomeSAindexStart[0])*(P->genomeSAindexNbases+1));  
    P->nSAi=P->genomeSAindexStart[P->genomeSAindexNbases];
    P->inOut->logMain << "Read from SAindex: genomeSAindexNbases=" << P->genomeSAindexNbases <<"  nSAi="<< P->nSAi <<endl <<flush;
    
    //search for the genome in shared memory, if found, shmLoad is false, shmID is the ID
    if (P->genomeLoad=="LoadAndKeep" || P->genomeLoad=="LoadAndRemove" || P->genomeLoad=="Remove" || P->genomeLoad=="LoadAndExit") {// find shared memory fragment
        shmID=shmget(shmKey,0,0);
        shmLoad=(shmID==-1);
        if (shmLoad && errno !=ENOENT) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR: problems with shared memory: error from shmget():" << strerror(errno) << "\n" <<flush;
            errOut << "SOLUTION: check shared memory settigns as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory\n" <<flush;     
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_SHM, *P);
        };
        if (!shmLoad) P->inOut->logMain <<"Found genome in shared memory\n"<<flush;
    };
    
    if (P->genomeLoad=="Remove") {//kill the genome and exit
        if (shmLoad) {//did not find genome in shared memory, nothing to kill
            ostringstream errOut;
            errOut << "EXITING: Did not find the genome in memory, did not remove any genomes from shared memory\n";
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_GENOME_FILES, *P);
        } else {
            struct shmid_ds *buf=NULL;
            int shmStatus=shmctl(shmID,IPC_RMID,buf);
            if (shmStatus==-1) {
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: problems with shared memory: error from shmctl() while trying to remove shared memory piece:" << strerror(errno) << "\n" <<flush;
                errOut << "SOLUTION: check shared memory settigns as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory\n" <<flush;     
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_SHM, *P);
            };            
            P->inOut->logMain <<"DONE: removed the genome from shared memory\n"<<flush;            
            exit(0);
        };
    } else if (P->genomeLoad=="NoSharedMemory" || shmLoad) {//find the size of the genome and SAs from files - the genome is not in shared memory
     
        GenomeIn.seekg (0, ios::end);
        P->nGenome=(uint) GenomeIn.tellg();
        GenomeIn.clear();        
        GenomeIn.seekg (0, ios::beg);
        P->inOut->logMain <<"Genome file size: "<<P->nGenome <<" bytes; state: good=" <<GenomeIn.good()\
                <<" eof="<<GenomeIn.eof()<<" fail="<<GenomeIn.fail()<<" bad="<<GenomeIn.bad()<<"\n"<<flush;

        SAin.seekg (0, ios::end);
        P->nSAbyte=(uint) (SAin.tellg());
        GenomeIn.clear();                
        SAin.seekg (0, ios::beg);
        P->inOut->logMain <<"SA file size: "<<P->nGenome <<" bytes; state: good=" <<SAin.good()\
                <<" eof="<<SAin.eof()<<" fail="<<SAin.fail()<<" bad="<<SAin.bad()<<"\n"<<flush;
        
    } else {//genome is in shared memory, attach it and record the sizes
    
        shmStart = (char*) shmat(shmID, NULL, 0);
        if (shmStart==((void *) -1)) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR: problems with shared memory: error from shmat() while trying to get address of the shared memory piece:" << strerror(errno) << "\n" <<flush;
            errOut << "SOLUTION: check shared memory settigns as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory\n" <<flush;     
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_SHM, *P);
        };          
        shmNG= (uint*) (shmStart+SHM_sizeG);
        shmNSA= (uint*) (shmStart+SHM_sizeSA);       
   
        uint iwait=0;
        while (*shmNG==0) {
            iwait++;
            P->inOut->logMain <<"Another job is still loading the genome, sleeping for 1 min\n" <<flush;
            sleep(60);                    
            if (iwait==100) {
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: waited too long for the other job to finish loading the genome" << strerror(errno) << "\n" <<flush;
                errOut << "SOLUTION: remove the shared memory chunk with ipcrm, or by running STAR with --genomeLoad Remove, and restart STAR" <<flush;     
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_GenomeLoadingWaitedTooLong, *P);                
            };
        };

        P->nGenome=*shmNG;
        P->nSAbyte=*shmNSA;
        
        P->inOut->logMain << "Using shared memory for genome. key=0x" <<hex<<shmKey<<dec<< ";   shmid="<<shmID<<endl<<flush;
    };

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
    
    
    /////////////////////////////////////// allocate arrays
    if (P->genomeLoad=="NoSharedMemory") {// simply allocate memory, do not use shared memory
        try {
            G1=new char[P->nGenome+L+L];        
            SA.allocateArray();
            if (P->annotScoreScale>0) sigG=new char[P->nGenome];
            SAi.allocateArray();
            P->inOut->logMain <<"Shared memory is not used for genomes. Allocated a private copy of the genome.\n"<<flush;
        } 
        catch (exception & exc) {
            ostringstream errOut;           
            errOut <<"EXITING: fatal error trying to allocate genome arrays, exception thrown: "<<exc.what()<<endl;
            errOut <<"Possible cause 1: not enough RAM. Check if you have enough RAM " << P->nGenome+L+L+SA.lengthByte+SAi.lengthByte+2000000000 << " bytes\n";
            errOut <<"Possible cause 2: not enough virtual memory allowed with ulimit. SOLUTION: run ulimit -v " <<  P->nGenome+L+L+SA.lengthByte+SAi.lengthByte+2000000000<<endl <<flush;
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_MemoryAllocation, *P);            
        };
    } else {//using shared memeory
        
        if (shmLoad) {//genome was not in shared memory: allocate shm
            shmSize=SA.lengthByte + P->nGenome+L+L+SHM_startG+8;
            shmSize+= SAi.lengthByte;                
            if (P->annotScoreScale>0) shmSize+=P->nGenome;
            shmID = shmget(shmKey, shmSize, IPC_CREAT | SHM_NORESERVE | 0666); //        shmID = shmget(shmKey, shmSize, IPC_CREAT | SHM_NORESERVE | SHM_HUGETLB | 0666);
            if (shmID < 0) {
                ostringstream errOut;
                errOut <<"EXITING: fatal error from shmget() trying to allocate shared memory piece: error type" << strerror(errno) <<"\n";
                errOut <<"Possible cause 1: not enough RAM. Check if you have enough RAM of at least" << shmSize+2000000000 << " bytes\n";
                errOut <<"Possible cause 2: not enough virtual memory allowed with ulimit. SOLUTION: run ulimit -v " <<  shmSize+2000000000 <<"\n";
                errOut <<"Possible cause 3: allowed shared memory size is not large enough. SOLUTIONS: (i) consult STAR manual on how to increase shared memory allocation; " \
                         "(ii) ask your system administrator to increase shared memory allocation; (iii) run STAR with --genomeLoad NoSharedMemory\n"<<flush;
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_MemoryAllocation, *P);            

            };
            shmStart = (char*) shmat(shmID, NULL, 0);
            shmNG= (uint*) (shmStart+SHM_sizeG);
            shmNSA= (uint*) (shmStart+SHM_sizeSA);                          
        };
        
        G1=shmStart+SHM_startG;
        SA.pointArray(G1+P->nGenome+L+L);
        char* shmNext=SA.charArray+P->nSAbyte;
        
        SAi.pointArray(shmNext);
        shmNext += SAi.lengthByte;

        if (P->annotScoreScale>0) {//optional allocation
            sigG = shmNext;
            shmNext += P->nGenome;
        };        
    };

    G=G1+L;

    if (P->genomeLoad=="NoSharedMemory" || shmLoad) {//load genome and SAs from files
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

    if (shmLoad && (P->genomeLoad=="LoadAndKeep" || P->genomeLoad=="LoadAndRemove" || P->genomeLoad=="LoadAndExit") ) {//record sizes. This marks the end of genome loading
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
        exit(0);
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


