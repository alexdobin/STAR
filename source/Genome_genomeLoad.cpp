#include <math.h>
#include "Genome.h"
#include "SuffixArrayFuns.h"
#include "PackedArray.h"
#include "ErrorWarning.h"
#include "streamFuns.h"
#include "SharedMemory.h"
#include "genomeScanFastaFiles.h"

//addresses with respect to shmStart of several genome values
#define SHM_sizeG 0
#define SHM_sizeSA 8
#define SHM_startG 16
// #define SHM_startSA 24
//
// //first available byt of the shm
// #define SHM_startSHM 32

void Genome::genomeLoad(){//allocate and load Genome

    time_t rawtime;
    time ( &rawtime );
    *(P.inOut->logStdOut) << timeMonthDayTime(rawtime) << " ..... loading genome\n" <<flush;

    uint *shmNG=NULL, *shmNSA=NULL;   //pointers to shm stored values , *shmSG, *shmSSA
    uint64 shmSize=0;//, shmStartG=0; shmStartSA=0;

    uint L=200,K=6;

    Parameters P1;

    //some initializations before reading the parameters
    GstrandBit=0;

    ifstream parFile((pGe.gDir+("/genomeParameters.txt")).c_str());
    if (parFile.good()) {
        P.inOut->logMain << "Reading genome generation parameters:\n";

        //read genome internal parameters
        while (parFile.good()) {
            string word1;
            parFile >> word1;
            if (word1=="###") {
                parFile >> word1;
                if (word1=="GstrandBit") {
                    uint gsb1=0;
                    parFile >> gsb1;
                    GstrandBit=(uint8) gsb1;
                    P.inOut->logMain << "### GstrandBit=" << (uint) GstrandBit <<"\n";
                } else {
                    P.inOut->logMain << "### " <<word1;
                    getline(parFile,word1);
                    P.inOut->logMain <<word1<<"\n";
                };
            };
        };
        parFile.clear();
        parFile.seekg(0,ios::beg);//rewind


        P1.inOut = P.inOut;
        P1.scanAllLines(parFile,3,-1);
        parFile.close();
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not open genome file "<< pGe.gDir+("/genomeParameters.txt") << endl;
        errOut << "SOLUTION: check that the path to genome files, specified in --genomeDir is correct and the files are present, and have user read permsissions\n" <<flush;
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };

    //check genome version
    if (P1.versionGenome.size()==0) {//
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: read no value for the versionGenome parameter from genomeParameters.txt file\n";
        errOut << "SOLUTION: please re-generate genome from scratch with the latest version of STAR\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    } else if (P1.versionGenome == P.versionGenome) {//
        P.inOut->logMain << "Genome version is compatible with current STAR\n";
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: Genome version: " << P1.versionGenome << " is INCOMPATIBLE with running STAR version: "<< STAR_VERSION <<"\n";
        errOut << "SOLUTION: please re-generate genome from scratch with running version of STAR, or with version: " << P.versionGenome <<"\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };

    //find chr starts from files
    chrInfoLoad();

    //check if sjdbInfo.txt exists => genome was generated with junctions
    bool sjdbInfoExists=false;
    struct stat sjdb1;
    if ( stat( (pGe.gDir+"/sjdbInfo.txt").c_str(), &sjdb1) == 0 ) {//file exists
        sjdbInfoExists=true;
    };

    if ( P.sjdbInsert.yes && sjdbInfoExists && P1.pGe.sjdbInsertSave=="") {
        //if sjdbInsert, and genome had junctions, and genome is old - it should be re-generated with new STAR
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: old Genome is INCOMPATIBLE with on the fly junction insertion\n";
        errOut << "SOLUTION: please re-generate genome from scratch with the latest version of STAR\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };

    //record required genome parameters in P
    pGe.gSAindexNbases=P1.pGe.gSAindexNbases;
    pGe.gChrBinNbits=P1.pGe.gChrBinNbits;
    genomeChrBinNbases=1LLU<<pGe.gChrBinNbits;
    pGe.gSAsparseD=P1.pGe.gSAsparseD;

    if (P1.pGe.gFileSizes.size()>0){//genomeFileSize was recorded in the genomeParameters file, copy the values to P
        pGe.gFileSizes = P1.pGe.gFileSizes;
    };

    if (P.parArray.at(pGe.sjdbOverhang_par)->inputLevel==0 && P1.pGe.sjdbOverhang>0) {
        //if --sjdbOverhang was not defined by user and it was defined >0 at the genome generation step, then use pGe.sjdbOverhang from the genome generation step
        pGe.sjdbOverhang=P1.pGe.sjdbOverhang;
        P.inOut->logMain << "--sjdbOverhang = " << pGe.sjdbOverhang << " taken from the generated genome\n";
    } else if (sjdbInfoExists && P.parArray.at(pGe.sjdbOverhang_par)->inputLevel>0 && pGe.sjdbOverhang!=P1.pGe.sjdbOverhang) {
        //if pGe.sjdbOverhang was defined at the genome generation step,the mapping step value has to agree with it
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: present --sjdbOverhang="<<pGe.sjdbOverhang << " is not equal to the value at the genome generation step ="<< P1.pGe.sjdbOverhang << "\n";
        errOut << "SOLUTION: \n" <<flush;
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };
    sjdbOverhang = pGe.sjdbOverhang;
    sjdbLength = pGe.sjdbOverhang==0 ? 0 : pGe.sjdbOverhang*2+1;
    
    pGe.gType=1;
    if (P1.pGe.gTypeString=="SuperTranscriptome")
        pGe.gType=101;

    P.inOut->logMain << "Started loading the genome: " << asctime (localtime ( &rawtime ))<<"\n"<<flush;

    ifstream GenomeIn, SAin, SAiIn;

    if (pGe.gFileSizes.size() < 2) {//no size info available
        pGe.gFileSizes.push_back(0);
        pGe.gFileSizes.push_back(0);
    };
    nGenome = OpenStream("Genome",GenomeIn,pGe.gFileSizes.at(0));
    nSAbyte = OpenStream("SA",SAin,pGe.gFileSizes.at(1));
    OpenStream("SAindex",SAiIn,1); //we do not need SAiIn siz, using a dummy value here to prevent from reading its size from the disk

    uint SAiInBytes=0;
    SAiInBytes += fstreamReadBig(SAiIn,(char*) &pGe.gSAindexNbases, sizeof(pGe.gSAindexNbases));
    genomeSAindexStart = new uint[pGe.gSAindexNbases+1];
    SAiInBytes += fstreamReadBig(SAiIn,(char*) genomeSAindexStart, sizeof(genomeSAindexStart[0])*(pGe.gSAindexNbases+1));
    nSAi=genomeSAindexStart[pGe.gSAindexNbases];
    P.inOut->logMain << "Read from SAindex: pGe.gSAindexNbases=" << pGe.gSAindexNbases <<"  nSAi="<< nSAi <<endl;

    /////////////////////////////////// at this point all array sizes should be known: calculate packed array lengths
    if (GstrandBit==0) {//not defined before
        GstrandBit = (uint) floor(log(nGenome)/log(2))+1;
        if (GstrandBit<32) 
            GstrandBit=32; //TODO: use simple access function for SA
    };

    GstrandMask = ~(1LLU<<GstrandBit);
    nSA=(nSAbyte*8)/(GstrandBit+1);
    SA.defineBits(GstrandBit+1,nSA);

    SAiMarkNbit=GstrandBit+1;
    SAiMarkAbsentBit=GstrandBit+2;

    SAiMarkNmaskC=1LLU << SAiMarkNbit;
    SAiMarkNmask=~SAiMarkNmaskC;
    SAiMarkAbsentMaskC=1LLU << SAiMarkAbsentBit;
    SAiMarkAbsentMask=~SAiMarkAbsentMaskC;

    SAi.defineBits(GstrandBit+3,nSAi);

    P.inOut->logMain << "nGenome=" << nGenome << ";  nSAbyte=" << nSAbyte <<endl<< flush;
    P.inOut->logMain <<"GstrandBit="<<int(GstrandBit)<<"   SA number of indices="<<nSA<<endl<<flush;

    shmSize=SA.lengthByte + nGenome+L+L+SHM_startG+8;
    shmSize+= SAi.lengthByte;

    if ((pGe.gLoad=="LoadAndKeep" ||
         pGe.gLoad=="LoadAndRemove" ||
         pGe.gLoad=="LoadAndExit" ||
         pGe.gLoad=="Remove") && sharedMemory == NULL) {

        bool unloadLast = pGe.gLoad=="LoadAndRemove";
        try {
            sharedMemory = new SharedMemory(shmKey, unloadLast);
            sharedMemory->SetErrorStream(P.inOut->logStdOut);

            if (!sharedMemory->NeedsAllocation())
            P.inOut->logMain <<"Found genome in shared memory\n"<<flush;

            if (pGe.gLoad=="Remove") {//kill the genome and exit
                if (sharedMemory->NeedsAllocation()) {//did not find genome in shared memory, nothing to kill
                ostringstream errOut;
                errOut << "EXITING: Did not find the genome in memory, did not remove any genomes from shared memory\n";
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
                } else {
                            sharedMemory->Clean();
                    P.inOut->logMain <<"DONE: removed the genome from shared memory\n"<<flush;
                            return;
                };
            };

            if (sharedMemory->NeedsAllocation()){
                P.inOut->logMain <<"Allocating shared memory for genome\n"<<flush;
                sharedMemory->Allocate(shmSize);
            };
        } catch (const SharedMemoryException & exc){
            HandleSharedMemoryException(exc, shmSize);
        };

        shmStart = (char*) sharedMemory->GetMapped();
        shmNG= (uint*) (shmStart+SHM_sizeG);
        shmNSA= (uint*) (shmStart+SHM_sizeSA);

        if (!sharedMemory->IsAllocator()) {
            // genome is in shared memory or being loaded
            // wait for the process that will populate it
            // and record the sizes
            uint iwait=0;
            while (*shmNG != nGenome) {
                iwait++;
                P.inOut->logMain <<"Another job is still loading the genome, sleeping for 1 min\n" <<flush;
                sleep(60);
                if (iwait==100) {
                    ostringstream errOut;
                    errOut << "EXITING because of FATAL ERROR: waited too long for the other job to finish loading the genome" << strerror(errno) << "\n" <<flush;
                    errOut << "SOLUTION: remove the shared memory chunk by running STAR with --genomeLoad Remove, and restart STAR" <<flush;
                    exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_LOADING_WAITED_TOO_LONG, P);
                };
            };
            if (nSAbyte!=*shmNSA){
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: the SA file size did not match what we found in shared memory" << "\n" << flush;
                errOut << "SOLUTION: remove the shared memory chunk by running STAR with --genomeLoad Remove, and restart STAR" << flush;
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
            };
            P.inOut->logMain << "Using shared memory for genome. key=0x" <<hex<<shmKey<<dec<< ";   shmid="<< sharedMemory->GetId() <<endl<<flush;
        };
        G1=shmStart+SHM_startG;
        SA.pointArray(G1+nGenome+L+L);
        char* shmNext=SA.charArray+nSAbyte;

        SAi.pointArray(shmNext);
        shmNext += SAi.lengthByte;
    } else if (pGe.gLoad=="NoSharedMemory") {// simply allocate memory, do not use shared memory
        genomeInsertL=0;
        genomeInsertChrIndFirst=nChrReal;
        if (pGe.gFastaFiles.at(0)!="-") {//will insert sequences in the genome, now estimate the extra size
           uint oldlen=chrStart.back();//record the old length
           genomeInsertL=genomeScanFastaFiles(P, G, false, *this)-oldlen;
        };

        try {
            if (P.sjdbInsert.pass1 || P.sjdbInsert.pass2) {
                //reserve extra memory for insertion at the 1st and/or 2nd step
                nGenomeInsert=nGenome+genomeInsertL;
                nSAinsert=nSA+2*genomeInsertL;

                nGenomePass1=nGenomeInsert;
                nSApass1=nSAinsert;
                if (P.sjdbInsert.pass1) {
                    nGenomePass1+=P.limitSjdbInsertNsj*sjdbLength;
                    nSApass1+=2*P.limitSjdbInsertNsj*sjdbLength;
                };

                nGenomePass2=nGenomePass1;
                nSApass2=nSApass1;
                if (P.sjdbInsert.pass2) {
                    nGenomePass2+=P.limitSjdbInsertNsj*sjdbLength;
                    nSApass2+=2*P.limitSjdbInsertNsj*sjdbLength;
                };

                G1=new char[nGenomePass2+L+L];

                SApass2.defineBits(GstrandBit+1,nSApass2);
                SApass2.allocateArray();

                SApass1.defineBits(GstrandBit+1,nSApass1);
                SApass1.pointArray(SApass2.charArray+SApass2.lengthByte-SApass1.lengthByte);

                SAinsert.defineBits(GstrandBit+1,nSAinsert);
                SAinsert.pointArray(SApass1.charArray+SApass1.lengthByte-SAinsert.lengthByte);

                SA.pointArray(SAinsert.charArray+SAinsert.lengthByte-SA.lengthByte);
            } else {//no sjdb insertions
                if (genomeInsertL==0) {// no sequence insertion, simple allocation
                    G1=new char[nGenome+L+L];
                    SA.allocateArray();
                } else {
                    G1=new char[nGenome+L+L+genomeInsertL];
                    SAinsert.defineBits(GstrandBit+1,nSA+2*genomeInsertL);//TODO: re-define GstrandBit if necessary
                    SAinsert.allocateArray();
                    SA.pointArray(SAinsert.charArray+SAinsert.lengthByte-SA.lengthByte);
                };
            };
            SAi.allocateArray();
            P.inOut->logMain <<"Shared memory is not used for genomes. Allocated a private copy of the genome.\n"<<flush;
        } catch (exception & exc) {
            ostringstream errOut;
            errOut <<"EXITING: fatal error trying to allocate genome arrays, exception thrown: "<<exc.what()<<endl;
            errOut <<"Possible cause 1: not enough RAM. Check if you have enough RAM " << nGenome+L+L+SA.lengthByte+SAi.lengthByte+2000000000 << " bytes\n";
            errOut <<"Possible cause 2: not enough virtual memory allowed with ulimit. SOLUTION: run ulimit -v " <<  nGenome+L+L+SA.lengthByte+SAi.lengthByte+2000000000<<endl <<flush;
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_MEMORY_ALLOCATION, P);
        };
    };

    G=G1+L;

    bool isAllocatorProcess = sharedMemory != NULL && sharedMemory->IsAllocator();

    if (pGe.gLoad=="NoSharedMemory" || isAllocatorProcess) {//load genome and SAs from files
        //load genome
        P.inOut->logMain <<"Genome file size: "<<nGenome <<" bytes; state: good=" <<GenomeIn.good()\
                <<" eof="<<GenomeIn.eof()<<" fail="<<GenomeIn.fail()<<" bad="<<GenomeIn.bad()<<"\n"<<flush;
        P.inOut->logMain <<"Loading Genome ... " << flush;
        uint64 genomeReadBytesN=fstreamReadBig(GenomeIn,G,nGenome);
        P.inOut->logMain <<"done! state: good=" <<GenomeIn.good()\
                         <<" eof="<<GenomeIn.eof()<<" fail="<<GenomeIn.fail()<<" bad="<<GenomeIn.bad()<<"; loaded "<<genomeReadBytesN<<" bytes\n" << flush;
        GenomeIn.close();

        for (uint ii=0;ii<L;ii++) {// attach a tail with the largest symbol
            G1[ii]=K-1;
            G[nGenome+ii]=K-1;
        };

        //load SAs
        P.inOut->logMain <<"SA file size: "<<SA.lengthByte <<" bytes; state: good=" <<SAin.good()\
                <<" eof="<<SAin.eof()<<" fail="<<SAin.fail()<<" bad="<<SAin.bad()<<"\n"<<flush;
        P.inOut->logMain <<"Loading SA ... " << flush;
        genomeReadBytesN=fstreamReadBig(SAin,SA.charArray, SA.lengthByte);
        P.inOut->logMain <<"done! state: good=" <<SAin.good()\
                <<" eof="<<SAin.eof()<<" fail="<<SAin.fail()<<" bad="<<SAin.bad()<<"; loaded "<<genomeReadBytesN<<" bytes\n" << flush;
        SAin.close();

        P.inOut->logMain <<"Loading SAindex ... " << flush;
        SAiInBytes +=fstreamReadBig(SAiIn,SAi.charArray, SAi.lengthByte);
        P.inOut->logMain <<"done: "<<SAiInBytes<<" bytes\n" << flush;
    };

    SAiIn.close();

    if ((pGe.gLoad=="LoadAndKeep" ||
         pGe.gLoad=="LoadAndRemove" ||
         pGe.gLoad=="LoadAndExit") && isAllocatorProcess )
    {
        //record sizes. This marks the end of genome loading
        *shmNG=nGenome;
        *shmNSA=nSAbyte;
    };

    time ( &rawtime );
    P.inOut->logMain << "Finished loading the genome: " << asctime (localtime ( &rawtime )) <<"\n"<<flush;

    #ifdef COMPILE_FOR_MAC
    {
        uint sum1=0;
        for (uint ii=0;ii<nGenome; ii++) sum1 +=  (uint) (unsigned char) G[ii];
        P.inOut->logMain << "Sum of all Genome bytes: " <<sum1 <<"\n"<<flush;
        sum1=0;
        for (uint ii=0;ii<SA.lengthByte; ii++) sum1 +=  (uint) (unsigned char) SA.charArray[ii];
        P.inOut->logMain << "Sum of all SA bytes: " <<sum1 <<"\n"<<flush;
        sum1=0;
        for (uint ii=0;ii<SAi.lengthByte; ii++) sum1 +=  (uint) (unsigned char) SAi.charArray[ii];
        P.inOut->logMain << "Sum of all SAi bytes: " <<sum1 <<"\n"<<flush;
    };
    #endif

    if (pGe.gLoad=="LoadAndExit") {
	uint shmSum=0;
	for (uint ii=0;ii<shmSize;ii++) shmSum+=shmStart[ii];
        P.inOut->logMain << "pGe.gLoad=LoadAndExit: completed, the genome is loaded and kept in RAM, EXITING now.\n"<<flush;
        return;
    };

    insertSequences();

    chrBinFill();

    //splice junctions database
    if (nGenome==chrStart[nChrReal]) {//no sjdb
        sjdbN=0;
        sjGstart=chrStart[nChrReal]+1; //not sure why I need that
    } else {//there are sjdb chromosomes
        ifstream sjdbInfo((pGe.gDir+"/sjdbInfo.txt").c_str());
        if (sjdbInfo.fail()) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL error, could not open file " << (pGe.gDir+"/sjdbInfo.txt") <<"\n";
            errOut << "SOLUTION: check that the path to genome files, specified in --genomeDir is correct and the files are present, and have user read permsissions\n" <<flush;
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };


        sjdbInfo >> sjdbN >> pGe.sjdbOverhang;
        P.inOut->logMain << "Processing splice junctions database sjdbN=" <<sjdbN<<",   pGe.sjdbOverhang=" <<pGe.sjdbOverhang <<" \n";

        sjChrStart=nChrReal;
        sjGstart=chrStart[sjChrStart];

        //fill the sj-db to genome translation array
        sjDstart=new uint [sjdbN];
        sjAstart=new uint [sjdbN];
        sjdbStart=new uint [sjdbN];
        sjdbEnd=new uint [sjdbN];

        sjdbMotif=new uint8 [sjdbN];
        sjdbShiftLeft=new uint8 [sjdbN];
        sjdbShiftRight=new uint8 [sjdbN];
        sjdbStrand=new uint8 [sjdbN];

        for (uint ii=0;ii<sjdbN;ii++) {//get the info about junctions from sjdbInfo.txt
            {
                uint16 d1,d2,d3,d4;
                sjdbInfo >> sjdbStart[ii] >> sjdbEnd[ii] >> d1 >> d2 >> d3 >> d4;
                sjdbMotif[ii]      = (uint8) d1;
                sjdbShiftLeft[ii]  = (uint8) d2;
                sjdbShiftRight[ii] = (uint8) d3;
                sjdbStrand[ii] = (uint8) d4;
            };
            sjDstart[ii]   = sjdbStart[ii]  - pGe.sjdbOverhang;
            sjAstart[ii]   = sjdbEnd[ii] + 1;
            if (sjdbMotif[ii]==0) {//shinon-canonical junctions back to their true coordinates
                sjDstart[ii] += sjdbShiftLeft[ii];
                sjAstart[ii] += sjdbShiftLeft[ii];
            };
        };
    };

    //check and redefine some parameters
    //max intron size
    if (P.alignIntronMax==0 && P.alignMatesGapMax==0) {
        P.inOut->logMain << "alignIntronMax=alignMatesGapMax=0, the max intron size will be approximately determined by (2^winBinNbits)*winAnchorDistNbins=" \
                << (1LLU<<P.winBinNbits)*P.winAnchorDistNbins <<endl;
    } else {
        //redefine winBinNbits
        P.winBinNbits = (uint) floor( log2( max( max(4LLU,P.alignIntronMax), (P.alignMatesGapMax==0 ? 1000LLU : P.alignMatesGapMax) ) /4 ) + 0.5);
        P.winBinNbits = max( P.winBinNbits, (uint) floor(log2(nGenome/40000+1)+0.5) );
        //ISSUE - to be fixed in STAR3: if alignIntronMax>0 but alignMatesGapMax==0, winBinNbits will be defined by alignIntronMax
        P.inOut->logMain << "To accommodate alignIntronMax="<<P.alignIntronMax<<" redefined winBinNbits="<< P.winBinNbits <<endl;
    };

    if (P.winBinNbits > pGe.gChrBinNbits) {
       P.inOut->logMain << "winBinNbits=" <<P.winBinNbits <<" > " << "pGe.gChrBinNbits=" << pGe.gChrBinNbits << "   redefining:\n";
       P.winBinNbits=pGe.gChrBinNbits;
       P.inOut->logMain << "winBinNbits=" <<P.winBinNbits <<endl;
    };


    if (P.alignIntronMax==0 && P.alignMatesGapMax==0) {
    } else {
        //redefine winFlankNbins,winAnchorDistNbins
        P.winFlankNbins=max(P.alignIntronMax,P.alignMatesGapMax)/(1LLU<<P.winBinNbits)+1;
        P.winAnchorDistNbins=2*P.winFlankNbins;
        P.inOut->logMain << "To accommodate alignIntronMax="<<P.alignIntronMax<<" and alignMatesGapMax="<<P.alignMatesGapMax<<\
                ", redefined winFlankNbins="<<P.winFlankNbins<<" and winAnchorDistNbins="<<P.winAnchorDistNbins<<endl;
    };

    P.winBinChrNbits=pGe.gChrBinNbits-P.winBinNbits;
    P.winBinN = nGenome/(1LLU << P.winBinNbits)+1;//this may be changed later
    
    if (pGe.gType==101) {//SuperTranscriptome
		superTr = new SuperTranscriptome (P);
        superTr->load(G, chrStart, chrLength);
        
        //genomeOut
        P.pGeOut.gDir=pGe.gDir+"/fullGenome/";
        genomeOut.convYes=true;
        genomeOut.gapsAreJunctions=true;
        genomeOut.convFile=pGe.gDir+"/fullGenome/conversionToFullGenome.tsv";
        
        genomeOut.g = new Genome(P,P.pGeOut);
        genomeOut.g->genomeOut=genomeOut;
        genomeOut.g->genomeOutLoad();
        
        //debug checks
//         for (auto &b : genomeOut.g->genomeOut.convBlocks) {
//             uint64 g1=*(genomeOut.g->G+b[2]);
//             uint64 g2=*(G+b[0]);
//             if (g1!=g2)
//                 cout <<b[0]<<" "<<b[1]<<" "<<b[2]<<endl;
//         };
    };
    
        
    if (P1.pGe.transform.typeString=="Haploid") {
        pGe.transform.type=1;
    } else if (P1.pGe.transform.typeString=="Diploid") {
        pGe.transform.type=2;
    } else {//TODO check for wrong values
        pGe.transform.type=0;
    };    
    
    if (pGe.transform.type != 0) {//convert genome coordinates
        //genomeOut
        P.pGeOut.gDir=pGe.gDir+"/normalGenome/";
        genomeOut.convYes=true;
        genomeOut.gapsAreJunctions=false;
        genomeOut.convFile=pGe.gDir+"/transformGenomeBlocks.tsv";

        
        genomeOut.g = new Genome(P,P.pGeOut);
        genomeOut.g->genomeOut=genomeOut;
        genomeOut.g->genomeOutLoad();
    };
    
    if (P.pGe.gLoad=="LoadAndExit" || P.pGe.gLoad=="Remove") {
        exit(0);
    };
};
