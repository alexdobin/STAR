#include "IncludeDefine.h"
#include "Parameters.h"
#include "SequenceFuns.h"
#include "Genome.h"
#include "ReadAlignChunk.h"
#include "ReadAlign.h"
#include "Stats.h"
#include "genomeGenerate.h"
#include "outputSJ.h"
#include "ThreadControl.h"
#include "GlobalVariables.cpp"
#include "TimeFunctions.h"
#include <sys/types.h>
#include <sys/stat.h>
#include "ErrorWarning.h"
#include "sysRemoveDir.h"
#include "BAMfunctions.h"
#include "Transcriptome.h"
#include "BAMbinSortByCoordinate.h"
#include "signalFromBAM.h"
#include "sjdbBuildIndex.h"
#include "mapThreadsSpawn.h"
#include "ErrorWarning.h"

#include "htslib/htslib/sam.h"
extern int bam_cat(int nfn, char * const *fn, const bam_hdr_t *h, const char* outbam);

int main(int argInN, char* argIn[]) {
   
    time(&g_statsAll.timeStart);
   
    Parameters *P = new Parameters; //all parameters
       
    P->inputParameters(argInN, argIn);
    
    *(P->inOut->logStdOut) << timeMonthDayTime(g_statsAll.timeStart) << " ..... Started STAR run\n" <<flush;           
    
    //generate genome
    if (P->runMode=="genomeGenerate") {
        genomeGenerate(P);
        (void) sysRemoveDir (P->outFileTmp);        
        P->inOut->logMain << "DONE: Genome generation, EXITING\n" << flush;
        exit(0);
    } else if (P->runMode!="alignReads") {
        P->inOut->logMain << "EXITING because of INPUT ERROR: unknown value of input parameter runMode=" <<P->runMode<<endl<<flush;
        exit(1);
    };
    
    Genome mainGenome (P);
    mainGenome.genomeLoad();
    //calculate genome-related parameters
    P->winBinN = P->nGenome/(1LLU << P->winBinNbits)+1;
    
    Transcriptome *mainTranscriptome=NULL;
    if (P->quantModeI>0) {//load transcriptome
        mainTranscriptome=new Transcriptome(P);
    };
/////////////////////////////////////////////////////////////////////////////////////////////////START
    if (P->runThreadN>1) {
        g_threadChunks.threadArray=new pthread_t[P->runThreadN];
        pthread_mutex_init(&g_threadChunks.mutexInRead, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutSAM, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutBAM1, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutUnmappedFastx, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutFilterBySJout, NULL);
        pthread_mutex_init(&g_threadChunks.mutexStats, NULL);
    };

    g_statsAll.progressReportHeader(P->inOut->logProgress);    
    
    if (P->twopass1readsN>0) {//2-pass
        //re-define P for the pass1
        
        Parameters *P1=new Parameters;
        *P1=*P;
        //turn off unnecessary calculations
        P1->outSAMtype[0]="None";
        P1->outSAMbool=false;
        P1->outBAMunsorted=false;
        P1->outBAMcoord=false;
    
        P1->chimSegmentMin=0;
        P1->quantModeI=0;
        P1->outFilterBySJoutStage=0;
        
        P1->outReadsUnmapped="None";
        
        P1->outFileNamePrefix=P->twopassDir;

        P1->readMapNumber=P->twopass1readsN;
//         P1->inOut->logMain.open((P1->outFileNamePrefix + "Log.out").c_str());

        g_statsAll.resetN();
        time(&g_statsAll.timeStartMap);
        P->inOut->logProgress << timeMonthDayTime(g_statsAll.timeStartMap) <<"\tStarted 1st pass mapping\n" <<flush;
        *P->inOut->logStdOut << timeMonthDayTime(g_statsAll.timeStartMap) << " ..... Started 1st pass mapping\n" <<flush;

        //run mapping for Pass1
        ReadAlignChunk *RAchunk1[P->runThreadN];        
        for (int ii=0;ii<P1->runThreadN;ii++) {
            RAchunk1[ii]=new ReadAlignChunk(P1, mainGenome, mainTranscriptome, ii);
        };    
        mapThreadsSpawn(P1, RAchunk1);
        outputSJ(RAchunk1,P1); //collapse and output junctions
//         for (int ii=0;ii<P1->runThreadN;ii++) {
//             delete [] RAchunk[ii];
//         };          
        
        time_t rawtime; time (&rawtime);
        P->inOut->logProgress << timeMonthDayTime(rawtime) <<"\tFinished 1st pass mapping\n";
        *P->inOut->logStdOut << timeMonthDayTime(rawtime) << " ..... Finished 1st pass mapping\n" <<flush;
        ofstream logFinal1 ( (P->twopassDir + "/Log.final.out").c_str());
        g_statsAll.reportFinal(logFinal1,P1);

        //re-build genome files
        
        sjdbBuildIndex (P, mainGenome.G, mainGenome.SA, mainGenome.SA2, mainGenome.SAi);
        time ( &rawtime ); 
        *P->inOut->logStdOut  << timeMonthDayTime(rawtime) << " ..... Finished inserting 1st pass junctions into genome" <<endl;
        //re-calculate genome-related parameters
        P->winBinN = P->nGenome/(1LLU << P->winBinNbits)+1;
        
        //reopen reads files
        P->closeReadsFiles();
        P->openReadsFiles();
    };

    //initialize Stats
    g_statsAll.resetN();
    g_statsAll.progressReportHeader(P->inOut->logProgress);    
    time(&g_statsAll.timeStartMap);
    *P->inOut->logStdOut << timeMonthDayTime(g_statsAll.timeStartMap) << " ..... Started mapping\n" <<flush;
    
    g_statsAll.timeLastReport=g_statsAll.timeStartMap;

    //open SAM/BAM files for output
    if (P->outSAMmode != "None") {//open SAM file and write header
        ostringstream samHeaderStream;

        if (P->outSAMheaderHD.at(0)!="-") {
            samHeaderStream << P->outSAMheaderHD.at(0);
            for (uint ii=1;ii<P->outSAMheaderHD.size(); ii++) {
                samHeaderStream << "\t" << P->outSAMheaderHD.at(ii);
            };
            samHeaderStream << "\n";
        } else {
            samHeaderStream << "@HD\tVN:1.4\n";
        };
        
        for (uint ii=0;ii<P->nChrReal;ii++) {
            samHeaderStream << "@SQ\tSN:"<< P->chrName.at(ii) <<"\tLN:"<<P->chrLength[ii]<<"\n";
        };

        if (P->outSAMheaderPG.at(0)!="-") {
            samHeaderStream << P->outSAMheaderPG.at(0);
            for (uint ii=1;ii<P->outSAMheaderPG.size(); ii++) {
                samHeaderStream << "\t" << P->outSAMheaderPG.at(ii);
            };
            samHeaderStream << "\n";
        };        
        
        samHeaderStream << "@PG\tID:STAR\tPN:STAR\tVN:" << STAR_VERSION <<"\tCL:" << P->commandLineFull <<"\n";
        
        if (P->outSAMheaderCommentFile!="-") {
            ifstream comstream (P->outSAMheaderCommentFile);
            while (comstream.good()) {
                string line1;
                getline(comstream,line1);
                if (line1.find_first_not_of(" \t\n\v\f\r")!=std::string::npos) {//skip blank lines
                    samHeaderStream << line1 <<"\n";
                };
            };
        };         
        

        for (uint32 ii=0;ii<P->outSAMattrRGlineSplit.size();ii++) {//@RG lines
            samHeaderStream << "@RG\t" << P->outSAMattrRGlineSplit.at(ii) <<"\n";
        };
 
        samHeaderStream <<  "@CO\t" <<"user command line: " << P->commandLine <<"\n";
        
        P->samHeader=samHeaderStream.str();
        
        if (P->outSAMbool) {//
            *P->inOut->outSAM << P->samHeader;
        };
        if (P->outBAMunsorted){
            outBAMwriteHeader(P->inOut->outBAMfileUnsorted,P->samHeader,P->chrName,P->chrLength);
        };
//             if (P->outBAMcoord){
//                 outBAMwriteHeader(P->inOut->outBAMfileCoord,P->samHeader,P->chrName,P->chrLength);            
//             };
        
        if ( (P->quantModeI & PAR_quantModeI_TranscritomeSAM) > 0) {
            samHeaderStream.str("");
            vector <uint> trlength;
            for (uint32 ii=0;ii<mainTranscriptome->trID.size();ii++) {
                uint32 iex1=mainTranscriptome->trExI[ii]+mainTranscriptome->trExN[ii]-1; //last exon of the transcript
                trlength.push_back(mainTranscriptome->exLenCum[iex1]+mainTranscriptome->exSE[2*iex1+1]-mainTranscriptome->exSE[2*iex1]+1);          
                samHeaderStream << "@SQ\tSN:"<< mainTranscriptome->trID.at(ii) <<"\tLN:"<<trlength.back()<<"\n";
            };            
            outBAMwriteHeader(P->inOut->outQuantBAMfile,samHeaderStream.str(),mainTranscriptome->trID,trlength);        
        };
        
    };
    
    if (P->chimSegmentMin>0) {
        P->inOut->outChimJunction.open((P->outFileNamePrefix + "Chimeric.out.junction").c_str());
        P->inOut->outChimSAM.open((P->outFileNamePrefix + "Chimeric.out.sam").c_str());
        P->inOut->outChimSAM << P->samHeader;
        pthread_mutex_init(&g_threadChunks.mutexOutChimSAM, NULL);   
        pthread_mutex_init(&g_threadChunks.mutexOutChimJunction, NULL);
    };
         
    // P->inOut->logMain << "mlock value="<<mlockall(MCL_CURRENT|MCL_FUTURE) <<"\n"<<flush;

    // prepare chunks and spawn mapping threads    
    ReadAlignChunk *RAchunk[P->runThreadN];
    for (int ii=0;ii<P->runThreadN;ii++) {
        RAchunk[ii]=new ReadAlignChunk(P, mainGenome, mainTranscriptome, ii);
    };    
    
    mapThreadsSpawn(P, RAchunk);
   
    if (P->outFilterBySJoutStage==1) {//completed stage 1, go to stage 2
        P->inOut->logMain << "Completed stage 1 mapping of outFilterBySJout mapping\n"<<flush;
        outputSJ(RAchunk,P);//collapse novel junctions
        P->readFilesIndex=-1;
        
        P->outFilterBySJoutStage=2;
        mapThreadsSpawn(P, RAchunk);
    };
    
    //close some BAM files
    if (P->inOut->outBAMfileUnsorted!=NULL) {
        bgzf_flush(P->inOut->outBAMfileUnsorted);
        bgzf_close(P->inOut->outBAMfileUnsorted);
    };
    if (P->inOut->outQuantBAMfile!=NULL) {
        bgzf_flush(P->inOut->outQuantBAMfile);
        bgzf_close(P->inOut->outQuantBAMfile);
    };      
    
    if (P->outBAMcoord && P->limitBAMsortRAM==0) {//make it equal ot the genome size
        P->limitBAMsortRAM=P->nGenome+mainGenome.SA.lengthByte+mainGenome.SAi.lengthByte;
    };
        
    //no need for genome anymore, free the memory
//     mainGenome.~Genome(); //need explicit call because of the delete P->inOut below
    mainGenome.freeMemory();
    
    if (P->runThreadN>1 && P->outSAMorder=="PairedKeepInputOrder") {//concatenate Aligned.* files
        RAchunk[0]->chunkFilesCat(P->inOut->outSAM, P->outFileTmp + "/Aligned.out.sam.chunk", g_threadChunks.chunkOutN);
    };    
    
    
    
    if (P->outBAMcoord) {//sort BAM if needed
        *P->inOut->logStdOut << timeMonthDayTime() << " ..... Started sorting BAM\n" <<flush;
        P->inOut->logMain << timeMonthDayTime() << " ..... Started sorting BAM\n" <<flush;
        
        //check max size needed for sorting
        uint maxMem=0;
        for (uint32 ibin=0; ibin<RAchunk[0]->chunkOutBAMcoord->nBins; ibin++) {
            uint binS=0;
            for (int it=0; it<P->runThreadN; it++) {//collect sizes from threads
                binS += RAchunk[it]->chunkOutBAMcoord->binTotalBytes[ibin];
            };        
            if (binS>maxMem) maxMem=binS;
        };
        P->inOut->logMain << "Max memory needed for sorting = "<<maxMem<<endl;
        if (maxMem>P->limitBAMsortRAM) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal ERROR: not enough memory for BAM sorting: \n";
            errOut <<"SOLUTION: re-run STAR with at least --limitBAMsortRAM " <<maxMem+1000000000;
            exitWithError(errOut.str(), std::cerr, P->inOut->logMain, EXIT_CODE_PARAMETER, *P);                                    
        };
        
        
        uint totalMem=0;
//         P->inOut->logMain << "Started sorting BAM ..." <<endl;
        #pragma omp parallel num_threads(P->runThreadN) 
        #pragma omp for schedule (dynamic,1)
        for (uint32 ibin=0; ibin<RAchunk[0]->chunkOutBAMcoord->nBins; ibin++) {
            
            uint binN=0, binS=0;
            for (int it=0; it<P->runThreadN; it++) {//collect sizes from threads
                binN += RAchunk[it]->chunkOutBAMcoord->binTotalN[ibin];
                binS += RAchunk[it]->chunkOutBAMcoord->binTotalBytes[ibin];
            };
            
            uint newMem=binS+binN*24;
            bool boolWait=true;
            while (boolWait) {
                #pragma omp critical
                if (totalMem+newMem < P->limitBAMsortRAM) {
                    boolWait=false;
                    totalMem+=newMem;
                };
                sleep(0.1);
            };
            BAMbinSortByCoordinate(ibin,binN,binS,P->runThreadN,P->outBAMsortTmpDir,P->inOut->outBAMfileCoord, P);
            #pragma omp critical
            totalMem-=newMem;//"release" RAM
        };
        //concatenate all BAM files, using bam_cat
        char **bamBinNames = new char* [RAchunk[0]->chunkOutBAMcoord->nBins];
        vector <string> bamBinNamesV;
        for (uint32 ibin=0; ibin<RAchunk[0]->chunkOutBAMcoord->nBins; ibin++) {
            
            bamBinNamesV.push_back(P->outBAMsortTmpDir+"/b"+to_string((uint) ibin));            
            struct stat buffer;
            if (stat (bamBinNamesV.back().c_str(), &buffer) != 0) {//check if file exists
                bamBinNamesV.pop_back();
            };
        };
        for (uint32 ibin=0; ibin<bamBinNamesV.size(); ibin++) {
                bamBinNames[ibin] = (char*) bamBinNamesV.at(ibin).c_str();
        };
        bam_cat(bamBinNamesV.size(), bamBinNames, 0, P->outBAMfileCoordName.c_str());
    };
    //wiggle output
    if (P->outWigFlags.yes) {
        *P->inOut->logStdOut << timeMonthDayTime() << " ..... Started wiggle output\n" <<flush;
        P->inOut->logMain << timeMonthDayTime() << " ..... Started wiggle output\n" <<flush;
        string wigOutFileNamePrefix=P->outFileNamePrefix + "Signal";
        signalFromBAM(P->outBAMfileCoordName, wigOutFileNamePrefix, *P);
    };
    
    //aggregate output (junctions, signal, etc)
    //collapse splice junctions from different threads/chunks, and output them
    outputSJ(RAchunk,P);
    
    g_statsAll.progressReport(P->inOut->logProgress);
    P->inOut->logProgress  << "ALL DONE!\n"<<flush;
    P->inOut->logFinal.open((P->outFileNamePrefix + "Log.final.out").c_str());
    g_statsAll.reportFinal(P->inOut->logFinal,P);
    *P->inOut->logStdOut << timeMonthDayTime(g_statsAll.timeFinish) << " ..... Finished successfully\n" <<flush;
    
    P->inOut->logMain  << "ALL DONE!\n"<<flush;
    sysRemoveDir (P->outFileTmp);
    
    P->closeReadsFiles();//this will kill the readFilesCommand processes if necessary
    delete P->inOut; //to close files
    delete P;
    
    return 0;    
};
