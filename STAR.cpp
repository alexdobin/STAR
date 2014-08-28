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
    
    
    {//2-pass
        time_t rawtime;
        time ( &rawtime );
        cout << timeMonthDayTime(rawtime) << "start sjdbBuild" <<endl;
        sjdbBuildIndex (P, mainGenome.G, mainGenome.SA);
        cout << timeMonthDayTime(rawtime) << "finished" <<endl;
    };
    //calculate genome-related parameters
    P->winBinN = P->nGenome/(1LLU << P->winBinNbits)+1;
    

    Transcriptome *mainTranscriptome=NULL;
    if (P->quantModeI>0) {//load transcriptome
        mainTranscriptome=new Transcriptome(P);
    };
/////////////////////////////////////////////////////////////////////////////////////////////////START
    
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
        
        samHeaderStream << "@PG\tID:STAR\tPN:STAR\tVN:" << SVN_VERSION_COMPILED <<"\tCL:" << P->commandLineFull <<"\n";
        
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
    
   
    ReadAlignChunk *RAchunk[P->runThreadN];
    for (int ii=0;ii<P->runThreadN;ii++) {
        RAchunk[ii]=new ReadAlignChunk(P, mainGenome, mainTranscriptome, ii);
        RAchunk[ii]->RA->iRead=0;
        RAchunk[ii]->iThread=ii;
    };
    
    if (P->runThreadN>1) {
        g_threadChunks.threadArray=new pthread_t[P->runThreadN];
        pthread_mutex_init(&g_threadChunks.mutexInRead, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutSAM, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutBAM1, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutUnmappedFastx, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutFilterBySJout, NULL);
        pthread_mutex_init(&g_threadChunks.mutexStats, NULL);
    };
    
    
    ///////////////////////////////////////////////////////////////////
    g_statsAll.progressReportHeader(P->inOut->logProgress);    
    time(&g_statsAll.timeStartMap);
    *P->inOut->logStdOut << timeMonthDayTime(g_statsAll.timeStartMap) << " ..... Started mapping\n" <<flush;
    
    g_statsAll.timeLastReport=g_statsAll.timeStartMap;
    
    for (int ithread=1;ithread<P->runThreadN;ithread++) {//spawn threads
        int threadStatus=pthread_create(&g_threadChunks.threadArray[ithread], NULL, &g_threadChunks.threadRAprocessChunks, (void *) RAchunk[ithread]);
        if (threadStatus>0) {//something went wrong with one of threads
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: phtread error while creating thread # " << ithread <<", error code: "<<threadStatus ;
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, 1, *P);
        };
        P->inOut->logMain << "Created thread # " <<ithread <<endl;
    };
    
    RAchunk[0]->processChunks(); //start main thread
    
    for (int ithread=1;ithread<P->runThreadN;ithread++) {//wait for all threads to complete
        int threadStatus = pthread_join(g_threadChunks.threadArray[ithread], NULL);
        if (threadStatus>0) {//something went wrong with one of threads
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: phtread error while joining thread # " << ithread <<", error code: "<<threadStatus ;
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, 1, *P);
        };
        P->inOut->logMain << "Joined thread # " <<ithread <<"\n";        
    };    
    
    if (P->outFilterBySJoutStage==1) {//completed stage 1, go to stage 2
        outputSJ(RAchunk,P);//collapse novel junctions
        P->readFilesIndex=-1;
        
        P->outFilterBySJoutStage=2;
        P->inOut->logMain << "starting stage 2: filtering BySJout\n";
        for (int ithread=1;ithread<P->runThreadN;ithread++) {//spawn threads
            int threadStatus = pthread_create(&g_threadChunks.threadArray[ithread], NULL, &g_threadChunks.threadRAprocessChunks, (void *) RAchunk[ithread]);
            if (threadStatus>0) {//something went wrong with one of threads
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: phtread error while creating thread # " << ithread <<", error code: "<<threadStatus ;
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, 1, *P);
            };
            P->inOut->logMain << "Created thread # " <<ithread <<"\n";
        };

        RAchunk[0]->processChunks(); //start main thread

        for (int ithread=1;ithread<P->runThreadN;ithread++) {//wait for all threads to complete
            int threadStatus = pthread_join(g_threadChunks.threadArray[ithread], NULL);
            if (threadStatus) {//something went wrong with one of threads
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: phtread error while joining thread # " << ithread <<", error code: "<<threadStatus ;
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, 1, *P);
            };
            P->inOut->logMain << "Joined thread # " <<ithread <<"\n";        
        };  
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
    mainGenome.~Genome(); //need explicit call because of the delete P->inOut below
    
    if (P->runThreadN>1 && P->outSAMorder=="PairedKeepInputOrder") {//concatenate Aligned.* files
        RAchunk[0]->chunkFilesCat(P->inOut->outSAM, P->outFileTmp + "/Aligned.out.sam.chunk", g_threadChunks.chunkOutN);
    };    
    
    
    
    if (P->outBAMcoord) {//sort BAM if needed
        *P->inOut->logStdOut << timeMonthDayTime() << " ..... Started sorting BAM\n" <<flush;
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
        signalFromBAM(P->outBAMfileCoordName, wigOutFileNamePrefix, P->outWigFlags.strand, P->outWigFlags.type, P->outWigReferencesPrefix);
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
    
    delete P->inOut; //to close files
    delete P;
    
    return 0;    
};
