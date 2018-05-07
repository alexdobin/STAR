#include <sys/types.h>
#include <sys/stat.h>

#include "IncludeDefine.h"
#include "Parameters.h"
#include "SequenceFuns.h"
#include "Genome.h"
#include "Chain.h"
#include "ReadAlignChunk.h"
#include "ReadAlign.h"
#include "Stats.h"
#include "genomeGenerate.h"
#include "outputSJ.h"
#include "ThreadControl.h"
#include "GlobalVariables.h"
#include "TimeFunctions.h"
#include "ErrorWarning.h"
#include "sysRemoveDir.h"
#include "BAMfunctions.h"
#include "Transcriptome.h"
#include "BAMbinSortByCoordinate.h"
#include "BAMbinSortUnmapped.h"
#include "signalFromBAM.h"
#include "mapThreadsSpawn.h"
#include "ErrorWarning.h"
#include "SjdbClass.h"
#include "sjdbInsertJunctions.h"
#include "Variation.h"
#include "bam_cat.h"

#include "htslib/htslib/sam.h"
#include "parametersDefault.xxd"



void usage() {
    printf("Usage: STAR  [options]... --genomeDir REFERENCE   --readFilesIn R1.fq R2.fq\n");
    printf("Spliced Transcripts Alignment to a Reference (c) Alexander Dobin, 2009-2015\n\n");

    std::cout << parametersDefault;

    printf("For more details see:\n");
    printf("<https://github.com/alexdobin/STAR>\n");
    printf("<https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>\n");

    exit(0);
}


int main(int argInN, char* argIn[]) {
    
    //debug
//     uint nn=4000000000;
//     uint* dummy=new uint[nn];
//     std::cout <<"Allocated"<<endl;
//     for (uint ii=0;ii<nn;ii++) {
//         dummy[ii]=14829735431805718528LLU;;
//     };
//     std::cout <<"Filled "<<dummy[0]<<endl;
//     delete[] dummy;     
//     std::cout <<"Deleted"<<endl;
    
    // If no argument is given, or the first argument is either '-h' or '--help', run usage()
    if((argInN == 1) || (argInN == 2 && (strcmp("-h",argIn[1]) == 0 || strcmp ("--help",argIn[1]) == 0 ))) {
        usage();
    }

    time(&g_statsAll.timeStart);

    Parameters P; //all parameters
    P.inputParameters(argInN, argIn);

    *(P.inOut->logStdOut) << timeMonthDayTime(g_statsAll.timeStart) << " ..... started STAR run\n" <<flush;

    //generate genome
    if (P.runMode=="alignReads") 
    {//continue
    } else if (P.runMode=="genomeGenerate") 
    {
        Genome mainGenome(P);
        mainGenome.genomeGenerate();
        (void) sysRemoveDir (P.outFileTmp);
        P.inOut->logMain << "DONE: Genome generation, EXITING\n" << flush;
        exit(0);
    } else if (P.runMode=="liftOver") 
    {        
        for (uint ii=0; ii<P.pGe.gChainFiles.size();ii++)
        {
            Chain chain(P,P.pGe.gChainFiles.at(ii));
            chain.liftOverGTF(P.pGe.sjdbGTFfile,P.outFileNamePrefix+"GTFliftOver_"+to_string(ii+1)+".gtf");
            P.inOut->logMain << "DONE: lift-over of GTF file, EXITING\n" << flush;
            exit(0);
        };
    } else {
        P.inOut->logMain << "EXITING because of INPUT ERROR: unknown value of input parameter runMode=" <<P.runMode<<endl<<flush;
        exit(1);        
    };

    Genome mainGenome (P);
    mainGenome.genomeLoad();
    
    
    if (P.pGe.gLoad=="LoadAndExit" || P.pGe.gLoad=="Remove")
    {
        return 0;
    };
    
    P.twoPass.pass2=false; //this is the 1st pass    
    

    //calculate genome-related parameters
    Transcriptome *mainTranscriptome=NULL;
    mainGenome.Var=new Variation(P, mainGenome.chrStart, mainGenome.chrNameIndex);


    if (P.pGe.gFastaFiles.at(0)!="-")
    {//insert sequences in the genome

    };

    SjdbClass sjdbLoci;

    if (P.sjdbInsert.pass1)
    {
        Genome mainGenome1=mainGenome;//not sure if I need to create the copy - mainGenome1 below should not be changed
        sjdbInsertJunctions(P, mainGenome, mainGenome1, sjdbLoci);
    };

/////////////////////////////////////////////////////////////////////////////////////////////////START
    if (P.runThreadN>1) {
        g_threadChunks.threadArray=new pthread_t[P.runThreadN];
        pthread_mutex_init(&g_threadChunks.mutexInRead, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutSAM, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutBAM1, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutUnmappedFastx, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutFilterBySJout, NULL);
        pthread_mutex_init(&g_threadChunks.mutexStats, NULL);
        pthread_mutex_init(&g_threadChunks.mutexBAMsortBins, NULL);
    };

    g_statsAll.progressReportHeader(P.inOut->logProgress);

    if (P.twoPass.yes) {//2-pass
        //re-define P for the pass1

        Genome mainGenome1=mainGenome;

        Parameters P1=P;
        //turn off unnecessary calculations
        P1.outSAMtype[0]="None";
        P1.outSAMbool=false;
        P1.outBAMunsorted=false;
        P1.outBAMcoord=false;

        P1.pCh.segmentMin=0;

        P1.quant.yes=false;
        P1.quant.trSAM.yes=false;
        P1.quant.geCount.yes=false;

        P1.outFilterBySJoutStage=0;

        P1.outReadsUnmapped="None";

        P1.outFileNamePrefix=P.twoPass.dir;

        P1.readMapNumber=P.twoPass.pass1readsN;
//         P1.inOut->logMain.open((P1.outFileNamePrefix + "Log.out").c_str());

        P1.wasp.outputMode="None"; //no WASP filtering on the 1st pass
        
        g_statsAll.resetN();
        time(&g_statsAll.timeStartMap);
        P.inOut->logProgress << timeMonthDayTime(g_statsAll.timeStartMap) <<"\tStarted 1st pass mapping\n" <<flush;
        *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeStartMap) << " ..... started 1st pass mapping\n" <<flush;

        //run mapping for Pass1
        ReadAlignChunk *RAchunk1[P.runThreadN];
        for (int ii=0;ii<P1.runThreadN;ii++) {
            RAchunk1[ii]=new ReadAlignChunk(P1, mainGenome, mainTranscriptome, ii);
        };
        mapThreadsSpawn(P1, RAchunk1);
        outputSJ(RAchunk1,P1); //collapse and output junctions
//         for (int ii=0;ii<P1.runThreadN;ii++) {
//             delete [] RAchunk[ii];
//         };

        time_t rawtime; time (&rawtime);
        P.inOut->logProgress << timeMonthDayTime(rawtime) <<"\tFinished 1st pass mapping\n";
        *P.inOut->logStdOut << timeMonthDayTime(rawtime) << " ..... finished 1st pass mapping\n" <<flush;
        ofstream logFinal1 ( (P.twoPass.dir + "/Log.final.out").c_str());
        g_statsAll.reportFinal(logFinal1);

        P.twoPass.pass2=true;//starting the 2nd pass
        P.twoPass.pass1sjFile=P.twoPass.dir+"/SJ.out.tab";

        sjdbInsertJunctions(P, mainGenome, mainGenome1, sjdbLoci);

        //reopen reads files
        P.closeReadsFiles();
        P.openReadsFiles();
    } else {//not 2-pass
        //nothing for now
    };

    if ( P.quant.yes ) {//load transcriptome
        mainTranscriptome=new Transcriptome(P);
    };

    //initialize Stats
    g_statsAll.resetN();
    time(&g_statsAll.timeStartMap);
    *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeStartMap) << " ..... started mapping\n" <<flush;

    g_statsAll.timeLastReport=g_statsAll.timeStartMap;

    //open SAM/BAM files for output
    if (P.outSAMmode != "None") {//open SAM file and write header
        ostringstream samHeaderStream;

        for (uint ii=0;ii<mainGenome.nChrReal;ii++) {
            samHeaderStream << "@SQ\tSN:"<< mainGenome.chrName.at(ii) <<"\tLN:"<<mainGenome.chrLength[ii]<<"\n";
        };

        mainGenome.chrNameAll=mainGenome.chrName;
        mainGenome.chrLengthAll=mainGenome.chrLength;
        {//add exra references
            ifstream extrastream (P.pGe.gDir + "/extraReferences.txt");
            while (extrastream.good()) {
                string line1;
                getline(extrastream,line1);
                istringstream stream1 (line1);
                string field1;
                stream1 >> field1;//should check for @SQ

                if (field1!="") {//skip blank lines
                    samHeaderStream << line1 <<"\n";

                    stream1 >> field1;
                    mainGenome.chrNameAll.push_back(field1.substr(3));
                    stream1 >> field1;
                    mainGenome.chrLengthAll.push_back((uint) stoll(field1.substr(3)));
                };
            };
            extrastream.close();
        };
        
        if (P.outSAMheaderPG.at(0)!="-") {
            samHeaderStream << P.outSAMheaderPG.at(0);
            for (uint ii=1;ii<P.outSAMheaderPG.size(); ii++) {
                samHeaderStream << "\t" << P.outSAMheaderPG.at(ii);
            };
            samHeaderStream << "\n";
        };

        samHeaderStream << "@PG\tID:STAR\tPN:STAR\tVN:" << STAR_VERSION <<"\tCL:" << P.commandLineFull <<"\n";

        if (P.outSAMheaderCommentFile!="-") {
            ifstream comstream (P.outSAMheaderCommentFile);
            while (comstream.good()) {
                string line1;
                getline(comstream,line1);
                if (line1.find_first_not_of(" \t\n\v\f\r")!=std::string::npos) {//skip blank lines
                    samHeaderStream << line1 <<"\n";
                };
            };
            comstream.close();
        };


        for (uint32 ii=0;ii<P.outSAMattrRGlineSplit.size();ii++) {//@RG lines
            samHeaderStream << "@RG\t" << P.outSAMattrRGlineSplit.at(ii) <<"\n";
        };


        samHeaderStream <<  "@CO\t" <<"user command line: " << P.commandLine <<"\n";
        
        samHeaderStream << P.samHeaderExtra;

        if (P.outSAMheaderHD.at(0)!="-") {
            P.samHeaderHD = P.outSAMheaderHD.at(0);
            for (uint ii=1;ii<P.outSAMheaderHD.size(); ii++) {
                P.samHeaderHD +="\t" + P.outSAMheaderHD.at(ii);
            };
        } else {
            P.samHeaderHD = "@HD\tVN:1.4";
        };


        P.samHeader=P.samHeaderHD+"\n"+samHeaderStream.str();
        //for the sorted BAM, need to add SO:cooridnate to the header line
        P.samHeaderSortedCoord=P.samHeaderHD + (P.outSAMheaderHD.size()==0 ? "" : "\tSO:coordinate") + "\n" + samHeaderStream.str();

        if (P.outSAMbool) {//
            *P.inOut->outSAM << P.samHeader;
        };
        if (P.outBAMunsorted){
            outBAMwriteHeader(P.inOut->outBAMfileUnsorted,P.samHeader,mainGenome.chrNameAll,mainGenome.chrLengthAll);
        };
//             if (P.outBAMcoord){
//                 outBAMwriteHeader(P.inOut->outBAMfileCoord,P.samHeader,mainGenome.chrName,mainGenome.chrLength);
//             };

        if ( P.quant.trSAM.yes ) {
            samHeaderStream.str("");
            vector <uint> trlength;
            for (uint32 ii=0;ii<mainTranscriptome->trID.size();ii++) {
                uint32 iex1=mainTranscriptome->trExI[ii]+mainTranscriptome->trExN[ii]-1; //last exon of the transcript
                trlength.push_back(mainTranscriptome->exLenCum[iex1]+mainTranscriptome->exSE[2*iex1+1]-mainTranscriptome->exSE[2*iex1]+1);
                samHeaderStream << "@SQ\tSN:"<< mainTranscriptome->trID.at(ii) <<"\tLN:"<<trlength.back()<<"\n";
            };
            for (uint32 ii=0;ii<P.outSAMattrRGlineSplit.size();ii++) {//@RG lines
                samHeaderStream << "@RG\t" << P.outSAMattrRGlineSplit.at(ii) <<"\n";
            };
            outBAMwriteHeader(P.inOut->outQuantBAMfile,samHeaderStream.str(),mainTranscriptome->trID,trlength);
        };

    };

    if (P.pCh.segmentMin>0) {
        if (P.pCh.out.junctions) {
            P.inOut->outChimJunction.open((P.outFileNamePrefix + "Chimeric.out.junction").c_str());
        };
        if (P.pCh.out.samOld) {
            P.inOut->outChimSAM.open((P.outFileNamePrefix + "Chimeric.out.sam").c_str());
            P.inOut->outChimSAM << P.samHeader;
        };
        pthread_mutex_init(&g_threadChunks.mutexOutChimSAM, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutChimJunction, NULL);
    };

    // P.inOut->logMain << "mlock value="<<mlockall(MCL_CURRENT|MCL_FUTURE) <<"\n"<<flush;

    // prepare chunks and spawn mapping threads
    ReadAlignChunk *RAchunk[P.runThreadN];
    for (int ii=0;ii<P.runThreadN;ii++) {
        RAchunk[ii]=new ReadAlignChunk(P, mainGenome, mainTranscriptome, ii);
    };

    mapThreadsSpawn(P, RAchunk);

    if (P.outFilterBySJoutStage==1) {//completed stage 1, go to stage 2
        P.inOut->logMain << "Completed stage 1 mapping of outFilterBySJout mapping\n"<<flush;
        outputSJ(RAchunk,P);//collapse novel junctions
        P.readFilesIndex=-1;


        P.outFilterBySJoutStage=2;
        if (P.outBAMcoord) {
            for (int it=0; it<P.runThreadN; it++) {//prepare the unmapped bin
                RAchunk[it]->chunkOutBAMcoord->coordUnmappedPrepareBySJout();
            };
        };

        mapThreadsSpawn(P, RAchunk);
    };

    //close some BAM files
    if (P.inOut->outBAMfileUnsorted!=NULL) {
        bgzf_flush(P.inOut->outBAMfileUnsorted);
        bgzf_close(P.inOut->outBAMfileUnsorted);
    };
    if (P.inOut->outQuantBAMfile!=NULL) {
        bgzf_flush(P.inOut->outQuantBAMfile);
        bgzf_close(P.inOut->outQuantBAMfile);
    };

    if (P.outBAMcoord && P.limitBAMsortRAM==0) {//make it equal ot the genome size
        P.limitBAMsortRAM=mainGenome.nGenome+mainGenome.SA.lengthByte+mainGenome.SAi.lengthByte;
    };

    //no need for genome anymore, free the memory
    mainGenome.freeMemory();

    if ( P.quant.geCount.yes )
    {//output gene quantifications
        for (int ichunk=1; ichunk<P.runThreadN; ichunk++)
        {//sum counts from all chunks into 0th chunk
            RAchunk[0]->chunkTr->quants->addQuants(*(RAchunk[ichunk]->chunkTr->quants));
        };
        RAchunk[0]->chunkTr->quantsOutput();
    };

    if (P.runThreadN>1 && P.outSAMorder=="PairedKeepInputOrder") {//concatenate Aligned.* files
        RAchunk[0]->chunkFilesCat(P.inOut->outSAM, P.outFileTmp + "/Aligned.out.sam.chunk", g_threadChunks.chunkOutN);
    };


    if (P.outBAMcoord) {//sort BAM if needed
        *P.inOut->logStdOut << timeMonthDayTime() << " ..... started sorting BAM\n" <<flush;
        P.inOut->logMain << timeMonthDayTime() << " ..... started sorting BAM\n" <<flush;
        uint32 nBins=P.outBAMcoordNbins;

        //check max size needed for sorting
        uint maxMem=0;
        for (uint32 ibin=0; ibin<nBins-1; ibin++) {//check all bins
            uint binS=0;
            for (int it=0; it<P.runThreadN; it++) {//collect sizes from threads
                binS += RAchunk[it]->chunkOutBAMcoord->binTotalBytes[ibin]+24*RAchunk[it]->chunkOutBAMcoord->binTotalN[ibin];
            };
            if (binS>maxMem) maxMem=binS;
        };
        P.inOut->logMain << "Max memory needed for sorting = "<<maxMem<<endl;
        if (maxMem>P.limitBAMsortRAM) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal ERROR: not enough memory for BAM sorting: \n";
            errOut <<"SOLUTION: re-run STAR with at least --limitBAMsortRAM " <<maxMem+1000000000;
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        };


        uint totalMem=0;
//         P.inOut->logMain << "Started sorting BAM ..." <<endl;
        #pragma omp parallel num_threads(P.outBAMsortingThreadNactual)
        #pragma omp for schedule (dynamic,1)
        for (uint32 ibin1=0; ibin1<nBins; ibin1++) {
            uint32 ibin=nBins-1-ibin1;//reverse order to start with the last bin - unmapped reads

            uint binN=0, binS=0;
            for (int it=0; it<P.runThreadN; it++) {//collect sizes from threads
                binN += RAchunk[it]->chunkOutBAMcoord->binTotalN[ibin];
                binS += RAchunk[it]->chunkOutBAMcoord->binTotalBytes[ibin];
            };

            if (binS==0) continue; //empty bin

            if (ibin == nBins-1) {//last bin for unmapped reads
                BAMbinSortUnmapped(ibin,P.runThreadN,P.outBAMsortTmpDir, P, mainGenome);
            } else {
            uint newMem=binS+binN*24;
            bool boolWait=true;
            while (boolWait) {
                #pragma omp critical
                if (totalMem+newMem < P.limitBAMsortRAM) {
                    boolWait=false;
                    totalMem+=newMem;
                };
                sleep(0.1);
            };
            BAMbinSortByCoordinate(ibin,binN,binS,P.runThreadN,P.outBAMsortTmpDir, P, mainGenome);
            #pragma omp critical
            totalMem-=newMem;//"release" RAM
        };
        };
        //concatenate all BAM files, using bam_cat
        char **bamBinNames = new char* [nBins];
        vector <string> bamBinNamesV;
        for (uint32 ibin=0; ibin<nBins; ibin++) {

            bamBinNamesV.push_back(P.outBAMsortTmpDir+"/b"+std::to_string((uint) ibin));
            struct stat buffer;
            if (stat (bamBinNamesV.back().c_str(), &buffer) != 0) {//check if file exists
                bamBinNamesV.pop_back();
            };
        };
        for (uint32 ibin=0; ibin<bamBinNamesV.size(); ibin++) {
                bamBinNames[ibin] = (char*) bamBinNamesV.at(ibin).c_str();
        };
        bam_cat(bamBinNamesV.size(), bamBinNames, 0, P.outBAMfileCoordName.c_str());
    };
    //wiggle output
    if (P.outWigFlags.yes) {
        *(P.inOut->logStdOut) << timeMonthDayTime() << " ..... started wiggle output\n" <<flush;
        P.inOut->logMain << timeMonthDayTime() << " ..... started wiggle output\n" <<flush;
        string wigOutFileNamePrefix=P.outFileNamePrefix + "Signal";
        signalFromBAM(P.outBAMfileCoordName, wigOutFileNamePrefix, P);
    };

    //aggregate output junctions
    //collapse splice junctions from different threads/chunks, and output them
    outputSJ(RAchunk,P);

    g_statsAll.progressReport(P.inOut->logProgress);
    P.inOut->logProgress  << "ALL DONE!\n"<<flush;
    P.inOut->logFinal.open((P.outFileNamePrefix + "Log.final.out").c_str());
    g_statsAll.reportFinal(P.inOut->logFinal);
    *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeFinish) << " ..... finished successfully\n" <<flush;

    P.inOut->logMain  << "ALL DONE!\n"<<flush;
    if (P.outTmpKeep=="None")
    {
        sysRemoveDir (P.outFileTmp);
    } else
    {
        //nothing to do
    };

    P.closeReadsFiles();//this will kill the readFilesCommand processes if necessary
    //mainGenome.~Genome(); //need explicit call because of the 'delete P.inOut' below, which will destroy P.inOut->logStdOut

    ///@todo create Destructor to destroy P.inOut
    delete P.inOut; //to close files

    return 0;
};
