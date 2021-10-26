#include "bamSortByCoordinate.h"
#include "BAMfunctions.h"
#include "BAMbinSortByCoordinate.h"
#include "BAMbinSortUnmapped.h"
#include "ErrorWarning.h"
#include "bam_cat.h"

void bamSortByCoordinate (Parameters &P, ReadAlignChunk **RAchunk, Genome &genome, Solo &solo) {
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

        uint64 unmappedReadsN = 0;
        for (int it=0; it<P.runThreadN; it++)
            unmappedReadsN += RAchunk[it]->chunkOutBAMcoord->binTotalN[nBins-1];

        P.inOut->logMain << "Max memory needed for sorting = "<<maxMem<<endl;
        if (maxMem>P.limitBAMsortRAM) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal ERROR: not enough memory for BAM sorting: \n";
            errOut <<"SOLUTION: re-run STAR with at least --limitBAMsortRAM " <<maxMem+1000000000;
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        } else if(maxMem==0 && unmappedReadsN==0) {//both mapped and unmapped reads are absent
            P.inOut->logMain << "WARNING: nothing to sort - no output alignments" <<endl;
            BGZF *bgzfOut;
            bgzfOut=bgzf_open(P.outBAMfileCoordName.c_str(),("w"+to_string((long long) P.outBAMcompression)).c_str());
            if (bgzfOut==NULL) {
                ostringstream errOut;
                errOut <<"EXITING because of fatal ERROR: could not open output bam file: " << P.outBAMfileCoordName << "\n";
                errOut <<"SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running STAR";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
            };
            outBAMwriteHeader(bgzfOut,P.samHeaderSortedCoord,genome.chrNameAll,genome.chrLengthAll);
            bgzf_close(bgzfOut);
        } else {//sort
            uint totalMem=0;
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
                    BAMbinSortUnmapped(ibin,P.runThreadN,P.outBAMsortTmpDir, P, genome, solo);
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
                    BAMbinSortByCoordinate(ibin,binN,binS,P.runThreadN,P.outBAMsortTmpDir, P, genome, solo);
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
    };    
};
