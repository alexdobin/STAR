#include <cmath>
#include <cstdlib>
#include <new>

#include "Genome.h"

#include "IncludeDefine.h"
#include "Parameters.h"
#include "SuffixArrayFuns.h"
#include "PackedArray.h"
#include "TimeFunctions.h"
#include "ErrorWarning.h"
#include "GTF.h"
#include "SjdbClass.h"
#include "sjdbLoadFromFiles.h"
#include "sjdbPrepare.h"
#include "genomeParametersWrite.h"
#include "sjdbInsertJunctions.h"
#include "genomeScanFastaFiles.h"
#include "genomeSAindex.h"

#include "serviceFuns.cpp"
#include "streamFuns.h"
#include "SequenceFuns.h"


char* globalG;
uint globalL;

static uint systemAvailableMemoryBytes()
{
#ifdef __linux__
    ifstream memInfo("/proc/meminfo");
    string field, unit;
    uint value=0;

    while (memInfo >> field >> value >> unit) {
        if (field=="MemAvailable:") {
            return value*1024LLU;
        };
    };
#endif
    return 0;
};


inline int funCompareSuffixesFromWord ( const void *a, const void *b, uint wordStart){

    uint *ga=(uint*)((globalG-7LLU)+(*((uint*)a)))-wordStart;
    uint *gb=(uint*)((globalG-7LLU)+(*((uint*)b)))-wordStart;

    uint jj=0;
    int  ii=0;
    uint va=0,vb=0;
    uint8 *va1, *vb1;

    while (jj+wordStart < globalL) {
        va=*(ga-jj);
        vb=*(gb-jj);

        #define has5(v) ((((v)^0x0505050505050505) - 0x0101010101010101) & ~((v)^0x0505050505050505) & 0x8080808080808080)

        if (has5(va) && has5(vb))
        {//there is 5 in the sequence - only compare bytes before 5
            va1=(uint8*) &va;
            vb1=(uint8*) &vb;
            for (ii=7;ii>=0;ii--)
            {
                if (va1[ii]>vb1[ii])
                {
                    return 1;
                } else if (va1[ii]<vb1[ii])
                {
                    return -1;
                } else if (va1[ii]==5)
                {//va=vb at the end of chr
                    if ( *((uint*)a) > *((uint*)b) )
                    {//anti-stable order,since indexes are sorted in the reverse order
                        return  -1;
                    } else
                    {//a cannot be equal to b
                        return 1;
                    };
                };
            };
        } else
        {//no 5, simple comparison
            if (va>vb)
            {
                return 1;
            } else if (va<vb)
            {
                return -1;
            };
        };
        jj++;
    };

    //suffixes are equal up to globalL, no simply compare the indexes
    if ( *((uint*)a) > *((uint*)b) )
    {//anti-stable order,since indexes are sorted in the reverse order
        return  -1;
    } else
    {//a cannot be equal to b
        return 1;
    };
};

inline int funCompareSuffixes ( const void *a, const void *b){
    return funCompareSuffixesFromWord(a,b,0);
};

inline int funCompareSuffixesSkipFirstWord ( const void *a, const void *b){
    return funCompareSuffixesFromWord(a,b,1);
};

inline uint funG2strLocus (uint SAstr, uint const N, char const GstrandBit, uint const GstrandMask) {
    bool strandG = (SAstr>>GstrandBit) == 0;
    SAstr &= GstrandMask;
    if ( !strandG ) SAstr += N;
    return SAstr;
};

inline uint funSAsortPrefix(char* G, uint ii, uint prefixLength) {
    uint prefix=0;
    for (uint jj=0; jj<prefixLength; jj++) {
        uint g=(uint) G[ii-jj];
        if (g>5) g=5;
        prefix=prefix*6+g;
    };
    return prefix;
};

inline uint funSAsortPrefixAtOffset(char* G, uint ii, uint offset, uint prefixLength) {
    uint prefix=0;
    for (uint jj=0; jj<prefixLength; jj++) {
        uint g=(uint) G[ii-offset-jj];
        if (g>5) g=5;
        prefix=prefix*6+g;
    };
    return prefix;
};

inline bool funSAsortPrefixFirstWordHas5(uint prefix, uint prefixLength) {
    uint wordBases=(uint) sizeof(uint);
    if (prefixLength>wordBases) {
        for (uint jj=0; jj<prefixLength-wordBases; jj++) {
            prefix /= 6;
        };
        prefixLength=wordBases;
    };
    for (uint jj=0; jj<prefixLength; jj++) {
        if (prefix % 6 == 5) {
            return true;
        };
        prefix /= 6;
    };
    return false;
};

void Genome::genomeGenerate() {

    //check parameters
	createDirectory(pGe.gDir, P.runDirPerm, "--genomeDir", P);

	{//move Log.out file into genome directory
		string logfn=pGe.gDir+"Log.out";
		if ( rename( P.outLogFileName.c_str(), logfn.c_str() ) ) {
			warningMessage("Could not move Log.out file from " + P.outLogFileName + " into " + logfn + ". Will keep " + P.outLogFileName +"\n", \
						   std::cerr, P.inOut->logMain, P);
		} else {
			P.outLogFileName=logfn;
		};
	};
    if (sjdbOverhang<=0 && (pGe.sjdbFileChrStartEnd.at(0)!="-" || pGe.sjdbGTFfile!="-")) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT PARAMETER ERROR: for generating genome with annotations (--sjdbFileChrStartEnd or --sjdbGTFfile options)\n";
        errOut << "you need to specify >0 --sjdbOverhang\n";
        errOut << "SOLUTION: re-run genome generation specifying non-zero --sjdbOverhang, which ideally should be equal to OneMateLength-1, or could be chosen generically as ~100\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    if (pGe.sjdbFileChrStartEnd.at(0)=="-" && pGe.sjdbGTFfile=="-") {
        if (P.parArray.at(P.pGe.sjdbOverhang_par)->inputLevel>0 && sjdbOverhang>0) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL INPUT PARAMETER ERROR: when generating genome without annotations (--sjdbFileChrStartEnd or --sjdbGTFfile options)\n";
            errOut << "do not specify >0 --sjdbOverhang\n";
            errOut << "SOLUTION: re-run genome generation without --sjdbOverhang option\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };
        sjdbOverhang=0;
    };
    
    //time
    time_t rawTime;
    string timeString;

    time(&rawTime);
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ... starting to generate Genome files\n" <<flush;
    *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... starting to generate Genome files\n" <<flush;

    //define some parameters from input parameters
    genomeChrBinNbases=1LLU << pGe.gChrBinNbits;

    nGenome = genomeScanFastaFiles(P,NULL,false,*this);//first scan the fasta file to find all the sizes
    genomeSequenceAllocate(nGenome, nG1alloc, G, G1);
    genomeScanFastaFiles(P,G,true,*this);    //load the genome sequence

    uint64 nGenomeTrue=0;
    for (auto &cl : chrLength)
    	nGenomeTrue += cl; //nGenomeTrue = sum of chr lengths

    P.inOut->logMain <<"Genome sequence total length = " << nGenomeTrue << "\n";
    P.inOut->logMain <<"Genome size with padding = "<< nGenome <<"\n";

    //consensusSequence(); //replace with consensus allele DEPRECATED
        
    SjdbClass sjdbLoci; //will be filled in transcriptGeneSJ below
    GTF mainGTF(*this, P, pGe.gDir, sjdbLoci); //this loads exonLoci and gene/transcript metadata only, sjdbLoci is not filled
    
    Genome::transformGenome(&mainGTF);
    
    mainGTF.superTranscript(); //this may change the genome into (Super)Transcriptome

    chrBinFill();//chrBin is first used in the transcriptGeneSJ below
    mainGTF.transcriptGeneSJ(pGe.gDir);
    
    sjdbLoadFromFiles(P,sjdbLoci);//this will not be transformed. TODO prevent this parameter combination

    if (pGe.gSAindexNbases > log2(nGenomeTrue)/2-1) {
        ostringstream warnOut; 
        warnOut << "--genomeSAindexNbases " << pGe.gSAindexNbases << " is too large for the genome size=" << nGenomeTrue;
        warnOut << ", which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases " << int(log2(nGenomeTrue)/2-1);
        warningMessage(warnOut.str(),P.inOut->logMain,std::cerr,P);
    };    
    
    //output genome metadata
    writeChrInfo(pGe.gDir);

    //preparing to generate SA
    for (uint ii=0;ii<nGenome;ii++) {//- strand
        //if (G[ii]>5)
        //    cerr << ii <<" "<< G[ii]<<"\n";
        G[2*nGenome-1-ii]=G[ii]<4 ? 3-G[ii] : G[ii];
    };   
    nSA=0;
    for (uint ii=0;ii<2*nGenome;ii+=pGe.gSAsparseD) {
        if (G[ii]<4) {
            nSA++;
        };
    };

    // GstrandBit
    GstrandBit = (char) (uint) floor(log(nGenome+P.limitSjdbInsertNsj*sjdbLength)/log(2))+1; //GstrandBit uses P.limitSjdbInsertNsj even if no insertion requested, in case it will be requested at the mapping stage
    if (GstrandBit<32) GstrandBit=32; //TODO: should not this be 31? Need to test for small genomes. TODO: use simple access function for SA
    P.inOut->logMain <<"Estimated genome size with padding and SJs: total=genome+SJ="<<nGenome+P.limitSjdbInsertNsj*sjdbLength<<" = "<<nGenome<<" + "<<P.limitSjdbInsertNsj*sjdbLength<<"\n";
    P.inOut->logMain << "GstrandBit=" << int(GstrandBit) <<"\n";
    GstrandMask = ~(1LLU<<GstrandBit);
    SA.defineBits(GstrandBit+1,nSA);

    if (P.sjdbInsert.yes) {//reserve space for junction insertion       
        SApass1.defineBits( GstrandBit+1, nSA+2*sjdbLength*min((uint64)sjdbLoci.chr.size(),P.limitSjdbInsertNsj) );//TODO: this allocation is wasteful, get a better estimate of the number of junctions
    } else {//same as SA
        SApass1.defineBits(GstrandBit+1,nSA);
    };

    P.inOut->logMain  << "Number of SA indices: "<< nSA << "\n"<<flush;

    //sort SA
    time ( &rawTime );
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ... starting to sort Suffix Array. This may take a long time...\n" <<flush;
    *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... starting to sort Suffix Array. This may take a long time...\n" <<flush;


//     if (false)
    {//sort SA chunks

        for (uint ii=0;ii<nGenome;ii++) {//re-fill the array backwards for sorting
            swap(G[2*nGenome-1-ii],G[ii]);
        };
        globalG=G;
        globalL=pGe.gSuffixLengthMax/sizeof(uint);
        // Count sort prefixes over STAR's 0..5 genome alphabet.
        uint indPrefLen=6;
        const uint64 saAvailableBytesForPrefix=P.limitGenomeGenerateRAM>nG1alloc ? P.limitGenomeGenerateRAM-nG1alloc : 0;
        if (P.runThreadN>=16 && nGenome>=50000000) {
            uint indPrefN8=1;
            for (uint ii=0;ii<8;ii++) indPrefN8 *= 6;
            const uint64 indPrefTempBytes8=(uint64) P.runThreadN*indPrefN8*2*sizeof(uint) + (indPrefN8+1)*sizeof(uint);
            if (indPrefTempBytes8 < saAvailableBytesForPrefix/4) {
                indPrefLen=8;
            };
        };
        uint indPrefN=1;
        for (uint ii=0;ii<indPrefLen;ii++) indPrefN *= 6;
        uint* indPrefCount = new uint [indPrefN];
        memset(indPrefCount,0,indPrefN*sizeof(indPrefCount[0]));
        const bool saSortProfile=getenv("STAR_PROFILE_SA_SORT")!=NULL;
        nSA=0;
        for (uint ii=0;ii<2*nGenome;ii+=pGe.gSAsparseD) {
            if (G[ii]<4) {
                uint p1=funSAsortPrefix(G,ii,indPrefLen);
                indPrefCount[p1]++;
                nSA++;
            };
        };

        uint indPrefEmptyN=0;
        uint indPrefSingletonN=0;
        uint indPrefNonEmptyN=0;
        uint indPrefSortableN=0;
        uint indPrefMaxCount=0;
        uint indPrefMaxIndex=0;
        if (saSortProfile) {
            for (uint ii=0; ii<indPrefN; ii++) {
                if (indPrefCount[ii]==0) {
                    indPrefEmptyN++;
                } else {
                    indPrefNonEmptyN++;
                    if (indPrefCount[ii]==1) {
                        indPrefSingletonN++;
                    } else {
                        indPrefSortableN++;
                    };
                    if (indPrefCount[ii]>indPrefMaxCount) {
                        indPrefMaxCount=indPrefCount[ii];
                        indPrefMaxIndex=ii;
                    };
                };
            };
        };

        uint saChunkSize=(P.limitGenomeGenerateRAM-nG1alloc)/8/P.runThreadN; //number of SA indexes per chunk
        saChunkSize=saChunkSize*6/10; //allow extra space for qsort
        //uint saChunkN=((nSA/saChunkSize+1)/P.runThreadN+1)*P.runThreadN;//ensure saChunkN is divisible by P.runThreadN
        //saChunkSize=nSA/saChunkN+100000;//final chunk size
        if (P.runThreadN>1) saChunkSize=min(saChunkSize,nSA/(P.runThreadN-1));
        uint saChunkN = nSA / saChunkSize + 1;//estimate
        uint* indPrefStart = new uint [saChunkN*2]; //start and stop, *2 just in case
        uint* indPrefChunkCount = new uint [saChunkN*2];
        indPrefStart[0]=0;
        saChunkN=0;//start counting chunks
        uint chunkSize1=indPrefCount[0];
        for (uint ii=1; ii<indPrefN; ii++) {
            chunkSize1 += indPrefCount[ii];
            if (chunkSize1 > saChunkSize) {
                saChunkN++;
                indPrefStart[saChunkN]=ii;
                indPrefChunkCount[saChunkN-1]=chunkSize1-indPrefCount[ii];
                chunkSize1=indPrefCount[ii];
            };
        };
        saChunkN++;
        indPrefStart[saChunkN]=indPrefN+1;
        indPrefChunkCount[saChunkN-1]=chunkSize1;

        uint* saChunkStart = new uint [saChunkN+1];
        saChunkStart[0]=0;
        for (uint iChunk=0; iChunk<saChunkN; iChunk++) {
            saChunkStart[iChunk+1]=saChunkStart[iChunk]+indPrefChunkCount[iChunk];
        };

        P.inOut->logMain  << "Number of chunks: " << saChunkN <<";   chunks size limit: " << saChunkSize*8 <<" bytes\n" <<flush;

        time ( &rawTime );
        P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ... sorting Suffix Array chunks...\n" <<flush;
        *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... sorting Suffix Array chunks...\n" <<flush;

        const uint genomeScanN=(2*nGenome+pGe.gSAsparseD-1)/pGe.gSAsparseD;
        const uint64 saAvailableBytes=P.limitGenomeGenerateRAM>nG1alloc ? P.limitGenomeGenerateRAM-nG1alloc : 0;
        const uint64 saAllChunkScatterBytes=nSA*sizeof(uint) + (uint64) P.runThreadN*indPrefN*2*sizeof(uint) + (indPrefN+1)*sizeof(uint) + (saChunkN+1)*sizeof(uint);
        const uint saChunkBatchSize=max(saChunkSize*(uint) P.runThreadN,saChunkSize);
        uint saChunkBatchN=0;
        if (P.runThreadN>=16 && saChunkN>1) {
            uint iChunkBatch=0;
            while (iChunkBatch<saChunkN) {
                uint iChunkBatchNext=iChunkBatch;
                uint saBatchN=0;
                while (iChunkBatchNext<saChunkN && (saBatchN==0 || saBatchN+indPrefChunkCount[iChunkBatchNext]<=saChunkBatchSize)) {
                    saBatchN += indPrefChunkCount[iChunkBatchNext];
                    iChunkBatchNext++;
                };
                saChunkBatchN++;
                iChunkBatch=iChunkBatchNext;
            };
        };

        // Batched fill trades one genome scan per chunk for two scans per batch.
        // It pays off once enough chunks are being sorted in parallel.
        const bool saChunkBatchFill=P.runThreadN>=16 && saChunkN>1;
        const uint64 saRetainedChunkBytes=nSA*sizeof(uint);
        const uint64 saRamPackPeakBytes=saRetainedChunkBytes+SApass1.lengthByte;
        const uint64 saRamPeakBytes=max(saAllChunkScatterBytes,saRamPackPeakBytes);
        const uint64 saRamHeadroomBytes=max(saRamPeakBytes/10,(uint64) 2000000000LLU);
        const uint64 saRamRequiredBytes=saRamPeakBytes+saRamHeadroomBytes;
        const uint64 systemAvailableBytes=systemAvailableMemoryBytes();
        const bool saRamLimitOK=saAvailableBytes>=saRamRequiredBytes;
        const bool saRamSystemOK=systemAvailableBytes>0 && systemAvailableBytes>=saRamRequiredBytes;
        bool saChunksInMemoryActive=saRamLimitOK && saRamSystemOK;
        P.inOut->logMain  << "SA chunk all-scatter available bytes: " << saAvailableBytes << "; required temporary bytes: " << saAllChunkScatterBytes << "\n" <<flush;
        P.inOut->logMain  << "SA chunk estimated batched-fill batches: " << saChunkBatchN << "\n" <<flush;
        P.inOut->logMain  << "SA chunk fill strategy: " << (saChunkBatchFill ? "batched" : "per-chunk") << "\n" <<flush;
        P.inOut->logMain  << "SA chunk sort prefix length: " << indPrefLen << "\n" <<flush;
        P.inOut->logMain  << "SA chunk sort granularity: " << (saChunkBatchFill ? "prefix-bin" : "chunk") << "\n" <<flush;
        P.inOut->logMain  << "SA chunk retained bytes: " << saRetainedChunkBytes << "; RAM peak bytes: " << saRamPeakBytes << "; RAM headroom bytes: " << saRamHeadroomBytes << "\n" <<flush;
        P.inOut->logMain  << "SA chunk system available bytes: " << systemAvailableBytes << "; limit available bytes: " << saAvailableBytes << "\n" <<flush;
        if (saSortProfile) {
            uint saChunkMinCount=indPrefChunkCount[0];
            uint saChunkMaxCount=indPrefChunkCount[0];
            uint saChunkMaxIndex=0;
            for (uint iChunk=1; iChunk<saChunkN; iChunk++) {
                if (indPrefChunkCount[iChunk]<saChunkMinCount) {
                    saChunkMinCount=indPrefChunkCount[iChunk];
                };
                if (indPrefChunkCount[iChunk]>saChunkMaxCount) {
                    saChunkMaxCount=indPrefChunkCount[iChunk];
                    saChunkMaxIndex=iChunk;
                };
            };
            P.inOut->logMain << "SA sort profile: enabled\n" <<flush;
            P.inOut->logMain << "SA sort profile prefix bins: total=" << indPrefN
                              << "; empty=" << indPrefEmptyN
                              << "; nonempty=" << indPrefNonEmptyN
                              << "; singleton=" << indPrefSingletonN
                              << "; sortable=" << indPrefSortableN
                              << "; max_bin=" << indPrefMaxCount
                              << "; max_prefix=" << indPrefMaxIndex << "\n" <<flush;
            P.inOut->logMain << "SA sort profile chunks: count=" << saChunkN
                              << "; min_indices=" << saChunkMinCount
                              << "; max_indices=" << saChunkMaxCount
                              << "; max_chunk=" << saChunkMaxIndex
                              << "; mean_indices=" << (saChunkN>0 ? nSA/saChunkN : 0) << "\n" <<flush;
        };

        uint* saChunksInMemory=NULL;
        if (saChunksInMemoryActive) {
            saChunksInMemory=new (nothrow) uint [nSA];
            if (saChunksInMemory==NULL) {
                saChunksInMemoryActive=false;
                P.inOut->logMain << "SA chunk storage fallback: disk-backed after RAM allocation failure\n" <<flush;
            };
        };
        P.inOut->logMain  << "SA chunk storage strategy: " << (saChunksInMemoryActive ? "RAM-retained" : "disk-backed") << "\n" <<flush;

        double saProfileFillCountSeconds=0.0;
        double saProfileFillScatterSeconds=0.0;
        double saProfileSortSeconds=0.0;
        double saProfileFinalizeSeconds=0.0;
        uint saProfileSkipFirstWordRangeN=0;
        uint saProfileSkipFirstWordIndexN=0;
        uint saProfileSubBinPrefixN=0;
        uint saProfileSubBinRangeN=0;
        uint saProfileSubBinIndexN=0;

        if (saChunkBatchFill) {
            uint iChunkBatch=0;
            while (iChunkBatch<saChunkN) {//fill a batch of chunks with two genome scans, then sort and write each chunk
                uint iChunkBatchNext=iChunkBatch;
                uint saBatchN=0;
                while (iChunkBatchNext<saChunkN && (saBatchN==0 || saBatchN+indPrefChunkCount[iChunkBatchNext]<=saChunkBatchSize)) {
                    saBatchN += indPrefChunkCount[iChunkBatchNext];
                    iChunkBatchNext++;
                };
                uint batchChunkN=iChunkBatchNext-iChunkBatch;
                uint* saBatch=saChunksInMemoryActive ? saChunksInMemory+saChunkStart[iChunkBatch] : new uint [saBatchN];
                uint* saBatchStart=new uint [batchChunkN+1];
                saBatchStart[0]=0;
                for (uint iChunk=0;iChunk<batchChunkN;iChunk++) {
                    saBatchStart[iChunk+1]=saBatchStart[iChunk]+indPrefChunkCount[iChunkBatch+iChunk];
                };

                uint batchPrefStart=indPrefStart[iChunkBatch];
                uint batchPrefEnd=min(indPrefStart[iChunkBatchNext],indPrefN);
                uint batchPrefN=batchPrefEnd-batchPrefStart;

                uint* saPrefStart=new uint [batchPrefN+1];
                saPrefStart[0]=0;
                for (uint iPref=0;iPref<batchPrefN;iPref++) {
                    saPrefStart[iPref+1]=saPrefStart[iPref]+indPrefCount[batchPrefStart+iPref];
                };
                if (saPrefStart[batchPrefN]!=saBatchN) {
                    ostringstream errOut;
                    errOut << "EXITING because of FATAL problem while generating the suffix array\n";
                    errOut << "The number of indices collected for suffix array batch " << iChunkBatch;
                    errOut << " = " << saPrefStart[batchPrefN] << " is not equal to expected " << saBatchN << "\n";
                    errOut << "SOLUTION: try to re-run suffix array generation, if it still does not work, report this problem to the author\n"<<flush;
                    exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                };

                uint* saThreadPrefCount=new uint [(uint) P.runThreadN*batchPrefN];
                memset(saThreadPrefCount,0,sizeof(saThreadPrefCount[0])*(uint) P.runThreadN*batchPrefN);

                double saProfileTime0=omp_get_wtime();
                #pragma omp parallel num_threads(P.runThreadN)
                {
                    int tid=omp_get_thread_num();
                    uint* threadPrefCount=saThreadPrefCount+(uint) tid*batchPrefN;
                    #pragma omp for schedule(static)
                    for (uint iScan=0;iScan<genomeScanN;iScan++) {
                        uint ii=iScan*pGe.gSAsparseD;
                        if (G[ii]<4) {
                            uint p1=funSAsortPrefix(G,ii,indPrefLen);
                            if (p1>=batchPrefStart && p1<batchPrefEnd) {
                                threadPrefCount[p1-batchPrefStart]++;
                            };
                        };
                    };
                };
                if (saSortProfile) {
                    saProfileFillCountSeconds += omp_get_wtime()-saProfileTime0;
                };

                uint* saThreadPrefStart=new uint [(uint) P.runThreadN*batchPrefN];
                for (uint iPref=0;iPref<batchPrefN;iPref++) {
                    uint prefCount=0;
                    for (int iThread=0;iThread<P.runThreadN;iThread++) {
                        saThreadPrefStart[(uint) iThread*batchPrefN+iPref]=saPrefStart[iPref]+prefCount;
                        prefCount += saThreadPrefCount[(uint) iThread*batchPrefN+iPref];
                    };
                    if (prefCount!=indPrefCount[batchPrefStart+iPref]) {
                        ostringstream errOut;
                        errOut << "EXITING because of FATAL problem while generating the suffix array\n";
                        errOut << "The number of indices collected for suffix array prefix " << batchPrefStart+iPref << " = " << prefCount;
                        errOut << " is not equal to expected " << indPrefCount[batchPrefStart+iPref] << "\n";
                        errOut << "SOLUTION: try to re-run suffix array generation, if it still does not work, report this problem to the author\n"<<flush;
                        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                    };
                };

                saProfileTime0=omp_get_wtime();
                #pragma omp parallel num_threads(P.runThreadN)
                {
                    int tid=omp_get_thread_num();
                    uint* threadPrefStart=saThreadPrefStart+(uint) tid*batchPrefN;
                    #pragma omp for schedule(static)
                    for (uint iScan=0;iScan<genomeScanN;iScan++) {
                        uint ii=iScan*pGe.gSAsparseD;
                        if (G[ii]<4) {
                            uint p1=funSAsortPrefix(G,ii,indPrefLen);
                            if (p1>=batchPrefStart && p1<batchPrefEnd) {
                                uint iPref=p1-batchPrefStart;
                                saBatch[threadPrefStart[iPref]]=ii;
                                threadPrefStart[iPref]++;
                            };
                        };
                    };
                };
                if (saSortProfile) {
                    saProfileFillScatterSeconds += omp_get_wtime()-saProfileTime0;
                };

                saProfileTime0=omp_get_wtime();
                const bool saCanSkipFirstWord=indPrefLen>=(uint) sizeof(uint) && globalL>1;
                const uint saSubBinExtraLen=3;
                const uint saSubBinN=216; // 6^3
                const uint saSubBinMinN=1000000;
                vector<uint*> saSortRange;
                vector<uint> saSortRangeN;
                vector<char> saSortRangeSkipFirstWord;
                saSortRange.reserve(batchPrefN);
                saSortRangeN.reserve(batchPrefN);
                saSortRangeSkipFirstWord.reserve(batchPrefN);
                for (uint iPref=0; iPref < batchPrefN; iPref++) {
                    uint* saPref=saBatch+saPrefStart[iPref];
                    uint saPrefN=saPrefStart[iPref+1]-saPrefStart[iPref];
                    if (saPrefN>1) {
                        bool skipFirstWord=saCanSkipFirstWord && !funSAsortPrefixFirstWordHas5(batchPrefStart+iPref,indPrefLen);
                        bool subBinned=false;
                        if (skipFirstWord && saPrefN>=saSubBinMinN) {
                            uint saSubBinCount[saSubBinN];
                            memset(saSubBinCount,0,sizeof(saSubBinCount));
                            bool subBinOK=true;
                            for (uint ii=0; ii<saPrefN; ii++) {
                                if (saPref[ii]<indPrefLen+saSubBinExtraLen-1) {
                                    subBinOK=false;
                                    break;
                                };
                                saSubBinCount[funSAsortPrefixAtOffset(G,saPref[ii],indPrefLen,saSubBinExtraLen)]++;
                            };

                            uint saSubBinNonEmptyN=0;
                            uint saSubBinStart[saSubBinN+1];
                            saSubBinStart[0]=0;
                            if (subBinOK) {
                                for (uint ii=0; ii<saSubBinN; ii++) {
                                    if (saSubBinCount[ii]>0) {
                                        saSubBinNonEmptyN++;
                                    };
                                    saSubBinStart[ii+1]=saSubBinStart[ii]+saSubBinCount[ii];
                                };
                            };

                            uint* saPrefTmp=NULL;
                            if (saSubBinNonEmptyN>1) {
                                saPrefTmp=new (nothrow) uint [saPrefN];
                            };
                            if (saPrefTmp!=NULL) {
                                uint saSubBinCursor[saSubBinN];
                                memcpy(saSubBinCursor,saSubBinStart,sizeof(saSubBinCursor));
                                for (uint ii=0; ii<saPrefN; ii++) {
                                    uint iSub=funSAsortPrefixAtOffset(G,saPref[ii],indPrefLen,saSubBinExtraLen);
                                    saPrefTmp[saSubBinCursor[iSub]]=saPref[ii];
                                    saSubBinCursor[iSub]++;
                                };
                                memcpy(saPref,saPrefTmp,saPrefN*sizeof(saPref[0]));
                                delete [] saPrefTmp;

                                uint saSubBinSortRangeN=0;
                                for (uint iSub=0; iSub<saSubBinN; iSub++) {
                                    uint saSubN=saSubBinStart[iSub+1]-saSubBinStart[iSub];
                                    if (saSubN>1) {
                                        saSortRange.push_back(saPref+saSubBinStart[iSub]);
                                        saSortRangeN.push_back(saSubN);
                                        saSortRangeSkipFirstWord.push_back(1);
                                        saSubBinSortRangeN++;
                                    };
                                };
                                if (saSortProfile) {
                                    saProfileSubBinPrefixN++;
                                    saProfileSubBinRangeN += saSubBinSortRangeN;
                                    saProfileSubBinIndexN += saPrefN;
                                };
                                subBinned=true;
                            };
                        };
                        if (!subBinned) {
                            saSortRange.push_back(saPref);
                            saSortRangeN.push_back(saPrefN);
                            saSortRangeSkipFirstWord.push_back(skipFirstWord ? 1 : 0);
                        };
                    };
                };

                #pragma omp parallel for num_threads(P.runThreadN) schedule(dynamic,16) reduction(+:saProfileSkipFirstWordRangeN,saProfileSkipFirstWordIndexN)
                for (int iRange=0; iRange < (int) saSortRangeN.size(); iRange++) {
                    uint* saRange=saSortRange[iRange];
                    uint saRangeN=saSortRangeN[iRange];
                    if (saSortRangeSkipFirstWord[iRange]) {
                        qsort(saRange,saRangeN,sizeof(saRange[0]),funCompareSuffixesSkipFirstWord);
                        if (saSortProfile) {
                            saProfileSkipFirstWordRangeN++;
                            saProfileSkipFirstWordIndexN += saRangeN;
                        };
                    } else {
                        qsort(saRange,saRangeN,sizeof(saRange[0]),funCompareSuffixes);
                    };
                };
                if (saSortProfile) {
                    saProfileSortSeconds += omp_get_wtime()-saProfileTime0;
                };

                saProfileTime0=omp_get_wtime();
                #pragma omp parallel for num_threads(P.runThreadN) schedule(dynamic,1)
                for (int iChunk=0; iChunk < (int) batchChunkN; iChunk++) {
                    uint iChunkAbs=iChunkBatch+(uint) iChunk;
                    uint* saChunk=saBatch+saBatchStart[iChunk];
                    for (uint ii=0;ii<indPrefChunkCount[iChunkAbs];ii++) {
                        saChunk[ii]=2*nGenome-1-saChunk[ii];
                    };
                    if (!saChunksInMemoryActive) {
                        string chunkFileName=pGe.gDir+"/SA_"+to_string(iChunkAbs);
                        ofstream & saChunkFile = ofstrOpen(chunkFileName,ERROR_OUT, P);
                        fstreamWriteBig(saChunkFile, (char*) saChunk, sizeof(saChunk[0])*indPrefChunkCount[iChunkAbs],chunkFileName,ERROR_OUT,P);
                        saChunkFile.close();
                    };
                };
                if (saSortProfile) {
                    saProfileFinalizeSeconds += omp_get_wtime()-saProfileTime0;
                };

                delete [] saThreadPrefStart;
                delete [] saThreadPrefCount;
                delete [] saPrefStart;
                delete [] saBatchStart;
                if (!saChunksInMemoryActive) {
                    delete [] saBatch;
                };
                iChunkBatch=iChunkBatchNext;
            };
        } else {
            #pragma omp parallel for num_threads(P.runThreadN) schedule(dynamic,1) reduction(+:saProfileFillScatterSeconds,saProfileSortSeconds,saProfileFinalizeSeconds)
            for (int iChunk=0; iChunk < (int) saChunkN; iChunk++) {//start the chunk cycle: sort each chunk with qsort and write to a file
                uint* saChunk=saChunksInMemoryActive ? saChunksInMemory+saChunkStart[iChunk] : new uint [indPrefChunkCount[iChunk]];//allocate local array for each chunk
                double saProfileTime0=omp_get_wtime();
                for (uint ii=0,jj=0;ii<2*nGenome;ii+=pGe.gSAsparseD) {//fill the chunk with SA indices
                    if (G[ii]<4) {
                        uint p1=funSAsortPrefix(G,ii,indPrefLen);
                        if (p1>=indPrefStart[iChunk] && p1<indPrefStart[iChunk+1]) {
                            saChunk[jj]=ii;
                            jj++;
                            if (jj==indPrefChunkCount[iChunk]) break;
                        };
                    };
                };
                if (saSortProfile) {
                    saProfileFillScatterSeconds += omp_get_wtime()-saProfileTime0;
                };

                //sort the chunk
                saProfileTime0=omp_get_wtime();
                if (indPrefChunkCount[iChunk]>1) {
                    qsort(saChunk,indPrefChunkCount[iChunk],sizeof(saChunk[0]),funCompareSuffixes);
                };
                if (saSortProfile) {
                    saProfileSortSeconds += omp_get_wtime()-saProfileTime0;
                };
                saProfileTime0=omp_get_wtime();
                for (uint ii=0;ii<indPrefChunkCount[iChunk];ii++) {
                    saChunk[ii]=2*nGenome-1-saChunk[ii];
                };
                //write files
                if (!saChunksInMemoryActive) {
                    string chunkFileName=pGe.gDir+"/SA_"+to_string( (uint) iChunk);
                    ofstream & saChunkFile = ofstrOpen(chunkFileName,ERROR_OUT, P);
                    fstreamWriteBig(saChunkFile, (char*) saChunk, sizeof(saChunk[0])*indPrefChunkCount[iChunk],chunkFileName,ERROR_OUT,P);
                    saChunkFile.close();
                    delete [] saChunk;
                    saChunk=NULL;
                };
                if (saSortProfile) {
                    saProfileFinalizeSeconds += omp_get_wtime()-saProfileTime0;
                };
            };
        };

        if (saSortProfile) {
            ostringstream saProfileOut;
            saProfileOut << fixed << setprecision(3);
            saProfileOut << "SA sort profile phases: fill_count_seconds=" << saProfileFillCountSeconds
                         << "; fill_scatter_seconds=" << saProfileFillScatterSeconds
                         << "; sort_seconds=" << saProfileSortSeconds
                         << "; finalize_seconds=" << saProfileFinalizeSeconds;
            P.inOut->logMain << saProfileOut.str() << "\n" <<flush;
            P.inOut->logMain << "SA sort profile comparator: skip_first_word_ranges=" << saProfileSkipFirstWordRangeN
                              << "; skip_first_word_indices=" << saProfileSkipFirstWordIndexN << "\n" <<flush;
            P.inOut->logMain << "SA sort profile sub-bins: prefixes=" << saProfileSubBinPrefixN
                              << "; ranges=" << saProfileSubBinRangeN
                              << "; indices=" << saProfileSubBinIndexN << "\n" <<flush;
        };

        time ( &rawTime );
        P.inOut->logMain     << timeMonthDayTime(rawTime) << (saChunksInMemoryActive ? " ... packing SA from RAM chunks...\n" : " ... loading chunks from disk, packing SA...\n") <<flush;
        *P.inOut->logStdOut  << timeMonthDayTime(rawTime) << (saChunksInMemoryActive ? " ... packing SA from RAM chunks...\n" : " ... loading chunks from disk, packing SA...\n") <<flush;

        //read chunks and pack into full SA
        SApass1.allocateArray();
        SA.pointArray(SApass1.charArray + SApass1.lengthByte-SA.lengthByte); //SA is shifted to have space for junction insertion
        uint N2bit= 1LLU << GstrandBit;
        uint packedInd=0;

        #ifdef genenomeGenerate_SA_textOutput
                ofstream SAtxtStream ((pGe.gDir + "/SAtxt").c_str());
        #endif

        if (saChunksInMemoryActive) {
            for (uint ii=0; ii<nSA; ii++) {
                SA.writePacked( ii, (saChunksInMemory[ii]<nGenome) ? saChunksInMemory[ii] : ( (saChunksInMemory[ii]-nGenome) | N2bit ) );

                #ifdef genenomeGenerate_SA_textOutput
                    SAtxtStream << saChunksInMemory[ii] << "\n";
                #endif
            };
            packedInd=nSA;
            delete [] saChunksInMemory;
            saChunksInMemory=NULL;
        } else {
            #define SA_CHUNK_BLOCK_SIZE 10000000
            uint* saIn=new uint[SA_CHUNK_BLOCK_SIZE]; //TODO make adjustable

            for (uint iChunk=0;iChunk<saChunkN;iChunk++) {//load files one by one and convert to packed
                ostringstream saChunkFileNameStream("");
                saChunkFileNameStream<< pGe.gDir << "/SA_" << iChunk;
                ifstream saChunkFile(saChunkFileNameStream.str().c_str());
                while (! saChunkFile.eof()) {//read blocks from each file
                    uint chunkBytesN=fstreamReadBig(saChunkFile,(char*) saIn,SA_CHUNK_BLOCK_SIZE*sizeof(saIn[0]));
                    for (uint ii=0;ii<chunkBytesN/sizeof(saIn[0]);ii++) {
                        SA.writePacked( packedInd+ii, (saIn[ii]<nGenome) ? saIn[ii] : ( (saIn[ii]-nGenome) | N2bit ) );

                        #ifdef genenomeGenerate_SA_textOutput
                            SAtxtStream << saIn[ii] << "\n";
                        #endif
                    };
                    packedInd += chunkBytesN/sizeof(saIn[0]);
                };
                saChunkFile.close();
                remove(saChunkFileNameStream.str().c_str());//remove the chunk file
            };
            delete [] saIn;
        };

        #ifdef genenomeGenerate_SA_textOutput
                SAtxtStream.close();
        #endif

        if (packedInd != nSA ) {//
            ostringstream errOut;
            errOut << "EXITING because of FATAL problem while generating the suffix array\n";
            errOut << "The number of indices read from chunks = "<<packedInd<<" is not equal to expected nSA="<<nSA<<"\n";
            errOut << "SOLUTION: try to re-run suffix array generation, if it still does not work, report this problem to the author\n"<<flush;
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };

        //DONE with suffix array generation

        for (uint ii=0;ii<nGenome;ii++) {//return to normal order for future use
            swap(G[2*nGenome-1-ii],G[ii]);
        };
        delete [] indPrefCount;
        delete [] indPrefStart;
        delete [] indPrefChunkCount;
        delete [] saChunkStart;
    };

    time ( &rawTime );
    timeString=asctime(localtime ( &rawTime ));
    timeString.erase(timeString.end()-1,timeString.end());
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ... finished generating suffix array\n" <<flush;
    *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... finished generating suffix array\n" <<flush;

    genomeSAindex(G, SA, P, SAi, *this);

    sjdbN=0;
    if (P.sjdbInsert.yes) {//insert junctions
        P.sjdbInsert.outDir=pGe.gDir;
        P.twoPass.pass2=false;

        Genome genome1(*this); //create copy here, *this will be changed by sjdbInsertJunctions
        sjdbInsertJunctions(P, *this, genome1, sjdbLoci);
    };

    pGe.gFileSizes.clear();
    pGe.gFileSizes.push_back(nGenome);
    pGe.gFileSizes.push_back(SA.lengthByte);

    //write genome parameters file
    genomeParametersWrite(pGe.gDir+("/genomeParameters.txt"), P, ERROR_OUT, *this);

    //write genome to disk
    time ( &rawTime );
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ... writing Genome to disk ...\n" <<flush;
    *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... writing Genome to disk ...\n" <<flush;

    writeGenomeSequence(pGe.gDir);

    //write SA
    time ( &rawTime );
    P.inOut->logMain  << "SA size in bytes: "<<SA.lengthByte << "\n"<<flush;

    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ... writing Suffix Array to disk ...\n" <<flush;
    *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... writing Suffix Array to disk ...\n" <<flush;

    ofstream & SAout = ofstrOpen(pGe.gDir+"/SA",ERROR_OUT, P);
    fstreamWriteBig(SAout,(char*) SA.charArray, (streamsize) SA.lengthByte,pGe.gDir+"/SA",ERROR_OUT,P);
    SAout.close();

    //write SAi
    time(&rawTime);
    P.inOut->logMain    << timeMonthDayTime(rawTime) <<" ... writing SAindex to disk\n" <<flush;
    *P.inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... writing SAindex to disk\n" <<flush;

    //write SAi to disk
    ofstream & SAiOut = ofstrOpen(pGe.gDir+"/SAindex",ERROR_OUT, P);

    fstreamWriteBig(SAiOut, (char*) &pGe.gSAindexNbases, sizeof(pGe.gSAindexNbases),pGe.gDir+"/SAindex",ERROR_OUT,P);
    fstreamWriteBig(SAiOut, (char*) genomeSAindexStart, sizeof(genomeSAindexStart[0])*(pGe.gSAindexNbases+1),pGe.gDir+"/SAindex",ERROR_OUT,P);
    fstreamWriteBig(SAiOut,  SAi.charArray, SAi.lengthByte,pGe.gDir+"/SAindex",ERROR_OUT,P);
    SAiOut.close();

    SApass1.deallocateArray();

    time(&rawTime);
    timeString=asctime(localtime ( &rawTime ));
    timeString.erase(timeString.end()-1,timeString.end());

    time(&rawTime);
    P.inOut->logMain    << timeMonthDayTime(rawTime) << " ..... finished successfully\n" <<flush;
    *P.inOut->logStdOut << timeMonthDayTime(rawTime) << " ..... finished successfully\n" <<flush;
};

void Genome::writeChrInfo(const string dirOut) 
{//write chr information
    ofstream & chrN = ofstrOpen(dirOut+"/chrName.txt",ERROR_OUT, P);
    ofstream & chrS = ofstrOpen(dirOut+"/chrStart.txt",ERROR_OUT, P);
    ofstream & chrL = ofstrOpen(dirOut+"/chrLength.txt",ERROR_OUT, P);
    ofstream & chrNL = ofstrOpen(dirOut+"/chrNameLength.txt",ERROR_OUT, P);

    for (uint ii=0;ii<nChrReal;ii++) {//output names, starts, lengths
        chrN<<chrName[ii]<<"\n";
        chrS<<chrStart[ii]<<"\n";
        chrL<<chrLength.at(ii)<<"\n";
        chrNL<<chrName[ii]<<"\t"<<chrLength.at(ii)<<"\n";
    };
    chrS<<chrStart[nChrReal]<<"\n";//size of the genome
    chrN.close();chrL.close();chrS.close(); chrNL.close();
};
void Genome::writeGenomeSequence(const string dirOut) 
{//write genome sequence
    ofstream &genomeOut = ofstrOpen(dirOut+"/Genome",ERROR_OUT, P);
    fstreamWriteBig(genomeOut,G,nGenome,dirOut+"/Genome",ERROR_OUT,P);
    genomeOut.close();
};
