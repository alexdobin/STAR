#include <cmath>

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


inline int funCompareSuffixes ( const void *a, const void *b){

    uint *ga=(uint*)((globalG-7LLU)+(*((uint*)a)));
    uint *gb=(uint*)((globalG-7LLU)+(*((uint*)b)));

    uint jj=0;
    int  ii=0;
    uint va=0,vb=0;
    uint8 *va1, *vb1;

    while (jj < globalL) {
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

inline uint funG2strLocus (uint SAstr, uint const N, char const GstrandBit, uint const GstrandMask) {
    bool strandG = (SAstr>>GstrandBit) == 0;
    SAstr &= GstrandMask;
    if ( !strandG ) SAstr += N;
    return SAstr;
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
        //count the number of indices with 4nt prefix
        uint indPrefN=1LLU << 16;
        uint* indPrefCount = new uint [indPrefN];
        memset(indPrefCount,0,indPrefN*sizeof(indPrefCount[0]));
        nSA=0;
        for (uint ii=0;ii<2*nGenome;ii+=pGe.gSAsparseD) {
            if (G[ii]<4) {
                uint p1=(G[ii]<<12) + (G[ii-1]<<8) + (G[ii-2]<<4) + G[ii-3];
                indPrefCount[p1]++;
                nSA++;
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

        P.inOut->logMain  << "Number of chunks: " << saChunkN <<";   chunks size limit: " << saChunkSize*8 <<" bytes\n" <<flush;

        time ( &rawTime );
        P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ... sorting Suffix Array chunks and saving them to disk...\n" <<flush;
        *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... sorting Suffix Array chunks and saving them to disk...\n" <<flush;

        #pragma omp parallel for num_threads(P.runThreadN) ordered schedule(dynamic,1)
        for (int iChunk=0; iChunk < (int) saChunkN; iChunk++) {//start the chunk cycle: sort each chunk with qsort and write to a file
            uint* saChunk=new uint [indPrefChunkCount[iChunk]];//allocate local array for each chunk
            for (uint ii=0,jj=0;ii<2*nGenome;ii+=pGe.gSAsparseD) {//fill the chunk with SA indices
                if (G[ii]<4) {
                    uint p1=(G[ii]<<12) + (G[ii-1]<<8) + (G[ii-2]<<4) + G[ii-3];
                    if (p1>=indPrefStart[iChunk] && p1<indPrefStart[iChunk+1]) {
                        saChunk[jj]=ii;
                        jj++;
                    };
                    //TODO: if (jj==indPrefChunkCount[iChunk]) break;
                };
            };


            //sort the chunk
            qsort(saChunk,indPrefChunkCount[iChunk],sizeof(saChunk[0]),funCompareSuffixes);
            for (uint ii=0;ii<indPrefChunkCount[iChunk];ii++) {
                saChunk[ii]=2*nGenome-1-saChunk[ii];
            };
            //write files
            string chunkFileName=pGe.gDir+"/SA_"+to_string( (uint) iChunk);
            ofstream & saChunkFile = ofstrOpen(chunkFileName,ERROR_OUT, P);
            fstreamWriteBig(saChunkFile, (char*) saChunk, sizeof(saChunk[0])*indPrefChunkCount[iChunk],chunkFileName,ERROR_OUT,P);
            saChunkFile.close();
            delete [] saChunk;
            saChunk=NULL;
        };

        time ( &rawTime );
        P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ... loading chunks from disk, packing SA...\n" <<flush;
        *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... loading chunks from disk, packing SA...\n" <<flush;

        //read chunks and pack into full SA
        SApass1.allocateArray();
        SA.pointArray(SApass1.charArray + SApass1.lengthByte-SA.lengthByte); //SA is shifted to have space for junction insertion
        uint N2bit= 1LLU << GstrandBit;
        uint packedInd=0;

        #define SA_CHUNK_BLOCK_SIZE 10000000
        uint* saIn=new uint[SA_CHUNK_BLOCK_SIZE]; //TODO make adjustable

        #ifdef genenomeGenerate_SA_textOutput
                ofstream SAtxtStream ((pGe.gDir + "/SAtxt").c_str());
        #endif

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

        #ifdef genenomeGenerate_SA_textOutput
                SAtxtStream.close();
        #endif
        delete [] saIn;

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
