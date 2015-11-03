#include "IncludeDefine.h"
#include "Parameters.h"
#include "SuffixArrayFuns.h"
#include "PackedArray.h"
#include <math.h>
#include "TimeFunctions.h"
#include "ErrorWarning.h"
#include "loadGTF.h"
#include "SjdbClass.h"
#include "sjdbLoadFromFiles.h"
#include "sjdbPrepare.h"
#include "genomeParametersWrite.h"
#include "sjdbInsertJunctions.h"
#include "genomeScanFastaFiles.h"
#include "genomeSAindex.h"

#include "serviceFuns.cpp"
#include "streamFuns.h"

char* globalG;
uint globalL;


inline int funCompareSuffixes ( const void *a, const void *b){
        
    uint *ga=(uint*)((globalG-7LLU)+(*((uint*)a)));
    uint *gb=(uint*)((globalG-7LLU)+(*((uint*)b)));

    uint jj=0;
    int  ii=0;
    uint va=0,vb=0;
    uint8 *va1, *vb1;

    bool aeqb=true;
    while (aeqb) {
        va=*(ga-jj);
        vb=*(gb-jj);
        va1=(uint8*) &va;
        vb1=(uint8*) &vb;
        for (ii=7;ii>=0;ii--)
        {
            if (va1[ii]!=vb1[ii] || va1[ii]==5)
            {
                aeqb=false;
                break;
            };
        };
        jj++;
    };

    if (va1[ii]>vb1[ii]) 
    {
        return 1;
    } else if (va1[ii]<vb1[ii]) 
    {
        return -1;
    } else 
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

// inline bool funCompareSuffixesBool ( const void *a, const void *b) 
// {
// 	uint jj=0LLU;
//         
// 	uint *ga=(uint*)((globalG-7LLU)+(*((uint*)a)));
// 	uint *gb=(uint*)((globalG-7LLU)+(*((uint*)b)));
//     uint va=0,vb=0;
// 
// 	while (va==vb && jj<globalL) {
// 		va=*(ga-jj);
// 		vb=*(gb-jj);
// 		jj++;
// 	};
// 
// 	if (va<vb) {
//         return true;
// 	} else {
// 		return false;
// 	};
// };


inline uint funG2strLocus (uint SAstr, uint const N, char const GstrandBit, uint const GstrandMask) {
    bool strandG = (SAstr>>GstrandBit) == 0;
    SAstr &= GstrandMask;
    if ( !strandG ) SAstr += N;
    return SAstr;
};

void genomeGenerate(Parameters *P) {
    
    //check parameters
    if (P->sjdbOverhang<=0 && (P->sjdbFileChrStartEnd.at(0)!="-" || P->sjdbGTFfile!="-")) 
    {
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT PARAMETER ERROR: for generating genome with annotations (--sjdbFileChrStartEnd or --sjdbGTFfile options)\n";
        errOut << "you need to specify >0 --sjdbOverhang\n";
        errOut << "SOLUTION: re-run genome generation specifying non-zero --sjdbOverhang, which ideally should be equal to OneMateLength-1, or could be chosen generically as ~100\n";        
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
    } 
    if (P->sjdbFileChrStartEnd.at(0)=="-" && P->sjdbGTFfile=="-") 
    {
        if (P->parArray.at(P->sjdbOverhang_par)->inputLevel>0 && P->sjdbOverhang>0)
        {
            ostringstream errOut;
            errOut << "EXITING because of FATAL INPUT PARAMETER ERROR: when generating genome without annotations (--sjdbFileChrStartEnd or --sjdbGTFfile options)\n";
            errOut << "do not specify >0 --sjdbOverhang\n";
            errOut << "SOLUTION: re-run genome generation without --sjdbOverhang option\n";        
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
        };
        P->sjdbOverhang=0;
    };
    
    //time
    time_t rawTime;
    string timeString;
    
    time(&rawTime);
    P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... Starting to generate Genome files\n" <<flush;
    *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... Starting to generate Genome files\n" <<flush;
    
    //define some parameters from input parameters
    P->genomeChrBinNbases=1LLU << P->genomeChrBinNbits;
    //write genome parameters file
    genomeParametersWrite(P->genomeDir+("/genomeParameters.txt"), P, ERROR_OUT);
    
    char *G=NULL, *G1=NULL;        
    uint nGenomeReal=genomeScanFastaFiles(P,G,false);//first scan the fasta file to find all the sizes  
    P->chrBinFill();

    uint L=10000;//maximum length of genome suffix    
    uint nG1alloc=(nGenomeReal + L)*2;
    G1=new char[nG1alloc];
    G=G1+L;
    
    memset(G1,GENOME_spacingChar,nG1alloc);//initialize to K-1 all bytes
 
    genomeScanFastaFiles(P,G,true);    //load the genome sequence   

    uint N = nGenomeReal;
    P->nGenome=N;
    uint N2 = N*2;     

    ofstream & chrN = ofstrOpen(P->genomeDir+"/chrName.txt",ERROR_OUT, P);
    ofstream & chrS = ofstrOpen(P->genomeDir+"/chrStart.txt",ERROR_OUT, P);
    ofstream & chrL = ofstrOpen(P->genomeDir+"/chrLength.txt",ERROR_OUT, P);
    ofstream & chrNL = ofstrOpen(P->genomeDir+"/chrNameLength.txt",ERROR_OUT, P);
    
    for (uint ii=0;ii<P->nChrReal;ii++) {//output names, starts, lengths
        chrN<<P->chrName[ii]<<"\n";
        chrS<<P->chrStart[ii]<<"\n";
        chrL<<P->chrLength.at(ii)<<"\n";
        chrNL<<P->chrName[ii]<<"\t"<<P->chrLength.at(ii)<<"\n";        
    };
    chrS<<P->chrStart[P->nChrReal]<<"\n";//size of the genome
    chrN.close();chrL.close();chrS.close(); chrNL.close();   
    
    if (P->limitGenomeGenerateRAM < (nG1alloc+nG1alloc/3)) {//allocate nG1alloc/3 for SA generation
        ostringstream errOut;                            
        errOut <<"EXITING because of FATAL PARAMETER ERROR: limitGenomeGenerateRAM="<< (P->limitGenomeGenerateRAM) <<"is too small for your genome\n";
        errOut <<"SOLUTION: please specify limitGenomeGenerateRAM not less than"<< nG1alloc+nG1alloc/3 <<" and make that much RAM available \n";
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
    };     
    
    //preparing to generate SA
    for (uint ii=0;ii<N;ii++) {//- strand
        G[N2-1-ii]=G[ii]<4 ? 3-G[ii] : G[ii];
    };      
    
    P->nSA=0;
    for (uint ii=0;ii<N2;ii+=P->genomeSAsparseD) {
        if (G[ii]<4) {
            P->nSA++;
        };
    };     
    
    P->GstrandBit = (uint) floor(log(N)/log(2))+1; 
    if (P->GstrandBit<32) P->GstrandBit=32; //TODO: use simple access function for SA
    
    P->GstrandMask = ~(1LLU<<P->GstrandBit);
    PackedArray SA1;//SA without sjdb
    SA1.defineBits(P->GstrandBit+1,P->nSA);
    PackedArray SA2;//SA with sjdb, reserve more space
    if (P->sjdbInsert.yes)
    {//reserve space for junction insertion
        SA2.defineBits(P->GstrandBit+1,P->nSA+2*P->limitSjdbInsertNsj*P->sjdbLength);//TODO: this allocation is wasteful, get a better estimate of the number of junctions
    } else
    {//same as SA1
        SA2.defineBits(P->GstrandBit+1,P->nSA);
    };
        
    P->inOut->logMain  << "Number of SA indices: "<< P->nSA << "\n"<<flush;    

    //sort SA
    time ( &rawTime );
    P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... starting to sort  Suffix Array. This may take a long time...\n" <<flush;   
    *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... starting to sort  Suffix Array. This may take a long time...\n" <<flush;
   

//     if (false)
    {//sort SA chunks
        
        for (uint ii=0;ii<N;ii++) {//re-fill the array backwards for sorting
            swap(G[N2-1-ii],G[ii]);
        };          
        globalG=G;
        globalL=L/sizeof(uint);
        //count the number of indices with 4nt prefix
        uint indPrefN=1LLU << 16;
        uint* indPrefCount = new uint [indPrefN];
        memset(indPrefCount,0,indPrefN*sizeof(indPrefCount[0]));
        P->nSA=0;
        for (uint ii=0;ii<N2;ii+=P->genomeSAsparseD) {
            if (G[ii]<4) {
                uint p1=(G[ii]<<12) + (G[ii-1]<<8) + (G[ii-2]<<4) + G[ii-3];
                indPrefCount[p1]++;
                P->nSA++;
            };
        };

        uint saChunkSize=(P->limitGenomeGenerateRAM-nG1alloc)/8/P->runThreadN; //number of SA indexes per chunk
        saChunkSize=saChunkSize*6/10; //allow extra space for qsort            
        //uint saChunkN=((P->nSA/saChunkSize+1)/P->runThreadN+1)*P->runThreadN;//ensure saChunkN is divisible by P->runThreadN
        //saChunkSize=P->nSA/saChunkN+100000;//final chunk size
        if (P->runThreadN>1) saChunkSize=min(saChunkSize,P->nSA/(P->runThreadN-1));

        uint saChunkN=P->nSA/saChunkSize;//estimate
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

        P->inOut->logMain  << "Number of chunks: " << saChunkN <<";   chunks size limit: " << saChunkSize*8 <<" bytes\n" <<flush;

        time ( &rawTime );
        P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... sorting Suffix Array chunks and saving them to disk...\n" <<flush;   
        *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... sorting Suffix Array chunks and saving them to disk...\n" <<flush;

        #pragma omp parallel for num_threads(P->runThreadN) ordered schedule(dynamic,1)
        for (int iChunk=0; iChunk < (int) saChunkN; iChunk++) {//start the chunk cycle: sort each chunk with qsort and write to a file
            uint* saChunk=new uint [indPrefChunkCount[iChunk]];//allocate local array for each chunk
            for (uint ii=0,jj=0;ii<N2;ii+=P->genomeSAsparseD) {//fill the chunk with SA indices
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
                saChunk[ii]=N2-1-saChunk[ii];
            };  
            //write files
            string chunkFileName=P->genomeDir+"/SA_"+to_string( (uint) iChunk);
            ofstream & saChunkFile = ofstrOpen(chunkFileName,ERROR_OUT, P);   
            fstreamWriteBig(saChunkFile, (char*) saChunk, sizeof(saChunk[0])*indPrefChunkCount[iChunk],chunkFileName,ERROR_OUT,P);
            saChunkFile.close();
            delete [] saChunk;
            saChunk=NULL;
        };

        time ( &rawTime );
        P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... loading chunks from disk, packing SA...\n" <<flush;   
        *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... loading chunks from disk, packing SA...\n" <<flush;    

        //read chunks and pack into full SA1
        SA2.allocateArray();
        SA1.pointArray(SA2.charArray + SA2.lengthByte-SA1.lengthByte); //SA1 is shifted to have space for junction insertion
        uint N2bit= 1LLU << P->GstrandBit;          
        uint packedInd=0;

        #define SA_CHUNK_BLOCK_SIZE 10000000
        uint* saIn=new uint[SA_CHUNK_BLOCK_SIZE]; //TODO make adjustable
        
        #ifdef genenomeGenerate_SA_textOutput
                ofstream SAtxtStream ((P->genomeDir + "/SAtxt").c_str());
        #endif

        for (uint iChunk=0;iChunk<saChunkN;iChunk++) {//load files one by one and convert to packed
            ostringstream saChunkFileNameStream("");
            saChunkFileNameStream<< P->genomeDir << "/SA_" << iChunk;
            ifstream saChunkFile(saChunkFileNameStream.str().c_str());
            while (! saChunkFile.eof()) {//read blocks from each file
                uint chunkBytesN=fstreamReadBig(saChunkFile,(char*) saIn,SA_CHUNK_BLOCK_SIZE*sizeof(saIn[0]));
                for (uint ii=0;ii<chunkBytesN/sizeof(saIn[0]);ii++) {
                    SA1.writePacked( packedInd+ii, (saIn[ii]<N) ? saIn[ii] : ( (saIn[ii]-N) | N2bit ) );
                    
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

        if (packedInd != P->nSA ) {//
            ostringstream errOut;                            
            errOut << "EXITING because of FATAL problem while generating the suffix array\n";
            errOut << "The number of indices read from chunks = "<<packedInd<<" is not equal to expected nSA="<<P->nSA<<"\n";
            errOut << "SOLUTION: try to re-run suffix array generation, if it still does not work, report this problem to the author\n"<<flush;
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
        };
        
        //DONE with suffix array generation
        
        for (uint ii=0;ii<N;ii++) {//return to normal order for future use
            swap(G[N2-1-ii],G[ii]);
        };         
        delete [] indPrefCount;
        delete [] indPrefStart;
        delete [] indPrefChunkCount;
    };    

    time ( &rawTime );
    timeString=asctime(localtime ( &rawTime ));
    timeString.erase(timeString.end()-1,timeString.end());
    P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... Finished generating suffix array\n" <<flush;  
    *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... Finished generating suffix array\n" <<flush;          

////////////////////////////////////////
//          SA index
//
//     PackedArray SAold;
// 
//     if (true)
//     {//testing: load SA from disk
//             //read chunks and pack into full SA1
//         
//         ifstream oldSAin("./DirTrue/SA");
//         oldSAin.seekg (0, ios::end);
//         P->nSAbyte=(uint) oldSAin.tellg();
//         oldSAin.clear();        
//         oldSAin.seekg (0, ios::beg);
// 
//         P->nSA=(P->nSAbyte*8)/(P->GstrandBit+1);
//         SAold.defineBits(P->GstrandBit+1,P->nSA);  
//         SAold.allocateArray();
//         
//         oldSAin.read(SAold.charArray,SAold.lengthByte);
//         oldSAin.close();
//         
//         SA1=SAold;
//         SA2=SAold;
//     };
    
    PackedArray SAip;
    genomeSAindex(G,SA1,P,SAip);

    if (P->sjdbFileChrStartEnd.at(0)!="-" || P->sjdbGTFfile!="-")
    {//insert junctions
        SjdbClass sjdbLoci;

        Genome mainGenome(P);
        mainGenome.G=G;
        mainGenome.SA=SA1;
        mainGenome.SApass1=SA2;
        mainGenome.SAi=SAip;
        P->sjdbInsert.outDir=P->genomeDir;
        P->sjdbN=0;//no junctions are loaded yet
        P->twoPass.pass2=false;
        
        Parameters *P1=new Parameters;
        *P1=*P;        
        
        sjdbInsertJunctions(P, P1, mainGenome, sjdbLoci);
        
        //write an extra 0 at the end of the array, filling the last bytes that otherwise are not accessible, but will be written to disk
        //this is - to avoid valgrind complaints. Note that SA2 is allocated with plenty of space to spare.
        SA1=mainGenome.SA;
        SA1.writePacked(P->nSA,0);
    };
    
    //write genome to disk
    time ( &rawTime );
    P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... writing Genome to disk ...\n" <<flush;   
    *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... writing Genome to disk ...\n" <<flush;   
    
    ofstream & genomeOut = ofstrOpen(P->genomeDir+"/Genome",ERROR_OUT, P);   
    fstreamWriteBig(genomeOut,G,P->nGenome,P->genomeDir+"/Genome",ERROR_OUT,P);
    genomeOut.close();  

    //write SA                
    time ( &rawTime );
    P->inOut->logMain  << "SA size in bytes: "<<SA1.lengthByte << "\n"<<flush;

    P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... writing Suffix Array to disk ...\n" <<flush;   
    *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... writing Suffix Array to disk ...\n" <<flush;   

    ofstream & SAout = ofstrOpen(P->genomeDir+"/SA",ERROR_OUT, P);   
    fstreamWriteBig(SAout,(char*) SA1.charArray, (streamsize) SA1.lengthByte,P->genomeDir+"/SA",ERROR_OUT,P);
    SAout.close();    
    
    //write SAi
    time(&rawTime);    
    P->inOut->logMain    << timeMonthDayTime(rawTime) <<" ... writing SAindex to disk\n" <<flush;   
    *P->inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... writing SAindex to disk\n" <<flush;   
    
    //write SAi to disk
    ofstream & SAiOut = ofstrOpen(P->genomeDir+"/SAindex",ERROR_OUT, P);

    fstreamWriteBig(SAiOut, (char*) &P->genomeSAindexNbases, sizeof(P->genomeSAindexNbases),P->genomeDir+"/SAindex",ERROR_OUT,P);
    fstreamWriteBig(SAiOut, (char*) P->genomeSAindexStart, sizeof(P->genomeSAindexStart[0])*(P->genomeSAindexNbases+1),P->genomeDir+"/SAindex",ERROR_OUT,P);        
    fstreamWriteBig(SAiOut,  SAip.charArray, SAip.lengthByte,P->genomeDir+"/SAindex",ERROR_OUT,P);
    SAiOut.close();    

    SA2.deallocateArray();

    time(&rawTime);
    timeString=asctime(localtime ( &rawTime ));
    timeString.erase(timeString.end()-1,timeString.end());
    
    time(&rawTime);        
    P->inOut->logMain    << timeMonthDayTime(rawTime) << " ..... Finished successfully\n" <<flush;    
    *P->inOut->logStdOut << timeMonthDayTime(rawTime) << " ..... Finished successfully\n" <<flush;
};
