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

#include "serviceFuns.cpp"
#include "streamFuns.h"

char* globalG;
uint globalL;


inline int funCompareSuffixes ( const void *a, const void *b){
	uint jj=0LLU;
        
	uint *ga=(uint*)((globalG-7LLU)+(*((uint*)a)));
	uint *gb=(uint*)((globalG-7LLU)+(*((uint*)b)));
    uint va=0,vb=0;

	while (va==vb && jj<globalL) {
		va=*(ga-jj);
		vb=*(gb-jj);
		jj++;
	};

	if (va>vb) {
		return 1;
	} else if (va==vb) {
		return 0;
	} else {
		return -1;
	};
};

inline bool funCompareSuffixesBool ( const void *a, const void *b) 
{
	uint jj=0LLU;
        
	uint *ga=(uint*)((globalG-7LLU)+(*((uint*)a)));
	uint *gb=(uint*)((globalG-7LLU)+(*((uint*)b)));
    uint va=0,vb=0;

	while (va==vb && jj<globalL) {
		va=*(ga-jj);
		vb=*(gb-jj);
		jj++;
	};

	if (va<vb) {
        return true;
	} else {
		return false;
	};
};


inline uint funG2strLocus (uint SAstr, uint const N, char const GstrandBit, uint const GstrandMask) {
    bool strandG = (SAstr>>GstrandBit) == 0;
    SAstr &= GstrandMask;
    if ( !strandG ) SAstr += N;
    return SAstr;
};

uint genomeScanFastaFiles (Parameters *P, char* G, bool flagRun) {//scans fasta files. flagRun=false: check and find full size, flaRun=true: collect all the data
    uint N=0; //total number of bases in the genome, including chr "spacers"
    ifstream fileIn;
    for (uint ii=0;ii<P->genomeFastaFiles.size();ii++) {//all the input files
        fileIn.open(P->genomeFastaFiles.at(ii).c_str());
        if ( !fileIn.good() ) {//
            ostringstream errOut;
            errOut << "EXITING because of INPUT ERROR: could not open genomeFastaFile: " <<P->genomeFastaFiles.at(ii) <<"\n";
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);            
        };
        while(!fileIn.eof()) {//read each file until eof
            string lineIn (4096,'.');
            getline(fileIn,lineIn);
            if (lineIn[0]=='>') {//new chromosome
                if (!flagRun) {
                    istringstream lineInStream (lineIn);
                    lineInStream.ignore(1,' ');
                    string chrName1;
                    lineInStream >> chrName1;
                    P->chrName.push_back(chrName1);
                };
                
                if (!flagRun && P->chrStart.size()>0) P->chrLength.push_back(N-P->chrStart.at(P->chrStart.size()-1)); //true length of the chr  
                
                if (N>0) {//pad the chromosomes to bins boudnaries
                    N = ( (N+1)/P->genomeChrBinNbases+1 )*P->genomeChrBinNbases;
                };

                if (!flagRun) {
                    P->chrStart.push_back(N);    
                    P->inOut->logMain << P->genomeFastaFiles.at(ii)<<" : chr # " << P->chrStart.size()-1 << "  \""<<P->chrName.at(P->chrStart.size()-1)<<"\" chrStart: "<<N<<"\n"<<flush;
                };
            } else {//char lines
                if (flagRun) lineIn.copy(G+N,lineIn.size(),0);
                N += lineIn.size();
            };
        };
        fileIn.close();        
    };
    
   
    if (!flagRun) P->chrLength.push_back(N-P->chrStart.at(P->chrStart.size()-1)); //true length of the chr  

    N = ( (N+1)/P->genomeChrBinNbases+1)*P->genomeChrBinNbases;
        
    if (!flagRun) { 
        P->nChrReal=P->chrStart.size();
        P->chrStart.push_back(N); //last chromosome end
        for (uint ii=0;ii<P->nChrReal;ii++) {
            P->chrNameIndex[P->chrName[ii]]=ii;
        };
    };
    
    return N;
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
    genomeParametersWrite(P->genomeDir+("/genomeParameters.txt"), P, "ERROR_00102");
    
    char *G=NULL, *G1=NULL;        
    uint nGenomeReal=genomeScanFastaFiles(P,G,false);//first scan the fasta file to find all the sizes  
    P->chrBinFill();

    uint L=10000;//maximum length of genome suffix    
    uint nG1alloc=(nGenomeReal + L)*2;
    G1=new char[nG1alloc];
    G=G1+L;
    
    memset(G1,GENOME_spacingChar,nG1alloc);//initialize to K-1 all bytes
 
    genomeScanFastaFiles(P,G,true);    //load the genome sequence
     
    //convert the genome to 0,1,2,3,4
    for (uint jj=0;jj<nGenomeReal;jj++) {
        switch (int(G[jj])){
            case(65): case(97):  G[jj]=char(0);break;//A
            case(67): case(99):  G[jj]=char(1);break;//C           
            case(71): case(103): G[jj]=char(2);break;//G                       
            case(84): case(116): G[jj]=char(3);break;//T                                
            case(78): case(110): G[jj]=char(4);break;//N
            case(48):            G[jj]=GENOME_spacingChar;break;//chromosomal breaks within the sequences
            default:              //anything else
                if (G[jj]!=GENOME_spacingChar) {
//                     P->inOut->logMain << "Unexpected character: char="<< G[jj] << "   int="<<int(G[jj])<<"   at " << jj << " , replacing with N\n";
                     G[jj]=char(4);                                 
                };
        };
    };    

    uint N = nGenomeReal;
    P->nGenome=N;
    uint N2 = N*2;     

    ofstream chrN,chrS,chrL,chrNL;
    
    ofstrOpen(P->genomeDir+"/chrName.txt","ERROR_00103", P, chrN);   
    ofstrOpen(P->genomeDir+"/chrStart.txt","ERROR_00103", P, chrS);   
    ofstrOpen(P->genomeDir+"/chrLength.txt","ERROR_00103", P, chrL);   
    ofstrOpen(P->genomeDir+"/chrNameLength.txt","ERROR_00103", P, chrNL);   
    
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
        
    P->nSAbyte=SA2.lengthByte;
    
    P->inOut->logMain  << "Number of SA indices: "<< P->nSA << "\n"<<flush;    

    //sort SA
    time ( &rawTime );
    P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... starting to sort  Suffix Array. This may take a long time...\n" <<flush;   
    *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... starting to sort  Suffix Array. This may take a long time...\n" <<flush;

    for (uint ii=0;ii<N;ii++) {//re-fill the array backwards for sorting
        swap(G[N2-1-ii],G[ii]);
    };          
    globalG=G;
    globalL=L/sizeof(uint);        

    {//not enough RAM, split into chunks          
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
        for (int iChunk=0; (uint)iChunk < saChunkN; iChunk++) {//start the chunk cycle: sort each chunk with qsort and write to a file
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
            ofstream saChunkFile;
            string chunkFileName=P->genomeDir+"/SA_"+to_string( (uint) iChunk);
            ofstrOpen(chunkFileName,"ERROR_00105", P, saChunkFile);   
            fstreamWriteBig(saChunkFile, (char*) saChunk, sizeof(saChunk[0])*indPrefChunkCount[iChunk],chunkFileName,"ERROR_00121",P);
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
    time(&rawTime);    
    timeString=asctime(localtime ( &rawTime ));
    timeString.erase(timeString.end()-1,timeString.end());
    P->inOut->logMain    << timeMonthDayTime(rawTime) <<" ... starting to generate Suffix Array index...\n" <<flush;   
    *P->inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... starting to generate Suffix Array index...\n" <<flush; 
        
    P->genomeSAindexStart = new uint [P->genomeSAindexNbases+1];
    P->genomeSAindexStart[0]=0;
    for (uint ii=1;ii<=P->genomeSAindexNbases;ii++) {//L-mer indices starts
        P->genomeSAindexStart[ii] = P->genomeSAindexStart[ii-1] + ( 1LLU<<(2*ii) );
    };
    P->nSAi = P->genomeSAindexStart[P->genomeSAindexNbases];
    
    uint* SAi=new uint[P->nSAi];
//     for (uint isa=0; isa<P->nSAi; isa++) {//initialize
//         SAi[isa]=P->nSA; //if the suffix is not found in the genome, it's location will be marked with this value
//     };
    
    uint* ind0=new uint[P->genomeSAindexNbases];
    uint* indSAlast=new uint[P->genomeSAindexNbases];

    for (uint ii=0; ii<P->genomeSAindexNbases; ii++) {
        ind0[ii]=-1;//this is needed in case "AAA...AAA",i.e. indPref=0 is not present in the genome for some lengths
        indSAlast[ii]=P->nSA;//that's probably not needed
    };

    P->SAiMarkNbit=P->GstrandBit+1;
    P->SAiMarkAbsentBit=P->GstrandBit+2;
    
    P->SAiMarkNmaskC=1LLU << P->SAiMarkNbit;
    P->SAiMarkNmask=~P->SAiMarkNmaskC;
    P->SAiMarkAbsentMaskC=1LLU << P->SAiMarkAbsentBit;
    P->SAiMarkAbsentMask=~P->SAiMarkAbsentMaskC;
       
    
    for (uint isa=0; isa<P->nSA; isa++) {//for all suffixes
        if (isa%100000000==0) P->inOut->logMain  << isa*100/P->nSA << "% " << flush;         
        
        uint SAstr=SA1[isa];
        bool dirG = (SAstr>>P->GstrandBit) == 0; //forward or reverse strand of the genome
        SAstr &= P->GstrandMask;
        if (!dirG) SAstr=P->nGenome-1-SAstr;

        uint indPref=0;
        for (uint iL=0; iL < P->genomeSAindexNbases; iL++) {//calculate index

            indPref <<= 2;
            
            uint g1= (uint) G[dirG ? SAstr+iL : SAstr-iL]; //reverese if (-) strand

            if (g1>3) {//if N, this suffix does not belong in SAi
                for (uint iL1=iL; iL1 < P->genomeSAindexNbases; iL1++) {
                    SAi[P->genomeSAindexStart[iL1]+ind0[iL1]] |= P->SAiMarkNmaskC;
                };
                break;
            };

            if (!dirG) g1=3-g1; //complement if (-) strand

            indPref += (uint) g1;
            
            if ( indPref > ind0[iL] || isa==0 ) {//new && good index, record it
                SAi[P->genomeSAindexStart[iL]+indPref]=isa;
                for (uint ii=ind0[iL]+1; ii<indPref; ii++) {//index is not present, record to the last present suffix
                    SAi[P->genomeSAindexStart[iL]+ii] = isa | P->SAiMarkAbsentMaskC; 
                };
                ind0[iL]=indPref;
//                 indSAlast[iL]=isa;
//             } else if (indPref==ind0[iL]) {
//                 indSAlast[iL]=isa;//last SA index with the same prefix 
            } else if ( indPref < ind0[iL] ) {
                ostringstream errOut;
                errOut << "BUG: next index is smaller than previous, EXITING\n" <<flush;
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
            };
        };
    };//for (uint isa=0; isa<P->nSA; isa++)
    P->inOut->logMain << " done\n"<<flush;
   
    //pack SAi
    PackedArray SAip;
    SAip.defineBits(P->GstrandBit+3,P->nSAi);//SAi uses an extra bit compared to SA because it needs to store values > nSA
    SAip.pointArray((char*) SAi);
    for (uint ii=0;ii<SAip.length;ii++) {
        SAip.writePacked(ii,SAi[ii]);
    };
    
    delete [] indSAlast;
    delete [] ind0;

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
        SA2.writePacked(P->nSA,0);
    };
    
    //write genome to disk
    time ( &rawTime );
    P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... writing Genome to disk ...\n" <<flush;   
    *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... writing Genome to disk ...\n" <<flush;   
    
    ofstream genomeOut;
    ofstrOpen(P->genomeDir+"/Genome","ERROR_00104", P, genomeOut);   
    fstreamWriteBig(genomeOut,G,P->nGenome,P->genomeDir+"/Genome","ERROR_00120",P);
    genomeOut.close();  

    //write SA                
    time ( &rawTime );
    P->inOut->logMain  << "SA size in bytes: "<< P->nSAbyte << "\n"<<flush;

    P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... writing Suffix Array to disk ...\n" <<flush;   
    *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... writing Suffix Array to disk ...\n" <<flush;   

    ofstream SAout;
    ofstrOpen(P->genomeDir+"/SA","ERROR_00106", P, SAout);   
    fstreamWriteBig(SAout,(char*) SA2.charArray, (streamsize) P->nSAbyte,P->genomeDir+"/SA","ERROR_00122",P);
    SAout.close();    
    
    //write SAi
    time(&rawTime);    
    P->inOut->logMain    << timeMonthDayTime(rawTime) <<" ... writing SAindex to disk\n" <<flush;   
    *P->inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... writing SAindex to disk\n" <<flush;   
    
    //write SAi to disk
    ofstream SAiOut;
    ofstrOpen(P->genomeDir+"/SAindex","ERROR_00107", P, SAiOut);   

    fstreamWriteBig(SAiOut, (char*) &P->genomeSAindexNbases, sizeof(P->genomeSAindexNbases),P->genomeDir+"/SAindex","ERROR_00123",P);
    fstreamWriteBig(SAiOut, (char*) P->genomeSAindexStart, sizeof(P->genomeSAindexStart[0])*(P->genomeSAindexNbases+1),P->genomeDir+"/SAindex","ERROR_00124",P);        
    fstreamWriteBig(SAiOut,  SAip.charArray, SAip.lengthByte,P->genomeDir+"/SAindex","ERROR_00125",P);
    SAiOut.close();    

    SA2.deallocateArray();
    delete [] SAi;    

    time(&rawTime);
    timeString=asctime(localtime ( &rawTime ));
    timeString.erase(timeString.end()-1,timeString.end());
    
    time(&rawTime);        
    P->inOut->logMain    << timeMonthDayTime(rawTime) << " ..... Finished successfully\n" <<flush;    
    *P->inOut->logStdOut << timeMonthDayTime(rawTime) << " ..... Finished successfully\n" <<flush;
};
