#include "ReadAlignChunk.h"
#include "GlobalVariables.h"
#include "ThreadControl.h"
#include "ErrorWarning.h"
#include "SequenceFuns.h"

void ReadAlignChunk::processChunks() {//read-map-write chunks
    noReadsLeft=false; //true if there no more reads left in the file
    bool newFile=false; //new file marker in the input stream
    while (!noReadsLeft) {//continue until the input EOF
            //////////////read a chunk from input files and store in memory
        if (P.outFilterBySJoutStage<2) {//read chunks from input file

            if (P.runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexInRead);

            uint chunkInSizeBytesTotal[2]={0,0};
            while (chunkInSizeBytesTotal[0] < P.chunkInSizeBytes && chunkInSizeBytesTotal[1] < P.chunkInSizeBytes && P.inOut->readIn[0].good() && P.inOut->readIn[1].good()) {
                char nextChar=P.inOut->readIn[0].peek();
                if (P.iReadAll==P.readMapNumber) {//do nto read any more reads
                    break;
                } else if (P.readFilesTypeN==10 && P.inOut->readIn[0].good() && P.outFilterBySJoutStage!=2) {//SAM input && not eof && not 2nd stage

                    string str1;
                    
                    if (nextChar=='@') {//with SAM input linest that start with @ are headers
                        getline(P.inOut->readIn[0], str1); //read line and skip it
                        continue;
                    };
                    
                    P.inOut->readIn[0] >> str1;
                    if (str1=="FILE") {
                        newFile=true;
                    } else {
                        P.iReadAll++; //increment read number

                        uint imate1=0;
                        for (uint imate=0;imate<P.readNmates;imate++) {
                            if (imate>0)
                                P.inOut->readIn[0] >> str1; //for imate=0 str1 was already read
                            uint flag;
                            P.inOut->readIn[0] >>flag; //read name and flag
                            char passFilterIllumina=(flag & 0x800 ? 'Y' : 'N');

                            if (imate==1) {//2nd line is always opposite of the 1st one
                                imate1=1-imate1;
                            } else if (P.readNmates==2 && (flag & 0x80)) {
                                imate1=1;
                            } else {
                                imate1=0;
                            };

                            //read ID or number
                            if (P.outSAMreadID=="Number") {
                                chunkInSizeBytesTotal[imate1] += sprintf(chunkIn[imate1] + chunkInSizeBytesTotal[imate1], "@%llu", P.iReadAll);
                            } else {
                                chunkInSizeBytesTotal[imate1] += sprintf(chunkIn[imate1] + chunkInSizeBytesTotal[imate1], "@%s", str1.c_str());
                            };

                            //iReadAll, passFilterIllumina, passFilterIllumina
                            chunkInSizeBytesTotal[imate1] += sprintf(chunkIn[imate1] + chunkInSizeBytesTotal[imate1], " %llu %c %i", P.iReadAll, passFilterIllumina, P.readFilesIndex);

                            for (int ii=3; ii<=9; ii++)
                                P.inOut->readIn[0] >> str1; //skip fields until sequence

                            string seq1,qual1;
                            P.inOut->readIn[0]  >> seq1 >> qual1;
                            if (flag & 0x10) {//sequence reverse-coomplemented
                                revComplementNucleotides(seq1);
                                reverse(qual1.begin(),qual1.end());
                            };
                            
                            getline(P.inOut->readIn[0],str1); //str1 is now all SAM attributes
                            chunkInSizeBytesTotal[imate1] += sprintf(chunkIn[imate1] + chunkInSizeBytesTotal[imate1], "%s\n%s\n+\n%s\n", str1.c_str(), seq1.c_str(), qual1.c_str());
                        };
                    };
                } else if (nextChar=='@') {//fastq, not multi-line
                    P.iReadAll++; //increment read number
                    for (uint imate=0; imate<P.readNmates; imate++) {//for all mates
                        uint32 iline=0;
                        if (P.outFilterBySJoutStage!=2) {//not the 2nd stage of the 2-stage mapping

                            //read or skip the 1st field of the read name line
                            if (P.outSAMreadID=="Number") {
                                chunkInSizeBytesTotal[imate] += sprintf(chunkIn[imate] + chunkInSizeBytesTotal[imate], "@%llu", P.iReadAll);
                                string dummy1;
                                P.inOut->readIn[imate] >> dummy1; //skip the first field of the read name
                            } else {
                                P.inOut->readIn[imate] >> (chunkIn[imate] + chunkInSizeBytesTotal[imate]);
                                chunkInSizeBytesTotal[imate] += strlen(chunkIn[imate] + chunkInSizeBytesTotal[imate]);
                            };

                            //read the second field of the read name line
                            char passFilterIllumina='N';
                            if (P.inOut->readIn[imate].peek()!='\n') {//2nd field exists
                                string field2;
                                P.inOut->readIn[imate] >> field2;
                                passFilterIllumina='N';
                                if (field2.length()>=3 && field2.at(2)=='Y') passFilterIllumina='Y';
                            };

                            //ignore the rest of the read name
                            P.inOut->readIn[imate].ignore(DEF_readNameSeqLengthMax,'\n');

                            chunkInSizeBytesTotal[imate] += sprintf(chunkIn[imate] + chunkInSizeBytesTotal[imate], " %llu %c %i \n", P.iReadAll, passFilterIllumina, P.readFilesIndex);

                            iline=1;
                        };
                        //do not need to worry about 2nd stage, that's read directly from the files
//                         else {//2nd stage of 2-stage mapping
//                         read index and file index are already recorded with the read name, simply copy it
//                         P.inOut->readIn[imate].getline(chunkIn[imate] + chunkInSizeBytesTotal[imate], DEF_readNameSeqLengthMax+1 );
//                         };
                        for (;iline<4;iline++) {//TODO ignore the 3rd line of fastq
                            P.inOut->readIn[imate].getline(chunkIn[imate] + chunkInSizeBytesTotal[imate], DEF_readNameSeqLengthMax+1 );
                            chunkInSizeBytesTotal[imate] += P.inOut->readIn[imate].gcount();
                            chunkIn[imate][chunkInSizeBytesTotal[imate]-1]='\n';
                        };
                    };
                } else if (nextChar=='>') {//fasta, can be multiline, which is converted to single line
                    P.iReadAll++; //increment read number
                    for (uint imate=0; imate<P.readNmates; imate++) {
                        if (P.outFilterBySJoutStage!=2) {//not the 2nd stage of the 2-stage mapping

                            if (P.outSAMreadID=="Number") {
                                chunkInSizeBytesTotal[imate] += sprintf(chunkIn[imate] + chunkInSizeBytesTotal[imate], ">%llu", P.iReadAll);
                            } else {
                                P.inOut->readIn[imate] >> (chunkIn[imate] + chunkInSizeBytesTotal[imate]);
                                chunkInSizeBytesTotal[imate] += strlen(chunkIn[imate] + chunkInSizeBytesTotal[imate]);
                            };

                            P.inOut->readIn[imate].ignore(DEF_readNameSeqLengthMax,'\n');

                            chunkInSizeBytesTotal[imate] += sprintf(chunkIn[imate] + chunkInSizeBytesTotal[imate], " %llu %c %i \n", P.iReadAll, 'N', P.readFilesIndex);


                        };
//                         else {//2nd stage of 2-stage mapping
//                         read index and file index are already recorded with the read name, simply copy it
//                         P.inOut->readIn[imate].getline(chunkIn[imate] + chunkInSizeBytesTotal[imate], DEF_readNameSeqLengthMax+1 );
//                         };
                        nextChar=P.inOut->readIn[imate].peek();
                        while (nextChar!='@' && nextChar!='>' && nextChar!=' ' && nextChar!='\n' && P.inOut->readIn[imate].good()) {//read multi-line fasta
                            P.inOut->readIn[imate].getline(chunkIn[imate] + chunkInSizeBytesTotal[imate], DEF_readSeqLengthMax + 1 );
                            if (P.inOut->readIn[imate].gcount()<2) break; //no more input
                            chunkInSizeBytesTotal[imate] += P.inOut->readIn[imate].gcount()-1;
                            nextChar=P.inOut->readIn[imate].peek();
                        };
                        chunkIn[imate][chunkInSizeBytesTotal[imate]]='\n';
                        chunkInSizeBytesTotal[imate] ++;
                    }; 
                } else if (nextChar==' ' || nextChar=='\n' || !P.inOut->readIn[0].good()) {//end of stream
                    P.inOut->logMain << "Thread #" <<iThread <<" end of input stream, nextChar="<<int(nextChar) <<endl;
                    break;
                } else {
                    string word1;
                    P.inOut->readIn[0] >> word1;
                    if (word1=="FILE") {//new file marker
                        newFile=true;
                    } else {//error
                        ostringstream errOut;
                        errOut << ERROR_OUT <<" EXITING because of FATAL ERROR in input reads: unknown file format: the read ID should start with @ or > \n";
                        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                    };
                };
                
                if (newFile) {
                        P.inOut->readIn[0] >> P.readFilesIndex;
                        pthread_mutex_lock(&g_threadChunks.mutexLogMain);
                        P.inOut->logMain << "Starting to map file # " << P.readFilesIndex<<"\n";
                        for (uint imate=0; imate<P.readFilesNames.size(); imate++) {
                            P.inOut->logMain << "mate " <<imate+1 <<":   "<<P.readFilesNames.at(imate).at(P.readFilesIndex) <<"\n";
                            P.inOut->readIn[imate].ignore(numeric_limits<streamsize>::max(),'\n');
                        };
                        P.inOut->logMain<<flush;
                        pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
                        newFile=false;
                };
            };
            //TODO: check here that both mates are zero or non-zero
            if (chunkInSizeBytesTotal[0]==0) {
                noReadsLeft=true; //true if there no more reads left in the file
                iChunkIn=g_threadChunks.chunkInN;//to keep things consistent
                g_threadChunks.chunkInN++;
            } else {
                noReadsLeft=false;
                iChunkIn=g_threadChunks.chunkInN;
                g_threadChunks.chunkInN++;
            };

            for (uint imate=0; imate<P.readNmates; imate++) chunkIn[imate][chunkInSizeBytesTotal[imate]]='\n';//extra empty line at the end of the chunks

            if (P.runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexInRead);

        } else {//read from one file per thread
            noReadsLeft=true;
            for (uint imate=0; imate<P.readNmates; imate++) {
                RA->chunkOutFilterBySJoutFiles[imate].flush();
                RA->chunkOutFilterBySJoutFiles[imate].seekg(0,ios::beg);
                RA->readInStream[imate]=& RA->chunkOutFilterBySJoutFiles[imate];
            };
        };

        mapChunk();

        if (iThread==0 && P.runThreadN>1 && P.outSAMorder=="PairedKeepInputOrder") {//concatenate Aligned.* files
            chunkFilesCat(P.inOut->outSAM, P.outFileTmp + "/Aligned.out.sam.chunk", g_threadChunks.chunkOutN);
        };

    };//cycle over input chunks

    if (P.outFilterBySJoutStage!=1 && RA->iRead>0) {//not the first stage of the 2-stage mapping
        if (P.outBAMunsorted) chunkOutBAMunsorted->unsortedFlush();
        if (P.outBAMcoord) chunkOutBAMcoord->coordFlush();
        if (chunkOutBAMquant!=NULL) chunkOutBAMquant->unsortedFlush();

        //the thread is finished mapping reads, concatenate the temp files into output files
        if (P.pCh.segmentMin>0) {
            chunkFstreamCat (RA->chunkOutChimSAM, P.inOut->outChimSAM, P.runThreadN>1, g_threadChunks.mutexOutChimSAM);
            chunkFstreamCat (*RA->chunkOutChimJunction, P.inOut->outChimJunction, P.runThreadN>1, g_threadChunks.mutexOutChimJunction);
        };
        if (P.outReadsUnmapped=="Fastx" ) {
            for (uint ii=0;ii<P.readNmates;ii++) {
                chunkFstreamCat (RA->chunkOutUnmappedReadsStream[ii],P.inOut->outUnmappedReadsStream[ii], P.runThreadN>1, g_threadChunks.mutexOutUnmappedFastx);
            };
        };
    };
    if (P.runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexLogMain);
    P.inOut->logMain << "Completed: thread #" <<iThread <<endl;
    if (P.runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
};

