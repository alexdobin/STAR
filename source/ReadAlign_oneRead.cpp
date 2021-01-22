#include "ReadAlign.h"
#include "readLoad.h"
#include "readBarcodeLoad.h"
#include "SequenceFuns.h"
#include "ErrorWarning.h"
#include "GlobalVariables.h"

int ReadAlign::oneRead() {//process one read: load, map, write

    //load read name, sequence, quality from the streams into internal arrays
    int readStatus[2];


    readStatus[0]=readLoad(*(readInStream[0]), P, readLength[0], readLengthOriginal[0], readNameMates[0], Read0[0], Read1[0], Qual0[0], Qual1[0], clipMates[0], iReadAll, readFilesIndex, readFilter, readNameExtra[0]);
    if (P.readNmates==2) {//read the 2nd mate
        readStatus[1]=readLoad(*(readInStream[1]), P, readLength[1], readLengthOriginal[1], readNameMates[1], Read0[1], Read1[0]+readLength[0]+1, Qual0[1], Qual1[0]+readLength[0]+1, clipMates[1], iReadAll, readFilesIndex, readFilter, readNameExtra[1]);

        if (readStatus[0]!=readStatus[1]) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR: Read1 and Read2 are not consistent, reached the end of the one before the other one\n";
            errOut << "SOLUTION: Check you your input files: they may be corrupted\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        } else if (readStatus[0]==-1) {//finished with the stream
            return -1;
        };

        //combine two reads together
        Lread=readLength[0]+readLength[1]+1;
        readLengthPairOriginal=readLengthOriginal[0]+readLengthOriginal[1]+1;
        if (Lread>DEF_readSeqLengthMax) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in reads input: Lread of the pair = " << Lread << "   while DEF_readSeqLengthMax=" << DEF_readSeqLengthMax <<endl;
            errOut << "Read Name="<<readNameMates[0]<<endl;
            errOut << "SOLUTION: increase DEF_readSeqLengthMax in IncludeDefine.h and re-compile STAR"<<endl<<flush;
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };

        Read1[0][readLength[0]]=MARK_FRAG_SPACER_BASE; //marker for spacer base
        Qual1[0][readLength[0]]=0;
        complementSeqNumbers(Read1[0]+readLength[0]+1,Read1[0]+readLength[0]+1,readLength[1]); //returns complement of Reads[ii]
        for (uint ii=0;ii<readLength[1]/2;ii++) {
            swap(Read1[0][Lread-ii-1],Read1[0][ii+readLength[0]+1]); //reverse complement
            swap(Qual1[0][Lread-ii-1],Qual1[0][ii+readLength[0]+1]); //reverse complement   ??? was Qualof the second mate populated
        };

    } else {//1 mate

        if (readStatus[0]==-1) {//finished with the stream
            return -1;
        };

        Lread=readLength[0];
        readLengthPairOriginal=readLengthOriginal[0];
        readLength[1]=0;

    };

    if (P.outFilterBySJoutStage != 2) {
        for (uint32 imate=0; imate<P.readNmates; imate++)
            g_statsAll.qualHistCalc(imate, Qual0[imate], readLengthOriginal[imate]);
    };
    
    loadBarcodeRead(P, readInStream, readBarcodeSeq, readBarcodeQual);
    
    
    readFileType=readStatus[0];

    complementSeqNumbers(Read1[0],Read1[1],Lread); //returns complement of Reads[ii]
    for (uint ii=0;ii<Lread;ii++) {//reverse
        Read1[2][Lread-ii-1]=Read1[1][ii];
        Qual1[1][Lread-ii-1]=Qual1[0][ii];
    };

    statsRA.readN++;
    statsRA.readBases += readLength[0]+readLength[1];

    //max number of mismatches allowed for this read
    outFilterMismatchNmaxTotal=min(P.outFilterMismatchNmax, (uint) (P.outFilterMismatchNoverReadLmax*(readLength[0]+readLength[1])));

    //map the read
    if (P.pGe.gType==101) {//SpliceGraph
        mapOneReadSpliceGraph();
    } else {//all other cases - standard alignment algorithm
        mapOneRead();
    };

    peOverlapMergeMap();
    
    multMapSelect();
    
    mappedFilter();  
    
    transformGenome();//for now genome transformation happens after multimapper selection, and mapping filter

    if (!peOv.yes) {//if the alignment was not mates merged - otherwise the chimeric detection was already done
        chimericDetection();
    };

    if (P.pCh.out.bam && chimRecord) {//chimeric alignment was recorded in main BAM files, and it contains the representative portion, so non-chimeric aligmnent is not output
        return 0;
    };

    waspMap();

    #ifdef OFF_BEFORE_OUTPUT
        #warning OFF_BEFORE_OUTPUT
        return 0;
    #endif

    //write out alignments
    outputAlignments();

    return 0;

};


