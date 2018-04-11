#include "ReadAlignChunk.h"
#include "Parameters.h"
#include "OutSJ.h"
#include <limits.h>
#include "ErrorWarning.h"

int compareUint(const void* i1, const void* i2) {//compare uint arrays
    uint s1=*( (uint*)i1 );
    uint s2=*( (uint*)i2 );

    if (s1>s2) {
        return 1;
    } else if (s1<s2) {
        return -1;
    } else {
        return 0;
    };
};

void outputSJ(ReadAlignChunk** RAchunk, Parameters& P) {//collapses junctions from all therads/chunks; outputs junctions to file

//     system("echo `date` ..... Writing splice junctions >> Log.timing.out");


    Junction oneSJ(RAchunk[0]->mapGen);
    char** sjChunks = new char* [P.runThreadN+1];
    #define OUTSJ_limitScale 5
    OutSJ allSJ (P.limitOutSJcollapsed*OUTSJ_limitScale,P,RAchunk[0]->mapGen);

    if (P.outFilterBySJoutStage!=1) {//chunkOutSJ
        for (int ic=0;ic<P.runThreadN;ic++) {//populate sjChunks with links to data
            sjChunks[ic]=RAchunk[ic]->chunkOutSJ->data;
            memset(sjChunks[ic]+RAchunk[ic]->chunkOutSJ->N*oneSJ.dataSize,255,oneSJ.dataSize);//mark the junction after last with big number
        };
    } else {//chunkOutSJ1
        for (int ic=0;ic<P.runThreadN;ic++) {//populate sjChunks with links to data
            sjChunks[ic]=RAchunk[ic]->chunkOutSJ1->data;
            memset(sjChunks[ic]+RAchunk[ic]->chunkOutSJ1->N*oneSJ.dataSize,255,oneSJ.dataSize);//mark the junction after last with big number
        };
    };

    while (true) {
        int icOut=-1;//chunk from which the junction is output
        for (int ic=0;ic<P.runThreadN;ic++) {//scan through all chunks, find the "smallest" junction
            if ( *(uint*)(sjChunks[ic])<ULONG_MAX && (icOut==-1 ||compareSJ((void*) sjChunks[ic], (void*) sjChunks[icOut])<0 ) ) {
                    icOut=ic;
                };
        };

        if (icOut<0) break; //no more junctions to output

        for (int ic=0;ic<P.runThreadN;ic++) {//scan through all chunks, find the junctions equal to icOut-junction
            if (ic!=icOut && compareSJ((void*) sjChunks[ic], (void*) sjChunks[icOut])==0) {
                oneSJ.collapseOneSJ(sjChunks[icOut],sjChunks[ic],P);//collapse ic-junction into icOut
                sjChunks[ic] += oneSJ.dataSize;//shift ic-chunk by one junction
            };
        };

        //write out the junction
        oneSJ.junctionPointer(sjChunks[icOut],0);//point to the icOut junction
        //filter the junction
        bool sjFilter;
        sjFilter=*oneSJ.annot>0 \
                || ( ( *oneSJ.countUnique>=(uint) P.outSJfilterCountUniqueMin[(*oneSJ.motif+1)/2] \
                    || (*oneSJ.countMultiple+*oneSJ.countUnique)>=(uint) P.outSJfilterCountTotalMin[(*oneSJ.motif+1)/2] )\
                && *oneSJ.overhangLeft >= (uint) P.outSJfilterOverhangMin[(*oneSJ.motif+1)/2] \
                && *oneSJ.overhangRight >= (uint) P.outSJfilterOverhangMin[(*oneSJ.motif+1)/2] \
                && ( (*oneSJ.countMultiple+*oneSJ.countUnique)>P.outSJfilterIntronMaxVsReadN.size() || *oneSJ.gap<=(uint) P.outSJfilterIntronMaxVsReadN[*oneSJ.countMultiple+*oneSJ.countUnique-1]) );

        if (sjFilter) {//record the junction in all SJ
            memcpy(allSJ.data+allSJ.N*oneSJ.dataSize,sjChunks[icOut],oneSJ.dataSize);
            allSJ.N++;
            if (allSJ.N == P.limitOutSJcollapsed*OUTSJ_limitScale ) {
                ostringstream errOut;
                errOut <<"EXITING because of fatal error: buffer size for SJ output is too small\n";
                errOut <<"Solution: increase input parameter --limitOutSJcollapsed\n";
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            };
        };

        sjChunks[icOut] += oneSJ.dataSize;//shift icOut-chunk by one junction
    };

    bool* sjFilter=new bool[allSJ.N];
    if (P.outFilterBySJoutStage!=2) {
        //filter non-canonical junctions that are close to canonical
        uint* sjA = new uint [allSJ.N*3];
        for (uint ii=0;ii<allSJ.N;ii++) {//scan through all junctions, filter by the donor ditance to a nearest donor, fill acceptor array
            oneSJ.junctionPointer(allSJ.data,ii);

            sjFilter[ii]=false;
            uint x1=0, x2=-1;
            if (ii>0)         x1=*( (uint*)(allSJ.data+(ii-1)*oneSJ.dataSize) ); //previous junction donor
            if (ii+1<allSJ.N) x2=*( (uint*)(allSJ.data+(ii+1)*oneSJ.dataSize) ); //next junction donor
            uint minDist=min(*oneSJ.start-x1, x2-*oneSJ.start);
            sjFilter[ii]= minDist >= (uint) P.outSJfilterDistToOtherSJmin[(*oneSJ.motif+1)/2];
            sjA[ii*3]=*oneSJ.start+(uint)*oneSJ.gap;//acceptor
            sjA[ii*3+1]=ii;

            if (*oneSJ.annot==0) {
                sjA[ii*3+2]=*oneSJ.motif;
            } else {
                sjA[ii*3+2]=SJ_MOTIF_SIZE+1;
            };

        };
        qsort((void*) sjA, allSJ.N, sizeof(uint)*3, compareUint);
        for (uint ii=0;ii<allSJ.N;ii++) {//
            if (sjA[ii*3+2]==SJ_MOTIF_SIZE+1) {//no filtering for annotated junctions
                sjFilter[sjA[ii*3+1]]=true;
            } else {
                uint x1=0, x2=-1;
                if (ii>0)         x1=sjA[ii*3-3]; //previous junction donor
                if (ii+1<allSJ.N) x2=sjA[ii*3+3]; //next junction donor
                uint minDist=min(sjA[ii*3]-x1, x2-sjA[ii*3]);
                sjFilter[sjA[ii*3+1]] = sjFilter[sjA[ii*3+1]] && ( minDist >= (uint) P.outSJfilterDistToOtherSJmin[(sjA[ii*3+2]+1)/2] );
            };
        };
    };

    //output junctions
    if (P.outFilterBySJoutStage!=1) {//output file
        ofstream outSJfileStream((P.outFileNamePrefix+"SJ.out.tab").c_str());
        for (uint ii=0;ii<allSJ.N;ii++) {//write to file
            if ( P.outFilterBySJoutStage==2 || sjFilter[ii]  ) {
                oneSJ.junctionPointer(allSJ.data,ii);
                oneSJ.outputStream(outSJfileStream);//write to file
            };
        };
        outSJfileStream.close();
    } else {//make sjNovel array in P
        P.sjNovelN=0;
        for (uint ii=0;ii<allSJ.N;ii++)
        {//count novel junctions
            if (sjFilter[ii])
            {//only those passing filter
                oneSJ.junctionPointer(allSJ.data,ii);
                if (*oneSJ.annot==0) P.sjNovelN++;
            };
        };
        P.sjNovelStart = new uint [P.sjNovelN];
        P.sjNovelEnd = new uint [P.sjNovelN];
        P.inOut->logMain <<"Detected " <<P.sjNovelN<<" novel junctions that passed filtering, will proceed to filter reads that contained unannotated junctions"<<endl;

        uint isj=0;
        for (uint ii=0;ii<allSJ.N;ii++) {//write to file
            if (sjFilter[ii]) {
                oneSJ.junctionPointer(allSJ.data,ii);
                if (*oneSJ.annot==0) {//unnnotated only
                    P.sjNovelStart[isj]=*oneSJ.start;
                    P.sjNovelEnd[isj]=*oneSJ.start+(uint)(*oneSJ.gap)-1;
                    isj++;
                };
            };
        };
    };
};
