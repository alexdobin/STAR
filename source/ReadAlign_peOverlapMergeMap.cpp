#include "ReadAlign.h"
#include "SequenceFuns.h"

void ReadAlign::peOverlapMergeMap() {

    if (!P.peOverlap.yes || P.readNmates!=2 ) {//no peOverlap //not readNends: this is alignment
        peOv.yes=false;
        return;
    };

    //debug
    //cout << ">" << readName+1;


    //merge PE mates into SE
    peMergeRA->copyRead(*this);
    peMergeRA->peMergeMates();
    peOv=peMergeRA->peOv;
    peOv.yes=false;

    if (peOv.nOv==0) {//check if mates can be merged, if not - return
        //cout <<"\n-1\n";
        return;
    };

    //change parameters for SE mapping
    //double P_alignSplicedMateMapLminOverLmate=P.alignSplicedMateMapLminOverLmate;
    //P.alignSplicedMateMapLminOverLmate=P.alignSplicedMateMapLminOverLmate*peMergeRA->readLength[0]/(readLength[0]+readLength[1]);

    //map SE
    peMergeRA->mapOneRead();
    if (peMergeRA->nW==0) { // || peMergeRA->trBest->maxScore+peOv.nOv < trBest->maxScore) {//no windows, score of the merged align is less. This is a preliminary check, more accurate check is done with alignment score calculated after transforming the SE back to PE
        //cout <<" -2\n";
        //for (uint ii=0;ii<peMergeRA->Lread;ii++) {
        //    cout <<P.genomeNumToNT[peMergeRA->Read1[0][ii]];
        //};
        //cout << "\n";
        return;
    };

    //convert best alignment SE to PE
    //trA=*trInit;
    //trA.peOverlapSEtoPE(peOv.nOv, *peMergeRA->trBest);
    //trA.alignScore(Read1,mapGen.G,P);
    //if (trA.maxScore<trBest->maxScore || trA.nMM > outFilterMismatchNmaxTotal) {//merged-mate SE alignment has lower score than the PE
    //    return;
    //};

    intScore peScore=trBest->maxScore;

    //convert SE to PE *this ReadAlign
    peMergeRA->peOv=peOv;
    peOverlapSEtoPE(*peMergeRA);

    //debug
    //if (oldScore>trBest->maxScore || trBest->maxScore<peMergeRA->trBest->maxScore)
    //    cout << readName << "   "<< oldScore << "   "<< peMergeRA->trBest->maxScore << "   "<<trBest->maxScore << endl;


    //chimeric detection for SE
    chimericDetectionPEmerged(*peMergeRA);

    //debug
    //cout << "\n";
    //for (uint ii=0;ii<peMergeRA->Lread;ii++) {
    //    cout <<P.genomeNumToNT[peMergeRA->Read1[0][ii]];
    //};
    //cout << "\n";

    //P.alignSplicedMateMapLminOverLmate=P_alignSplicedMateMapLminOverLmate;

    if (peScore<=trBest->maxScore || chimRecord) {//otherwise peOv.yes=false
        peOv.yes=true;
    };

    return;
};

void ReadAlign::peMergeMates() {

    uint s1=localSearchNisMM(Read1[0],readLength[0],Read1[0]+readLength[0]+1,readLength[1],P.peOverlap.MMp);
    uint s0=localSearchNisMM(Read1[0]+readLength[0]+1,readLength[1],Read1[0],readLength[0],P.peOverlap.MMp);

    uint o1=min(readLength[1],readLength[0]-s1);
    uint o0=min(readLength[0],readLength[1]-s0);

    peOv.nOv=max(o0,o1);

    if (peOv.nOv<P.peOverlap.NbasesMin) {//overlap is smaller than minimum allowed
        peOv.nOv=0;
        return;
    };

    if (o1>=o0) {
        peOv.mateStart[0]=0;
        peOv.mateStart[1]=s1;
        if (o1<readLength[1]) {//otherwise, if o1==readLength[1], read2 is entirely contained in read1
            //move unoverlapped portion of read2 to the end of read1
            memmove(Read1[0]+readLength[0], Read1[0]+readLength[0]+1+o1, readLength[1]-o1);
        };
    } else {
        peOv.mateStart[1]=0;
        peOv.mateStart[0]=s0;
        memmove(Read1[0]+Lread, Read1[0], readLength[0]);//temp move 0
        memmove(Read1[0], Read1[0]+readLength[0]+1, readLength[1]); //move 1 into 0
        if (o0<readLength[0]) {
            memmove(Read1[0]+readLength[1], Read1[0]+Lread+o0, readLength[0]-o0); //move 0 into 1
        };
    };

    //uint nMM=0;
    //for (uint ii=peOv.ovS; ii<readLength[0]; ii++) {//check for MM in the overlap area
    //    if (Read1[0][ii]!=Read1[0][ii-peOv.ovS+readLength[0]+1]) {
    //        Read1[0][ii]=4; //replace mismatched base with N
    //        ++nMM;
    //    };
    //};


    Lread=Lread-peOv.nOv-1;
    readLength[0]=Lread;
    readLength[1]=0;
    readLengthOriginal[0]=Lread;
    readLengthOriginal[1]=0;
    readNmates=1; //not readNends: this is alignment

    //fill Read1[1,2]
    complementSeqNumbers(Read1[0],Read1[1],Lread); //returns complement of Reads[ii]
    for (uint ii=0;ii<Lread;ii++) {//reverse
        Read1[2][Lread-ii-1]=Read1[1][ii];
    };

    return;
};

void Transcript::peOverlapSEtoPE(uint* mateStart, const Transcript &t) {//convert alignment from merged-SE to PE

    uint mLen[2];
    mLen[0]=readLength[t.Str];
    mLen[1]=readLength[1-t.Str];

    uint mSta2[2];
    mSta2[0]=0;//mates starts in the PE read
    mSta2[1]=mLen[0]+1;

    uint mSta[2];
    mSta[0]=mateStart[0];//mates starts in the merged SE read
    mSta[1]=mateStart[1];
    if (t.Str==1) {
	for (uint ii=0;ii<2;ii++) {
            mSta[ii]=t.Lread-readLength[ii]-mSta[ii];
        };
        swap(mSta[0],mSta[1]);
    };

    uint mEnd[2];
    mEnd[0]=mSta[0]+mLen[0];
    mEnd[1]=mSta[1]+mLen[1];

//     uint iex=0;
//     for ( ; iex<t.nExons; iex++) {//first, cycle through the exons from mate1
//         if (t.exons[iex][EX_R] >= mEnd[0] || t.exons[iex][EX_R]+t.exons[iex][EX_L] < mSta[0]) {//this exon is only in mate2, break this cycle
//             break;
//         };
//         //record these exons for mate1
//
//         exons[iex][EX_iFrag]=t.Str;
//         exons[iex][EX_sjA]=t.exons[iex][EX_sjA];
//         canonSJ[iex]=t.canonSJ[iex];
//         sjAnnot[iex]=t.sjAnnot[iex];
//         sjStr[iex]=t.sjStr[iex];
//         shiftSJ[iex][0]=t.shiftSJ[iex][0];
//         shiftSJ[iex][1]=t.shiftSJ[iex][1];
//
//         exons[iex][EX_R]=t.exons[iex][EX_R]-mSta[0];
//         exons[iex][EX_G]=t.exons[iex][EX_G];
//         if (t.exons[iex][EX_R]+t.exons[iex][EX_L] < mEnd[0]) {//exon is fully in mate1
//             exons[iex][EX_L]=t.exons[iex][EX_L];
//         } else {
//             exons[iex][EX_L]=mEnd[0]-t.exons[iex][EX_R];
//         };
//     };

    nExons=0;
    for (uint imate=0; imate<2; imate++) {//cycle over mate 1,2
        for (uint iex=0; iex<t.nExons; iex++) {//cycle through the exons
            if (t.exons[iex][EX_R] >= mEnd[imate] || t.exons[iex][EX_R]+t.exons[iex][EX_L] <= mSta[imate]) {//this exon is only in mate2, do not record here
                continue;
            };

            exons[nExons][EX_iFrag]=(imate==0 ? t.Str : 1-t.Str);
            exons[nExons][EX_sjA]=t.exons[iex][EX_sjA];
            if (iex<t.nExons-1) {
                canonSJ[nExons]=t.canonSJ[iex];
                sjAnnot[nExons]=t.sjAnnot[iex];
                sjStr[nExons]=t.sjStr[iex];
                shiftSJ[nExons][0]=t.shiftSJ[iex][0];
                shiftSJ[nExons][1]=t.shiftSJ[iex][1];
            };
            //record these exons for mate2
            if (t.exons[iex][EX_R]>=mSta[imate]) {//exon left is inside the mate
                exons[nExons][EX_G]=t.exons[iex][EX_G];
                exons[nExons][EX_L]=t.exons[iex][EX_L];
                exons[nExons][EX_R]=t.exons[iex][EX_R]-mSta[imate]+mSta2[imate];
            } else {//need to split the exon
                exons[nExons][EX_R]=mSta2[imate];//exon starts at the mate start
                uint delta=mSta[imate]-t.exons[iex][EX_R]; //shorten exon by this length
                exons[nExons][EX_L]=t.exons[iex][EX_L]-delta;
                exons[nExons][EX_G]=t.exons[iex][EX_G]+delta;
            };

            if (t.exons[iex][EX_R]+t.exons[iex][EX_L] > mEnd[imate]) {//exon right is to the left of the mate end, shorten the exon
                exons[nExons][EX_L]-=t.exons[iex][EX_R]+t.exons[iex][EX_L]-mEnd[imate];
            };

            ++nExons;
            if (nExons>MAX_N_EXONS) {//cannot transform this alignment to PE, too many exons
                maxScore=0;
                nExons=0;
                return;
            };
        };
        canonSJ[nExons-1]=-3; //marks "junction" between mates
        sjAnnot[nExons-1]=0;
        sjStr[nExons-1]=0;
        shiftSJ[nExons-1][0]=0;
        shiftSJ[nExons-1][1]=0;
    };

    //copy scalar variables
    for (uint ii=0;ii<3;ii++) {
        intronMotifs[ii]=t.intronMotifs[ii];
    };
    sjMotifStrand=t.sjMotifStrand;
    //iFrag; //do not need it
    Chr=t.Chr;
    Str=t.Str;
    roStr=t.roStr;
    gStart=t.gStart;
    gLength=t.gLength;
    cStart=t.cStart;

    rLength=0;
    for (uint iex=0;iex<nExons;iex++) {//caclulate total mapped length
        rLength += exons[iex][EX_L];
    };
    mappedLength=rLength ;
    rStart = exons[0][EX_R];
    roStart = (roStr == 0) ? rStart : Lread - rStart - rLength;

    //extendL; //do not need

    nGap=t.nGap;
    lGap=t.lGap;
    nDel=t.nDel;
    nIns=t.nIns;
    lDel=t.nDel;
    lIns=t.lIns;

    nUnique=t.nUnique;
    nAnchor=t.nAnchor;

    return;
};

void ReadAlign::peOverlapSEtoPE(ReadAlign &seRA) {//ReAdAlign: convert SE to PE and copy

    //nW=seRA.nW;
    //memcpy((void*) nWinTr, (void*) seRA.nWinTr, nW*sizeof(*nWinTr));

    uint trNtotal=0;
    intScore bestScore=-10*Lread;
    trBest=trArray;//just to initialize - to the 0th spot in the trArray
    
    uint64 iW1=0;
    for (uint iW=0; iW<seRA.nW; iW++) {//scan windows
        trAll[iW1]=trArrayPointer+trNtotal;
        uint64 iTr1=0;
        for (uint iTr=0; iTr<seRA.nWinTr[iW]; iTr++) {//scan transcripts
            *trAll[iW1][iTr1]=*trInit;
            
            trAll[iW1][iTr1]->peOverlapSEtoPE(peOv.mateStart, *seRA.trAll[iW][iTr]);
            if (trAll[iW1][iTr1]->nExons==0)
                continue; //conversion did not work
            
            trAll[iW1][iTr1]->alignScore(Read1,mapGen.G,P);
            if (trAll[iW1][iTr1]->maxScore > trAll[iW1][0]->maxScore) {
                swap(trAll[iW][iTr1],trAll[iW][0]);
            };
            
            ++iTr1;
            ++trNtotal;            
        };
        if (iTr1>0) {//if conversion worked for at least one align
            nWinTr[iW1]=iTr1;
            if (trAll[iW1][0]->maxScore > bestScore) {
                trBest=trAll[iW1][0];
                bestScore=trBest->maxScore;
            };            
            ++iW1;
        };
    };

    nW=iW1;
    return;
};

void ReadAlign::peOverlapChimericSEtoPE(const Transcript *seTrIn1, const Transcript *seTrIn2, Transcript *peTrOut1, Transcript *peTrOut2) {

    //convert merged into PE
    Transcript tempTrChim[2];
    tempTrChim[0]=*trInit;
    tempTrChim[1]=*trInit;
    tempTrChim[0].peOverlapSEtoPE(peOv.mateStart,*seTrIn1);
    tempTrChim[1].peOverlapSEtoPE(peOv.mateStart,*seTrIn2);

    uint segLen[2][2]; //segment length [tempTrChim][mate]
    uint segEx[2];//last exon of the mate0 [tempTrChim]
    uint i1=0,i2=0; //indices of mate to eliminate i1=[tempTrChim], i2=[mate]
    uint posOfJunctionInRead=0; //position of chimeric junction in read
    for (uint ii=0; ii<2; ii++) {
        segLen[ii][0]=0;
        segLen[ii][1]=0;
        for (uint iex=0; iex<tempTrChim[ii].nExons; iex++) {
            if (tempTrChim[ii].exons[iex][EX_iFrag]==tempTrChim[ii].exons[0][EX_iFrag]) {
                segLen[ii][0]+=tempTrChim[ii].exons[iex][EX_L];
                segEx[ii]=iex;
            } else {
                segLen[ii][1]+=tempTrChim[ii].exons[iex][EX_L];
            };
        };

        //find mate with shortest mapped segment
        //in case of tie, use longest mate as tie-breaker (where posOfJunctionInRead is highest)
        uint readLen0=readLengthOriginal[tempTrChim[ii].exons[0][EX_iFrag]];
        uint readLen1=readLengthOriginal[1-tempTrChim[ii].exons[0][EX_iFrag]];
        for (uint jj=0; jj<2; jj++) {
            uint curPosOfJunctionInRead = tempTrChim[ii].exons[jj][EX_R]>readLen0 ? readLen0+readLen1+1-tempTrChim[ii].exons[jj][EX_R] : tempTrChim[ii].exons[jj][EX_R];
            if (segLen[ii][jj]<segLen[i1][i2] || (segLen[ii][jj]==segLen[i1][i2] && curPosOfJunctionInRead>posOfJunctionInRead)) {
                posOfJunctionInRead=curPosOfJunctionInRead;
                i1=ii;//tempTrChim of the shortest segment length
                i2=jj;//mate of the shortest segment length
            };
        };
    };

    if (i2==1) {//eliminate mate1: simply cut the exons that belong to mate1
        tempTrChim[i1].nExons=segEx[i1]+1;
    } else {//eliminate mate 0: shift mate1 exon to the beginning
        for (uint iex=0; iex<tempTrChim[i1].nExons; iex++) {
            uint iex1=iex+segEx[i1]+1;
            for (uint ii=0; ii<EX_SIZE; ii++) {
                tempTrChim[i1].exons[iex][ii]=tempTrChim[i1].exons[iex1][ii];
            };
            tempTrChim[i1].canonSJ[iex]=tempTrChim[i1].canonSJ[iex1];
            tempTrChim[i1].sjAnnot[iex]=tempTrChim[i1].sjAnnot[iex1];
            tempTrChim[i1].sjStr[iex]=tempTrChim[i1].sjStr[iex1];
            tempTrChim[i1].shiftSJ[iex][0]=tempTrChim[i1].shiftSJ[iex1][0];
            tempTrChim[i1].shiftSJ[iex][1]=tempTrChim[i1].shiftSJ[iex1][1];
        };
        tempTrChim[i1].nExons=tempTrChim[i1].nExons-segEx[i1]-1;
    };

    *peTrOut1=tempTrChim[0];
    *peTrOut2=tempTrChim[1];
    return;
};

