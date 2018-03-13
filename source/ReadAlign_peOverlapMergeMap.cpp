#include "ReadAlign.h"
#include "SequenceFuns.h"

void ReadAlign::peOverlapMergeMap() {
    
    peOv.yes=false;
    
    if (!P.peOverlap.yes || P.readNmates!=2 ) {//no peOverlap
        return;
    };

    //debug
    //cout << ">" << readName+1;

        
    //merge PE mates into SE
    peMergeRA->copyRead(*this);
    peMergeRA->peMergeMates();
    peOv=peMergeRA->peOv;


    
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
    
    peMergeRA->peOv=peOv;
    //SE alignment is better, copy it to *this ReadAlign
    peOverlapSEtoPE(*peMergeRA);
    
    //chimeric detection for SE
    chimericDetectionPEmerged(*peMergeRA);

    //debug    
    //cout << "\n";
    //for (uint ii=0;ii<peMergeRA->Lread;ii++) {
    //    cout <<P.genomeNumToNT[peMergeRA->Read1[0][ii]];
    //};
    //cout << "\n";

    //P.alignSplicedMateMapLminOverLmate=P_alignSplicedMateMapLminOverLmate;
    peOv.yes=true;
    return;
};

void ReadAlign::peMergeMates() {
    
    peOv.ovS=localSearch(Read1[0],readLength[0],Read1[0]+readLength[0]+1,readLength[1],P.peOverlap.MMp);
    peOv.nOv=readLength[0]-peOv.ovS;
    
    if (peOv.nOv<P.peOverlap.NbasesMin) {//overlap is smaller than minimum allowed
        peOv.nOv=0;
        return;
    };
 
    //uint nMM=0;
    //for (uint ii=peOv.ovS; ii<readLength[0]; ii++) {//check for MM in the overlap area
    //    if (Read1[0][ii]!=Read1[0][ii-peOv.ovS+readLength[0]+1]) {
    //        Read1[0][ii]=4; //replace mismatched base with N
    //        ++nMM;
    //    };
    //};

    memmove(Read1[0]+readLength[0], Read1[0]+readLength[0]+1+peOv.nOv, readLength[1]-peOv.nOv);

    Lread=Lread-peOv.nOv-1;
    readLength[0]=Lread;    
    readLength[1]=0;
    readNmates=1;
    
    //fill Read1[1,2]
    complementSeqNumbers(Read1[0],Read1[1],Lread); //returns complement of Reads[ii]
    for (uint ii=0;ii<Lread;ii++) {//reverse
        Read1[2][Lread-ii-1]=Read1[1][ii];
        if (Read1[1][ii]<4) {
            Qual1[0][ii]=1;
            Qual1[1][Lread-ii-1]=1;
        } else {
            Qual1[0][ii]=0;
            Qual1[1][Lread-ii-1]=0;
        };

    };    

    return;
};

void Transcript::peOverlapSEtoPE(uint nOv, Transcript &t) {//convert alignment from merged-SE to PE
    
    uint len1=readLength[t.Str];//this transcrip=trInit which was initilized for original PE read
    uint start2=len1-nOv;//first base of the 2nd mate
    
    uint iex=0;
    for ( ; iex<t.nExons; iex++) {//first, cycle through the exons from mate1 
        if (t.exons[iex][EX_R] >= len1) {//this exon is only in mate2, break this cycle
            break;
        };
        //record these exons for mate1

        exons[iex][EX_iFrag]=t.Str;
        exons[iex][EX_sjA]=t.exons[iex][EX_sjA];
        canonSJ[iex]=t.canonSJ[iex];
        sjAnnot[iex]=t.sjAnnot[iex];
        sjStr[iex]=t.sjStr[iex];
        shiftSJ[iex][0]=t.shiftSJ[iex][0];
        shiftSJ[iex][1]=t.shiftSJ[iex][1];

        exons[iex][EX_R]=t.exons[iex][EX_R];
        exons[iex][EX_G]=t.exons[iex][EX_G];        
        if (t.exons[iex][EX_R]+t.exons[iex][EX_L] < len1) {//exon is fully in mate1
            exons[iex][EX_L]=t.exons[iex][EX_L];
        } else {
            exons[iex][EX_L]=len1-t.exons[iex][EX_R];
        };
    };
        
    nExons=iex;
    canonSJ[nExons-1]=-3; //marks "junction" between mates
    
    for (iex=0; iex<t.nExons; iex++) {//cycle through the exons from mate2
        if (t.exons[iex][EX_R]+t.exons[iex][EX_L] <= start2) {//this exon is only in mate2, do not record here
            continue;
        };
        
        exons[nExons][EX_iFrag]=1-t.Str;
        exons[nExons][EX_sjA]=t.exons[iex][EX_sjA];
        canonSJ[nExons]=t.canonSJ[iex];
        sjAnnot[nExons]=t.sjAnnot[iex];
        sjStr[nExons]=t.sjStr[iex];
        shiftSJ[nExons][0]=t.shiftSJ[iex][0];
        shiftSJ[nExons][1]=t.shiftSJ[iex][1];
        
        //record these exons for mate2
        if (t.exons[iex][EX_R]>=start2) {//exon is fully in mate2
            exons[nExons][EX_G]=t.exons[iex][EX_G];  
            exons[nExons][EX_L]=t.exons[iex][EX_L];            
            exons[nExons][EX_R]=t.exons[iex][EX_R]+nOv+1;
        } else {//need to split the exon 
            exons[nExons][EX_G]=t.exons[iex][EX_G]+start2-t.exons[iex][EX_R];  
            exons[nExons][EX_L]=t.exons[iex][EX_L]-start2+t.exons[iex][EX_R];            
            exons[nExons][EX_R]=len1+1;            
        };
        ++nExons;
    };
    
    
    
    //copy scalar variables
    for (uint ii=0;ii<4;ii++) {
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
  
    nW=seRA.nW;
    memcpy((void*) nWinTr, (void*) seRA.nWinTr, nW*sizeof(*nWinTr));
    
    uint trNtotal=0;
    intScore bestScore=0;
    trBest=NULL;
    for (uint iW=0; iW<nW; iW++) {//scan windows
        trAll[iW]=trArrayPointer+trNtotal;        
        for (uint iTr=0; iTr<nWinTr[iW]; iTr++) {//scan transcripts
            ++trNtotal;    
            *trAll[iW][iTr]=*trInit;
            trAll[iW][iTr]->peOverlapSEtoPE(peOv.nOv, *seRA.trAll[iW][iTr]);
            trAll[iW][iTr]->alignScore(Read1,mapGen.G,P);   
            if (trAll[iW][iTr]->maxScore > trAll[iW][0]->maxScore) {
                swap(trAll[iW][iTr],trAll[iW][0]);
            };
        };
        if (trAll[iW][0]->maxScore>bestScore) {
            trBest=trAll[iW][0];
            bestScore=trBest->maxScore;
        };        
    };   
    
    return;
};
