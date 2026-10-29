#include "Transcript.h"

Transcript::Transcript()
{
    reset();
};

Transcript::Transcript(const Transcript &t)
    : cigar(t.cigar), varInd(t.varInd), varGenCoord(t.varGenCoord), varReadCoord(t.varReadCoord), varAllele(t.varAllele), alignGenes(t.alignGenes)
{
    copyRowsScalars(t);
};

Transcript& Transcript::operator=(const Transcript &t)
{
    if (this != &t) {
        copyRowsScalars(t);
        cigar=t.cigar;
        varInd=t.varInd;
        varGenCoord=t.varGenCoord;
        varReadCoord=t.varReadCoord;
        varAllele=t.varAllele;
        alignGenes=t.alignGenes;
    };
    return *this;
};

void Transcript::copyRowsScalars(const Transcript &t)
{//the per-exon arrays are copied for rows 0..nExons (row nExons is read by a few callers), not all MAX_N_EXONS
 //rows past those keep the destination's old values. Every scalar member must be listed here.
    uint nRows = t.nExons<MAX_N_EXONS ? t.nExons+1 : MAX_N_EXONS;
    const uint kRows = MAX_N_EXONS<4 ? MAX_N_EXONS : 4; //fixed-size copy of the first rows, the rest only when live
    memcpy(exons, t.exons, kRows*sizeof(exons[0]));
    memcpy(shiftSJ, t.shiftSJ, kRows*sizeof(shiftSJ[0]));
    memcpy(canonSJ, t.canonSJ, kRows*sizeof(canonSJ[0]));
    memcpy(sjAnnot, t.sjAnnot, kRows*sizeof(sjAnnot[0]));
    memcpy(sjStr, t.sjStr, kRows*sizeof(sjStr[0]));
    if (nRows>kRows) {
        memcpy(exons+kRows, t.exons+kRows, (nRows-kRows)*sizeof(exons[0]));
        memcpy(shiftSJ+kRows, t.shiftSJ+kRows, (nRows-kRows)*sizeof(shiftSJ[0]));
        memcpy(canonSJ+kRows, t.canonSJ+kRows, (nRows-kRows)*sizeof(canonSJ[0]));
        memcpy(sjAnnot+kRows, t.sjAnnot+kRows, (nRows-kRows)*sizeof(sjAnnot[0]));
        memcpy(sjStr+kRows, t.sjStr+kRows, (nRows-kRows)*sizeof(sjStr[0]));
    };

    intronMotifs[0]=t.intronMotifs[0]; intronMotifs[1]=t.intronMotifs[1]; intronMotifs[2]=t.intronMotifs[2];
    sjMotifStrand=t.sjMotifStrand;
    sjYes=t.sjYes;
    nExons=t.nExons;
    readLengthOriginal=t.readLengthOriginal; readLength=t.readLength;
    Lread=t.Lread; readLengthPairOriginal=t.readLengthPairOriginal;
    iRead=t.iRead;
    readNmates=t.readNmates;
    readName=t.readName;
    iFrag=t.iFrag;
    rStart=t.rStart; roStart=t.roStart; rLength=t.rLength; gStart=t.gStart; gLength=t.gLength; cStart=t.cStart;
    Chr=t.Chr; Str=t.Str; roStr=t.roStr;
    haploType=t.haploType;
    primaryFlag=t.primaryFlag;
    nMatch=t.nMatch;
    nMM=t.nMM;
    mappedLength=t.mappedLength;
    extendL=t.extendL;
    maxScore=t.maxScore;
    nGap=t.nGap; lGap=t.lGap;
    nDel=t.nDel;
    nIns=t.nIns;
    lDel=t.lDel;
    lIns=t.lIns;
    nUnique=t.nUnique; nAnchor=t.nAnchor;
};

void Transcript::reset() {
    extendL=0;

//     for (uint ii=0;ii<4;ii++) {
//         polyXlength[ii]=0;
//         polyXnMM[ii]=0;
//     };
    primaryFlag=false;

    rStart=0; roStart=0; rLength=0; gStart=0; gLength=0; //read and genomic coordinates

    maxScore=0;
    nMatch=0;
    nMM=0;

    nGap=0; lGap=0; lDel=0; lIns=0; nDel=0; nIns=0;

    nUnique=nAnchor=0;
};

void Transcript::add(Transcript *trIn) {
    maxScore+=trIn->maxScore;
    nMatch+=trIn->nMatch;
    nMM+=trIn->nMM;
    nGap+=trIn->nGap; lGap+=trIn->lGap;
    lDel+=trIn->lDel; nDel+=trIn->nDel;
    lIns+=trIn->lIns; nIns+=trIn->nIns;
    nUnique+=trIn->nUnique;
};

void Transcript::extractSpliceJunctions(vector<array<uint64,2>> &sjOut, bool &annotYes)
{
    annotYes=true;
    for (uint64 iex=0; iex<nExons-1; iex++) {//record all junctions
        if (canonSJ[iex]>=0) {//only record junctions, not indels or mate gap
            array<uint64,2> sj;
            sj[0]=exons[iex][EX_G]+exons[iex][EX_L];//start
            sj[1]=exons[iex+1][EX_G] - sj[0]; //gap
            sjOut.push_back(sj);
            if (sjAnnot[iex]==0)
                annotYes=false;//if one of the SJs is unannoated, annotYes=false
        };
    };
};

uint64 Transcript::chrStartLengthExtended()
{
    uint64 start1  = cStart - exons[0][EX_R];
    uint64 length1 = exons[nExons-1][EX_G] + Lread - exons[nExons-1][EX_R] - exons[0][EX_G] + exons[0][EX_R];
    
    return (start1 << 32) | length1;
};
