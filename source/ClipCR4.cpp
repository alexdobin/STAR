#include "ClipCR4.h"

ClipCR4::ClipCR4()
{
    dbN = 64;
    
    scoreMatrix = {
        1, -2, -2, -2, -2,
       -2,  1, -2, -2, -2,
       -2, -2,  1, -2, -2,
       -2, -2, -2,  1, -2,
       -2, -2, -2,  -2, 0
    };
    
    readLen = 91;
    alphabetLength = 5;
    gapOpen = 2;
    gapExt = 2;
    
    
    dbSeqArr = new uint8_t[dbN*readLen];
    dbSeqs = new uint8_t*[dbN];
    for (int id=0; id<dbN; id++) {
        dbSeqs[id]=dbSeqArr+id*readLen;
    };
    
    dbSeqsLen = new int[dbN];
    for (int id=0; id<dbN; id++) {
        dbSeqsLen[id] = readLen; //rL;
    };  
    
    storeClip.resize(dbN);
    opalRes.resize(dbN);
    opalResP = new OpalSearchResult*[dbN];
    
    for (int i = 0; i < dbN; i++) {
        opalResP[i]=&opalRes[i];
    };     

};


void ClipCR4::opalFillOneSeq(uint32 idb, char *seq, uint32 seqL)
{
    uint32 minLen = min(seqL, readLen);
    
    memcpy(dbSeqs[idb], seq, minLen);
    
    for (uint32 ib=0; ib<minLen; ib++) {//TODO: vectorize with multiplications, i.e. dbSeq= (seq=='C')+(seq=='G')*2+(seq=='T')*3+(seq=='N')*4
        switch (dbSeqs[idb][ib]) {
            case 'A':
                dbSeqs[idb][ib]=0;
                break;
            case 'C':
                dbSeqs[idb][ib]=1;
                break;        
            case 'G':
                dbSeqs[idb][ib]=2;
                break;   
            case 'T':
                dbSeqs[idb][ib]=3;
                break;   
            default:
                dbSeqs[idb][ib]=4;
        };
    };
    
    if (seqL<readLen) //fill the rest of the sequence with Ns
        memset(dbSeqs[idb]+seqL, 4, readLen-seqL);
};

void ClipCR4::opalAlign(uint8_t *query, uint32 queryLen, int dbN1)
{
    for (int idb = 0; idb < dbN1; idb++) {
        opalInitSearchResult(opalResP[idb]);
    };    
    opalSearchDatabase(query, queryLen, dbSeqs, dbN1, dbSeqsLen,
                        gapOpen, gapExt, scoreMatrix.data(), alphabetLength, opalResP,
                        OPAL_SEARCH_SCORE_END, OPAL_MODE_OV, OPAL_OVERFLOW_BUCKETS);//, OPAL_OVERFLOW_SIMPLE);//
};

uint32 ClipCR4::polyTail3p(char *seq, uint32 seqLen)
{//clip polyA tail
    //hardcoded for CR4 trimming
    if (seqLen<20)
        return 0; //do not trim reads that are too short

    uint32_t ib1=seqLen-1;
    int32_t score=0, score1=0;
    for (uint32_t ib=1; ib<=seqLen; ib++) {
        if ( seq[seqLen-ib] == 0 ) {//A-tail
            score += 1; //+1 for matches
            if (score*10 >= (int)ib*7) {//score>=0.7*clipL
                ib1 = ib;
                score1 = score;
            };
        } else {
            score -= 2; //-2 for mismatches
            if ( ib-score > 27 ) 
                break; //score drop is too big
        };
    };
    if ( score1<20 )
        ib1=0; //score is too small
        
    return ib1;
};


