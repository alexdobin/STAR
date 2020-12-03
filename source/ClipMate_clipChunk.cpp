#include "ParametersClip.h"
#include "Parameters.h"
#include "SequenceFuns.h"

inline char* findChar(char *arr, char c);

void ClipMate::clipChunk(char *chArr, uint64 chSize)
{//clipping adapters from a chunk of reads
    if (type != 10) //type=10 for CellRanger4 5' clipping.
        return;
    
    char *chA1 = chArr;
    bool chGood = true; //=true after teh end of the chunk
    while (chGood) {//cycle over all
        
        int dbN1 = cr4->dbN; //maybe changed to a smaller value in the loop
        int idb = 0;
        for ( ; idb<cr4->dbN; idb++) {
            chA1 = findChar(chA1, '\n')+1; //skip read name
            
            char *chA2 = findChar(chA1, '\n');
            uint32 rL = (uint32) (chA2-chA1);
            //debug
            string tmp1(chA1, 91);
            
            
            
            cr4->opalFillOneSeq(idb, chA1, rL);
            
            cr4->storeClip[idb] = (uint8*) (chA2+1);//store the position of "+" character - we will record the clipped length there
            
            //before the next one    
            chA1 = chA2 + 3 + rL + 1; //start of the next read: skip \n+\n, quality=read lengt, \n
            if (chA1 > chArr+chSize) {
                chGood = false;
                dbN1 = idb+1;
                break;
            };
        };

        cr4->opalAlign((uint8_t*) adSeqNum.data(), adSeqNum.size(), dbN1);

        for (int idb=0; idb<dbN1; idb++) {//store results
            int L = cr4->opalRes[idb].endLocationTarget+1;
            int S = cr4->opalRes[idb].score;
            
            bool L0 = S<20 || (S==20 && L>26) || (S==21 && L>30);
            
            *cr4->storeClip[idb] = (uint8) (L0 ? 0 : L);
        };
    };
    
};

inline char* findChar(char *arr, char c)
{//find character in a character array. No check for out of boundary.
    char* cArr=arr;
    while (*cArr != c)
        ++cArr;
    
    return cArr;
};