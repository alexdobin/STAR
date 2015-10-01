#include "funCompareUintAndSuffixes.h"
#define uint64 unsigned long long

char* globalGenomeArray;

int funCompareUintAndSuffixes ( const void *a, const void *b){
    uint64* va= ((uint64*) a);
    uint64* vb= ((uint64*) b);

    if (va[0]>vb[0]) {
            return 1;
        } else if (va[0]<vb[0]) {
            return -1;
        } else {//compare suffixes
            char* ga=globalGenomeArray + *( ((uint64*) a)+1);
            char* gb=globalGenomeArray + *( ((uint64*) b)+1);
            int ig=0;
            while (true) {
                if (ga[ig]>gb[ig]) 
                {// second condition: reached the end of ga, it's >= than any character, but = does not matter
                    return 1;
                } else if (ga[ig]<gb[ig]) 
                {
                    return -1;
                } else if (ga[ig]==5) 
                {//reached the end of chr, now simply compare the indexes for stable search
                    if (va[1]>vb[1])
                    {
                        return 1;
                    } else
                    {//va cannot be equal to vb
                        return -1;
                    };
                } else 
                {//continue
                    ig++;
                };
            };
        };
};
