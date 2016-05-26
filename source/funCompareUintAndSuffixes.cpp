#include "funCompareUintAndSuffixes.h"

char* g_funCompareUintAndSuffixes_G;
uint64_t g_funCompareUintAndSuffixes_L;

int funCompareUintAndSuffixes ( const void *a, const void *b){
    uint64_t* va= ((uint64_t*) a);
    uint64_t* vb= ((uint64_t*) b);

    if (va[0]>vb[0]) {
            return 1;
        } else if (va[0]<vb[0]) {
            return -1;
        } else {//compare suffixes
            char* ga=g_funCompareUintAndSuffixes_G + *( ((uint64_t*) a)+1);
            char* gb=g_funCompareUintAndSuffixes_G + *( ((uint64_t*) b)+1);
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
