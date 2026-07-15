#include "funCompareUintAndSuffixesMemcmp.h"
#include "IncludeDefine.h"

char* g_funCompareUintAndSuffixesMemcmp_G;
uint64_t g_funCompareUintAndSuffixesMemcmp_N;
uint64_t g_funCompareUintAndSuffixesMemcmp_L;

int funCompareUintAndSuffixesMemcmp ( const void *a, const void *b)
{
    const uint64_t* va=static_cast<const uint64_t*>(a);
    const uint64_t* vb=static_cast<const uint64_t*>(b);

    if (va[0]>vb[0])
    {
        return 1;
    } else if (va[0]<vb[0])
    {
        return -1;
    } else
    {//compare suffixes
        const char* ga=g_funCompareUintAndSuffixesMemcmp_G+va[1];
        const char* gb=g_funCompareUintAndSuffixesMemcmp_G+vb[1];
        uint64_t compareLength=0;
        if (va[1]<g_funCompareUintAndSuffixesMemcmp_N && vb[1]<g_funCompareUintAndSuffixesMemcmp_N)
        {
            compareLength=min(g_funCompareUintAndSuffixesMemcmp_L,
                              min(g_funCompareUintAndSuffixesMemcmp_N-va[1]+1,
                                  g_funCompareUintAndSuffixesMemcmp_N-vb[1]+1));
            const char* sentinel=static_cast<const char*>(memchr(ga, GENOME_spacingChar, compareLength));
            if (sentinel!=NULL)
            {
                compareLength=sentinel-ga+1;
            };
        };

        int comp=memcmp(ga, gb, compareLength);
        if (comp!=0)
        {
            return comp;
        };

        if (va[1]>vb[1])
        {
            return 1;
        } else if (va[1]<vb[1])
        {
            return -1;
        };
        return 0;
    };
};
