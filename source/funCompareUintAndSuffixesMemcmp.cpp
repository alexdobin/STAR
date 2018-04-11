#include "funCompareUintAndSuffixesMemcmp.h"
#include <string.h> //for memcmp

char* g_funCompareUintAndSuffixesMemcmp_G;
uint64_t g_funCompareUintAndSuffixesMemcmp_L;

int funCompareUintAndSuffixesMemcmp ( const void *a, const void *b)
{   
    uint64_t* va= ((uint64_t*) a);
    uint64_t* vb= ((uint64_t*) b);

    if (va[0]>vb[0]) 
    {
        return 1;
    } else if (va[0]<vb[0]) 
    {
        return -1;
    } else 
    {//compare suffixes
//         char *p5=(char*) memchr(g_funCompareUintAndSuffixesMemcmp_G+va[1],5,g_funCompareUintAndSuffixesMemcmp_L); //first encounter of char=5
//         char *p5=g_funCompareUintAndSuffixesMemcmp_G+va[1]+g_funCompareUintAndSuffixesMemcmp_L;
        //compare suffixes but only until first char=5
//         int comp=memcmp(g_funCompareUintAndSuffixesMemcmp_G+va[1],g_funCompareUintAndSuffixesMemcmp_G+vb[1],p5+1-(g_funCompareUintAndSuffixesMemcmp_G+va[1]));
        int comp=memcmp(g_funCompareUintAndSuffixesMemcmp_G+va[1],g_funCompareUintAndSuffixesMemcmp_G+vb[1],g_funCompareUintAndSuffixesMemcmp_L);
        
        if (comp==0)
        {
            comp=va[1]>vb[1] ? 1 : -1;
        };
//         int comp=va[1]>vb[1] ? 1 : -1;
        return comp;
    };
};
