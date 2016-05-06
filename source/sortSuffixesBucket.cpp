#include "sortSuffixesBucket.h"
#include <string.h> //for memset
#include <iostream> //fro cout - debug

#define GENOME_charN 6
//max number of characters in the text (genome), presently 0-5

void sortSuffixesBucket(char *G, void *ind, int indN, int indSkip)
{
    
    //boundaries of unsorted intervals
    int *uB1=new int[indN/2];
    int *uB2=new int[indN/2];
    int *uuB1=new int[indN/2];
    int *uuB2=new int[indN/2];
    
    int uN;
    
    if (false) //TODO implement flag to switch this option
    {//no initial bucketing
        uB1[0]=0;
        uB2[0]=indN;
        uN=1;//number of unsorted intervals
    } else
    {//initial bucketing according to the ind
        uN=0;
        uint64_t iprev=*(uint64_t*)(ind+0*indSkip);
        int un=1;
        for (int id=1; id<indN; id++)
        {
            if (*(uint64_t*)(ind+id*indSkip)==iprev)
            {
                un++;
            } else
            {
                if (un>1)
                {
                    uB1[uN]=id-un;
                    uB2[uN]=id;
                    uN++;
                    un=1;
                };
                iprev=*(uint64_t*)(ind+id*indSkip);
            };
        };
    };
    
    char *ind1=new char[indN*indSkip];//array to store sorted indices    
    int charShift=0;//character position to be sorted
    int charStart[GENOME_charN];
    int charCount[GENOME_charN];
    
    while (uN>0)
    {
        int uuN=0;
        for (int iu=0; iu<uN; iu++)
        {//sort iu-th interval
            memset( (void*) charCount, 0, sizeof(charCount));

            for (int id=uB1[iu]; id<uB2[iu]; id++)
            {//count the number of chars in each bucket
                charCount[G[*(uint64_t*)(ind+id*indSkip+8) + charShift]]++;//TODO can template uint64_t
            };

            charStart[0]=uB1[iu];

            bool needToSort=false;
            
            for (int ic=0;ic<GENOME_charN-1;ic++) //last char is the chromosome separator: 5=end of chromosome, if all chars are 5, no sorting required, since the indexes are presumed pre-sorted
            {//set new unsort bounds and charStart
                charStart[ic+1]=charStart[ic]+charCount[ic];
                if (charCount[ic]>0 && charCount[ic]<(uB2[iu]-uB1[iu]) ) needToSort=true; //need to sort only if at least one char > 0 and more than one chars in this bucket
                if (charCount[ic]>1)
                {
                    uuB1[uuN]=charStart[ic];
                    uuB2[uuN]=charStart[ic+1];                  
                    uuN++;
                };
            };

            if (needToSort)
            {//otherwise it's either all sorted or all unsorted, no need to move indexes in either case
                for (int id=uB1[iu]; id<uB2[iu]; id++)
                {//bucket-sort all indexes 
                    char c=G[*(uint64_t*)(ind+id*indSkip+8) + charShift];
                    memcpy(ind1+charStart[c]*indSkip, ind+id*indSkip, indSkip);
                    charStart[c]++;
                };
                memcpy(ind+uB1[iu]*indSkip, ind1+uB1[iu]*indSkip, (uB2[iu]-uB1[iu])*indSkip);
            };
        };
        //going to the next cycle
        charShift++;
        uN=uuN;
        int *p;
        p=uB1; uB1=uuB1; uuB1=p;
        p=uB2; uB2=uuB2; uuB2=p;
        std::cout << charShift <<"   " << uN <<"\n";
    };
};