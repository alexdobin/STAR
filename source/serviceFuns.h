#ifndef CODE_serviceFuns
#define CODE_serviceFuns

#include "IncludeDefine.h"

#include <string>
#include <sstream>
#include <vector>
int splitString(const std::string &s, char delim, std::vector<std::string> &elems);

int funCompareUint2 (const void *a, const void *b);
int funCompareUint1 (const void *a, const void *b);

//////////////////////////// TEMPLATES
template <class arrayType, int arraySize>
inline int funCompareArrays (const void *a, const void *b) {
    arrayType* va= (uint*) a;
    arrayType* vb= (uint*) b;

    for (int ii=0;ii<arraySize;ii++) {
        if (va[ii]>vb[ii]) {
            return 1;
        } else if (va[ii]<vb[ii]) {
            return -1;
        };
    };

    return 0;        

};

template <class argType>
inline uint32 binarySearch1(argType x, argType *X, uint32 N) {
    //binary search in the sorted list
    //check the boundaries first
    if (x>X[N-1] || x<X[0]) return -1;
    
    uint32 i1=0, i2=N-1, i3=N/2;
    while (i2>i1+1) {//binary search
        i3=(i1+i2)/2;
        if (X[i3]>x) {
            i2=i3;
        } else {
            i1=i3;
        };
    };
    
    while (i1<N-1 && x==X[i1+1]) ++i1; //go forward to check for equals
    return i1;
};

template <class argType>
inline int32 binarySearch1a(argType x, argType *X, int32 N) {
    //binary search in the sorted list
    //check the boundaries first
    //1a is different when x is larger than the last element X[N-1], it returns N-1
    
    if (x>X[N-1]) {
        return N-1;
    } else if (x<X[0]) {
        return -1;
    };
    
    int32 i1=0, i2=N-1, i3=N/2;
    while (i2>i1+1) {//binary search
        i3=(i1+i2)/2;
        if (X[i3]>x) {
            i2=i3;
        } else {
            i1=i3;
        };
    };
    
    while (i1<N-1 && x==X[i1+1]) ++i1; //go forward to check for equals
    return i1;
};

template <class argType>
inline int32 binarySearch1b(argType x, argType *X, int32 N) 
{
    //binary search in the sorted list
    //check the boundaries first
    //1b returns the first X element that is >= x
    //X are all distinct
    //if x>X[N-1], -1 is returned
    
    if (x>X[N-1]) {
        return -1;
    } else if (x<=X[0]) {
        return 0;
    };
    
    int32 i1=0, i2=N-1, i3=N/2;
    while (i2>i1+1) {//binary search
        i3=(i1+i2)/2;
        if (X[i3]>=x) {
            i2=i3;
        } else {
            i1=i3;
        };
    };
    
    return i2;
};

#endif