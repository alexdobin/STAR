#ifndef DEF_serviceFuns
#define DEF_serviceFuns

#include "IncludeDefine.h"

template <class T>
    T sum1D(T* a, uint N) {
        T s=0;
        for (uint ii=0;ii<N;ii++) s+=a[ii];
        return s;
    };

inline int funCompareUint1 (const void *a, const void *b) {
    uint va= *((uint*) a);
    uint vb= *((uint*) b);

    if (va>vb) {
		return 1;
	} else if (va==vb) {
		return 0;        
	} else {
		return -1;
	};
};
    
inline int funCompareUint2 (const void *a, const void *b) {
    uint va= *((uint*) a);
    uint vb= *((uint*) b);
    uint va1=*(((uint*) a)+1);
    uint vb1=*(((uint*) b)+1);

    if (va>vb) {
		return 1;
	} else if (va==vb && va1>vb1) {
		return 1;
	} else if (va==vb && va1==vb1) {
		return 0;        
	} else {
		return -1;
	};
};

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


#endif

