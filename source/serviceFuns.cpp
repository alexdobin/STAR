#ifndef DEF_serviceFuns
#define DEF_serviceFuns

#include "IncludeDefine.h"

template <class T>
    T sum1D(T* a, uint N) {
        T s=0;
        for (uint ii=0;ii<N;ii++) s+=a[ii];
        return s;
    };

template <class numT> int funCompareNumbers (const void *a, const void *b) {
    numT va= *((numT*) a);
    numT vb= *((numT*) b);

    if (va>vb) {
        return 1;
    } else if (va==vb) {
        return 0;
    } else {
        return -1;
    };
};

template <class numT> int funCompareNumbersReverse (const void *a, const void *b) {
    numT va= *((numT*) a);
    numT vb= *((numT*) b);

    if (va>vb) {
        return -1;
    } else if (va==vb) {
        return 0;
    } else {
        return 1;
    };
};

template <class numT, int Shift> int funCompareNumbersReverseShift (const void *a, const void *b) {
    numT va= *((numT*) a + Shift);
    numT vb= *((numT*) b + Shift);

    if (va>vb) {
        return -1;
    } else if (va==vb) {
        return 0;
    } else {
        return 1;
    };
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
    arrayType* va= (arrayType*) a;
    arrayType* vb= (arrayType*) b;

    for (int ii=0;ii<arraySize;ii++) {
        if (va[ii]>vb[ii]) {
            return 1;
        } else if (va[ii]<vb[ii]) {
            return -1;
        };
    };

    return 0;

};

template <class arrayType, int arraySize>
inline int funCompareArraysReverse (const void *a, const void *b) {
    arrayType* va= (arrayType*) a;
    arrayType* vb= (arrayType*) b;

    for (int ii=0;ii<arraySize;ii++) {
        if (va[ii]>vb[ii]) {
            return -1;
        } else if (va[ii]<vb[ii]) {
            return 1;
        };
    };

    return 0;

};

template <class arrayType, int arraySize, int Shift>
inline int funCompareArraysShift (const void *a, const void *b) {
    arrayType* va= ((arrayType*) a) + Shift;
    arrayType* vb= ((arrayType*) b) + Shift;

    for (int ii=0;ii<arraySize;ii++) {
        if (va[ii]>vb[ii]) {
            return 1;
        } else if (va[ii]<vb[ii]) {
            return -1;
        };
    };

    return 0;

};

template <class Type>
inline int funCompareTypeSecondFirst (const void *a, const void *b) {
    Type va= *( ((Type*) a) + 1 );
    Type vb= *( ((Type*) b) + 1 );
    Type va1= *( ((Type*) a) + 0 );
    Type vb1= *( ((Type*) b) + 0 );

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

template <class Type, int Shift>
inline int funCompareTypeShift (const void *a, const void *b) {
    Type va= *( ((Type*) a)+Shift );
    Type vb= *( ((Type*) b)+Shift );

    if (va>vb) {
        return 1;
    } else if (va==vb) {
        return 0;
    } else {
        return -1;
    };

};

inline int splitString(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    int maxL=0;
    elems.clear();
    while (std::getline(ss, item, delim)) {
        maxL=max(maxL, (int)item.size());
        elems.push_back(item);
    };
    return maxL;//returns max string size
};

/*
 * std::vector<std::string>  splitString(const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> elems;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    };
    return elems;
};
*/

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
inline bool binarySearch_leLeft(argType x, argType *X, uint32 N, uint32 & i1) {
    //binary search in the sorted list
    //Retrun false if x is outside of X[0], X[N-1]
    //Return i1 = index of the element which is <= x, the leftmost element if equal.
    
    if (x>X[N-1] || x<X[0]) 
        return false;

    i1=0;
    uint32 i2=N-1, i3=N/2;
    while (i2>i1+1) {//binary search
        i3=(i1+i2)/2;
        if (X[i3]>x) {
            i2=i3;
        } else {// X[i3] <= x
            i1=i3;
        };
    };

    while (i1>0 && x==X[i1-1]) 
        --i1; //go left to check for equals
        
    return true;
};


template <class argType>
inline int32 binarySearch1a(argType x, argType *X, int32 N) {
    //binary search in the sorted list
    //check the boundaries first
    //1a returns the last X element that is <= x

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

    while (i1<N-1 && x==X[i1+1]) 
        ++i1; //go forward to check for equals
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

template <class argType>
inline int64 binarySearchExact(argType x, argType *X, uint64 N) {
    //binary search in the sorted list
    //check the boundaries first
    //returns -1 if no match found
    //if X are not all distinct, no guarantee which element is returned

    if (x>X[N-1] || x<X[0])
        return -1;


    int32 i1=0, i2=N-1, i3=N/2;
    while (i2>i1+1) {//binary search
        i3=(i1+i2)/2;
        if (X[i3]>=x) {
            i2=i3;
        } else {
            i1=i3;
        };
    };

    if (x==X[i2]) {
        return i2;
    } else if (x==X[i1]) {
        return i1;
    } else {
        return -1;
    };
};

// template <class argType>
// inline uint64 binarySearchStride(argType *X, uint64 N, artgType x, iStart, stride )
// {
//     //binary search in the sorted list wit stride
//     //1b returns the first X element that is >= x
//     //X are all distinct
//     //if x>X[N-1], -1 is returned
// 
//     if (x>X[N-1]) {
//         return -1;
//     } else if (x<=X[0]) {
//         return 0;
//     };
// 
//     int32 i1=0, i2=N-1, i3=N/2;
//     while (i2>i1+1) {//binary search
//         i3=(i1+i2)/2;
//         if (X[i3]>=x) {
//             i2=i3;
//         } else {
//             i1=i3;
//         };
//     };
// 
//     return i2;
// };


#endif
