#include "binarySearch2.h"

int binarySearch2(uint x, uint y, uint *X, uint *Y, int N) {
    //binary search in the sorted list to find the junction
    //check the boundaries first
    if (N==0 || x>X[N-1] || x<X[0]) return -1;

    int i1=0, i2=N-1, i3=N/2;
    while (i2>i1+1) {//binary search
        i3=(i1+i2)/2;
        if (X[i3]>x) {
            i2=i3;
        } else {
            i1=i3;
        };
    };

    if (x==X[i1]) {
        i3=i1;
    } else if (x==X[i2]) {
        i3=i2;
    } else {
        return -1;
    };

    for (int jj=i3;jj>=0;jj--) {//go back
        if (x!=X[jj]) {
            break;//next try forward
        } else if (y==Y[jj]) {
            return jj;
        };
    };

    for (int jj=i3;jj<N;jj++) {//go forward
        if (x!=X[jj]) {
            return -1;
        } else if (y==Y[jj]) {
            return jj;
        };
    };

    return -2; //this will happen if jj went past N
};
