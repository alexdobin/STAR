#ifndef H_FastResetVector
#define H_FastResetVector

#include <vector>

// Implementation of a vector that can be reset to a default value in ~O(1).
// It is more efficient than an ordinary vector if:
// - the values are sparse, such that only a few elements need to be deleted.
// - there are not too many accesses to elements, since upon each access
//   it must be checked if the element is stale.

template <class T> class FastResetVector {

    private:
        vector<T> data; // contains actual data to be stored
        T defaultValue; // all elements of `data` are initialized with this value
        unsigned int incarnation = 0; // increasing the incarnation invalidates all elements of `data` (=reset)
        vector<unsigned int> lastUpdate; // for each element in `data`, keep track of the incarnation that it was last updated

    public:
        FastResetVector(const size_t s, const T& d): data(s,d), defaultValue(d), lastUpdate(s) {}

        inline T& operator[](const size_t i) { // whenever an element is accessed, check if it's stale
            if (incarnation != lastUpdate[i]) { // is it stale?
                data[i] = defaultValue; // reset to defaut value
                lastUpdate[i] = incarnation; // mark as fresh for this incarnation
            };
            return data[i];
        }

        void reset() {
            incarnation++; // we can invalidate `data` simply through "reincarnation"
            if (incarnation == 0) // only when there is an integer overflow, we need to reinitialize
                fill(data.begin(), data.end(), defaultValue);
        }

};

#endif

