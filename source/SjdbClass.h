#ifndef DEF_SjdbClass
#define DEF_SjdbClass

#include "IncludeDefine.h"
#include <set>

class SjdbClass {
public:
    vector <string> chr;
    vector <uint64> start,end;
    vector <char> str;
    vector <uint8> priority;
    
    vector<set<uint64>> gene;
};

#endif
