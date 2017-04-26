#ifndef DEF_Chain
#define DEF_Chain

#include "IncludeDefine.h"
#include "Parameters.h"
#include "ErrorWarning.h"

class OneChain 
{
    public:
        uint bN;
        string chr1,chr2;//1/2 (old/new) chr names
        vector <uint> bStart1, bStart2, bLen; //blocks starts in 1/2, lengths
};

class Chain {
    public:
// //         uint bN;//number of blocks
// //         vector <uint> bStart1, bStart2, bLen; //blocks starts in 1/2, lengths
        
        Chain(Parameters &Pin, string chainFileNameIn);
        void liftOverGTF(string gtfFileName, string outFileName);
    private:
        Parameters &P;
        string chainFileName;
        void chainLoad();
        std::map <string,OneChain> chrChains;
};

#endif