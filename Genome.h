#ifndef GENOME_DEF
#define GENOME_DEF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "PackedArray.h"
class Genome {
    public:
        char *G, *sigG;
        PackedArray SA;
        PackedArray SAi;
        void genomeLoad();
        
        Genome (Parameters* Pin ) : P(Pin) {};
        Genome () {};
        ~Genome();
    private:
        char *G1; //pointer -200 of G
        Parameters* P;
        int shmID;
};
#endif
