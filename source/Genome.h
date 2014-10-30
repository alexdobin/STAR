#ifndef GENOME_DEF
#define GENOME_DEF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "PackedArray.h"
class Genome {
    public:
        char *G, *sigG;
        PackedArray SA,SA2;
        PackedArray SAi;
        
        Genome (Parameters* Pin ) : P(Pin) {};
        Genome () {};//empty constructor
        ~Genome ();
        void freeMemory();
        void genomeLoad();

    private:
        char *G1; //pointer -200 of G
        Parameters* P;
        int shmID;
};
#endif
