#ifndef GENOME_DEF
#define GENOME_DEF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "PackedArray.h"
#include "SharedMemory.h"

class Genome {
    public:
        char *G, *sigG;
        PackedArray SA,SA2;
        PackedArray SAi;
        
        Genome (Parameters* Pin );

        Genome () {};//empty constructor
        void freeMemory();
        void genomeLoad();

    private:
	Parameters* P;
    key_t shmKey;  
    char *shmStart;
    char *G1; //pointer -200 of G
    SharedMemory * sharedMemory = NULL;
    uint OpenStream(string name, ifstream & stream);
};
#endif
