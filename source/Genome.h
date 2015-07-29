#ifndef GENOME_DEF
#define GENOME_DEF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "PackedArray.h"

#ifdef _WIN32
#include "SharedMemorySegment.h"
#else
#include "SharedMemory.h"
#endif

#include "CrossPlatform.h"

class Genome {
    public:
        char *G, *sigG;
        PackedArray SA,SApass1,SApass2;
        PackedArray SAi;
        
        uint nGenomePass1, nGenomePass2, nSApass1, nSApass2; 

        Genome (Parameters* Pin );
        Genome () {};//empty constructor
        ~Genome();
        
        void freeMemory();
        void genomeLoad();

    private:
    Parameters* P;
    key_t shmKey;  
    char *shmStart;
    char *G1; //pointer -200 of G

#if !defined(_WIN32) && defined(USE_UNIX_SHM)
    SharedMemory * sharedMemory;
#else
	SharedMemorySegment* sharedMemory; 
#endif

    uint OpenStream(string name, ifstream & stream);
    void HandleSharedMemoryException(const SharedMemoryException & exc, uint64 shmSize);
};
#endif
