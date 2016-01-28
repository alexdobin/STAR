#ifndef GENOME_DEF
#define GENOME_DEF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "PackedArray.h"

#if !defined(_WIN32) && defined(USE_UNIX_SHM)
#include "SharedMemory.h"
#else
#include "SharedMemorySegment.h"
#endif

#include "CrossPlatform.h"

class Genome {
    public:
        char *G, *sigG;
        PackedArray SA,SAinsert,SApass1,SApass2;
        PackedArray SAi;
        
        uint nGenomeInsert, nGenomePass1, nGenomePass2, nSAinsert, nSApass1, nSApass2; 

        Genome (Parameters* Pin );
        Genome () {};//empty constructor
        ~Genome();
        
        void freeMemory();
        void genomeLoad();

        void insertSequences();


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
