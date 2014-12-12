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
        
        Genome (Parameters* Pin );

        Genome () {};//empty constructor
        ~Genome ();
        void freeMemory();
        void genomeLoad();

    private:
	Parameters* P;
        key_t shmKey;  
        char *shmStart;
        char *G1; //pointer -200 of G
        int shmID;

        bool GetSharedObjectByKey(key_t shmKey, int * shmID);
        int CreateSharedObject(key_t shmKey, uint64 shmSize);
        int SharedObjectsUseCount(int shmID);
        void * MapSharedObjectToMemory(int shmID);
        const char * GetPosixObjectKey(key_t shmKey);
        struct stat GetSharedObjectInfo(int shmID);
        void RemoveSharedObject(int shmID, void * * ptr, key_t shmKey);
};
#endif
