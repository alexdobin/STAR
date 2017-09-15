#ifndef GENOME_DEF
#define GENOME_DEF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "PackedArray.h"
#include "SharedMemory.h"
#include "Variation.h"

class Genome {
    public:
        char *G, *sigG;
        PackedArray SA,SAinsert,SApass1,SApass2;
        PackedArray SAi;
        Variation *Var;
        
        uint nGenomeInsert, nGenomePass1, nGenomePass2, nSAinsert, nSApass1, nSApass2;
        
        ParametersGenome &pGe;
        
        //chr parameters
        vector <uint> chrStart, chrLength, chrLengthAll;        
        uint genomeChrBinNbases, chrBinN, *chrBin;
        vector <string> chrName, chrNameAll;
        map <string,uint> chrNameIndex;
        
        uint *genomeSAindexStart;//starts of the L-mer indices in the SAindex, 1<=L<=pGe.gSAindexNbases

        uint nGenome, nSA, nSAbyte, nChrReal;//genome length, SA length, # of chromosomes, vector of chromosome start loci
        uint nGenome2, nSA2, nSAbyte2, nChrReal2; //same for the 2nd pass
        uint nSAi; //size of the SAindex
        unsigned char GstrandBit, SAiMarkNbit, SAiMarkAbsentBit; //SA index bit for strand information
        uint GstrandMask, SAiMarkAbsentMask, SAiMarkAbsentMaskC, SAiMarkNmask, SAiMarkNmaskC;//maske to remove strand bit from SA index, to remove mark from SAi index
   
        //SJ database parameters
        uint sjdbOverhang, sjdbLength; //length of the donor/acceptor, length of the sj "chromosome" =2*pGe.sjdbOverhang+1 including spacer
        uint sjChrStart,sjdbN; //first sj-db chr
        uint sjGstart; //start of the sj-db genome sequence
        uint *sjDstart,*sjAstart,*sjStr, *sjdbStart, *sjdbEnd; //sjdb loci
        uint8 *sjdbMotif; //motifs of annotated junctions
        uint8 *sjdbShiftLeft, *sjdbShiftRight; //shifts of junctions
        uint8 *sjdbStrand; //junctions strand, not used yet        
        
       //sequence insert parameters
        uint genomeInsertL; //total length of the sequence to be inserted on the fly
        uint genomeInsertChrIndFirst; //index of the first inserted chromosome

        Genome (Parameters &Pin );
        ~Genome();

        void freeMemory();
        void genomeLoad();
        void chrBinFill();
        void chrInfoLoad();
        
        void insertSequences();

        void genomeGenerate();
        
    private:
        Parameters &P;
        key_t shmKey;
        char *shmStart;
        char *G1; //pointer -200 of G
        SharedMemory * sharedMemory;
        uint OpenStream(string name, ifstream & stream, uint size);
        void HandleSharedMemoryException(const SharedMemoryException & exc, uint64 shmSize);

};
#endif
