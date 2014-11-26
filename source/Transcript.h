#ifndef TRANSCRIPT_DEF
#define TRANSCRIPT_DEF

#include "IncludeDefine.h"
#include "Parameters.h"

class Transcript {
    public:
        //uint rMM[DEF_readSeqLengthMax]; //read-space MM coordinates. TODO: remove? reduce the size, do not save the whole array
//         uint gMap[DEF_readSeqLengthMax]; //map read bases to genome. TODO: allocate with smaller size
        uint exons[MAX_N_EXONS][EX_SIZE]; //coordinates of all exons: r-start, g-start, length
        uint shiftSJ[MAX_N_EXONS][2]; //shift of the SJ coordinates due to genomic micro-repeats
        int canonSJ[MAX_N_EXONS]; //canonicity of each junction
        uint8 sjAnnot[MAX_N_EXONS]; //anotated or not
        uint8 sjStr[MAX_N_EXONS]; //anotated or not
        
        uint intronMotifs[3];
        uint8 sjMotifStrand;
        
        uint nExons; //number of exons in the read transcript
        
        uint Lread, readLengthPairOriginal;
        uint iRead; //read identifier
        
        int iFrag; //frag number of the transcript, if the the transcript contains only one frag
        
        //loci
        uint rStart, roStart, rLength, gStart, gLength, cStart; //read, original read, and genomic start/length, chromosome start
        uint Chr,Str,roStr; //chromosome and strand and original read Strand
        
        bool primaryFlag;
       
//         uint polyXlength[4], polyXnMM[4] ;
        uint nWAmax; //number of genes, max number of alignments per gene, total number of alignments for this reads
        
        uint nMatch;//min number of matches    
        uint nMM;//max number of mismatches        
        uint mappedLength; //total mapped length, sum of lengths of all blocks(exons)
        
        uint extendL; //extension length        
        intScore maxScore; //maximum Score
        intScore nextTrScore; //next after maximum Tr Score
        

        uint nGap, lGap; //number of genomic gaps (>alignIntronMin) and their total length
        uint nDel; //number of genomic deletions (ie genomic gaps)
        uint nIns; //number of (ie read gaps)        
        uint lDel; //total genomic deletion length
        uint lIns; //total genomic insertion length
        
        uint nUnique, nAnchor; //number of unique pieces in the alignment, number of anchor pieces in the alignment   
        
        Transcript(); //resets to 0
        void reset(); //reset to 0
        void resetMapG(); // reset map to 0
        void resetMapG(uint); // reset map to 0 for Lread bases
        void add(Transcript*); // add
        void alignScore(char **Read1, char *G, Parameters *P);
};

#endif
