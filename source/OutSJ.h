#ifndef CODE_OutSJ
#define CODE_OutSJ

#include "Parameters.h"
#include "Genome.h"

class Junction {//one junction
public:
    const static uint startP=0;
    const static uint gapP=startP+sizeof(uint);
    const static uint strandP=gapP+sizeof(uint32);
    const static uint motifP=strandP+sizeof(char);
    const static uint annotP=motifP+sizeof(char);
    const static uint countUniqueP=annotP+sizeof(char);
    const static uint countMultipleP=countUniqueP+sizeof(uint32);
    const static uint overhangLeftP=countMultipleP+sizeof(uint32);
    const static uint overhangRightP=overhangLeftP+sizeof(uint16);

    uint *start;
    uint32 *gap;
    char *strand, *motif, *annot;
    uint32 *countUnique, *countMultiple;
    uint16 *overhangLeft, *overhangRight;

    const static uint dataSize=overhangRightP+sizeof(uint16);

    Junction(Genome &genomeIn);
    void junctionPointer(char* sjPoint, uint isj);
    void outputStream(ostream &outStream);
    void collapseOneSJ(char* isj1P, char* isjP, Parameters& P);
    
private:
    Genome &mapGen;
};

class OutSJ {

public:
    //all junctions
    char* data; //sj array[Njunctions][dataSize]
    uint N; //number of junctions stored
    Junction oneSJ;

    OutSJ(uint nSJmax, Parameters &Pin, Genome &genomeIn);
    void collapseSJ();//collapse the junctions in data
//     int compareSJ(void* i1, void* i2);

private:
    Parameters &P;
    Genome &mapGen;  
};

int compareSJ(const void* i1, const void* i2); //external functions

#endif

