#ifndef OUT_SJ_DEF
#define OUT_SJ_DEF

#include "Parameters.h"

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
    
    void junctionPointer(char* sjPoint, uint isj);
    void outputStream(ostream &outStream, Parameters* P);
    void collapseOneSJ(char* isj1P, char* isjP, Parameters* P);
};

class OutSJ {

public:
    //all junctions
    char* data; //sj array[Njunctions][dataSize]
    
    uint N; //number of junctions stored
    
    Junction oneSJ;
    
    Parameters *P;
       
    OutSJ(uint nSJmax, Parameters *P);
        
    void collapseSJ();//collapse the junctions in data
//     int compareSJ(void* i1, void* i2);
    
};

int compareSJ(const void* i1, const void* i2); //external functions

#endif

