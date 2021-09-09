#include "Transcriptome.h"
#include "serviceFuns.cpp"
#include "ReadAnnotations.h"

void Transcriptome::geneFullAlignOverlap_CR(uint nA, Transcript **aAll, int32 strandType, ReadAnnotations &readAnnot)
{
    readAnnot.geneFull_Ex50pAS={};
        
    for (uint32 iA=0; iA<nA; iA++) {//only includes aligns that are entirely inside genes (?)
        Transcript &a = *aAll[iA];
            
        uint64 aS = a.exons[0][EX_G]; //align start
        uint64 aE = a.exons[a.nExons-1][EX_G] + a.exons[a.nExons-1][EX_L]-1; //align end
        //TODO: for paired end, need to look fo max min
            
        int64 gi1=binarySearch1a<uint64>(aS, geneFull.s, (int32) nGe); //search align-start against gene-starts. Find last gene which still starts to the left of align-start
                
        while (gi1>=0 && geneFull.eMax[gi1]>=aE) {//these genes may overlap this block
            if (geneFull.e[gi1]>=aE) {//this gene contains the block:  gene-end is to the right of block start
                int32 str1 = geneFull.str[gi1]==1 ? a.Str : 1-a.Str;
                if (strandType==-1 || strandType==str1)  {
                    readAnnot.geneFull_Ex50pAS.insert(geneFull.g[gi1]);
                    readAnnot.geneFull_Ex50pAS_Tr=iA;
                };
            };
            --gi1;// go to the previous gene
        };
    };
};


