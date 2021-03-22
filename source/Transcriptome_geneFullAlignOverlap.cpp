#include "Transcriptome.h"
#include "serviceFuns.cpp"
#include "ReadAnnotations.h"

void Transcriptome::geneFullAlignOverlap(uint nA, Transcript **aAll, int32 strandType, ReadAnnotations &readAnnot)
{
    readAnnot.geneFull={};
        
    for (uint32 iA=0; iA<nA; iA++) {
        Transcript &a = *aAll[iA];//one unique alignment only
        
        /*
        for (int64 ib=a.nExons-1; ib>=0; ib--) {//scan through all blocks of the alignments
            
            uint64 be1=a.exons[ib][EX_G]+a.exons[ib][EX_L]-1;//end of the block
            int64 gi1=binarySearch1a<uint64>(be1, geneFull.s, (int32) nGe); //search block-end against gene-starts. Find last gene which still starts to the left of block-end
                
            while (gi1>=0 && geneFull.eMax[gi1]>=a.exons[ib][EX_G]) {//these genes may overlap this block
                if (geneFull.e[gi1]>=a.exons[ib][EX_G]) {//this gene overlaps the block:  gene-end is to the right of block start
                    int32 str1 = geneFull.str[gi1]==1 ? a.Str : 1-a.Str;
                    if (strandType==-1 || strandType==str1)  {
                        readAnnot.geneFull.insert(geneFull.g[gi1]);
                        readAnnot.geneFullTr=iA;
                    };
                };
                --gi1;// go to the previous gene
            };
        };
        */
        
        uint64 aS = a.exons[0][EX_G]; //align start
        uint64 aE = aS + a.exons[a.nExons-1][EX_L]-1; //align end
        //TODO: for paired end, need to look fo max min
        
        int64 gi1=binarySearch1a<uint64>(aS, geneFull.s, (int32) nGe); //search align-start against gene-starts. Find last gene which still starts to the left of align-start
                
        while (gi1>=0 && geneFull.eMax[gi1]>=aE) {//these genes may overlap this block
            if (geneFull.e[gi1]>=aE) {//this gene contains the block:  gene-end is to the right of block start
                int32 str1 = geneFull.str[gi1]==1 ? a.Str : 1-a.Str;
                if (strandType==-1 || strandType==str1)  {
                    readAnnot.geneFull.insert(geneFull.g[gi1]);
                    readAnnot.geneFullTr=iA;
                };
            };
            --gi1;// go to the previous gene
        };
    };
};


