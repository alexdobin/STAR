#include "Transcriptome.h"
#include "serviceFuns.cpp"

void Transcriptome::geneFullAlignOverlap(uint nA, Transcript **aAll, int32 strandType, set<uint32> &geneOverlap)
{
     for (uint32 iA=0; iA<nA; iA++) {
         Transcript &a = *aAll[iA];//one unique alignment only

         int64 gi1=-1;

         for (int64 ib=a.nExons-1; ib>=0; ib--) {//scan through all blocks of the alignments

             uint64 be1=a.exons[ib][EX_G]+a.exons[ib][EX_L]-1;//end of the block
             gi1=binarySearch1a<uint64>(be1, geneFull.s, (int32) nGe);

             while (gi1>=0 && geneFull.eMax[gi1]>=a.exons[ib][EX_G]) {//these exons may overlap this block
                 if (geneFull.e[gi1]>=a.exons[ib][EX_G]) {//this gene overlaps the block
                     int32 str1 = geneFull.str[gi1]==1 ? a.Str : 1-a.Str;
                     if (strandType==-1 || strandType==str1)
                         geneOverlap.insert(geneFull.g[gi1]);
                 };
                 --gi1;// go to the previous gene
             };
         };
     };
};


