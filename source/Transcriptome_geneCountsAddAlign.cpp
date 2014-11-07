#include "Transcriptome.h"
#include "serviceFuns.cpp"

void Transcriptome::geneCountsAddAlign(uint nA, Transcript **aAll) {//add alignments from one read to gene counts
     if (nA>1) {
         quants->geneCounts.cMulti++;
     } else {
         Transcript& a=*aAll[0];//the unique alignment

         int32 gene1=-1;
         for (int ib=a.nExons-1; ib>=0; ib--) {//scan through all blocks of the alignments

             uint64 g1=a.exons[ib][EX_G]+a.exons[ib][EX_L]-1;//end of the block
             //search end of the block among the starts of exons
             int64 e1=-1;
             if (ib==a.nExons-1) {//binary search for the first time
                 e1=binarySearch1a<uint64>(g1, exG.s, (int32) exG.nEx);
             } else {//simple backwards scan 
                 while (e1>=0 && exG.s[e1]>a.exons[ib][EX_G]) {//stop when exon start is less than block end
                     --e1;
                 };
             };

             while (e1>=0 && exG.e[e1]>=a.exons[ib][EX_G]) {//these exons overlap this block
                 if (gene1==-1) {//first gene overlapping this read
                     gene1=exG.g[e1];
                 } else if (gene1!=exG.g[e1]) {//another gene overlaps this read
                     quants->geneCounts.cAmbig++;
                     return;
                 };//otherwise it's the same gene
                 --e1;
             };
         };

         if (gene1!=-1) {
              quants->geneCounts.uStr[a.Str][gene1]++;
         }; 


     };
};


