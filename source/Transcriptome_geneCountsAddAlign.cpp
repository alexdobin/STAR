#include "Transcriptome.h"
#include "serviceFuns.cpp"

void Transcriptome::geneCountsAddAlign(uint nA, Transcript **aAll) {//add alignments from one read to gene counts
     if (nA>1) {
         quants->geneCounts.cMulti++;
     } else {
         Transcript& a=*aAll[0];//one unique alignment only

         vector<int32> gene1(quants->geneCounts.nType,-1);

         int64 e1=-1;

         for (int ib=a.nExons-1; ib>=0; ib--)
         {//scan through all blocks of the alignments

             uint64 g1=a.exons[ib][EX_G]+a.exons[ib][EX_L]-1;//end of the block

//              if ((uint)ib==a.nExons-1)
//              {//binary search for the first time: end of the block among the starts of exons
                 e1=binarySearch1a<uint64>(g1, exG.s, (int32) exG.nEx);
//              } else
//              {//simple backwards scan
//                  while (e1>=0 && exG.s[e1]>g1)
//                  {//stop when exon start is less than block end
//                      --e1;
//                  };
//              };

             while (e1>=0 && exG.eMax[e1]>=a.exons[ib][EX_G])
             {//these exons may overlap this block
                 if (exG.e[e1]>=a.exons[ib][EX_G])
                 {//this exon overlaps the block
                     uint str1=(uint)exG.str[e1]-1;
                     for (int itype=0; itype<quants->geneCounts.nType; itype++)
                     {
                         //str1<2 (i.e. strand=0) requirement means that genes w/o strand will accept reads from both strands
                         if ( itype==1 && a.Str!=str1 && str1<2) continue; //same strand
                         if ( itype==2 && a.Str==str1 && str1<2) continue; //reverse strand

                         if (gene1.at(itype)==-1)
                         {//first gene overlapping this read
                             gene1[itype]=exG.g[e1];
                         } else if (gene1.at(itype)==-2)
                         {
                             continue;//this align was already founf to be ambig for this strand
                         } else if (gene1.at(itype)!=(int32)exG.g[e1]) {//another gene overlaps this read
                             gene1[itype]=-2;//mark ambiguous
                         };//otherwise it's the same gene
                     };
                 };

                 --e1;// go to the previous exon
             };
         };

        for (int itype=0; itype<quants->geneCounts.nType; itype++)
        {
            if (gene1.at(itype)==-1)
            {
                quants->geneCounts.cNone[itype]++;
            } else if (gene1.at(itype)==-2)
            {
                quants->geneCounts.cAmbig[itype]++;
            } else
            {
                quants->geneCounts.gCount[itype][gene1.at(itype)]++;
            };
         };
     };
};


