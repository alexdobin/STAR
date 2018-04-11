#include "Transcript.h"
#include "serviceFuns.cpp"

int Transcript::variationAdjust(const Genome &mapGen, char *R)
{
    Variation &Var=*mapGen.Var;
    
    if (!Var.yes)
    {//no variation
        return 0;
    };
    
    int dScore=0;//change in the score
    uint nMM1=0;
    
    //for each block, check whether it overlaps one or more SNPs
    for (uint ie=0; ie<nExons; ie++)
    {
        //binary search to find nearest SNP
        int32 isnp=binarySearch1b <uint> (exons[ie][EX_G], Var.snp.loci, Var.snp.N);
        if (isnp>=0)
        {
            while ((uint)isnp<Var.snp.N && exons[ie][EX_G]+exons[ie][EX_L]>Var.snp.loci[isnp])
            {//these SNPs overlap the block
                varInd.push_back(isnp); //record snp index
                varGenCoord.push_back(Var.snp.loci[isnp]-mapGen.chrStart[Chr]);
                
                varReadCoord.push_back(exons[ie][EX_R]+Var.snp.loci[isnp]-exons[ie][EX_G]);
                char ntR=R[varReadCoord.back()];//nt of the read in the SNP position, already trnasformed to + genome strand
                
                uint8 igt;
                if (ntR>3) {
                    igt=4;
                } else {
                    for (igt=1; igt<3; igt++) {//1st or 2nd allele, =3 of none
                        if (Var.snp.nt[isnp][igt]==ntR) {
                            break;
                        };
                    };
                };

                //if (ntR == Var.snp.nt[isnp][0])
                //{//mark snp that agrees with the reference
                //    igt*=10;
                //};

                varAllele.push_back(igt);
                
                if (igt<3 && ntR != Var.snp.nt[isnp][0])
                {//non-reference allele, correct nMM and score
                    ++nMM1;
                };

                ++isnp;
            };
        };
    };
    
    #define VAR_noScoreCorrection
    #ifndef VAR_noScoreCorrection
    if (nMM1>0)
    {//one or more mismtaches need to be corrected
        uint nMMold=nMM;       
        alignScore(Read1, G, P);
        nMM-=nMM1;
        nMatch+=nMM1;
        dScore=2*(nMMold-nMM);//score only changes if the number of mismatches is reduced after SNP adjustment
    };
    #else
        //#warning VAR_noScoreCorrection set: no variation score correction
    #endif
    
    return dScore;
};
