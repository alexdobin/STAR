#ifndef CODE_Variation
#define CODE_Variation

#include "IncludeDefine.h"
#include "Parameters.h"
#include <array>

// struct SNPnt
// {
//     char ref;
//     char a1;
//     char a2;
// };

class SNP
{
public:
    uint32 N; //number of snps
    uint *loci; //snp coordinate
//     SNPnt* nt; //reference and alternative bases
//     char **nt; //reference and alternative bases
//     char *nt1; //1D array to store nt
    vector<array<char,3>> nt;//reference and alternative bases
};

class Variation
{
public:
    //methods
    Variation (Parameters* Pin); //create transcriptome structure, load and initialize parameters
    void loadVCF(string fileIn); //load VCF file
    
    //variables
    bool yes;    
    SNP snp;
    
    Parameters* P;
    
private:
    string vcfFile;
    string varOutFileName;
    ofstream varOutStream;//output file for variations
};

#endif