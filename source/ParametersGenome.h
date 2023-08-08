#ifndef CODE_ParametersGenome
#define CODE_ParametersGenome

#include "IncludeDefine.h"
#include <unordered_set>

class Parameters;

class ParametersGenome {//"constant" genome parameters - user input
public:
    string gDir;
    string gLoad;
    
    uint32 gType;//type code
    string gTypeString;
    
    vector <string> gFastaFiles;
    vector <string> gChainFiles;
    //string gConsensusFile; DEPRECATED

    struct {
        int32 type;
        string typeString;
        string vcfFile;
        vector<string> output; //which output to transform
        bool outYes, outSAM, outSJ, outQuant;
    } transform;
    
    uint gSAindexNbases;//length of the SA pre-index strings
    uint gChrBinNbits;
    uint gSAsparseD;//SA sparsity
    uint gSuffixLengthMax;//maximum length of the suffixes, has to be longer than read length
    vector <uint> gFileSizes;//size of the genome files

    vector <string> sjdbFileChrStartEnd;
    string sjdbGTFfile;
    string sjdbGTFchrPrefix;
    
    string sjdbGTFfeatureExon;
    string sjdbGTFtagExonParentTranscript;
    string sjdbGTFtagExonParentGene;
    vector<string> sjdbGTFtagExonParentGeneName;
    vector<string> sjdbGTFtagExonParentGeneType;
    
    string sjdbInsertSave;
    uint sjdbOverhang;
    int sjdbOverhang_par;
    int sjdbScore;

    struct {
        vector<string> mitoStrings;
        unordered_set<uint64> mito;
    } chrSet;

    
    void initialize(Parameters *Pin);

private:
    Parameters *pP;
    
};

#endif
