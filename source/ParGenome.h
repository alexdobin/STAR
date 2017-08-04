#ifndef DEF_ParGenome
#define DEF_ParGenome

class ParGenome {//"constant" genome parameters
    public:
        vector <uint> chrStart, chrLength, chrLengthAll;
        string genomeDir,genomeLoad;
        vector <string> genomeFastaFiles, genomeChainFiles;
        uint genomeSAsparseD;//sparsity=distance between indices
        uint genomeInsertL; //total length of the sequence to be inserted on the fly
        uint genomeInsertChrIndFirst; //index of the first inserted chromosome
        uint genomeSuffixLengthMax; //maximum length of the suffixes, has to be longer than read length
        vector <uint> genomeFileSizes; //size of the genome files

        uint genomeChrBinNbits, genomeChrBinNbases, chrBinN, *chrBin;
        
        uint genomeSAindexNbases; //length of the SA pre-index strings
        uint *genomeSAindexStart;//starts of the L-mer indices in the SAindex, 1<=L<=genomeSAindexNbases

       //SJ database parameters
        vector <string> sjdbFileChrStartEnd;
        string sjdbGTFfile, sjdbGTFchrPrefix, sjdbGTFfeatureExon, sjdbGTFtagExonParentTranscript, sjdbGTFtagExonParentGene;
        uint sjdbOverhang,sjdbLength; //length of the donor/acceptor, length of the sj "chromosome" =2*sjdbOverhang+1 including spacer
        int sjdbOverhang_par;
        int sjdbScore;

        uint sjChrStart,sjdbN; //first sj-db chr
        uint sjGstart; //start of the sj-db genome sequence
        uint *sjDstart,*sjAstart,*sjStr, *sjdbStart, *sjdbEnd; //sjdb loci
        uint8 *sjdbMotif; //motifs of annotated junctions
        uint8 *sjdbShiftLeft, *sjdbShiftRight; //shifts of junctions
        uint8 *sjdbStrand; //junctions strand, not used yet

        //Genome parameters
        uint nGenome, nSA, nSAbyte, nChrReal;//genome length, SA length, # of chromosomes, vector of chromosome start loci
        uint nGenome2, nSA2, nSAbyte2, nChrReal2; //same for the 2nd pass
        uint nSAi; //size of the SAindex
        vector <string> chrName, chrNameAll;
        map <string,uint> chrNameIndex;
        unsigned char GstrandBit, SAiMarkNbit, SAiMarkAbsentBit; //SA index bit for strand information
        uint GstrandMask, SAiMarkAbsentMask, SAiMarkAbsentMaskC, SAiMarkNmask, SAiMarkNmaskC;//maske to remove strand bit from SA index, to remove mark from SAi index


};

#endif