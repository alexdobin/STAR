#ifndef CODE_Transcriptome
#define CODE_Transcriptome

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"

class Transcriptome {
public:
    
    Parameters* Ptr; //transcriptomic parameters (i.e. chrName,...), to be used with RAtr for output
    Transcriptome (Parameters* Pin); //create transcriptome structure, load and initialize parameters
    uint32 quantAlign (Transcript &aG, Transcript *aTall);//transform coordinates for all aligns from genomic in RA to transcriptomic in RAtr
    vector <string> trID; //transcript IDs
    uint *trS, *trE, *trEmax; //transcripts start,end,end-max
   
    uint32 nTr, nEx; //number of transcript/exons
    uint16 *trExN; //number of exons per transcript
    uint32 *trExI; //index of the first exon for each transcript in exSE
    uint8 *trStr; //transcript strand
    uint32 *exSE; //exons start/end
    uint32 *exLenCum; //cumulative length of previous exons

private:
    Parameters* P; //normal "genomic" parameters
    
};

#endif
