/* 
 * inserts sequences into the SA and SAi indices
 */
#include "Genome.h"
#include "genomeScanFastaFiles.h"
#include "insertSeqSA.h"
#include "TimeFunctions.h"

void Genome::insertSequences()
{
if (P->genomeFastaFiles.at(0)!="-")
{
    time_t rawtime;
    time ( &rawtime );
    P->inOut->logMain  << timeMonthDayTime(rawtime) << " ..... Inserting extra sequences into genome indexes" <<endl;       
    //move the junctions to free up space for seqs
    // chrStart/Name/Length nChrReal include the extra sequences
    // nGenome is the old, small genome size
    uint sjdblen=P->nGenome-(P->chrStart.back()-P->genomeInsertL);//length of sjdb sequences
    memmove(G+P->chrStart.back(),G+P->chrStart.back()-P->genomeInsertL,sjdblen);
    memset(G+P->chrStart.back()-P->genomeInsertL, GENOME_spacingChar, P->genomeInsertL);//fill empty space with spacing characters
        
    genomeScanFastaFiles(P, G+P->chrStart.back()-P->genomeInsertL, true); //read the seqs from file(s) into the free space
    uint64 nGenomeOld=P->nGenome;
    P->nGenome=P->chrStart.back()+sjdblen; 
    //insert new sequences into the SA
    insertSeqSA(SA, SAinsert, SAi, G, G+P->chrStart.back()-P->genomeInsertL, nGenomeOld-sjdblen, P->genomeInsertL, sjdblen, P);

    //insert new sequences into the SAi
    //update P
    //save the genome if necessary
};
};
