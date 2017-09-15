/*
 * inserts sequences into the SA and SAi indices
 */
#include "Genome.h"
#include "genomeScanFastaFiles.h"
#include "insertSeqSA.h"
#include "TimeFunctions.h"

void Genome::insertSequences()
{
if (pGe.gFastaFiles.at(0)!="-")
{
    time_t rawtime;
    time ( &rawtime );
    P.inOut->logMain  << timeMonthDayTime(rawtime) << " ..... inserting extra sequences into genome indexes" <<endl;
    //move the junctions to free up space for seqs
    // chrStart/Name/Length nChrReal include the extra sequences
    // nGenome is the old, small genome size
    uint sjdblen=nGenome-(chrStart.back()-genomeInsertL);//length of sjdb sequences
    memmove(G+chrStart.back(),G+chrStart.back()-genomeInsertL,sjdblen);
    memset(G+chrStart.back()-genomeInsertL, GENOME_spacingChar, genomeInsertL);//fill empty space with spacing characters

    genomeScanFastaFiles(P, G+chrStart.back()-genomeInsertL, true, *this); //read the seqs from file(s) into the free space
    uint64 nGenomeOld=nGenome;
    nGenome=chrStart.back()+sjdblen;
    //insert new sequences into the SA
    insertSeqSA(SA, SAinsert, SAi, G, G+chrStart.back()-genomeInsertL, nGenomeOld-sjdblen, genomeInsertL, sjdblen, P, *this);

    //insert new sequences into the SAi
    //update P
    //save the genome if necessary
};
};
