#include "Genome.h"
#include "TimeFunctions.h"
#include "genomeParametersWrite.h"
#include "streamFuns.h"
#include "GenomeInsertAnnotations.h"

namespace {
void writeAnnotationSidecars(Genome &genome, const string &dirOut)
{
    if (!(genome.P.runMode=="genomeInsert" && genome.P.sjdbInsert.pass1 && genome.pGe.sjdbGTFfile!="-")) {
        genomeInsertCopyAnnotationSidecars(genome.pGe.gDir, dirOut);
        return;
    };

    genomeInsertMergeAnnotationSidecars(genome.pGe.gDir, genome.P.sjdbInsert.outDir, dirOut, genome.P, true);
}
}

void Genome::writeGenomeIndex(const string dirOut)
{
    string dirOut1=dirOut;
    if (dirOut1.back()=='/') {
        dirOut1.erase(dirOut1.end()-1);
    };

    createDirectory(dirOut1+"/", P.runDirPerm, "--genomeInsertOutDir", P);

    writeAnnotationSidecars(*this, dirOut1);

    writeChrInfo(dirOut1);

    pGe.gFileSizes.clear();
    pGe.gFileSizes.push_back(nGenome);
    pGe.gFileSizes.push_back(SA.lengthByte);

    // Record both the loaded base FASTA files and the inserted FASTA files in
    // genomeParameters.txt, while leaving the in-memory parameters unchanged.
    vector<string> genomeFastaFilesSaved=pGe.gFastaFiles;
    if (genomeInsertL>0 && genomeFastaFilesSaved.at(0)!="-" && genomeFastaFilesLoaded.size()>0) {
        pGe.gFastaFiles=genomeFastaFilesLoaded;
        if (pGe.gFastaFiles.size()==1 && pGe.gFastaFiles.at(0)=="-") {
            pGe.gFastaFiles.clear();
        };
        pGe.gFastaFiles.insert(pGe.gFastaFiles.end(), genomeFastaFilesSaved.begin(), genomeFastaFilesSaved.end());
    };

    genomeParametersWrite(dirOut1+"/genomeParameters.txt", P, ERROR_OUT, *this);
    pGe.gFastaFiles=genomeFastaFilesSaved;

    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain    << timeMonthDayTime(rawTime) <<" ... writing Genome to disk ...\n" <<flush;
    *P.inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... writing Genome to disk ...\n" <<flush;
    writeGenomeSequence(dirOut1);

    time(&rawTime);
    P.inOut->logMain    << timeMonthDayTime(rawTime) <<" ... writing Suffix Array to disk ...\n" <<flush;
    *P.inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... writing Suffix Array to disk ...\n" <<flush;
    ofstream &saOut = ofstrOpen(dirOut1+"/SA", ERROR_OUT, P);
    fstreamWriteBig(saOut, (char*) SA.charArray, (streamsize) SA.lengthByte, dirOut1+"/SA", ERROR_OUT, P);
    saOut.close();

    time(&rawTime);
    P.inOut->logMain    << timeMonthDayTime(rawTime) <<" ... writing SAindex to disk\n" <<flush;
    *P.inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... writing SAindex to disk\n" <<flush;
    ofstream &saIndexOut = ofstrOpen(dirOut1+"/SAindex", ERROR_OUT, P);
    fstreamWriteBig(saIndexOut, (char*) &pGe.gSAindexNbases, sizeof(pGe.gSAindexNbases), dirOut1+"/SAindex", ERROR_OUT, P);
    fstreamWriteBig(saIndexOut, (char*) genomeSAindexStart, sizeof(genomeSAindexStart[0])*(pGe.gSAindexNbases+1), dirOut1+"/SAindex", ERROR_OUT, P);
    fstreamWriteBig(saIndexOut, SAi.charArray, SAi.lengthByte, dirOut1+"/SAindex", ERROR_OUT, P);
    saIndexOut.close();

    time(&rawTime);
    P.inOut->logMain    << timeMonthDayTime(rawTime) << " ..... finished writing updated genome index\n" <<flush;
    *P.inOut->logStdOut << timeMonthDayTime(rawTime) << " ..... finished writing updated genome index\n" <<flush;
};
