#include "twoPassRunPass1.h"
#include "mapThreadsSpawn.h"
#include "ReadAlignChunk.h"
#include "Stats.h"
#include "GlobalVariables.h"
#include "outputSJ.h"
#include "sjdbInsertJunctions.h"

void twoPassRunPass1(Parameters &P, Genome &genomeMain, Transcriptome *transcriptomeMain, SjdbClass &sjdbLoci) //2-pass: 1st pass
{
    if (!P.twoPass.yes)
        return;

    //re-define P and genomeMain for the pass1
    Genome genomeMain1=genomeMain;

    Parameters P1=P;
    //turn off unnecessary calculations
    P1.outSAMtype[0]="None";
    P1.outSAMbool=false;
    P1.outBAMunsorted=false;
    P1.outBAMcoord=false;

    P1.pCh.segmentMin=0;

    P1.quant.yes=false;
    P1.quant.trSAM.yes=false;
    P1.quant.trSAM.bamYes=false;
    P1.quant.geneFull.yes=false;
    P1.quant.geCount.yes=false;
    P1.quant.gene.yes=false;

    P1.outSAMunmapped.within=false;

    P1.outFilterBySJoutStage=0;

    P1.outReadsUnmapped="None";

    P1.outFileNamePrefix=P.twoPass.dir;

    P1.readMapNumber=min(P.twoPass.pass1readsN, P.readMapNumber);
//         P1.inOut->logMain.open((P1.outFileNamePrefix + "Log.out").c_str());

    P1.wasp.outputMode="None"; //no WASP filtering on the 1st pass
    P1.wasp.yes = false;
    P1.wasp.SAMtag = false;
    
    P1.pSolo.type=P1.pSolo.SoloTypes::None; //no solo in the first pass
    P1.pSolo.yes = false;
    
    g_statsAll.resetN();
    time(&g_statsAll.timeStartMap);
    P.inOut->logProgress << timeMonthDayTime(g_statsAll.timeStartMap) <<"\tStarted 1st pass mapping\n" <<flush;
    *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeStartMap) << " ..... started 1st pass mapping\n" <<flush;

    //no genome conversion for 1st pass
    P1.pGe.transform.outYes = false;
    P1.pGe.transform.outQuant = false;
    P1.pGe.transform.outSAM = false;
    P1.pGe.transform.outSJ = false;
    
    auto convYes = genomeMain.genomeOut.convYes;
    auto gOut = genomeMain.genomeOut.g;

    genomeMain.genomeOut.convYes = false;
    genomeMain.genomeOut.g = &genomeMain;

    //run mapping for Pass1
    ReadAlignChunk *RAchunk1[P.runThreadN];
    for (int ii=0;ii<P1.runThreadN;ii++) {
        RAchunk1[ii]=new ReadAlignChunk(P1, genomeMain, transcriptomeMain, ii);
    };
    mapThreadsSpawn(P1, RAchunk1);
    outputSJ(RAchunk1,P1); //collapse and output junctions
//         for (int ii=0;ii<P1.runThreadN;ii++) {
//             delete [] RAchunk[ii];
//         };

    //back to requested genome conversions
    genomeMain.genomeOut.convYes = convYes;
    genomeMain.genomeOut.g = gOut;

    time_t rawtime; time (&rawtime);
    P.inOut->logProgress << timeMonthDayTime(rawtime) <<"\tFinished 1st pass mapping\n";
    *P.inOut->logStdOut << timeMonthDayTime(rawtime) << " ..... finished 1st pass mapping\n" <<flush;
    ofstream logFinal1 ( (P.twoPass.dir + "/Log.final.out").c_str());
    g_statsAll.reportFinal(logFinal1);

    P.twoPass.pass2=true;//starting the 2nd pass
    P.twoPass.pass1sjFile=P.twoPass.dir+"/SJ.out.tab";

    sjdbInsertJunctions(P, genomeMain, genomeMain1, sjdbLoci);

    //reopen reads files
    P.closeReadsFiles();
    P.openReadsFiles();
};