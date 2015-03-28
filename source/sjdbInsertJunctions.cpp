#include "sjdbInsertJunctions.h"
#include "sjdbLoadFromStream.h"
#include "sjdbLoadFromFiles.h"
#include "sjdbPrepare.h"
#include "ErrorWarning.h"
// #include "streamFuns.h"
#include "SjdbClass.h"
#include "loadGTF.h"
#include "sjdbBuildIndex.h"
#include "streamFuns.h"

void sjdbInsertJunctions(Parameters *P, Parameters *P1, Genome &genome) {
        
    SjdbClass sjdbLoci;        
    time_t rawtime;

    //load 1st pass junctions
    if (P->twoPass.pass1sjFile.size()>0)
    {
        ifstream sjdbStreamIn ( P->twoPass.pass1sjFile.c_str() );   
        if (sjdbStreamIn.fail()) {
            ostringstream errOut;
            errOut << "FATAL INPUT error, could not open input file with junctions from the 1st pass=" << P->twoPass.pass1sjFile <<"\n";
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
        };
        sjdbLoadFromStream(sjdbStreamIn, sjdbLoci);
        time ( &rawtime );
        P->inOut->logMain << timeMonthDayTime(rawtime) << "   Loaded database junctions from the 1st pass file: " << P->twoPass.pass1sjFile <<": "<<sjdbLoci.chr.size()<<" total junctions\n\n";
    };
    
    //load from junction files
    if (P->sjdbFileChrStartEnd.at(0)!="-")
    {
        sjdbLoadFromFiles(P,sjdbLoci);
        P->inOut->logMain << timeMonthDayTime(rawtime) << "   Loaded database junctions from the sjdbFileChrStartEnd file(s), " << sjdbLoci.chr.size()<<" total junctions\n\n";        
    };
    
    //load from GTF
    if (P->sjdbGTFfile!="-")
    {
        loadGTF(sjdbLoci, P, P->genomeDirOut);
        P->inOut->logMain << timeMonthDayTime(rawtime) << "   Loaded database junctions from the GTF file: " << P->sjdbGTFfile<<": "<<sjdbLoci.chr.size()<<" total junctions\n\n";
    };
    
    //load from the already generated genome
    if (P->sjdbN>0)
    {
        ifstream sjdbStreamIn;
        ifstrOpen(P->genomeDir+"/sjdbList.out.tab", "ERROR_012003", "SOLUTION: re-generate the genome in genomeDir=" + P->genomeDir, P, sjdbStreamIn);
        sjdbLoadFromStream(sjdbStreamIn, sjdbLoci);
        time ( &rawtime );
        P->inOut->logMain << timeMonthDayTime(rawtime) << "   Loaded database junctions from the generated genome " << P->genomeDir+"/sjdbList.out.tab" <<": "<<sjdbLoci.chr.size()<<" total junctions\n\n";
    };

    sjdbPrepare (sjdbLoci, P, genome.G, P->chrStart[P->nChrReal], P->twoPass.dir);//P->nGenome - change when replacing junctions
    time ( &rawtime );
    P->inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished preparing junctions" <<endl;

    //insert junctions into the genome and SA and SAi
    sjdbBuildIndex (P, P1, genome.G, genome.SA, genome.SA2, genome.SAi);
    time ( &rawtime ); 
    *P->inOut->logStdOut  << timeMonthDayTime(rawtime) << " ..... Finished inserting 1st pass junctions into genome" <<endl;
    //re-calculate genome-related parameters
    P->winBinN = P->nGenome/(1LLU << P->winBinNbits)+1;
};