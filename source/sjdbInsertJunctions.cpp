#include "sjdbInsertJunctions.h"
#include "sjdbLoadFromStream.h"
#include "sjdbLoadFromFiles.h"
#include "sjdbPrepare.h"
#include "ErrorWarning.h"
// #include "streamFuns.h"
#include "SjdbClass.h"
#include "loadGTF.h"
#include "sjdbBuildIndex.h"

void sjdbInsertJunctions(Parameters *P, Genome &genome) {
        
    SjdbClass sjdbLoci;        
    time_t rawtime;

    //load 1st pass junctions
    if (P->twopassSJpass1file.size()>0)
    {
        ifstream sjdbStreamIn ( P->twopassSJpass1file.c_str() );   
        if (sjdbStreamIn.fail()) {
            ostringstream errOut;
            errOut << "FATAL INPUT error, could not open input file with junctions from the 1st pass=" << P->twopassSJpass1file <<"\n";
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
        };
        sjdbLoadFromStream(sjdbStreamIn, sjdbLoci);
        time ( &rawtime );
        P->inOut->logMain << timeMonthDayTime(rawtime) << "   Loaded database junctions from the 1st pass file: " << P->twopassSJpass1file <<": "<<sjdbLoci.chr.size()<<" total junctions\n\n";
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

    sjdbPrepare (sjdbLoci, P, genome.G, P->nGenome, P->twopassDir);//P->nGenome - change when replacing junctions
    time ( &rawtime );
    P->inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished preparing junctions" <<endl;

    //insert junctions into the genome and SA and SAi
    sjdbBuildIndex (P, genome.G, genome.SA, genome.SA2, genome.SAi);
    time ( &rawtime ); 
    *P->inOut->logStdOut  << timeMonthDayTime(rawtime) << " ..... Finished inserting 1st pass junctions into genome" <<endl;
    //re-calculate genome-related parameters
    P->winBinN = P->nGenome/(1LLU << P->winBinNbits)+1;
};