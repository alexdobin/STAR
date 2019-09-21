#include "sjdbLoadFromFiles.h"
#include "sjdbLoadFromStream.h"
#include "ErrorWarning.h"
#include "TimeFunctions.h"

void sjdbLoadFromFiles(Parameters &P, SjdbClass &sjdbLoci) {

    if (P.pGe.sjdbFileChrStartEnd.at(0)!="-") {
        for (uint ifile=0;ifile<P.pGe.sjdbFileChrStartEnd.size(); ifile++) {
            ifstream sjdbStreamIn ( P.pGe.sjdbFileChrStartEnd.at(ifile).c_str() );
            if (sjdbStreamIn.fail()) {
                ostringstream errOut;
                errOut << "FATAL INPUT error, could not open input file pGe.sjdbFileChrStartEnd=" << P.pGe.sjdbFileChrStartEnd.at(ifile) <<"\n";
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            };

            sjdbLoadFromStream(sjdbStreamIn, sjdbLoci);

            sjdbLoci.priority.resize(sjdbLoci.chr.size(),10);
            
            time_t rawtime;
            time ( &rawtime );
            P.inOut->logMain << timeMonthDayTime(rawtime) << "   Loaded database junctions from the pGe.sjdbFileChrStartEnd file(s), total number of junctions:" << sjdbLoci.chr.size()<<"\n\n";
        };
    }; //if (P.pGe.sjdbFileChrStartEnd!="-")

};