#include "sjdbLoadFromFiles.h"
#include "sjdbLoadFromStream.h"
#include "ErrorWarning.h"

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

            P.inOut->logMain << "Loaded database junctions from file: " << P.pGe.sjdbFileChrStartEnd.at(ifile) <<", total number of junctions: "<<sjdbLoci.chr.size()<<" junctions\n\n";
        };
    }; //if (P.pGe.sjdbFileChrStartEnd!="-")

};