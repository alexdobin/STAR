#include "Genome.h"
#include "ErrorWarning.h"

void Genome::validateGenomeInsertAnnotations()
{
    if (pGe.sjdbGTFfile=="-") {
        return;
    };

    if (sjdbOverhang==0) {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: --runMode genomeInsert with --sjdbGTFfile requires --sjdbOverhang > 0\n";
        errOut << "SOLUTION: use the same --sjdbOverhang value as the base index, or regenerate the base index with splice junction annotations\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    };

    if (chrNameIndex.size()==0) {
        for (uint64 ii=0; ii<nChrReal; ii++) {
            chrNameIndex[chrName[ii]]=ii;
        };
    };

    ifstream gtfIn(pGe.sjdbGTFfile.c_str());
    if (gtfIn.fail()) {
        ostringstream errOut;
        errOut << "FATAL error, could not open file pGe.sjdbGTFfile=" << pGe.sjdbGTFfile <<"\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    uint64 exonN=0;
    uint64 lineN=0;
    while (gtfIn.good()) {
        string oneLine, chr1, ddd2, featureType;
        getline(gtfIn, oneLine);
        lineN++;

        istringstream oneLineStream(oneLine);
        oneLineStream >> chr1 >> ddd2 >> featureType;
        if (chr1.size()==0 || chr1.substr(0,1)=="#" || featureType!=pGe.sjdbGTFfeatureExon) {
            continue;
        };

        exonN++;
        if (pGe.sjdbGTFchrPrefix!="-") {
            chr1=pGe.sjdbGTFchrPrefix + chr1;
        };

        auto chrIt=chrNameIndex.find(chr1);
        if (chrIt==chrNameIndex.end()) {
            ostringstream errOut;
            errOut << "EXITING because of fatal INPUT FILE error: --runMode genomeInsert could not find GTF chromosome " << chr1 << " in the expanded genome\n";
            errOut << "Offending file: " << pGe.sjdbGTFfile << ", line " << lineN << "\n";
            errOut << "SOLUTION: make sure --sjdbGTFfile contains annotations for sequences supplied with --genomeFastaFiles\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };

        if (chrIt->second<genomeInsertChrIndFirst) {
            ostringstream errOut;
            errOut << "EXITING because of fatal INPUT FILE error: --runMode genomeInsert expects --sjdbGTFfile to contain only annotations for inserted sequences\n";
            errOut << "Found exon on base genome chromosome " << chr1 << " in file " << pGe.sjdbGTFfile << ", line " << lineN << "\n";
            errOut << "SOLUTION: provide a GTF for the added sequences only; existing annotations from --genomeDir are preserved automatically\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };
    };

    if (exonN==0) {
        ostringstream errOut;
        errOut << "Fatal INPUT FILE error, no exon lines in the GTF file: " << pGe.sjdbGTFfile <<"\n";
        errOut << "Solution: check the formatting of the GTF file, or use --sjdbGTFfeatureExon if exons are marked with a different word.\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };
};
