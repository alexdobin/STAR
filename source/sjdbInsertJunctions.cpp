#include "sjdbInsertJunctions.h"
#include "sjdbLoadFromStream.h"
#include "sjdbLoadFromFiles.h"
#include "sjdbPrepare.h"
#include "ErrorWarning.h"
#include "loadGTF.h"
#include "sjdbBuildIndex.h"
#include "streamFuns.h"
#include "genomeParametersWrite.h"

void sjdbInsertJunctions(Parameters & P, Genome & mapGen, Genome & mapGen1, SjdbClass & sjdbLoci)
{
    time_t rawtime;

    if (mapGen.sjdbN>0 && sjdbLoci.chr.size()==0)
    {//load from the saved genome, only if the loading did not happen already (if sjdb insertion happens at the 1st pass, sjdbLoci will be populated
        ifstream & sjdbStreamIn = ifstrOpen(P.pGe.gDir+"/sjdbList.out.tab", ERROR_OUT, "SOLUTION: re-generate the genome in pGe.gDir=" + P.pGe.gDir, P);
        sjdbLoadFromStream(sjdbStreamIn, sjdbLoci);
        sjdbLoci.priority.resize(sjdbLoci.chr.size(),30);
        time ( &rawtime );
        P.inOut->logMain << timeMonthDayTime(rawtime) << "   Loaded database junctions from the generated genome " << P.pGe.gDir+"/sjdbList.out.tab" <<": "<<sjdbLoci.chr.size()<<" total junctions\n\n";
    };

    if (P.twoPass.pass2)
    {//load 1st pass new junctions
     //sjdbLoci already contains the junctions from before 1st pass
        ifstream sjdbStreamIn ( P.twoPass.pass1sjFile.c_str() );
        if (sjdbStreamIn.fail()) {
            ostringstream errOut;
            errOut << "FATAL INPUT error, could not open input file with junctions from the 1st pass=" << P.twoPass.pass1sjFile <<"\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };
        sjdbLoadFromStream(sjdbStreamIn, sjdbLoci);
        sjdbLoci.priority.resize(sjdbLoci.chr.size(),0);
        time ( &rawtime );
        P.inOut->logMain << timeMonthDayTime(rawtime) << "   Loaded database junctions from the 1st pass file: " << P.twoPass.pass1sjFile <<": "<<sjdbLoci.chr.size()<<" total junctions\n\n";
    } else
    {//loading junctions from GTF or tab or from the saved genome is only allowed at the 1st pass
     //at the 2nd pass these are already in the sjdbLoci

        if (P.pGe.sjdbFileChrStartEnd.at(0)!="-")
        {//load from junction files
            sjdbLoadFromFiles(P,sjdbLoci);
            sjdbLoci.priority.resize(sjdbLoci.chr.size(),10);
            time ( &rawtime );
            P.inOut->logMain << timeMonthDayTime(rawtime) << "   Loaded database junctions from the pGe.sjdbFileChrStartEnd file(s), " << sjdbLoci.chr.size()<<" total junctions\n\n";
        };

        if (P.pGe.sjdbGTFfile!="-")
        {//load from GTF
            loadGTF(sjdbLoci, P, P.sjdbInsert.outDir, mapGen);
            sjdbLoci.priority.resize(sjdbLoci.chr.size(),20);
            time ( &rawtime );
            P.inOut->logMain << timeMonthDayTime(rawtime) << "   Loaded database junctions from the GTF file: " << P.pGe.sjdbGTFfile<<": "<<sjdbLoci.chr.size()<<" total junctions\n\n";
        };
    };

    char *Gsj=new char [2*mapGen.sjdbLength*sjdbLoci.chr.size()*(P.var.yes ? 2:1)+1];//array to store junction sequences, will be filled in sjdbPrepare
    sjdbPrepare (sjdbLoci, P, mapGen.chrStart[mapGen.nChrReal], P.sjdbInsert.outDir, mapGen, Gsj);//mapGen.nGenome - change when replacing junctions
    time ( &rawtime );
    P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished preparing junctions" <<endl;

    if (mapGen.sjdbN>P.limitSjdbInsertNsj)
    {
        ostringstream errOut;
        errOut << "Fatal LIMIT error: the number of junctions to be inserted on the fly ="<<mapGen.sjdbN<<" is larger than the limitSjdbInsertNsj="<<P.limitSjdbInsertNsj<<"\n";                errOut << "Fatal LIMIT error: the number of junctions to be inserted on the fly ="<<mapGen.sjdbN<<" is larger than the limitSjdbInsertNsj="<<P.limitSjdbInsertNsj<<"\n";
        errOut << "SOLUTION: re-run with at least --limitSjdbInsertNsj "<<mapGen.sjdbN<<"\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };
    //insert junctions into the genome and SA and SAi
    sjdbBuildIndex (P, Gsj, mapGen.G, mapGen.SA, (P.twoPass.pass2 ? mapGen.SApass2 : mapGen.SApass1), mapGen.SAi, mapGen, mapGen1);
    delete [] Gsj; //junction sequences have been added to G
    time ( &rawtime );
    P.inOut->logMain     << timeMonthDayTime(rawtime) << " ..... finished inserting junctions into genome" <<endl;

    if (P.pGe.sjdbInsertSave=="All")
    {//save and copy all genome files into sjdbInsert.outDir, except those created above
        if (P.pGe.gDir != P.sjdbInsert.outDir)
        {
            copyFile(P.pGe.gDir+"/chrName.txt", P.sjdbInsert.outDir+"/chrName.txt");
            copyFile(P.pGe.gDir+"/chrStart.txt", P.sjdbInsert.outDir+"/chrStart.txt");
            copyFile(P.pGe.gDir+"/chrNameLength.txt", P.sjdbInsert.outDir+"/chrNameLength.txt");
            copyFile(P.pGe.gDir+"/chrLength.txt", P.sjdbInsert.outDir+"/chrLength.txt");
        };

        genomeParametersWrite(P.sjdbInsert.outDir+("/genomeParameters.txt"), P, ERROR_OUT, mapGen);

        ofstream & genomeOut = ofstrOpen(P.sjdbInsert.outDir+"/Genome",ERROR_OUT, P);
        fstreamWriteBig(genomeOut,mapGen.G,mapGen.nGenome,P.sjdbInsert.outDir+"/Genome",ERROR_OUT,P);
        genomeOut.close();

        ofstream & saOut = ofstrOpen(P.sjdbInsert.outDir+"/SA",ERROR_OUT, P);
        fstreamWriteBig(saOut,(char*) mapGen.SA.charArray, (streamsize) mapGen.SA.lengthByte, P.sjdbInsert.outDir+"/SA",ERROR_OUT,P);
        saOut.close();

        ofstream & saIndexOut = ofstrOpen(P.sjdbInsert.outDir+"/SAindex",ERROR_OUT, P);
        fstreamWriteBig(saIndexOut, (char*) &P.pGe.gSAindexNbases, sizeof(P.pGe.gSAindexNbases),P.sjdbInsert.outDir+"/SAindex",ERROR_OUT,P);
        fstreamWriteBig(saIndexOut, (char*) mapGen.genomeSAindexStart, sizeof(mapGen.genomeSAindexStart[0])*(P.pGe.gSAindexNbases+1),P.sjdbInsert.outDir+"/SAindex",ERROR_OUT,P);
        fstreamWriteBig(saIndexOut,  mapGen.SAi.charArray, mapGen.SAi.lengthByte,P.sjdbInsert.outDir+"/SAindex",ERROR_OUT,P);
        saIndexOut.close();
    };

    //re-calculate genome-related parameters
    P.winBinN = mapGen.nGenome/(1LLU << P.winBinNbits)+1;
};
